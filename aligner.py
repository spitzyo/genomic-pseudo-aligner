####################         IMPORTS          ####################
from collections import Counter, defaultdict
from typing import List, Dict, Any, Tuple
from kmers import Kmer
from extcoverage import Coverage
from extvartrack import VariantTracker

####################         CONSTANTS        ####################
MAX_MISMATCH_THRESHOLD = 0.1  # 10% maximum allowed error rate
_COMPLEMENT = str.maketrans("ACGTN", "TGCAN")

####################         FUNCTIONS        ####################

def reverse_complement(sequence: str) -> str:
    """Returns the reverse complement of a DNA sequence (N stays N)."""
    return sequence.translate(_COMPLEMENT)[::-1]

def extract_kmers(sequence, quality, k,
                  min_kmer_quality=None) -> List[Tuple[int, str]]:
    """Returns (offset in sequence, k-mer) pairs. k-mers containing N are
    skipped, and so are k-mers whose mean quality is below min_kmer_quality
    (EXTQUALITY). The offset is kept so skipped k-mers don't shift the
    positions of the ones after them."""
    quality_ok = len(quality) == len(sequence)
    kmers = []
    for i in range(len(sequence) - k + 1):
        k_mer = sequence[i:i + k]
        if "N" in k_mer:
            continue
        if min_kmer_quality is not None:
            # a read whose quality and sequence lengths differ has quality 0
            kmer_quality = (sum(quality[i:i + k]) / k) if quality_ok else 0.0
            if kmer_quality < min_kmer_quality:
                continue
        kmers.append((i, k_mer))
    return kmers

def aligned_span(offset, read_len, genome_len) -> Tuple[int, int]:
    """Returns the [first, last) read indices that fall inside the genome
    when the read starts at reference position offset. The offset may be
    negative, and the read may run past the genome end."""
    first = max(0, -offset)
    last = min(read_len, genome_len - offset)
    return first, max(first, last)

####################         CLASSES          ####################

class Read:
    """
    This class holds all data related to a read: its identifier,
    full DNA sequence and a quality string. It also contains the mapping
    status of every read, and the genomes that it was mapped to (if any).
    """

    def __init__(self, identifier, sequence, quality_str, status="unmapped"):
        """
        Initialize a Read object.
        :param identifier: as given in the header of the input file.
        :param sequence: a full DNA sequence.
        :param quality_str: converted by the Phred33 format.
        :param status: the status of the read (default is "unmapped").
        """
        self.identifier = identifier
        self.sequence = sequence
        self.quality = quality_str

        # Status (unmapped/ambiguous/unique) would be decided upon execution:
        self.status = status
        self.mapped_genomes = [] # list of mapped Reference's identifiers
        # mapped genome id -> 0-based reference position of the read's first
        # base. Negative when the read hangs off the start of the genome.
        self.positions = {}
        # '+' if the read matched the reference as given, '-' if its reverse
        # complement did, None while unmapped/filtered.
        self.strand = None

    def oriented(self, strand):
        """Returns (sequence, quality) on the given strand. For '-' this is
        the reverse complement with the qualities reversed, so index i of
        both still describes the same base."""
        if strand == '-':
            return reverse_complement(self.sequence), self.quality[::-1]
        return self.sequence, self.quality

    def get_mean_quality(self) -> float:
        """This method returns the mean (ave) quality value for the read."""
        if len(self.quality) != len(self.sequence):
            return 0.0 # if quality score and seq length do not match
        return sum(self.quality) / len(self.quality)

    def get_kmer_quality(self, start_pos, k) -> float:
        """This method returns the mean (average) quality value for the kmer
        that begins in the start_pos, with a given k length."""
        if len(self.quality) != len(self.sequence):
            return 0.0
        quality_of_kmer = self.quality[start_pos:start_pos + k]
        return sum(quality_of_kmer) / len(quality_of_kmer)

    def get_kmers(self, k, min_kmer_quality=None) -> List[Kmer]:
        """This function extracts k-mers from this Read and returns a list.
        EXTQUALITY: kmers may be filtered by their quality scores."""
        return [Kmer(k_mer) for _, k_mer in
                extract_kmers(self.sequence, self.quality, k,
                              min_kmer_quality)]

class Aligner:
    """
    This class executes the pseudo-alignment algorithm, while using
    the KmerCollection to align every k-mer Read against the Reference k-mers.
    The final alignment is determined by the num of matches between the k-mers.
    """

    def __init__(self, kmer_collection):
        self._kmer_collection = kmer_collection # kmers of all ref genomes
        # alignment_summary connected between a genome id -> mapping stats:
        self.alignment_summary = {}

        # These attributes are relevant for extensions (not the basic algo)
        self._coverage = None # EXTCOVERAGE
        self._quality_stats = {'filtered_quality_reads': 0,
                                'filtered_quality_kmers': 0,
                                'filtered_hr_kmers': 0} # EXTQUALITY
        self._variant_tracker = None # EXTVARTRACK

    def enable_coverage(self, genome_lens: Dict[str, int], min_coverage = 1):
        """This method enables coverage: EXTCOVERAGE, introduces a new attr."""
        self._coverage = Coverage(genome_lens, min_coverage)

    def _apply_quality(self, sequence, quality, k, min_kmer_quality,
                       max_genomes) -> Tuple[List[Tuple[int, str]], int, int]:
        """This EXTQUALITY helper method extracts the (offset, k-mer) pairs of
        one read orientation, applying the k-mer filters given by the user.
        Returns (k-mers, num low quality k-mers, num high-repeat k-mers)."""
        kmers = extract_kmers(sequence, quality, k, min_kmer_quality)
        low_quality = 0
        if min_kmer_quality: # keeping track of the filtered kmers
            # calcs the maximal number of kmers that can be made from Read:
            max_kmer_num = max(0, len(sequence) - k + 1)
            low_quality = max_kmer_num - len(kmers)
            # len(kmers) represents the num of kmers after filtering with MKQ
        high_repeat = 0
        if max_genomes:
            filtered_kmers = []
            for offset, kmer_seq in kmers:
                ref_kmer = self._kmer_collection.get_kmer(kmer_seq)
                if ref_kmer is not None:
                    num_genomes = len(ref_kmer.get_genomes())
                    if num_genomes <= max_genomes: # filtering
                        filtered_kmers.append((offset, kmer_seq))
            high_repeat = len(kmers) - len(filtered_kmers)
            kmers = filtered_kmers
        return kmers, low_quality, high_repeat

    def enable_vartracking(self, min_quality, min_coverage, min_af=None):
        """This method enables variant tracking: EXTVARTRACK,
        as it introduces a new attribute to be used."""
        self._variant_tracker = VariantTracker(min_quality, min_coverage,
                                               min_af)

    def align_read(self, read, k, unique_diff=1, total_diff=1,
                   min_read_quality=None, min_kmer_quality=None,
                   max_genomes=None) -> None:
        """
        This method aligns a single Read to the ref genomes k-mers.
        :param read: a Read instance (to be aligned).
        :param k: length of k-mer.
        :param unique_diff: threshold for distinguishing unique vs ambiguous.
        :param total_diff: max diff in total k-mers to maintain unique status.
        :param min_read_quality: filter reads below mean quality (EXTQUALITY).
        :param min_kmer_quality: filter k-mers below quality (EXTQUALITY).
        :param max_genomes: filter k-mers from too many genomes (EXTQUALITY).
        :return: should set the status of the Read as an attr and return None.
        """
        read.mapped_genomes = [] # reset, in case the read is aligned again
        read.positions = {}
        read.strand = None

        # EXTQUALITY: read-level filter
        if min_read_quality and read.get_mean_quality() < min_read_quality:
            self._quality_stats['filtered_quality_reads'] += 1
            read.status = 'filtered'
            return None # and continue to the next Read

        # Strand: count k-mer hits for the read as given ('+') and for its
        # reverse complement ('-'). Keep the orientation with more specific
        # hits, then more total hits; a full tie keeps '+'.
        best = None
        for strand in ('+', '-'):
            sequence, quality = read.oriented(strand)
            kmers, low_quality, high_repeat = self._apply_quality(
                sequence, quality, k, min_kmer_quality, max_genomes)
            counts = self._count_kmers_with_pos(kmers)
            score = (sum(counts[0].values()), sum(counts[1].values()))
            if best is None or score > best[0]:
                best = (score, strand, sequence, quality, kmers,
                        low_quality, high_repeat, counts)
        (_, strand, sequence, quality, kmers, low_quality, high_repeat,
         (specific_counts, mapped_genomes, offset_votes)) = best

        self._quality_stats['filtered_quality_kmers'] += low_quality
        self._quality_stats['filtered_hr_kmers'] += high_repeat
        if not kmers:
            read.status = 'filtered'
            return None

        # Basic algorithm:
        # No specific k-mers -> unmapped
        if not specific_counts:
            # this case is also relevant for a read length < k size
            read.status = 'unmapped'
            self._accumulate_status('unmapped', [])
            return None

        # Compare specific k-mer counts
        sorted_specific = sorted(specific_counts.items(),
                                 key=lambda x: x[1], reverse=True)
        top_genome = sorted_specific[0][0]
        top_specific_count = sorted_specific[0][1]

        if len(sorted_specific) > 1:
            second_specific_count = sorted_specific[1][1]
        else:
            second_specific_count = 0

        if top_specific_count - second_specific_count >= unique_diff:
            # Validation with total k-mer counts
            mapped_count = mapped_genomes[top_genome]
            max_count = max(mapped_genomes.values()) if mapped_genomes else 0

            # Check if the difference is too large
            if max_count - mapped_count > total_diff:
                read.status = 'ambiguous'
                # Include all genomes with high enough counts
                qualified_genomes = []
                for gen_id, count in mapped_genomes.items():
                    if count >= mapped_count:
                        qualified_genomes.append(gen_id)
                read.mapped_genomes = qualified_genomes
            else:
                read.status = 'unique'
                read.mapped_genomes = [top_genome]
        else:
            read.status = 'ambiguous'
            read.mapped_genomes = list(specific_counts.keys())

        verified_genomes = []
        for genome_id in read.mapped_genomes:
            genome = self._kmer_collection.get_genome(genome_id)
            votes = offset_votes.get(genome_id)
            if genome is None or not votes:
                continue
            start_pos, mismatch_rate = self._choose_offset(sequence, genome,
                                                           votes)

            # Filter using the CONSTANT
            if mismatch_rate <= MAX_MISMATCH_THRESHOLD:
                verified_genomes.append(genome_id)
                read.positions[genome_id] = start_pos
                # Only call variants if verified and unique
                if self._variant_tracker and read.status == 'unique':
                    self._process_variants(sequence, quality, genome_id,
                                           start_pos)

        if verified_genomes:
            read.mapped_genomes = verified_genomes
            read.strand = strand
        else:
            # If it failed the Hamming check, mark as unmapped
            read.status = 'unmapped'
            read.mapped_genomes = []
            read.positions = {}
            self._accumulate_status('unmapped', [])
            return

        # update statistics
        self._accumulate_status(read.status, read.mapped_genomes)

        # EXTCOVERAGE: add coverage
        if self._coverage:
            # add the read coverage attr to the Coverage inst (of this Aligner)
            if read.status != 'unmapped': # either unique or ambiguous
                for gen_id in read.mapped_genomes:
                    # every reference base the read covers at its offset,
                    # clipped to the genome
                    start_pos = read.positions[gen_id]
                    genome = self._kmer_collection.get_genome(gen_id)
                    first, last = aligned_span(start_pos, len(read.sequence),
                                               genome.total_bases)
                    self._coverage.add_read_coverage(
                        read, gen_id,
                        range(start_pos + first, start_pos + last))

    def _count_kmers_with_pos(self, kmers: List[Tuple[int, str]]):
        """This method counts specific and total kmers for every genome.
        It also counts, per genome, how many k-mer hits vote for each
        read start position (reference pos of the k-mer minus its offset
        in the read)."""

        mapped_genomes = {}  # For total k-mer counts
        specific_counts = {}  # For specific k-mer counts
        offset_votes = defaultdict(Counter)  # genome id -> {offset: votes}

        for i, kmer_seq in kmers:
            ref_kmer = self._kmer_collection.get_kmer(kmer_seq)
            if not ref_kmer:
                continue

            # Handle specific k-mers
            if ref_kmer.is_specific():
                genome = ref_kmer.get_genomes()[0]
                genome_id = genome.identifier
                if genome_id not in specific_counts:
                    specific_counts[genome_id] = 0
                specific_counts[genome_id] += 1

            # Count all k-mers
            for genome in ref_kmer.get_genomes():
                genome_id = genome.identifier
                if genome_id not in mapped_genomes:
                    mapped_genomes[genome_id] = 0
                mapped_genomes[genome_id] += 1

                for pos in ref_kmer.get_positions(genome):
                    offset_votes[genome_id][pos - i] += 1

        return specific_counts, mapped_genomes, offset_votes

    def _choose_offset(self, sequence, genome, votes) -> Tuple[int, float]:
        """Picks the read start position with the most k-mer votes. Ties are
        broken by the lower mismatch rate, then by the smaller position.
        Returns (position, mismatch rate)."""
        top_votes = max(votes.values())
        best = None
        for offset in sorted(o for o, n in votes.items() if n == top_votes):
            mismatches, overlap = self._calculate_mismatches(sequence, genome,
                                                             offset)
            rate = mismatches / overlap if overlap else 1.0
            if best is None or rate < best[1]:
                best = (offset, rate)
        return best

    def _accumulate_status(self, status, mapped_genomes: List[str]) -> None:
        """
        This helper method accumulates the alignment status into the attribute.
        :param status: the status of the Read (either unambiguous or unique).
        :param mapped_genomes: the mapped genomes of a Read.
        """
        for genome_identifier in mapped_genomes:
            # if genome is not in summary:
            if genome_identifier not in self.alignment_summary:
                self.alignment_summary[genome_identifier] = {'unique': 0,
                                                               'ambiguous': 0,
                                                               'unmapped': 0,}
            self.alignment_summary[genome_identifier][status] += 1

    def align_reads(self, reads, k, unique_diff=1, total_diff=1) -> None:
        """
        Aligns a Reads list to  k-mer collection (of ref genomes).
        Params are the same as in prev method align_reads().
        """
        for read in reads:
            self.align_read(read, k, unique_diff, total_diff)
        return None

    def add_genome_to_summary(self, genome_id: str) -> None:
        """This method adds a genome to alignment summary even when unmapped"""
        if genome_id not in self.alignment_summary:
            self.alignment_summary[genome_id] = {'unique': 0,
                                                   'ambiguous':0,
                                                   'unmapped': 0}

    # EXTCOVERAGE
    def dump_coverage(self, selected_genomes: List[str] = None,
                      full_coverage: bool = False) -> [Dict[str, int],
                                                       None]:
        """This COVERAGE EXTENSION method calls the dump coverage
        method of Coverage class to dump it to the console."""
        if self._coverage is None:
            return  # no coverage if the attr is None
        return self._coverage.dump_coverage(selected_genomes,
                                             full_coverage, False)

    # EXTQUALITY
    def get_quality_stats(self) -> dict:
        """This QUALITY EXTENSION method returns Reads' quality stats."""
        return self._quality_stats if self._quality_stats else {}
        # in case there are no quality statistics, return empty dict

    # EXTVARTRACK
    def _process_variants(self, read_seq, read_quality, genome_id,
                          start_pos):
        """This method performs variant detection after the alignment,
        by directly comparing the read to the reference at start_pos.
        read_seq and read_quality are in reference orientation (reverse
        complemented / reversed for minus-strand reads).
        Bases hanging off either end of the genome are ignored. Every
        processed read counts toward the depth of the positions it covers."""

        genome = self._kmer_collection.get_genome(genome_id)
        if not genome:
            return

        first, last = aligned_span(start_pos, len(read_seq),
                                   genome.total_bases)
        self._variant_tracker.add_aligned_read(
            genome_id, start_pos + first, read_seq[first:last],
            read_quality[first:last],
            genome.seq[start_pos + first:start_pos + last])

    def dump_variants(self, selected_genomes=None) -> [Dict[str, Any],
                                                       None]:
        """This VARIANT TRACKING EXTENSION method calls the dump variants
        method of the relevant class to return the final VarTrack data."""
        if self._variant_tracker is None:
            return
        return self._variant_tracker.dump_variants(selected_genomes)

    @staticmethod
    def _calculate_mismatches(sequence, genome,
                              start_pos) -> Tuple[int, int]:
        """Helper to calculate Hamming Distance between the read and the
        genome at start_pos, over the bases that overlap the genome.
        Returns (mismatches, number of overlapping bases)."""
        first, last = aligned_span(start_pos, len(sequence),
                                   genome.total_bases)
        ref_seq = genome.seq[start_pos + first:start_pos + last]
        mismatches = sum(1 for r, g in zip(sequence[first:last], ref_seq)
                         if r != g)
        return mismatches, last - first
