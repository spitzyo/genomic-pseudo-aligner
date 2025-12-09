####################         IMPORTS          ####################
from collections import defaultdict
from typing import List, Dict, Any
from kmers import Kmer
from extcoverage import Coverage
from extvartrack import VariantTracker


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

    def get_mean_quality(self) -> float:
        """This method returns the mean (ave) quality value for the read."""
        if len(self.quality) != len(self.sequence):
            return 0.0 # if quality score and seq length do not match
        return sum(self.quality) / len(self.quality)

    def get_kmer_quality(self, start_pos, k) -> float:
        """This method returns the mean (average) quality value for the kmer
        that begins in the start_pos, with a given k length."""
        if len(self.quality) != len(self.sequence):
            return 0.0 # quality score and seq len do not match
        quality_of_kmer = self.quality[start_pos:start_pos + k]
        return sum(quality_of_kmer) / len(quality_of_kmer)

    def get_kmers(self, k, min_kmer_quality=None) -> List[Kmer]:
        """This function extracts k-mers from this Read and returns a list.
        EXTQUALITY: kmers may be filtered by their quality scores."""
        kmers = []
        for i in range(len(self.sequence)-k+1):
            k_mer = self.sequence[i:i+k]
            if "N" in k_mer:
                continue # do not add kmers with "N"

            # EXTQUALITY
            if min_kmer_quality is not None:
                kmer_quality = self.get_kmer_quality(i, k)
                if kmer_quality < min_kmer_quality:
                    continue # do not add a kmer below threshold

            kmer_instance = Kmer(k_mer)
            kmers.append(kmer_instance)
        return kmers

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

    def _apply_quality(self, read, k, min_read_quality,
                        min_kmer_quality, max_genomes) -> [None, List]:
        """This EXTQUALITY helper method applies quality filters (as given
        by user) and returns only the filters kmers to be aligned."""
        if min_read_quality and read.get_mean_quality() < min_read_quality:
            self._quality_stats['filtered_quality_reads'] += 1
            return None # keeping track in stats of filtered reads
        kmers = read.get_kmers(k, min_kmer_quality) # get all filtered kmers
        if min_kmer_quality: # keeping track of the filtered kmers
            # calcs the maximal number of kmers that can be made from Read:
            max_kmer_num = len(read.sequence)-k+1
            self._quality_stats['filtered_quality_kmers'] += \
                (max_kmer_num - len(kmers))
            # len(kmers) represents the num of kmers after filtering with MKQ
        if max_genomes:
            num_of_kmers = len(kmers)
            filtered_kmers = []
            for kmer in kmers:
                kmer_seq = kmer.sequence
                ref_kmer = self._kmer_collection.get_kmer(kmer_seq)
                if ref_kmer is not None:
                    num_genomes = len(ref_kmer.get_genomes())
                    if num_genomes <= max_genomes: # filtering
                        filtered_kmers.append(kmer) # add to filtered kmers
            self._quality_stats['filtered_hr_kmers'] += \
                (num_of_kmers - len(filtered_kmers)) # update filtering stats
            if not filtered_kmers:
                return None # and Read would be marked as unmapped
            kmers = filtered_kmers # update the kmers after filter
        return kmers # after three possible filtering processes

    def enable_vartracking(self, min_quality, min_coverage):
        """This method enables variant tracking: EXTVARTRACK,
        as it introduces a new attribute to be used."""
        self._variant_tracker = VariantTracker(min_quality, min_coverage)

    def align_read(self, read, k, m=1, p=1, min_read_quality=None,
                   min_kmer_quality=None, max_genomes=None) -> None:
        """
        This method aligns a single Read to the ref genomes k-mers.
        :param read: a Read instance (to be aligned).
        :param k: length of k-mer.
        :param m: max diff in the num of k-mers for ambiguous mappings.
        :param p: max diff in total k-mers to maintain unique status.
        :param min_read_quality: filter reads below mean quality (EXTQUALITY).
        :param min_kmer_quality: filter k-mers below quality (EXTQUALITY).
        :param max_genomes: filter k-mers from too many genomes (EXTQUALITY).
        :return: should set the status of the Read as an attr and return None.
        """

        # EXTQUALITY
        if any([min_read_quality, min_kmer_quality, max_genomes]):
            kmers = self._apply_quality(read, k, min_read_quality,
                                         min_kmer_quality, max_genomes)
            if kmers is None: # Read filtered by EXTQUALITY
                read.status = 'filtered'
                return None # and continue to the next Read
        else: # if no quality filtering is required
            kmers = read.get_kmers(k)

        if not kmers: # EXTQUALITY: no kmers left after filtering
            read.status = 'filtered'
            return None # finish with this read

        specific_counts, mapped_genomes, genome_positions =\
            self._count_kmers_with_pos(kmers, k)  # helper func

        # EXTVARTRACK: process vars if tracker is on
        if self._variant_tracker and genome_positions:
            for genome_id, positions in genome_positions.items():
                if positions:
                    self._process_variants(read, genome_id, positions)

        # Basic algorithm:
        # Step 2.1: No specific k-mers -> unmapped
        if not specific_counts:
            # this case is also relevant for a read length < k size
            read.status = 'unmapped'
            self._accumulate_status('unmapped', [])
            return None

        # Step 2.2: Compare specific k-mer counts
        sorted_specific = sorted(specific_counts.items(),
                                 key=lambda x: x[1], reverse=True)
        top_genome = sorted_specific[0][0]
        top_specific_count = sorted_specific[0][1]

        if len(sorted_specific) > 1:
            second_specific_count = sorted_specific[1][1]
        else:
            second_specific_count = 0

        if top_specific_count - second_specific_count >= m:
            # Step 2.3: Validation with total k-mer counts
            mapped_count = mapped_genomes[top_genome]
            max_count = max(mapped_genomes.values()) if mapped_genomes else 0

            # Step 2.3.2.4: Check if the difference is too large
            if max_count - mapped_count > p:
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

        # update statistics
        self._accumulate_status(read.status, read.mapped_genomes)

        # EXTCOVERAGE: add coverage
        if self._coverage:
            # add the read coverage attr to the Coverage inst (of this Aligner)
            if read.status != 'unmapped': # either unique or ambiguous
                for gen_id in read.mapped_genomes:
                    self._coverage.add_read_coverage(read, gen_id,
                                                      genome_positions[gen_id])

    def _count_kmers_with_pos(self, kmers: List[Kmer], k):
        """This method counts specific and total kmers for every genome,
        and tracks positions of any matched k-mers here (genome_positions)."""

        mapped_genomes = {}  # For total k-mer counts
        specific_counts = {}  # For specific k-mer counts
        genome_positions = defaultdict(set) # default dict with an empty set

        for kmer in kmers:
            ref_kmer = self._kmer_collection.get_kmer(kmer.sequence)
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

                # getting all kmer positions in this genome:
                for pos in ref_kmer.get_positions(genome):
                    for i in range(k): # adding all pos for a k-long kmer
                        genome_positions[genome_id].add(pos + i)
        return specific_counts, mapped_genomes, genome_positions

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
            self.alignment_summary[genome_identifier][status] += 1 # update

    def align_reads(self, reads, k, m=1, p=1) -> None:
        """
        This method aligns a Reads list to  k-mer collection (of ref genomes).
        Params are the same as in prev method align_reads().
        """
        for read in reads:
            self.align_read(read, k, m, p)
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
    def _process_variants(self, read, genome_id, positions):
        """This method performs variant detection after the alignment,
        by directly comparing the sequences."""
        read_seq = read.sequence
        read_quality = read.quality

        # get the genome instance:
        genome = None
        for g in self._kmer_collection.get_all_genomes():
            if g.identifier == genome_id:
                genome = g
                break

        if not genome or not positions:
            return

        start_pos = min(positions)
        end_pos = min(start_pos + len(read_seq), genome.total_bases)
        ref_seq = genome.seq[start_pos:end_pos] # get ref sequence

        for i, (read_base, ref_base) in enumerate(zip(read_seq[:
        len(ref_seq)], ref_seq)):
            # iterating through bases simultaneously of ref and read seqs.
            if read_base != ref_base:  # diff between ref and read - found var
                pos = start_pos + i
                quality = read_quality[i]
                self._variant_tracker.add_direct_variant(genome_id, pos,
                                                          ref_base,
                                                          read_base,
                                                          quality)

    def dump_variants(self, selected_genomes=None) -> [Dict[str, Any],
                                                       None]:
        """This VARIANT TRACKING EXTENSION method calls the dump variants
        method of the relevant class to return the final VarTrack data."""
        if self._variant_tracker is None:
            return
        return self._variant_tracker.dump_variants(selected_genomes)