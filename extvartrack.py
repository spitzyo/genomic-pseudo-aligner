####################         IMPORTS          ####################
from bisect import bisect_right
from collections import defaultdict
from typing import Dict, List, Tuple

####################         CLASSES          ####################

class Variant:
    """
    This class represents one alternate allele at one reference position,
    detected during the pseudo-alignment algorithm processing. Different
    alternate bases at the same position are separate Variant instances.
    """

    def __init__(self, pos: int, ref_base: str, alt_base: str, quality: float):
        """
        :param pos: position in the reference genome's sequence.
        :param ref_base: a single Reference nucleotide/base.
        :param alt_base: an alternative base seen in the Read
        :param quality: the quality score of the first base seen with this
        variant (used as quality_score until supporting reads are added).
        """
        self._position = pos
        self._ref_base = ref_base
        self._alt_base = alt_base
        self._first_quality = quality
        self._quality_sum = 0.0
        self.alt_count = 0 # num of reads showing this alternate base
        self.ref_count = 0 # num of reads showing the reference base
        self.depth = 0 # num of uniquely mapped reads covering this position

    @property
    def position(self) -> int:
        return self._position

    @property
    def coverage(self) -> int:
        """Older name for depth."""
        return self.depth

    @coverage.setter
    def coverage(self, value: int) -> None:
        self.depth = value

    @property
    def quality_score(self) -> float:
        """Mean Phred quality of the bases supporting this variant."""
        if self.alt_count == 0:
            return self._first_quality
        return self._quality_sum / self.alt_count

    @property
    def allele_frequency(self) -> float:
        """Fraction of the reads covering the position that show this alt."""
        return self.alt_count / self.depth if self.depth else 0.0

    def add_supporting_read(self, quality: float) -> None:
        """Records one more read showing this alternate base."""
        self.alt_count += 1
        self._quality_sum += quality

    def update_variant_counts(self, is_variant: bool=True):
        """Counts one more read covering this position; it supports the
        variant if is_variant, the reference otherwise."""
        self.depth += 1
        if is_variant:
            self.alt_count += 1
            self._quality_sum += self._first_quality
        else:
            self.ref_count += 1

    def get_bases(self) -> Tuple[str, str]:
        """Returns a tuple of (reference_base, alternate_base)."""
        return self._ref_base, self._alt_base


class VariantTracker:
    """
    This class actually performs the variant tracking process.
    An instance of this class can be stored as an Aligner instance's attribute.
    It tracks, per genome, the alternate bases seen at each position and the
    span of every uniquely mapped read, so that the depth at a position counts
    all reads covering it, not only the ones carrying a variant.
    """

    def __init__(self, min_quality=None, min_cov=None, min_af=None):
        """
        :param min_quality: the minimal quality score in order to be variant.
        Default determined by the Q20 quality score (1% error probability)
        :param min_cov: minimal depth (reads covering the position) to report
        a variant. Default determined to be the standard min coverage.
        :param min_af: minimal allele frequency (alt_count / depth) to report
        a variant. Default 0 (no filter).
        """
        self._min_quality = min_quality if min_quality is not None else 20.0
        self._min_coverage = min_cov if min_cov is not None else 10
        self._min_af = min_af if min_af is not None else 0.0
        # genome_id -> {position -> {alt_base -> Variant}}
        self._variants = {}
        # genome_id -> {position -> reads not showing the reference base},
        # counted at any quality, so ref_count = depth - non_ref
        self._non_ref = {}
        # genome_id -> start / end (exclusive) positions of the reads
        self._read_starts = defaultdict(list)
        self._read_ends = defaultdict(list)
        self._unsorted = set() # genomes whose start/end lists need sorting

    ####################      RECORDING READS     ####################

    def _add_read_span(self, genome_id, start, end) -> None:
        self._read_starts[genome_id].append(start)
        self._read_ends[genome_id].append(end)
        self._unsorted.add(genome_id)

    def _add_non_ref_base(self, genome_id, position, ref_base, alt_base,
                          quality) -> None:
        """Counts a read base that differs from the reference, and records
        it as a candidate variant if it passes the quality threshold."""
        counts = self._non_ref.setdefault(genome_id, {})
        counts[position] = counts.get(position, 0) + 1
        if quality < self._min_quality:
            return
        if alt_base == 'N' or ref_base == 'N':
            return # an unknown base is not a variant
        alleles = self._variants.setdefault(genome_id, {}).setdefault(
            position, {})
        if alt_base not in alleles:
            alleles[alt_base] = Variant(position, ref_base, alt_base, quality)
        alleles[alt_base].add_supporting_read(quality)

    def add_aligned_read(self, genome_id, ref_start, read_bases,
                         read_qualities, ref_bases) -> None:
        """
        Records one uniquely mapped read. The three sequences must be
        aligned with each other and clipped to the genome, in reference
        orientation: read_bases[i] sits on reference position ref_start + i.
        """
        if not read_bases:
            return
        self._add_read_span(genome_id, ref_start, ref_start + len(read_bases))
        for i, (read_base, ref_base) in enumerate(zip(read_bases, ref_bases)):
            if read_base != ref_base:
                self._add_non_ref_base(genome_id, ref_start + i, ref_base,
                                       read_base, read_qualities[i])

    def add_direct_variant(self, genome_id, position, ref_base, alt_base,
                           quality) -> None:
        """Records a single read that covers position and shows alt_base
        there. Counts toward depth whatever its quality."""
        self._add_read_span(genome_id, position, position + 1)
        self._add_non_ref_base(genome_id, position, ref_base, alt_base,
                               quality)

    ####################         QUERYING         ####################

    def get_depth(self, genome_id, position) -> int:
        """Number of recorded reads that cover position."""
        starts = self._read_starts.get(genome_id)
        if not starts:
            return 0
        ends = self._read_ends[genome_id]
        if genome_id in self._unsorted:
            starts.sort()
            ends.sort()
            self._unsorted.discard(genome_id)
        # reads starting at or before position, minus those already ended
        return bisect_right(starts, position) - bisect_right(ends, position)

    def _candidates(self, genome_id) -> List[Variant]:
        """All candidate variants of a genome, with depth and ref_count
        filled in, sorted by (position, alternate base)."""
        candidates = []
        for pos, alleles in sorted(self._variants.get(genome_id, {}).items()):
            depth = self.get_depth(genome_id, pos)
            ref_count = depth - self._non_ref[genome_id].get(pos, 0)
            for alt_base in sorted(alleles):
                variant = alleles[alt_base]
                variant.depth = depth
                variant.ref_count = ref_count
                candidates.append(variant)
        return candidates

    def _passes(self, variant) -> bool:
        return (variant.depth >= self._min_coverage and
                variant.allele_frequency >= self._min_af)

    def get_variant_calls(self, genome_id) -> List[Variant]:
        """Returns every variant (one per position and alternate base) that
        meets the depth and allele frequency thresholds."""
        return [v for v in self._candidates(genome_id) if self._passes(v)]

    def get_variants(self, genome_id) -> Dict[int, Variant]:
        """Returns {position: Variant} for variants meeting the thresholds.
        Where several alternate bases pass at one position, this holds the
        one with the most reads; get_variant_calls() returns all of them.
        Read-only: calling it any number of times gives the same result."""
        variants = {}
        for variant in self.get_variant_calls(genome_id):
            current = variants.get(variant.position)
            if current is None or variant.alt_count > current.alt_count:
                variants[variant.position] = variant
        return variants

    def get_stats(self, genome_id) -> Dict[str, int]:
        """Counts this genome's candidate variants: total_variants, how many
        pass the filters (reported_variants) and how many fail them
        (filtered_variants). Computed from scratch on every call."""
        candidates = self._candidates(genome_id)
        reported = sum(1 for v in candidates if self._passes(v))
        return {'total_variants': len(candidates),
                'reported_variants': reported,
                'filtered_variants': len(candidates) - reported}

    def dump_variants(self, selected_genomes=None) -> Dict:
        """Formats the variant data and statistics for JSON output.
        Variants[genome][position] is a list with one entry per alternate
        base at that position."""

        output = {"Variants": {}, "Statistics": {}}
        if selected_genomes is None:
            genomes = list(self._variants.keys())
        else:
            genomes = selected_genomes

        for genome_id in genomes:
            if genome_id not in self._variants:
                continue # skip if this genome has no variants
            output['Statistics'][genome_id] = self.get_stats(genome_id)
            variants = self.get_variant_calls(genome_id)
            if not variants:
                continue
            genome_variants = {}
            for variant in variants:
                ref_base, alt_base = variant.get_bases()
                genome_variants.setdefault(str(variant.position), []).append(
                    {'reference': ref_base,
                     'alternate': alt_base,
                     'depth': variant.depth,
                     'ref_count': variant.ref_count,
                     'alt_count': variant.alt_count,
                     'allele_frequency': round(variant.allele_frequency, 4),
                     'quality_score': round(variant.quality_score, 2)})
            output['Variants'][genome_id] = genome_variants
        return output
