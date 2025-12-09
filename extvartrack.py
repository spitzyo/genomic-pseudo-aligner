####################         IMPORTS          ####################
from typing import Dict, Tuple

####################         CLASSES          ####################

class Variant:
    """
    This class represents a single variant, that would be detected during
    the pseudo-alignment algorithm processing.
    """

    def __init__(self, pos: int, ref_base: str, alt_base: str, quality: float):
        """
        :param pos: position in the reference genome's sequence.
        :param ref_base: a single Reference nucleotide/base.
        :param alt_base: an alternative base seen in the Read
        :param quality: the quality score of this variant base
        """
        self.__position = pos
        self.__ref_base = ref_base
        self.__alt_base = alt_base
        self.quality_score = quality
        self.__supporting_reads = 0 # num of reads showing with spec variant
        self.coverage = 0 # num of reads covering this base position

    def update_variant_counts(self, is_variant: bool=True):
        """Updates coverage and support counts for this variant position."""
        self.coverage += 1
        if is_variant:
            self.__supporting_reads += 1

    def get_bases(self) -> Tuple[str, str]:
        """Returns a tuple of (reference_base, alternate_base)."""
        return self.__ref_base, self.__alt_base


class VariantTracker:
    """
    This class actually performs the variant tracking process.
    An instance of this class can be stored as an Aligner instance's attribute.
    It actually tracks variants for every genome, and calculates filtering,
    while basing on quality scores (as provided) and coverage thresholds.
    """

    def __init__(self, min_quality=None, min_cov=None):
        """
        :param min_quality: the minimal quality score in order to be variant.
        Default determined by the Q20 quality score (1% error probability)
        :param min_cov: minimal position coverage to call variants.
        We need coverage to hit the min so we would have a large enough sample.
        Default determined to be the standard min coverage for variants.
        """
        self.__min_quality = min_quality if min_quality is not None else 20.0
        self.__min_coverage = min_cov if min_cov is not None else 10
        self.__variants = {} # genome_id -> {position -> Variant instance}
        self.__stats = {}    # genome_id -> statistics_dict

    def __init_stats(self, genome_id: str):
        """Initializes the statistics dictionary for a specific genome."""
        if genome_id not in self.__stats:
            self.__stats[genome_id] = {'total_variants': 0,
                                       'filtered_variants': 0}

    def add_direct_variant(self, genome_id, position, ref_base, alt_base,
                           quality) -> None:
        """Records a potential variant if it meets quality thresholds."""
        if quality < self.__min_quality:
            return # do not process low quality reads

        if genome_id not in self.__variants:
            self.__variants[genome_id] = {}
            self.__init_stats(genome_id)

        if position not in self.__variants[genome_id]:
            self.__variants[genome_id][position] = Variant(position, ref_base,
                                                           alt_base, quality)
        self.__variants[genome_id][position].update_variant_counts()

    def get_variants(self, genome_id) -> Dict[int, Variant]:
        """Returns variants that meet the minimum coverage threshold."""
        if genome_id not in self.__variants:
            return {}

        filtered_variants = {}
        for pos, variant in self.__variants[genome_id].items():
            total_coverage = variant.coverage
            if total_coverage >= self.__min_coverage:
                filtered_variants[pos] = variant
                self.__stats[genome_id]['total_variants'] += 1
            else:
                self.__stats[genome_id]['total_variants'] += 1
        return filtered_variants

    def dump_variants(self, selected_genomes=None) -> Dict:
        """Formats the variant data and statistics for JSON output."""

        output = {"Variants": {}, "Statistics": self.__stats}
        if selected_genomes is None:
            genomes = list(self.__variants.keys())
        else:
            genomes = selected_genomes

        for genome_id in genomes:
            if genome_id not in self.__variants:
                continue # skip if this genome has no variants
            variants = self.get_variants(genome_id)
            if not variants:
                continue
            genome_variants = {}
            for pos, variant in variants.items():
                ref_base, alt_base = variant.get_bases()
                genome_variants[str(pos)] = \
                    {'reference': ref_base,
                    'alternate': alt_base,
                    'quality_score': variant.quality_score,
                    'coverage': variant.coverage}
            output['Variants'][genome_id] = genome_variants
        return output