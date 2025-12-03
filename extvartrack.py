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
        self.__quality_score = quality
        self.__supporting_reads = 0 # num of reads showing with spec variant
        self.__coverage = 0 # num of reads covering this base position

    def update_variant_counts(self, is_variant: bool=True):
        """This method updates counts for this specific variant position.
        The is_variant flag differentiates between whether this read supports
        the variant or whether it supports the Reference genome."""
        self.__coverage += 1
        if is_variant:
            self.__supporting_reads += 1

    def get_bases(self) -> Tuple[str, str]:
        """This method returns a tuple that consists of the base in the
        reference genome, and the variation / alternative base."""
        return self.__ref_base, self.__alt_base

    def get_quality_score(self) -> float:
        return self.__quality_score

    def get_coverage(self) -> float:
        return self.__coverage


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
        """This private method inits the stats for a certain genome,
        and stores it in the stats attribute. It is called by add_variant()."""

        if genome_id not in self.__stats:
            self.__stats[genome_id] = {'total_variants': 0,
                                       'filtered_variants': 0}
            # just inits the place in the dict with empty values (0 for each)

    def add_direct_variant(self, genome_id, position, ref_base, alt_base,
                           quality) -> None:
        """
        This method adds a potential variant to tracker, by looking at
        the scenario in which a kmer does not match the Reference genome
        perfectly, and then we check each pos to find the mismatch's pos.
        """
        if quality < self.__min_quality:
            return # do not process low quality reads

        if genome_id not in self.__variants:
            self.__variants[genome_id] = {}
            self.__init_stats(genome_id)

        if position not in self.__variants[genome_id]:
            self.__variants[genome_id][position] = Variant(position, ref_base,
                                                           alt_base, quality)
        self.__variants[genome_id][position].update_variant_counts() # count

    def get_variants(self, genome_id) -> Dict[int, Variant]:
        """
        This method returns variants after filtering (by quality thresholds).
        It uses our given min_quality and min_coverage.
        """
        if genome_id not in self.__variants:
            return {} # in case no such genome to check variants for

        filtered_variants = {}
        for pos, variant in self.__variants[genome_id].items():
            total_coverage = variant.get_coverage()
            if total_coverage >= self.__min_coverage:
                filtered_variants[pos] = variant
                self.__stats[genome_id]['total_variants'] += 1
            else:
                self.__stats[genome_id]['total_variants'] += 1 # count in stats
        return filtered_variants

    def dump_variants(self, selected_genomes=None) -> Dict:
        """
        This method dumps the final variant information in a JSON format.
        :param selected_genomes: list of genomes to dump (if needed)
        """

        output = {"Variants": {}, "Statistics": self.__stats}
        if selected_genomes is None:
            genomes = list(self.__variants.keys())
        else: # in case only some genomes need to be dumped
            genomes = selected_genomes

        for genome_id in genomes:
            if genome_id not in self.__variants:
                continue # skip if this genome has no variants
            variants = self.get_variants(genome_id)
            if not variants:
                continue # if there are no variants to process
            genome_variants = {}
            for pos, variant in variants.items(): # iterating through the vars
                ref_base, alt_base = variant.get_bases()
                genome_variants[str(pos)] = \
                    {'reference': ref_base,
                    'alternate': alt_base,
                    'quality_score': variant.get_quality_score(),
                    'coverage': variant.get_coverage()}
            output['Variants'][genome_id] = genome_variants # for e/ genome
        return output # to be dumped or save in a JSON format