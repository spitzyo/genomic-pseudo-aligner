####################         IMPORTS          ####################
from typing import Dict, List, Set
from collections import defaultdict
import json


####################         CLASSES          ####################

class Coverage:
    """
    Tracks how many reads (both unique and ambiguous) "cover" each position of
    the reference genomes (or some).
    """
    def __init__(self, genome_lens: Dict[str, int], min_coverage: int = 1):
        """
        :param genome_lens: a dictionary between genomes and their lengths.
        :param min_coverage: a given param to filter position below it.
        """
        self._unique_counts = defaultdict(lambda: defaultdict(int))
        self._ambig_counts = defaultdict(lambda: defaultdict(int))
        self._genome_lens = genome_lens
        self._min_coverage = min_coverage

    def add_read_coverage(self, read, genome_id: str,
                          positions: Set[int]) -> None:
        """Records coverage for a specific read."""
        if not positions:
            return
        if genome_id not in self._genome_lens:
            return

        if read.status == 'unique':
            coverage_dict = self._unique_counts
        else: # ambiguous
            coverage_dict = self._ambig_counts
        for pos in positions:
            coverage_dict[genome_id][pos] += 1

    def _count_covered_bases(self, coverage_list) -> int:
        """Counts positions meeting the minimum coverage threshold."""
        counter = 0
        for coverage in coverage_list:
            if coverage >= self._min_coverage:
                counter += 1
        return counter

    @staticmethod # this helper method doesn't use the Coverage inst by itself
    def _calculate_mean_coverage(coverage_list, genome_len) -> float:
        """Calculates the mean coverage, rounded to one decimal place."""
        total_coverage = 0
        for coverage in coverage_list:
            total_coverage += coverage
        mean_coverage = total_coverage / genome_len
        return round(mean_coverage, 1)

    def get_coverage_stats(self, selected_genomes: List[str]=None) -> Dict:
        """Generates stats summary and coverage of genome positions."""
        if selected_genomes is None: # for full coverage -> get all genomes
            selected_genomes = list(self._genome_lens.keys())
        stats_output = {"Coverage": {}, "Details": {}}

        for genome_id in selected_genomes:
            genome_len = self._genome_lens[genome_id]

            unique_coverage = [0]*genome_len
            ambiguous_coverage = [0]*genome_len
            for pos, count in self._unique_counts[genome_id].items():
                unique_coverage[pos] = count
            for pos, count in self._ambig_counts[genome_id].items():
                ambiguous_coverage[pos] = count

            covered_unique = self._count_covered_bases(unique_coverage)
            covered_ambiguous = self._count_covered_bases(ambiguous_coverage)
            mean_unique = self._calculate_mean_coverage(unique_coverage,
                                                         genome_len)
            mean_ambiguous = self._calculate_mean_coverage(ambiguous_coverage,
                                                            genome_len)

            stats_output["Coverage"][genome_id] = {
                "covered_bases_unique": covered_unique,
                "covered_bases_ambiguous": covered_ambiguous,
                "mean_coverage_unique": mean_unique,
                "mean_coverage_ambiguous": mean_ambiguous
            }

            stats_output["Details"][genome_id] = {
                "unique_cov": unique_coverage,
                "ambiguous_cov": ambiguous_coverage
            }

        return stats_output

    def dump_coverage(self, selected_genomes, full_coverage=True,
                      print_output=True):
        """Outputs coverage statistics to the console in JSON format."""
        stats = self.get_coverage_stats(selected_genomes)
        if full_coverage is False:
            if "Details" in stats:
                del stats["Details"]
        if print_output:
            print(json.dumps(stats, indent=4))
        return stats