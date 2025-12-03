####################         IMPORTS          ####################
from typing import Dict, List, Set
from collections import defaultdict
import json


####################         CLASSES          ####################

class Coverage:
    """
    This new class helps with the new EXTCOVERAGE extension of the program,
    in which we aim to track how many reads (both unique and ambiguous)
    "cover" each position of the reference genomes (or some of them).
    """
    def __init__(self, genome_lens: Dict[str, int], min_coverage: int = 1):
        """
        :param genome_lens: a dictionary between genomes and their lengths.
        :param min_coverage: a given param to filter position below it.
        """
        # initiating two distinct counters of genome->position->count:
        self.__unique_counts = defaultdict(lambda: defaultdict(int))
        self.__ambig_counts = defaultdict(lambda: defaultdict(int))
        self.__genome_lens = genome_lens
        self.__min_coverage = min_coverage

    def add_read_coverage(self, read, genome_id: str,
                          positions: Set[int]) -> None:
        """This method adds coverage for a given read.
        It avoids double counting the positions."""
        if not positions:
            return # can happen with empty positions (in low quality reads)
        if genome_id not in self.__genome_lens:
            return # another validity check

        if read.get_status() == 'unique':
            coverage_dict = self.__unique_counts
        else: # hence read is ambiguous
            coverage_dict = self.__ambig_counts
        for pos in positions:
            coverage_dict[genome_id][pos] += 1 # in either of the attrs

    def __count_covered_bases(self, coverage_list) -> int:
        """A helper method for get_coverage_stats that counts all the positions
        that have coverage greater or equal to the min_coverage threshold."""
        counter = 0
        for coverage in coverage_list:
            if coverage >= self.__min_coverage:
                counter += 1
        return counter

    @staticmethod # this helper method doesn't use the Coverage inst by itself
    def __calculate_mean_coverage(coverage_list, genome_len) -> float:
        """A helper method for get_coverage_stats that calculates the mean
         coverage and rounds it to one decimal place after the dot."""
        total_coverage = 0
        for coverage in coverage_list:
            total_coverage += coverage
        mean_coverage = total_coverage / genome_len
        return round(mean_coverage, 1) # round to one place after the dot

    def get_coverage_stats(self, selected_genomes: List[str]=None) -> Dict:
        """Method generates stats summary and coverage of genome positions."""
        if selected_genomes is None: # for full coverage -> get all genomes
            selected_genomes = list(self.__genome_lens.keys())
        stats_output = {"Coverage": {}, "Details": {}} # inits output structure

        for genome_id in selected_genomes:
            genome_len = self.__genome_lens[genome_id]

            unique_coverage = [0]*genome_len
            ambiguous_coverage = [0]*genome_len # inits lists for coverage
            for pos, count in self.__unique_counts[genome_id].items():
                unique_coverage[pos] = count
            for pos, count in self.__ambig_counts[genome_id].items():
                ambiguous_coverage[pos] = count # filling the lists with counts

            covered_unique = self.__count_covered_bases(unique_coverage)
            covered_ambiguous = self.__count_covered_bases(ambiguous_coverage)
            mean_unique = self.__calculate_mean_coverage(unique_coverage,
                                                         genome_len)
            mean_ambiguous = self.__calculate_mean_coverage(ambiguous_coverage,
                                                            genome_len)

            # final output:
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
        """This method dumps the coverage stats to console in a JSON format.
        The method would be called by a method in Aligner class.
        The statistics can be either printed or returned."""
        stats = self.get_coverage_stats(selected_genomes)
        if full_coverage is False:
            if "Details" in stats:
                del stats["Details"] # removed details key if selected genomes
        if print_output:
            print(json.dumps(stats, indent=4)) # dump to console
        return stats