####################         IMPORTS          ####################
from collections import defaultdict
from typing import List, Set
import sys


####################         CLASSES          ####################

class Reference:
    """
    This class holds all data related to a reference genome,
    usually imported from a FASTA file.
    The data includes: identifier, the sequence itself, dict of kmers
    and positions and the genome length.
    """

    def __init__(self, identifier: str, sequence: str):
        """
        Initialize a Reference genome object.
        :param identifier: a unique identifier for the reference genome.
        :param sequence: the DNA sequence of the reference genome (string).
        """
        self.identifier = identifier
        self.seq = sequence
        self.__ref_kmers = {} # dict to hold kmers and their positions
        self.total_bases = len(sequence) # num of bases == seq length

    def add_ref_kmers(self, k, kmer_collection) -> None:
        """
        This function adds k-mers of this certain reference genome
        to a given storage.
        :param k: length of every k-mer.
        :param kmer_collection: storage of k-mers (an KmerStorage instance).
        """
        if not isinstance(k, int) or k <= 0: # k validity
            print(f"Error: k size has to be a positive integer, got {k}.")
            sys.exit(1) # exit if invalid
        if k > self.total_bases: # handling invalid k values
            return None # do not process kmers if size is longer than seq

        kmer = self.seq[0:k] # inits a first kmer
        if "N" not in kmer: # N exclusion (N is a wildcard)
            kmer_collection.add_kmers([kmer], self, [0])
        for i in range(1, self.total_bases - k + 1):
            kmer = kmer[1:] + self.seq[i+k-1] # sliding by one base to rt
            if "N" not in kmer: kmer_collection.add_kmers([kmer], self, [i])


class Kmer:
    """
    This class represents a k-length sequence (k-mer), and has attributes
    of the sequence itself, the genomes in which it appears, and the starter
    positions within the genomes in which it appears.
    """

    def __init__(self, sequence: str):
        """
        Initialize a k-mer object.
        :param sequence: a string of the DNA sequence.
        """
        self.sequence = sequence
        self.__genome_pos = defaultdict(set)
        # a default dictionary of genomes and sets of positions
        self.__specific = None
        # defines if this k-mer is specific to one genome

    def add_position(self, genome: Reference, pos) -> None:
        """The method receives a ref genome and its position, and updates it
         into the kmer's genome pos attribute, and updates specificity."""
        self.__genome_pos[genome].add(pos) # puts it in the default dict
        self.update_specificity() # update it anytime we add a pos

    def update_specificity(self) -> None:
        """This method updates the specificity of a k-mer."""
        self.__specific = len(self.__genome_pos) == 1

    def is_specific(self) -> bool:
        if self.__specific is None:
            self.update_specificity()
        return self.__specific

    def get_genomes(self) -> List[Reference]:
        """This method returns the genomes in which the kmer appears."""
        return list(self.__genome_pos.keys())

    def get_positions(self, genome) -> list:
        """This method returns the positions of this kmer in a genome.
        It returns an empty list if nothing is found."""
        return self.__genome_pos.get(genome, [])


class KmerCollection:
    """
    This class stores and manages a collection of k-mers from diff genomes.
    """

    def __init__(self):
        """
        Initialize a collection of k-mers.
        Keys are the sequences (strings) and each value is
        a relevant Kmer instance.
        """
        self.__kmers = {} # initializes a dictionary
        self.__genome_index = defaultdict(set) # index for genome lookups
        self.__genome_order = [] # to retain order in the FASTA file

    def add_kmers(self, kmer_seqs, genome, positions) -> None:
        """
        This method adds or updates multiple k-mers to the storage.
        :param kmer_seqs: a list of k-mer sequences.
        :param genome: a reference genome (of the k-mer).
        :param positions: the positions within the genome of k-mer.
        """
        if (genome not in self.__genome_index and genome
                not in self.__genome_order):
            self.__genome_order.append(genome) # to track the genomes order
        for kmer_seq, position in zip(kmer_seqs, positions):
            if kmer_seq not in self.__kmers:
                self.__kmers[kmer_seq] = Kmer(kmer_seq)
            kmer = self.__kmers[kmer_seq]
            kmer.add_position(genome, position) # adding the new pos to dict
            self.__genome_index[genome].add(kmer_seq) # adding the index

    def get_kmer(self, kmer_seq) -> Kmer:
        """This method returns a Kmer instance of a given sequence (string)"""
        return self.__kmers.get(kmer_seq)

    def get_all_genomes(self) -> Set[Reference]:
        """This method returns the genomes in which at least 1 kmer appears."""
        genomes = set() # to avoid repetitiveness
        for kmer in self.__kmers.values():
            genomes.update(kmer.get_genomes())
        return genomes

    def get_all_kmers(self):
        """This method returns all the kmer instances."""
        return self.__kmers.values()

    def get_kmers_for_genome(self, genome: Reference) -> Set[str]:
        """This getter method returns all kmers of a genome."""
        return self.__genome_index[genome]

    def get_ordered_genomes(self) -> List[Reference]:
        return self.__genome_order