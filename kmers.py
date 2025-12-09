####################         IMPORTS          ####################
from collections import defaultdict
from typing import List, Set


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
        self._ref_kmers = {} # dict to hold kmers and their positions
        self.total_bases = len(sequence) # num of bases == seq length

    def add_ref_kmers(self, k, kmer_collection) -> None:
        """
        This function adds k-mers of this certain reference genome
        to a given storage.
        :param k: length of every k-mer.
        :param kmer_collection: storage of k-mers (an KmerStorage instance).
        """
        if not isinstance(k, int) or k <= 0: # k validity
            raise ValueError(f"k size must be a positive integer, got {k}.")
        if k > self.total_bases: # handling invalid k values
            return None # do not process kmers if size is longer than seq

        kmer = self.seq[0:k]
        if "N" not in kmer:
            kmer_collection.add_kmers([kmer], self, [0])
        for i in range(1, self.total_bases - k + 1):
            kmer = kmer[1:] + self.seq[i+k-1]
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
        self._genome_pos = defaultdict(set)
        self._specific = None # if k-mer is specific to one genome

    def add_position(self, genome: Reference, pos) -> None:
        """The method receives a ref genome and its position, and updates it
         into the kmer's genome pos attribute, and updates specificity."""
        self._genome_pos[genome].add(pos)
        self.update_specificity()

    def update_specificity(self) -> None:
        """Updates the specificity status of kmer."""
        self._specific = len(self._genome_pos) == 1

    def is_specific(self) -> bool:
        if self._specific is None:
            self.update_specificity()
        return self._specific

    def get_genomes(self) -> List[Reference]:
        """This method returns the genomes in which the kmer appears."""
        return list(self._genome_pos.keys())

    def get_positions(self, genome) -> list:
        """Returns the positions of this kmer in a genome or an empty list."""
        return self._genome_pos.get(genome, [])


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
        self._kmers = {}
        self._genome_index = defaultdict(set)
        self._genome_order = [] # to retain order in the original FASTA file

    def add_kmers(self, kmer_seqs, genome, positions) -> None:
        """
        This method adds or updates multiple k-mers to the storage.
        :param kmer_seqs: a list of k-mer sequences.
        :param genome: a reference genome (of the k-mer).
        :param positions: the positions within the genome of k-mer.
        """
        if (genome not in self._genome_index and genome
                not in self._genome_order):
            self._genome_order.append(genome)
        for kmer_seq, position in zip(kmer_seqs, positions):
            if kmer_seq not in self._kmers:
                self._kmers[kmer_seq] = Kmer(kmer_seq)
            kmer = self._kmers[kmer_seq]
            kmer.add_position(genome, position)
            self._genome_index[genome].add(kmer_seq)

    def get_kmer(self, kmer_seq) -> Kmer:
        """This method returns a Kmer instance of a given sequence (string)"""
        return self._kmers.get(kmer_seq)

    def get_all_genomes(self) -> Set[Reference]:
        """This method returns the genomes in which at least 1 kmer appears."""
        genomes = set()
        for kmer in self._kmers.values():
            genomes.update(kmer.get_genomes())
        return genomes

    def get_all_kmers(self):
        """This method returns all the kmer instances."""
        return self._kmers.values()

    def get_kmers_for_genome(self, genome: Reference) -> Set[str]:
        """Returns all k-mer sequences associated with the given genome."""
        return self._genome_index[genome]

    def get_ordered_genomes(self) -> List[Reference]:
        return self._genome_order