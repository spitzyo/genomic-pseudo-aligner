####################         IMPORTS          ####################
from kmers import Reference, Kmer, KmerCollection

####################          TESTS           ####################

# Reference class:
def test_reference_getters():
    reference = Reference("genome1", "ACGTACGT")
    assert reference.identifier == "genome1"
    assert reference.seq == "ACGTACGT"
    assert reference.total_bases == 8

def test_reference_add_kmers():
    reference = Reference("genome1", "ACGTACGT")
    collection = KmerCollection()
    reference.add_ref_kmers(3, collection) # create kmers with size 3
    assert collection.get_kmer("ACG") is not None # assert ACG is there
    assert collection.get_kmer("CGT") is not None # same for CGT
    assert collection.get_kmer("AAA") is None  # shouldn't exist

def test_reference_with_n():
    reference = Reference("genome1", "ACNGT")
    collection = KmerCollection()
    reference.add_ref_kmers(2, collection)
    assert collection.get_kmer("AC") is not None # AC should be there
    assert collection.get_kmer("CN") is None  # N should be excluded
    assert collection.get_kmer("GT") is not None # GT should exist

# Kmer class:
def test_kmer_object():
    kmer = Kmer("ACG")
    assert kmer.sequence == "ACG"
    assert len(kmer.get_genomes()) == 0 # it shouldn't have any linked genomes

def test_kmer_add_position():
    kmer = Kmer("ACG")
    reference1 = Reference("genome1", "ACGTACGT")
    kmer.add_position(reference1, 0)
    kmer.add_position(reference1, 4)
    assert len(kmer.get_genomes()) == 1 # only one genome is attached to kmer
    assert reference1 in kmer.get_genomes()
    assert 0 in kmer.get_positions(reference1) # ACG is in position 0
    assert 4 in kmer.get_positions(reference1) # and also in position 4

def test_kmer_specificity():
    kmer = Kmer("ACG")
    reference1 = Reference("genome1", "ACGTACGT")
    reference2 = Reference("genome2", "ACGTGCA")

    kmer.add_position(reference1, 0)
    assert kmer.is_specific() is True # check if ACG is only in one genome
    kmer.add_position(reference2, 0)
    assert kmer.is_specific() is False # check if ACG is connected to two

# KmerCollection:
def test_kmer_collection():
    collection = KmerCollection() # this is an empty collection
    assert len(list(collection.get_all_kmers())) == 0
    assert len(collection.get_all_genomes()) == 0

def test_add_kmers_to_collection():
    collection = KmerCollection()
    reference1 = Reference("genome1", "ACGTACGT")
    collection.add_kmers(["ACG", "CGT"], reference1, [0, 1])
    assert collection.get_kmer("ACG") is not None
    assert collection.get_kmer("CGT") is not None
    assert collection.get_kmer("AAA") is None # this kmer shouldn't exist

def test_kmers_for_genome():
    collection = KmerCollection()
    reference1 = Reference("genome1", "ACGTACGT")
    reference2 = Reference("genome2", "TGCATGCA")
    collection.add_kmers(["ACG", "CGT"], reference1, [0, 1])
    collection.add_kmers(["TGC", "GCA"], reference2, [0, 1])
    genome1_kmers = collection.get_kmers_for_genome(reference1)
    assert "ACG" in genome1_kmers
    assert "CGT" in genome1_kmers
    assert "TGC" not in genome1_kmers # it should only be a part of genome2