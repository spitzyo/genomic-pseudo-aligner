####################         IMPORTS          ####################
from aligner import Read, Aligner
from kmers import Reference, KmerCollection

####################          TESTS           ####################

def test_unmapped_read():
    genome1 = Reference("genome1", "ACGT")
    genome2 = Reference("genome2", "TGCA")
    kmercollection = KmerCollection()
    kmercollection.add_kmers(["ACG"], genome1, [0])
    kmercollection.add_kmers(["TGC"], genome2, [0])

    aligner = Aligner(kmercollection)
    read = Read("r1", "TAAT", "FFFF")
    aligner.align_read(read, k=3)
    assert read.status == "unmapped"
    # it's clear that neither 'ACG' nor 'TGA' would match with the kmers
    # in the read that are ['TAA', 'AAT'] -> should be unmapped.

def test_unique_read():
    genome1 = Reference("genome1", "ACGT")
    genome2 = Reference("genome2", "TGCA")
    kmercollection = KmerCollection()
    kmercollection.add_kmers(["ACG"], genome1, [0])
    kmercollection.add_kmers(["TGC"], genome2, [0])

    aligner = Aligner(kmercollection)
    read = Read("r2", "ACGT", "FFFF")
    aligner.align_read(read, k=3)
    assert read.status == "unique"
    assert read.mapped_genomes == ["genome1"]
    # same ref genomes like the prev test, but now the read has the following
    # kmers: ['ACG', 'GCT'] so it should match with genome1 uniquely.


def test_ambiguous_read():
    genome1 = Reference("genome1", "ACGT")
    genome2 = Reference("genome2", "CGTA")
    kmercollection = KmerCollection()
    kmercollection.add_kmers(["ACG"], genome1, [0])
    kmercollection.add_kmers(["CGT"], genome2, [0])

    aligner = Aligner(kmercollection)
    read = Read("r3", "ACGTA", "FFFFF")
    # this Read contains kmers from both genomes

    # make it ambiguous with m=2 parameter passed:
    aligner.align_read(read, k=3, m=2)
    assert read.status == "ambiguous"
    assert set(read.mapped_genomes) == {"genome1", "genome2"}


def test_read_quality():
    read = Read("test_id", "ACGTACGT", [20, 25, 30, 35, 20, 25, 30, 35])
    assert read.get_mean_quality() > 0
    assert read.get_kmer_quality(0, 4) > 0