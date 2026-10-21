####################         IMPORTS          ####################
from genomic_pseudo_aligner.aligner import Read, Aligner
from genomic_pseudo_aligner.kmers import Reference, KmerCollection

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
    # Neither 'ACG' nor 'TGC' would match with the kmers
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
    aligner.align_read(read, k=3, unique_diff=2)
    assert read.status == "ambiguous"
    assert set(read.mapped_genomes) == {"genome1", "genome2"}


def test_read_quality():
    read = Read("test_id", "ACGTACGT", [20, 25, 30, 35, 20, 25, 30, 35])
    assert read.get_mean_quality() > 0
    assert read.get_kmer_quality(0, 4) > 0

####################   OFFSETS AND OVERHANGS   ####################
import random


def _random_genome(length, k, seed=0):
    """Random sequence whose k-mers all occur once (keeps tests exact)."""
    rng = random.Random(seed)
    while True:
        seq = "".join(rng.choice("ACGT") for _ in range(length))
        kmers = [seq[i:i + k] for i in range(length - k + 1)]
        if len(set(kmers)) == len(kmers):
            return seq


def _aligner_for(genome_seq, k, identifier="g1"):
    genome = Reference(identifier, genome_seq)
    collection = KmerCollection()
    genome.add_ref_kmers(k, collection)
    return Aligner(collection)


def test_offset_with_most_votes_wins():
    k = 8
    base = _random_genome(60, k, seed=1)
    read_seq = base[30:50]
    # Copy the read's first k-mer to position 5 as well: that k-mer now votes
    # for start positions 5 and 30, every other k-mer votes only for 30.
    genome_seq = base[:5] + read_seq[:k] + base[5 + k:]
    assert genome_seq[30:50] == read_seq

    aligner = _aligner_for(genome_seq, k)
    read = Read("r", read_seq, [30] * len(read_seq))
    aligner.align_read(read, k)
    assert read.status == "unique"
    assert read.positions == {"g1": 30}


def test_offset_ignores_skipped_kmers():
    """k-mers skipped because of an N must not shift the later offsets."""
    k = 8
    genome_seq = _random_genome(60, k, seed=2)
    read_seq = genome_seq[30:50]
    read_seq = read_seq[:3] + "N" + read_seq[4:]
    aligner = _aligner_for(genome_seq, k)
    read = Read("r", read_seq, [30] * len(read_seq))
    aligner.align_read(read, k)
    assert read.status == "unique"
    assert read.positions == {"g1": 30}


def test_read_overhanging_genome_start():
    k = 8
    genome_seq = _random_genome(60, k, seed=3)
    read_seq = "TTTTT" + genome_seq[:15]  # 5 bases before the genome start
    aligner = _aligner_for(genome_seq, k)
    aligner.enable_coverage({"g1": len(genome_seq)})
    aligner.enable_vartracking(min_quality=0, min_coverage=1)
    read = Read("r", read_seq, [30] * len(read_seq))
    aligner.align_read(read, k)

    assert read.status == "unique"
    assert read.positions == {"g1": -5}
    unique_cov = aligner._coverage.get_coverage_stats()["Details"]["g1"][
        "unique_cov"]
    assert unique_cov == [1] * 15 + [0] * 45
    assert aligner.dump_variants()["Variants"] == {}


def test_read_overhanging_genome_end():
    k = 8
    genome_seq = _random_genome(60, k, seed=4)
    read_seq = genome_seq[45:] + "TTTTT"  # 5 bases past the genome end
    aligner = _aligner_for(genome_seq, k)
    aligner.enable_coverage({"g1": len(genome_seq)})
    aligner.enable_vartracking(min_quality=0, min_coverage=1)
    read = Read("r", read_seq, [30] * len(read_seq))
    aligner.align_read(read, k)

    assert read.status == "unique"
    assert read.positions == {"g1": 45}
    unique_cov = aligner._coverage.get_coverage_stats()["Details"]["g1"][
        "unique_cov"]
    assert unique_cov == [0] * 45 + [1] * 15
    assert aligner.dump_variants()["Variants"] == {}


####################          STRANDS          ####################
from genomic_pseudo_aligner.aligner import reverse_complement


def test_reverse_complement():
    assert reverse_complement("AACGTN") == "NACGTT"


def test_forward_and_reverse_complement_reads_both_map_unique():
    k = 8
    genome_seq = _random_genome(80, k, seed=5)
    aligner = _aligner_for(genome_seq, k)
    forward = Read("fwd", genome_seq[20:50], [30] * 30)
    reverse = Read("rev", reverse_complement(genome_seq[20:50]), [30] * 30)
    aligner.align_read(forward, k)
    aligner.align_read(reverse, k)

    assert (forward.status, forward.strand) == ("unique", "+")
    assert (reverse.status, reverse.strand) == ("unique", "-")
    assert forward.positions == reverse.positions == {"g1": 20}


def test_unmapped_read_has_no_strand():
    k = 8
    aligner = _aligner_for(_random_genome(80, k, seed=5), k)
    read = Read("r", "A" * 30, [30] * 30)
    aligner.align_read(read, k)
    assert read.status == "unmapped"
    assert read.strand is None


def test_minus_strand_variant_reported_in_reference_orientation():
    k = 8
    genome_seq = _random_genome(80, k, seed=6)
    variant_pos = 35  # read covers 20..49, so this is read index 15
    ref_base = genome_seq[variant_pos]
    alt_base = "A" if ref_base != "A" else "C"
    fragment = (genome_seq[20:variant_pos] + alt_base
                + genome_seq[variant_pos + 1:50])
    # Forward-orientation qualities: only the variant base is high quality.
    fwd_quality = [5] * 30
    fwd_quality[variant_pos - 20] = 35
    # The sequencer reported the minus strand: reverse complement the bases
    # and reverse the qualities.
    read = Read("rev", reverse_complement(fragment), fwd_quality[::-1])

    aligner = _aligner_for(genome_seq, k)
    aligner.enable_vartracking(min_quality=20, min_coverage=1)
    aligner.align_read(read, k)

    assert (read.status, read.strand) == ("unique", "-")
    variants = aligner._variant_tracker.get_variants("g1")
    assert list(variants) == [variant_pos]
    assert variants[variant_pos].get_bases() == (ref_base, alt_base)
    assert variants[variant_pos].quality_score == 35
