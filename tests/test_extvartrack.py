####################         IMPORTS          ####################
from extvartrack import Variant, VariantTracker

####################          TESTS           ####################

def test_variant_object():
    my_variant = Variant(10, "A", "C", 25.0)
    assert my_variant.get_bases() == ("A", "C")
    assert my_variant.coverage == 0

    # Add some reads: one that supports the variant and other that doesn't
    my_variant.update_variant_counts() # default is True (supports variant)
    my_variant.update_variant_counts(False) # if bool is false = supports ref.
    assert my_variant.coverage == 2 # check how many reads cover position


def test_varianttracker():
    # Set min_cov=1 so a single variant will pass the filter
    vartrack = VariantTracker(min_cov=1)

    # Add a SNP: (T) -> (A) at position 14
    vartrack.add_direct_variant("bacteria", 14, "T", "A", 30)
    variants = vartrack.get_variants("bacteria")
    assert len(variants) == 1 # check if the variant we added exists

    pos = 14
    assert pos in variants
    snp = variants[pos]
    assert snp.get_bases() == ("T", "A") # check if we can use the variant
    # getter method to get a tuple with the original base and its alternative.


def test_coverage_filter():
    variantracker = VariantTracker(min_cov=3)

    # Add variant once - shouldn't appear in results
    variantracker.add_direct_variant("sample1", 7, "T", "C", 40)
    assert variantracker.get_variants("sample1") == {}

    # Add twice more to pass filter (total coverage minimum is set to 3)
    variantracker.add_direct_variant("sample1", 7, "T", "C", 40)
    variantracker.add_direct_variant("sample1", 7, "T", "C", 40)
    assert len(variantracker.get_variants("sample1")) == 1 # only one should be

def test_variant_stats_count_filtered_and_are_idempotent():
    tracker = VariantTracker(min_cov=2)
    tracker.add_direct_variant("g1", 5, "A", "C", 30)  # coverage 1: filtered
    tracker.add_direct_variant("g1", 9, "G", "T", 30)
    tracker.add_direct_variant("g1", 9, "G", "T", 30)  # coverage 2: reported

    expected = {"g1": {"total_variants": 2, "reported_variants": 1,
                       "filtered_variants": 1}}
    tracker.get_variants("g1")  # must not change the stats
    first = tracker.dump_variants()
    second = tracker.dump_variants()
    assert first["Statistics"] == expected
    assert second["Statistics"] == expected
    assert first == second


####################   DEPTH AND ALLELE COUNTS   ####################
import random
import pytest
from aligner import Aligner, Read, reverse_complement
from kmers import Reference, KmerCollection
from main import readargs

K = 8
VARIANT_POS = 40


def _random_genome(length, k, seed=0):
    """Random sequence whose k-mers all occur once (keeps tests exact)."""
    rng = random.Random(seed)
    while True:
        seq = "".join(rng.choice("ACGT") for _ in range(length))
        kmers = [seq[i:i + k] for i in range(length - k + 1)]
        if len(set(kmers)) == len(kmers):
            return seq


def _align(reads_alt_bases, min_cov=1, min_af=None, seed=7):
    """Aligns one 30 bp read per entry of reads_alt_bases, all covering
    VARIANT_POS; None means the read shows the reference base there.
    Reads alternate between the + and - strand and start at different
    positions. Returns (reference base, tracker)."""
    genome_seq = _random_genome(100, K, seed)
    collection = KmerCollection()
    Reference("g1", genome_seq).add_ref_kmers(K, collection)
    aligner = Aligner(collection)
    aligner.enable_vartracking(min_quality=20, min_coverage=min_cov,
                               min_af=min_af)
    for n, alt in enumerate(reads_alt_bases):
        start = VARIANT_POS - 10 - n  # read covers VARIANT_POS at index 10+n
        bases = list(genome_seq[start:start + 30])
        if alt is not None:
            bases[VARIANT_POS - start] = alt
        seq = "".join(bases)
        if n % 2:
            seq = reverse_complement(seq)
        read = Read(f"r{n}", seq, [30] * 30)
        aligner.align_read(read, K)
        assert read.status == "unique"
    return genome_seq[VARIANT_POS], aligner._variant_tracker


def _other_bases(ref_base):
    return [b for b in "ACGT" if b != ref_base]


def test_depth_counts_reference_reads_too():
    """4 reference reads + 1 alt read: depth 5, alt_count 1, AF 0.2."""
    ref_base = _align([])[0]
    alt = _other_bases(ref_base)[0]
    _, tracker = _align([None, None, alt, None, None])

    records = tracker.dump_variants()["Variants"]["g1"][str(VARIANT_POS)]
    assert records == [{"reference": ref_base, "alternate": alt,
                        "depth": 5, "ref_count": 4, "alt_count": 1,
                        "allele_frequency": 0.2, "quality_score": 30.0}]
    assert tracker.get_depth("g1", VARIANT_POS) == 5


def test_different_alts_at_one_position_are_kept_separate():
    ref_base = _align([])[0]
    alt1, alt2 = _other_bases(ref_base)[:2]
    _, tracker = _align([alt1, None, alt2, alt1])

    records = tracker.dump_variants()["Variants"]["g1"][str(VARIANT_POS)]
    by_alt = {r["alternate"]: r for r in records}
    assert set(by_alt) == {alt1, alt2}
    assert by_alt[alt1]["alt_count"] == 2
    assert by_alt[alt2]["alt_count"] == 1
    for record in records:
        assert record["reference"] == ref_base
        assert record["depth"] == 4
        assert record["ref_count"] == 1
    assert by_alt[alt1]["allele_frequency"] == 0.5
    assert by_alt[alt2]["allele_frequency"] == 0.25
    assert tracker.get_stats("g1") == {"total_variants": 2,
                                       "reported_variants": 2,
                                       "filtered_variants": 0}


def test_min_coverage_filters_on_depth():
    ref_base = _align([])[0]
    alt = _other_bases(ref_base)[0]
    reads = [None, None, alt, None, None]
    assert _align(reads, min_cov=5)[1].get_variants("g1")  # depth 5 passes
    assert not _align(reads, min_cov=6)[1].get_variants("g1")


def test_min_allele_frequency_filter():
    ref_base = _align([])[0]
    alt = _other_bases(ref_base)[0]
    reads = [None, None, alt, None, None]  # AF 0.2
    assert _align(reads, min_af=0.2)[1].get_variants("g1")
    _, tracker = _align(reads, min_af=0.25)
    assert not tracker.get_variants("g1")
    assert tracker.get_stats("g1")["filtered_variants"] == 1


def test_low_quality_alt_counts_toward_depth_not_alt_count():
    tracker = VariantTracker(min_quality=20, min_cov=1)
    tracker.add_aligned_read("g1", 10, "AAC", [30, 30, 30], "AAA")
    tracker.add_aligned_read("g1", 10, "AAC", [30, 30, 5], "AAA")
    tracker.add_aligned_read("g1", 11, "AA", [30, 30], "AA")
    [variant] = tracker.get_variant_calls("g1")
    assert variant.position == 12
    assert (variant.depth, variant.alt_count, variant.ref_count) == (3, 1, 1)


def test_min_variant_af_argument():
    args = readargs(["-t", "align", "--min-variant-af", "0.3"])
    assert args.min_variant_af == 0.3
    with pytest.raises(SystemExit):
        readargs(["-t", "align", "--min-variant-af", "30"])
