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