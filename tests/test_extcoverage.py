####################         IMPORTS          ####################
from extcoverage import Coverage
from aligner import Read

####################          TESTS           ####################

def test_coverage():
    genome_lens = {"genome1": 100, "genome2": 200}
    coverage = Coverage(genome_lens) # inits a new instance to test
    read = Read("read1", "ACGT", [30, 30, 30, 30])
    read.status = "unique"
    coverage.add_read_coverage(read, "genome1", {5, 6, 7})

    stats = coverage.get_coverage_stats() # get coverage statistics and check:
    assert stats["Coverage"]["genome1"]["covered_bases_unique"] == 3
    assert stats["Coverage"]["genome1"]["mean_coverage_unique"] == 0.0


def test_coverage_threshold():
    coverage = Coverage({"genome1": 20}, min_coverage=2)
    read = Read("read1", "ACGT", [30, 30, 30, 30])
    read.status = "unique"

    # Add coverage for same position twice
    coverage.add_read_coverage(read, "genome1", {5})
    coverage.add_read_coverage(read, "genome1", {5})
    coverage.add_read_coverage(read, "genome1", {10})
    stats = coverage.get_coverage_stats() # create statistics
    assert stats["Coverage"]["genome1"]["covered_bases_unique"] == 1


def test_coverage_details():
    coverage = Coverage({"genome1": 10})

    # Add unique and ambiguous read coverage
    unique_read = Read("read1", "ACGT", [30, 30, 30, 30])
    unique_read.status = "unique"
    coverage.add_read_coverage(unique_read, "genome1", {2, 3})

    ambig_read = Read("read2", "TGCA", [30, 30, 30, 30])
    ambig_read.status = "ambiguous"
    coverage.add_read_coverage(ambig_read, "genome1", {4, 5})

    stats = coverage.get_coverage_stats() # create stats & check final details
    assert stats["Details"]["genome1"]["unique_cov"][2] == 1
    assert stats["Details"]["genome1"]["ambiguous_cov"][4] == 1


def test_full_coverage():
    coverage = Coverage({"genome1": 10})
    unique_read = Read("read1", "ACGT", [30, 30, 30, 30])
    unique_read.status = "unique"
    coverage.add_read_coverage(unique_read, "genome1", {2, 3})
    result = coverage.dump_coverage(None, True, False)
    assert "Details" in result