####################         IMPORTS          ####################

import pytest
from genomic_pseudo_aligner.kmers import Reference, KmerCollection
import gzip
from genomic_pseudo_aligner.helper_functions import (
    import_fasta, import_fastq, load_kdb_file, save_kmer_collection,
    is_pos_int, open_gz_file, AppErr)

####################          TESTS           ####################

def test_open_gz_file_nogzip(tmp_path):
    txt_file = tmp_path / "test.txt"
    txt_file.write_text("test data")
    with open_gz_file(str(txt_file), "r") as f:
        content = f.read() # read a regular file using the gzip opener
    assert content == "test data"

def test_open_gz_file_gzip(tmp_path):
    gz_file = tmp_path / "test.gz"
    with gzip.open(str(gz_file), 'wt') as f:
        f.write("compressed data") # create a gzip
    with open_gz_file(str(gz_file), "r") as f:
        content = f.read() # read the gzip
    assert content == "compressed data"


def test_import_fasta(tmp_path):
    fasta_file = tmp_path / "test.fa"
    fasta_file.write_text(">Seq1\nACGT\n>Seq2\nTGCA\n")
    references = import_fasta(str(fasta_file))

    assert len(references) == 2 # make sure two genomes were imported
    assert references[0].identifier == "Seq1"
    assert references[0].seq == "ACGT"
    assert references[1].identifier == "Seq2"
    assert references[1].seq == "TGCA"

def test_import_fasta_invalid(tmp_path):
    """This tests checks if import_fasta correctly handles
    sequences with invalid nucleotides (called 'Invalid')."""
    fasta_file = tmp_path / "test.fa"
    fasta_file.write_text(">Valid\nACGT\n>Invalid\nILOVEINTRO\n")
    references = import_fasta(str(fasta_file))
    assert len(references) == 1 # assure only one valid genome was imported
    assert references[0].identifier == "Valid"

def test_import_fasta_soft_masked(tmp_path):
    """Soft-masked (lowercase) bases, as in UCSC/Ensembl references, are
    accepted and uppercased."""
    fasta_file = tmp_path / "test.fa"
    fasta_file.write_text(">chr1\nACGTacgtNNnn\nttagGCAT\n>chr2\nacgt\n")
    references = import_fasta(str(fasta_file))
    assert [r.identifier for r in references] == ["chr1", "chr2"]
    assert references[0].seq == "ACGTACGTNNNNTTAGGCAT"
    assert references[1].seq == "ACGT"

def test_import_fasta_lowercase_invalid_still_rejected(tmp_path, capsys):
    fasta_file = tmp_path / "test.fa"
    fasta_file.write_text(">Valid\nacgt\n>Invalid\nacgt\nacxt\n")
    references = import_fasta(str(fasta_file))
    assert [r.identifier for r in references] == ["Valid"]
    assert "skipped 'Invalid': line 5" in capsys.readouterr().err


def test_save_load_kdb(tmp_path):
    kmer_collection = KmerCollection()
    ref = Reference("test_ref", "ACGTACGT")
    kmer_collection.add_kmers(["ACG", "CGT"], ref, [0, 1])

    kdb_file = tmp_path / "test.kdb" # save those kmers to a temporary kdb file
    save_kmer_collection(kmer_collection, str(kdb_file))
    loaded_collection = load_kdb_file(str(kdb_file)) # try load the saved kdb
    loaded_kmer = loaded_collection.get_kmer("ACG")

    assert loaded_kmer is not None # if import is successful and data retained
    assert loaded_kmer.sequence == "ACG"
    assert len(loaded_kmer.get_genomes()) == 1
    assert loaded_kmer.get_genomes()[0].identifier == "test_ref"


def test_import_fastq(tmp_path):
    fastq_file = tmp_path / "test.fq"
    fastq_file.write_text("@Read1\nACGT\n+\nY!@#\n@Read2\nTGCA\n+\nQWER\n")
    reads = import_fastq(str(fastq_file))

    assert len(reads) == 2 # check if the two reads were imported from fq file
    assert reads[0].identifier == "Read1"
    assert reads[0].sequence == "ACGT"
    assert reads[1].identifier == "Read2"
    assert reads[1].sequence == "TGCA"

def test_import_fastq_mixed_case(tmp_path):
    fastq_file = tmp_path / "test.fq"
    fastq_file.write_text("@Read1\nacgtN\n+\nIIIII\n@Read2\nTgCan\n+\nIIIII\n")
    reads = import_fastq(str(fastq_file))
    assert [r.sequence for r in reads] == ["ACGTN", "TGCAN"]

def test_import_fastq_lowercase_invalid_still_rejected(tmp_path):
    fastq_file = tmp_path / "test.fq"
    fastq_file.write_text("@Valid\nacgt\n+\nIIII\n@Invalid\nacxt\n+\nIIII\n")
    reads = import_fastq(str(fastq_file))
    assert [r.identifier for r in reads] == ["Valid"]

def test_import_fastq_quality(tmp_path):
    fastq_file = tmp_path / "test.fq"
    fastq_file.write_text("@Valid\nACGT\n+\nIIII\n@Invalid\nACGT\n+\nIII\n")
    reads = import_fastq(str(fastq_file))

    # only one read should be imported because len(ACGT) != len(III):
    assert len(reads) == 1
    assert reads[0].identifier == "Valid"

def test_is_pos_int():
    assert is_pos_int("!!!") is False
    assert is_pos_int(5) is True
    assert is_pos_int(0) is False
    assert is_pos_int(-3) is False
    assert is_pos_int(5.1) is False
    assert is_pos_int("5") is False

def test_missing_file():
    """Test that missing files raise AppErr instead of crashing."""
    with pytest.raises(AppErr):
        import_fasta("ghost_file.fa")

def test_kdb_from_before_the_package_move_gives_clear_error(tmp_path):
    """.kdb files saved before 0.3.0 point at the old top-level 'kmers'
    module. Loading one says to rebuild it."""
    old_kdb = tmp_path / "old.kdb"
    # protocol 2 pickle of kmers.KmerCollection(), as the old layout wrote
    with gzip.open(str(old_kdb), "wb") as f:
        f.write(b"\x80\x02ckmers\nKmerCollection\n)\x81.")
    with pytest.raises(AppErr, match="older version.*rebuild it"):
        load_kdb_file(str(old_kdb))
