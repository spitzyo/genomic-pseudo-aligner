####################         IMPORTS          ####################
import os
import sys
from kmers import Reference, KmerCollection
from aligner import Read, Aligner
from helper_functions import save_kmer_collection
from main import readargs
from basic_tasks import (build_reference, dump_reference, reference_task,
                         dumpref_task, create_kmer_details,
                         create_genome_summary, dump_alignment)

####################          TESTS           ####################

def test_build_reference(tmp_path):
    fasta_file = tmp_path / "test.fa"
    fasta_file.write_text(">Seq1\nACGTACGT\n>Seq2\nTGCATGCA\n")
    collection = KmerCollection()
    build_reference(str(fasta_file), 3, collection) # using k_size: 3

    assert len(collection.get_all_genomes()) == 2 # checks if all were loaded
    assert collection.get_kmer("ACG") is not None
    assert collection.get_kmer("TGC") is not None
    assert collection.get_kmer("AAA") is None # shouldn't exist in collection


def test_create_kmer_details():
    collection = KmerCollection() # inits an empty kmer collection instance
    ref = Reference("genome1", "ACGTACGT")
    collection.add_kmers(["ACG", "CGT"], ref, [0, 1])
    details = create_kmer_details(collection)

    assert "ACG" in details # check if ACG exists in the function's output
    assert "genome1" in details["ACG"]
    assert details["ACG"]["genome1"] == [0] # ACG is in position 0 in details


def test_create_genome_summary():
    collection = KmerCollection() # an empty kmer_collection instance
    ref1 = Reference("genome1", "ACGTACGT")
    ref2 = Reference("genome2", "AACGAAAA")
    collection.add_kmers(["ACG", "CGT"], ref1, [0, 1])
    collection.add_kmers(["AAC", "ACG"], ref2, [0, 1])
    # "ACG" exists in both genomes, and both have on unique kmer in addition

    summary = create_genome_summary(collection)
    assert summary["genome1"]["total_bases"] == 8 # count genome1 length
    assert summary["genome1"]["unique_kmers"] == 1 # only GCT is unique
    assert summary["genome2"]["unique_kmers"] == 1  # only AAC is unique
    assert summary["genome2"]["multi_mapping_kmers"] == 1 # "ACG" is shared
    # No soft-masking info was set on these References → expect zero values
    assert summary["genome1"]["soft_masked_bases"] == 0
    assert summary["genome1"]["soft_masked_fraction"] == 0.0

def test_create_genome_summary_with_masking():
    """Genome summary must reflect soft_masked_bases when present."""
    # 2 out of 8 bases are soft-masked (positions 0 and 1)
    ref = Reference("masked_genome", "ACGTACGT", masked_bases={0, 1})
    collection = KmerCollection()
    collection.add_kmers(["ACG"], ref, [0])
    summary = create_genome_summary(collection)
    assert summary["masked_genome"]["soft_masked_bases"] == 2
    assert summary["masked_genome"]["soft_masked_fraction"] == 0.25


def test_dump_alignment(tmp_path):
    read = Read("test_read", "ACGT", "FFFF")
    read.status = "unique"
    read.mapped_genomes = ["genome1"]
    output_file = tmp_path / "alignment_dump.json"
    dump_alignment(reads=[read], output_file=str(output_file))
    assert os.path.exists(output_file) # just check if the file is created


def test_dump_alignment_soft_masking_stats(tmp_path):
    """dump_alignment must report total_read_bases and soft-masking fractions."""
    import json
    # Read with 2 out of 4 bases soft-masked
    r1 = Read("r1", "ACGT", [40, 40, 40, 40],
              soft_masked_positions={0, 1})
    r1.status = "unique"
    r1.mapped_genomes = ["g1"]
    # Read with no soft-masking
    r2 = Read("r2", "TGCA", [40, 40, 40, 40])
    r2.status = "unmapped"
    output_file = tmp_path / "dump.json"
    dump_alignment(reads=[r1, r2], output_file=str(output_file))
    with open(output_file) as fh:
        result = json.load(fh)
    stats = result["Statistics"]
    assert stats["total_read_bases"] == 8        # 4 + 4
    assert stats["soft_masked_read_bases"] == 2  # only r1 had masking
    assert stats["soft_masked_read_fraction"] == 0.25  # 2 / 8


def test_dump_alignment_with_aligner(tmp_path):
    read = Read("test_read", "ACGT", "FFFF")
    read.status = "unique"
    read.mapped_genomes = ["genome1"]
    collection = KmerCollection()
    aligner = Aligner(collection)
    output_file = tmp_path / "alignment_dump.json"
    dump_alignment(reads=[read], aligner=aligner, output_file=str(output_file))
    assert os.path.exists(output_file) # same here, check if file exists
    fasta_file = tmp_path / "test.fa"
    fasta_file.write_text(">Seq1\nACGTACGT\n")
    kdb_file = tmp_path / "test.kdb"

    args = readargs(["-t", "reference",
                     "-g", str(fasta_file),
                     "-r", str(kdb_file),
                     "-k", "3"]) # get arguments through readargs
    reference_task(args)
    assert os.path.exists(kdb_file) # check if reference kdb file was created


def test_dumpref_task(tmp_path, capsys):
    collection = KmerCollection()
    ref = Reference("test_ref", "ACGTACGT")
    collection.add_kmers(["ACG", "CGT"], ref, [0, 1])
    kdb_file = tmp_path / "test.kdb"
    output_file = tmp_path / "dump.json"
    save_kmer_collection(collection, str(kdb_file))  # save kmers to kdb file
    args = readargs(["-t", "dumpref", "-r", str(kdb_file)])  # get arguments
    dumpref_task(args)
    captured = capsys.readouterr()  # Capture the output
    content = captured.out

    assert "Kmers" in content # checking the reference is in correct format
    assert "Summary" in content # same format check (general)
    assert "ACG" in content # checking for the specific kmer that should be
    assert "test_ref" in content # checking for the specific reference genome


def test_dump_reference_to_file(tmp_path):
    collection = KmerCollection()
    ref = Reference("test_ref", "ACGTACGT")
    collection.add_kmers(["ACG", "CGT"], ref, [0, 1])
    output_file = tmp_path / "dump.json"
    dump_reference(collection, is_file=False, output_file=str(output_file))

    # Just check if file has been created
    assert os.path.exists(output_file)


def test_dump_alignment(tmp_path):
    read = Read("test_read", "ACGT", "FFFF")
    read.status = "unique"
    read.mapped_genomes = ["genome1"]
    output_file = tmp_path / "alignment_dump.json"
    dump_alignment(reads=[read], output_file=str(output_file))
    assert os.path.exists(output_file) # just check if the file is created


def test_dump_alignment_with_aligner(tmp_path):
    read = Read("test_read", "ACGT", "FFFF")
    read.status = "unique"
    read.mapped_genomes = ["genome1"]
    collection = KmerCollection()
    aligner = Aligner(collection)
    output_file = tmp_path / "alignment_dump.json"
    dump_alignment(reads=[read], aligner=aligner, output_file=str(output_file))
    assert os.path.exists(output_file) # same here, check if file exists