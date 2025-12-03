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


def test_reference_task(tmp_path):
    fasta_file = tmp_path / "test.fa"
    fasta_file.write_text(">Seq1\nACGTACGT\n")
    kdb_file = tmp_path / "test.kdb"

    args = readargs(["-t", "reference",
                     "-g", str(fasta_file),
                     "-r", str(kdb_file),
                     "-k", "3"]) # get arguments through readargs
    reference_task(args)
    assert os.path.exists(kdb_file) # check if reference kdb file was created


def test_dumpref_task(tmp_path):
    collection = KmerCollection()
    ref = Reference("test_ref", "ACGTACGT")
    collection.add_kmers(["ACG", "CGT"], ref, [0, 1])
    kdb_file = tmp_path / "test.kdb"
    output_file = tmp_path / "dump.json"
    save_kmer_collection(collection, str(kdb_file)) # save kmers to kdb file
    args = readargs(["-t", "dumpref", "-r", str(kdb_file)]) # get arguments
    dumpref_task(args) # run the function
    original_stdout = sys.stdout # saves whatever was printed to console as var
    with open(output_file, 'w') as f:
        sys.stdout = f # writes the dumped data into the file
        dumpref_task(args) # runs the function again, but info will save to f
    sys.stdout = original_stdout

    with open(output_file, 'r') as f:
        content = f.read() # reading from the file

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
    read.set_status("unique")
    read.set_mapped_genomes(["genome1"])
    output_file = tmp_path / "alignment_dump.json"
    dump_alignment(reads=[read], output_file=str(output_file))
    assert os.path.exists(output_file) # just check if the file is created


def test_dump_alignment_with_aligner(tmp_path):
    read = Read("test_read", "ACGT", "FFFF")
    read.set_status("unique")
    read.set_mapped_genomes(["genome1"])
    collection = KmerCollection()
    aligner = Aligner(collection)
    output_file = tmp_path / "alignment_dump.json"
    dump_alignment(reads=[read], aligner=aligner, output_file=str(output_file))
    assert os.path.exists(output_file) # same here, check if file exists