####################         IMPORTS          ####################
import os
import sys
import json
import pytest
from genomic_pseudo_aligner.kmers import Reference, KmerCollection
from genomic_pseudo_aligner.aligner import Read, Aligner
from genomic_pseudo_aligner.helper_functions import (save_kmer_collection,
                                                      UsageErr)
from genomic_pseudo_aligner.cli import readargs
from genomic_pseudo_aligner.tasks import (build_reference, dump_reference,
                                          reference_task, dumpref_task,
                                          create_kmer_details,
                                          create_genome_summary,
                                          dump_alignment)

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

def test_dump_alignment_reports_quality_stats_and_all_genomes(capsys):
    read = Read("test_read", "ACGT", "FFFF")
    read.status = "unique"
    read.mapped_genomes = ["genome1"]
    aligner = Aligner(KmerCollection())
    aligner.add_genome_to_summary("genome1")
    aligner.add_genome_to_summary("genome2")  # no reads mapped here
    aligner._quality_stats['filtered_quality_reads'] = 3
    dump_alignment(reads=[read], aligner=aligner, quality_filtering=True)

    result = json.loads(capsys.readouterr().out)
    assert result["Statistics"]["filtered_quality_reads"] == 3
    assert result["Summary"]["genome2"] == {"unique_reads": 0,
                                            "ambiguous_reads": 0}


def test_threshold_arguments_are_integers():
    args = readargs(["-t", "align", "-m", "3", "-p", "2"])
    assert (args.unique_threshold, args.ambiguous_threshold) == (3, 2)


####################   DUMPALIGN LOADS THE REFERENCE ONCE   ####################
from genomic_pseudo_aligner import tasks
from genomic_pseudo_aligner.tasks import (
    dumpalign_task, align_task, load_or_build_reference)

READS_FQ = "@r1\nACGTACGTAC\n+\nIIIIIIIIII\n"


def _count_calls(monkeypatch, name):
    """Wraps tasks.<name> so each call is recorded."""
    calls = []
    original = getattr(tasks, name)
    def wrapper(*args, **kwargs):
        calls.append(args)
        return original(*args, **kwargs)
    monkeypatch.setattr(tasks, name, wrapper)
    return calls


def _inputs(tmp_path):
    fasta = tmp_path / "ref.fa"
    fasta.write_text(">from_fasta\nACGTACGTACGT\n")
    reads = tmp_path / "reads.fq"
    reads.write_text(READS_FQ)
    return str(fasta), str(reads)


def test_dumpalign_builds_fasta_reference_once(tmp_path, monkeypatch, capsys):
    fasta, reads = _inputs(tmp_path)
    builds = _count_calls(monkeypatch, "build_reference")
    dumpalign_task(readargs(["-t", "dumpalign", "-g", fasta, "-k", "3",
                             "--reads", reads, "--coverage"]))
    out = capsys.readouterr().out
    assert len(builds) == 1
    assert "unique_mapped_reads" in out and "Coverage" in out


def test_dumpalign_loads_kdb_once(tmp_path, monkeypatch, capsys):
    fasta, reads = _inputs(tmp_path)
    kdb = str(tmp_path / "ref.kdb")
    reference_task(readargs(["-t", "reference", "-g", fasta, "-k", "3",
                             "-r", kdb]))
    loads = _count_calls(monkeypatch, "load_kdb_file")
    dumpalign_task(readargs(["-t", "dumpalign", "-r", kdb,
                             "--reads", reads]))
    assert len(loads) == 1
    assert "unique_mapped_reads" in capsys.readouterr().out


def test_kdb_wins_over_fasta_when_both_given(tmp_path, monkeypatch, capsys):
    fasta, reads = _inputs(tmp_path)
    other_fasta = tmp_path / "other.fa"
    other_fasta.write_text(">from_kdb\nACGTACGTACGT\n")
    kdb = str(tmp_path / "ref.kdb")
    reference_task(readargs(["-t", "reference", "-g", str(other_fasta),
                             "-k", "3", "-r", kdb]))
    builds = _count_calls(monkeypatch, "build_reference")

    collection, k = load_or_build_reference(
        readargs(["-t", "align", "-g", fasta, "-k", "5", "-r", kdb]))
    assert k == 3  # from the .kdb, not -k
    assert [g.identifier for g in collection.get_ordered_genomes()] == [
        "from_kdb"]

    dumpalign_task(readargs(["-t", "dumpalign", "-g", fasta, "-k", "3",
                             "-r", kdb, "--reads", reads]))
    out = capsys.readouterr().out
    assert builds == []
    assert "from_kdb" in out and "from_fasta" not in out


def test_dumpalign_without_reads_raises_usage_error(tmp_path, capsys):
    fasta, _ = _inputs(tmp_path)
    with pytest.raises(UsageErr, match="dumpalign task needs"):
        dumpalign_task(readargs(["-t", "dumpalign", "-g", fasta, "-k",
                                 "3"]))
    assert capsys.readouterr().out == ""


def test_missing_reference_raises_usage_error(tmp_path, capsys):
    _, reads = _inputs(tmp_path)
    with pytest.raises(UsageErr, match="a reference is needed"):
        align_task(readargs(["-t", "align", "--reads", reads]))
    assert capsys.readouterr().out == ""


def test_align_saves_aln_that_dumpalign_reads(tmp_path, capsys):
    fasta, reads = _inputs(tmp_path)
    aln = str(tmp_path / "out.aln")
    result = align_task(readargs(["-t", "align", "-g", fasta, "-k", "3",
                                  "--reads", reads, "-a", aln]))
    reads_out, aligner = result
    assert reads_out[0].status == "unique"
    capsys.readouterr()
    dumpalign_task(readargs(["-t", "dumpalign", "-a", aln]))
    assert '"unique_mapped_reads": 1' in capsys.readouterr().out


def test_dumpalign_lists_genomes_in_fasta_order(tmp_path, capsys):
    """Coverage lists genomes in FASTA order. So does Summary for genomes
    without mapped reads (here: all of them); genomes with mapped reads
    are listed first, in the order their reads appear."""
    import random
    rng = random.Random(0)
    names = [f"genome{i}" for i in (5, 1, 4, 0, 3, 2)]
    fasta = tmp_path / "ref.fa"
    fasta.write_text("".join(
        f">{name}\n{''.join(rng.choice('ACGT') for _ in range(40))}\n"
        for name in names))
    reads = tmp_path / "reads.fq"
    reads.write_text(READS_FQ)
    dumpalign_task(readargs(["-t", "dumpalign", "-g", str(fasta), "-k", "8",
                             "--reads", str(reads), "--coverage"]))
    decoder = json.JSONDecoder()
    text = capsys.readouterr().out
    alignment, end = decoder.raw_decode(text)
    coverage, _ = decoder.raw_decode(text[end:].lstrip())
    assert list(alignment["Summary"]) == names
    assert list(coverage["Coverage"]) == names
