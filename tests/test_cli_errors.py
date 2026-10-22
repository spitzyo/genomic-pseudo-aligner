"""Exit status and stderr of the command-line interface on failures.

Every failure must exit non-zero with one clear message on stderr, and
leave stdout empty (stdout is for results only).
"""
####################         IMPORTS          ####################
import json
import os
import subprocess
import sys

import pytest

from genomic_pseudo_aligner.cli import main

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ERROR = "genomic-pseudo-aligner: error: "
WARNING = "genomic-pseudo-aligner: warning: "

GENOME = "ACGTACGTTAGCATCGATCGGATCATCGACTAGCATCAGCTAGCATCAG"
READ = GENOME[5:25]


def run_cli(capsys, *argv):
    """Runs the CLI in-process. Returns (exit status, stdout, stderr)."""
    try:
        main([str(a) for a in argv])
        code = 0
    except SystemExit as e:
        code = e.code if isinstance(e.code, int) else 1
    out, err = capsys.readouterr()
    return code, out, err


def assert_fails(result, code, *expected_in_stderr):
    """Checks the exit status, empty stdout, and exactly one error line on
    stderr. Other stderr lines may only be warnings, argparse's usage
    block, or the --help hint after a usage error."""
    status, out, err = result
    assert status == code, err
    assert out == ""
    lines = err.splitlines()
    error_lines = [line for line in lines if line.startswith(ERROR)]
    assert len(error_lines) == 1, err
    for line in lines:
        assert (line.startswith((ERROR, WARNING, "usage:", " ", "Run '"))
                ), f"unexpected stderr line: {line!r}"
    for text in expected_in_stderr:
        assert text in err, err


@pytest.fixture
def files(tmp_path):
    """A small valid FASTA and FASTQ, plus helpers for writing bad ones."""
    (tmp_path / "ref.fa").write_text(f">g1\n{GENOME}\n")
    (tmp_path / "reads.fq").write_text(f"@r1\n{READ}\n+\n{'I' * len(READ)}\n")
    return tmp_path


####################      SUCCESS PATH       ####################

def test_success_writes_only_json_to_stdout(files, capsys):
    status, out, err = run_cli(capsys, "-t", "dumpalign",
                               "-g", files / "ref.fa", "-k", 5,
                               "--reads", files / "reads.fq")
    assert (status, err) == (0, "")
    assert json.loads(out)["Statistics"]["unique_mapped_reads"] == 1


def test_align_without_output_file_warns_on_stderr(files, capsys):
    status, out, err = run_cli(capsys, "-t", "align", "-g", files / "ref.fa",
                               "-k", 5, "--reads", files / "reads.fq")
    assert (status, out) == (0, "")
    assert err.startswith(WARNING) and "no -a given" in err


####################      MISSING FILES      ####################

@pytest.mark.parametrize("argv, missing", [
    (["-t", "reference", "-g", "{d}/missing.fa", "-k", 5, "-r", "{d}/o.kdb"],
     "missing.fa"),
    (["-t", "dumpalign", "-g", "{d}/ref.fa", "-k", 5,
      "--reads", "{d}/missing.fq"], "missing.fq"),
    (["-t", "dumpalign", "-r", "{d}/missing.kdb", "--reads", "{d}/reads.fq"],
     "missing.kdb"),
    (["-t", "dumpref", "-r", "{d}/missing.kdb"], "missing.kdb"),
    (["-t", "dumpalign", "-a", "{d}/missing.aln"], "missing.aln"),
])
def test_missing_file(files, capsys, argv, missing):
    argv = [str(a).format(d=files) for a in argv]
    assert_fails(run_cli(capsys, *argv), 1, missing,
                 "No such file or directory")


def test_missing_aln_is_not_reported_as_a_write_error(files, capsys):
    _, _, err = run_cli(capsys, "-t", "dumpalign", "-a", files / "x.aln")
    assert "writ" not in err


def test_unwritable_output_file(files, capsys):
    assert_fails(run_cli(capsys, "-t", "reference", "-g", files / "ref.fa",
                         "-k", 5, "-r", files / "no_such_dir" / "o.kdb"),
                 1, "o.kdb", "No such file or directory")


####################    MISSING ARGUMENTS    ####################

@pytest.mark.parametrize("argv, expected", [
    (["-t", "bogus"], "invalid choice: 'bogus'"),
    ([], "the following arguments are required: -t"),
    (["-t", "reference", "-g", "{d}/ref.fa", "-k", 5],
     "the reference task needs -g (input FASTA), -k and -r"),
    (["-t", "reference", "-g", "{d}/ref.fa", "-k", 5, "-r", "{d}/out.txt"],
     "must end in .kdb"),
    (["-t", "dumpalign", "-g", "{d}/ref.fa", "-k", 5],
     "dumpalign task needs either -a"),
    (["-t", "align", "--reads", "{d}/reads.fq"], "a reference is needed"),
    (["-t", "align", "-g", "{d}/ref.fa", "--reads", "{d}/reads.fq"],
     "a reference is needed"),
    (["-t", "align", "-g", "{d}/ref.fa", "-k", 5], "needs --reads"),
    (["-t", "dumpref"], "a reference is needed"),
])
def test_invalid_arguments_exit_2(files, capsys, argv, expected):
    argv = [str(a).format(d=files) for a in argv]
    assert_fails(run_cli(capsys, *argv), 2, expected)


####################     REAL PROCESSES      ####################

@pytest.mark.parametrize("entry_point", [
    ["main.py"], ["-m", "genomic_pseudo_aligner"]])
def test_process_exit_status_and_streams(files, entry_point):
    """The same checks through a real process, for both entry points."""
    result = subprocess.run(
        [sys.executable, *entry_point, "-t", "reference",
         "-g", str(files / "missing.fa"), "-k", "5",
         "-r", str(files / "o.kdb")],
        cwd=REPO_ROOT, capture_output=True, text=True)
    assert result.returncode == 1
    assert result.stdout == ""
    assert result.stderr.startswith(ERROR)
    assert "missing.fa: No such file or directory" in result.stderr


####################  INVALID ARGUMENT VALUES  ####################

@pytest.mark.parametrize("option, value, expected", [
    ("-k", "0", "must be a positive integer, got '0'"),
    ("-k", "abc", "must be a positive integer, got 'abc'"),
    ("-k", "2.5", "must be a positive integer, got '2.5'"),
    ("-m", "-1", "must be a positive integer"),
    ("-p", "0", "must be a positive integer"),
    ("--max-genomes", "0", "must be a positive integer"),
    ("--min-coverage", "0", "must be a positive integer"),
    ("--min-read-quality", "-5", "must be an integer >= 0"),
    ("--min-kmer-quality", "x", "must be an integer >= 0"),
    ("--min-variant-coverage", "-1", "must be an integer >= 0"),
    ("--min-variant-quality", "-0.5", "must be a number >= 0"),
    ("--min-variant-af", "2", "must be a number between 0 and 1"),
    ("--min-variant-af", "high", "must be a number between 0 and 1"),
])
def test_invalid_option_value_exit_2(files, capsys, option, value, expected):
    result = run_cli(capsys, "-t", "dumpalign", "-g", files / "ref.fa",
                     "-k", 5, "--reads", files / "reads.fq", option, value)
    assert_fails(result, 2, option, expected)


def test_unknown_genome_fails_before_any_output(files, capsys):
    assert_fails(run_cli(capsys, "-t", "dumpalign", "-g", files / "ref.fa",
                         "-k", 5, "--reads", files / "reads.fq",
                         "--coverage", "--genomes", "g1,nope"),
                 2, "--genomes: not in the reference: nope",
                 "the reference has: g1")


def test_genomes_without_coverage_or_variants_warns(files, capsys):
    status, out, err = run_cli(capsys, "-t", "dumpalign",
                               "-g", files / "ref.fa", "-k", 5,
                               "--reads", files / "reads.fq",
                               "--genomes", "g1")
    assert status == 0 and json.loads(out)
    assert err.startswith(WARNING) and "--genomes has no effect" in err


####################        BAD FASTA        ####################

def _reference(capsys, fasta, k=5):
    return run_cli(capsys, "-t", "reference", "-g", fasta, "-k", k,
                   "-r", fasta.parent / "out.kdb")


@pytest.mark.parametrize("name, content, expected", [
    ("notes.txt", f">g1\n{GENOME}\n", "is not a FASTA file"),
    ("ref.fa.zip", f">g1\n{GENOME}\n", "is not a FASTA file"),
])
def test_fasta_with_wrong_extension_exit_2(files, capsys, name, content,
                                           expected):
    (files / name).write_text(content)
    assert_fails(_reference(capsys, files / name), 2, name, expected)


@pytest.mark.parametrize("content, expected", [
    ("", "has no usable FASTA records (the file is empty)"),
    ("\n\n", "has no usable FASTA records (the file is empty)"),
    (f"{GENOME}\n", "line 1 comes before the first '>' header"),
    (f"@r1\n{READ}\n+\n{'I' * len(READ)}\n",  # a FASTQ named .fa
     "line 1 comes before the first '>' header"),
    (">bad\nILOVEINTRO\n", "has no usable FASTA records"),
    (">g1\n>g2\n", "has no usable FASTA records"),
])
def test_unusable_fasta_exit_1(files, capsys, content, expected):
    (files / "bad.fa").write_text(content)
    assert_fails(_reference(capsys, files / "bad.fa"), 1, "bad.fa", expected)
    assert not (files / "out.kdb").exists()  # no empty reference written


def test_fasta_genomes_all_shorter_than_k_exit_1(files, capsys):
    (files / "short.fa").write_text(">tiny\nACGT\n")
    status, out, err = _reference(capsys, files / "short.fa")
    assert (status, out) == (1, "")
    assert WARNING + "genome 'tiny' (4 bases) is shorter than k=5" in err
    assert ERROR + "no k-mers could be built" in err


def test_fasta_bad_records_are_skipped_with_warnings(files, capsys):
    (files / "mixed.fa").write_text(
        f">good\n{GENOME}\n>bad\nACGT\nACXT\n>\n{GENOME}\n>empty\n")
    status, out, err = run_cli(capsys, "-t", "dumpref",
                               "-g", files / "mixed.fa", "-k", 5)
    assert status == 0
    assert list(json.loads(out)["Summary"]) == ["good"]
    assert err.count(WARNING) == 3, err
    assert "skipped 'bad': line 5 has characters other than A, C, G, T, N" \
        in err
    assert "skipped the record at line 6: empty '>' header" in err
    assert "skipped 'empty': no sequence" in err


def test_soft_masked_fasta_and_mixed_case_reads_align_like_uppercase(
        files, capsys):
    """A soft-masked reference and mixed-case reads give the same result
    as the all-uppercase inputs, with no warnings."""
    masked = GENOME[:10].lower() + GENOME[10:30] + GENOME[30:].lower()
    (files / "masked.fa").write_text(f">g1\n{masked}\n")
    mixed = READ[:5].lower() + READ[5:12] + READ[12:].lower()
    (files / "mixed.fq").write_text(f"@r1\n{mixed}\n+\n{'I' * len(READ)}\n")

    expected = _dumpalign(capsys, files, files / "reads.fq")
    status, out, err = run_cli(capsys, "-t", "dumpalign",
                               "-g", files / "masked.fa", "-k", 5,
                               "--reads", files / "mixed.fq")
    assert (status, err) == (0, "")
    assert json.loads(out) == json.loads(expected[1])
    assert json.loads(out)["Statistics"]["unique_mapped_reads"] == 1


####################        BAD FASTQ        ####################

def _dumpalign(capsys, files, reads):
    return run_cli(capsys, "-t", "dumpalign", "-g", files / "ref.fa",
                   "-k", 5, "--reads", reads)


def test_fastq_with_wrong_extension_exit_2(files, capsys):
    (files / "reads.txt").write_text("@r1\nACGT\n+\nIIII\n")
    assert_fails(_dumpalign(capsys, files, files / "reads.txt"), 2,
                 "reads.txt", "is not a FASTQ file")


@pytest.mark.parametrize("content, expected", [
    ("", "has no valid FASTQ reads (the file is empty)"),
    ("@r1\nACGT\n+\nII\n", "line 1: the quality and sequence lengths differ"),
    ("r1\nACGT\n+\nIIII\n", "line 1: the header doesn't start with '@'"),
    ("@r1\nACGT\nIIII\n", "line 1: the third line doesn't start with '+'"),
    ("@r1\nACXT\n+\nIIII\n", "line 1: the sequence has characters other"),
    ("@r1\n\n+\n\n", "line 1: the sequence is empty"),
    (f">g1\n{GENOME}\n", "line 1: the header doesn't start with '@'"),
])
def test_unusable_fastq_exit_1(files, capsys, content, expected):
    (files / "bad.fq").write_text(content)
    status, out, err = _dumpalign(capsys, files, files / "bad.fq")
    assert (status, out) == (1, "")
    assert ERROR + f"{files / 'bad.fq'} has no valid FASTQ reads" in err
    assert expected in err


def test_fastq_bad_records_are_skipped_with_one_warning(files, capsys):
    good = f"@r1\n{READ}\n+\n{'I' * len(READ)}\n"
    (files / "mixed.fq").write_text(
        good + "@bad1\nACGT\n+\nII\n" + "\n" + good.replace("r1", "r2")
        + "@bad2\nACXT\n+\nIIII\n")
    status, out, err = _dumpalign(capsys, files, files / "mixed.fq")
    assert status == 0
    assert json.loads(out)["Statistics"]["unique_mapped_reads"] == 2
    assert err.count(WARNING) == 1, err
    assert "skipped 2 malformed record(s); first at line 5: the quality " \
        "and sequence lengths differ" in err


####################   UNREADABLE .kdb / .aln   ####################
import gzip
import pickle
from genomic_pseudo_aligner.kmers import KmerCollection


@pytest.fixture
def saved(files, capsys):
    """A valid .kdb and .aln built from the small inputs."""
    assert run_cli(capsys, "-t", "reference", "-g", files / "ref.fa", "-k", 5,
                   "-r", files / "good.kdb")[0] == 0
    assert run_cli(capsys, "-t", "align", "-g", files / "ref.fa", "-k", 5,
                   "--reads", files / "reads.fq",
                   "-a", files / "good.aln")[0] == 0
    capsys.readouterr()
    return files


def _write_gzip(path, data):
    with gzip.open(str(path), "wb") as f:
        f.write(data)


BAD_KDB = {
    "garbage": (lambda p: p.write_bytes(b"not a database"),
                "is not a file saved by genomic-pseudo-aligner"),
    "not_a_pickle": (lambda p: _write_gzip(p, b"hello"),
                     "is corrupt or incomplete"),
    "truncated": (lambda p: p.write_bytes(
                      (p.parent / "good.kdb").read_bytes()[:30]),
                  "is corrupt or incomplete"),
    "alignment": (lambda p: p.write_bytes(
                      (p.parent / "good.aln").read_bytes()),
                  "does not contain a k-mer reference (it holds a tuple)"),
    "empty": (lambda p: _write_gzip(p, pickle.dumps(KmerCollection())),
              "contains no k-mers"),
    "old_version": (lambda p: _write_gzip(
                        p, b"\x80\x02ckmers\nKmerCollection\n)\x81."),
                    "was saved by an older version"),
    "directory": (lambda p: p.mkdir(), "Is a directory"),
}


@pytest.mark.parametrize("kind", list(BAD_KDB))
def test_unreadable_kdb_exit_1(saved, capsys, kind):
    make, expected = BAD_KDB[kind]
    kdb = saved / f"{kind}.kdb"
    make(kdb)
    for argv in (["-t", "dumpalign", "-r", kdb, "--reads",
                  saved / "reads.fq"],
                 ["-t", "dumpref", "-r", kdb]):
        assert_fails(run_cli(capsys, *argv), 1, f"{kind}.kdb", expected)


def test_kdb_with_wrong_extension_exit_2(saved, capsys):
    (saved / "good.db").write_bytes((saved / "good.kdb").read_bytes())
    assert_fails(run_cli(capsys, "-t", "dumpref", "-r", saved / "good.db"),
                 2, "good.db is not a .kdb file")


@pytest.mark.skipif(hasattr(os, "geteuid") and os.geteuid() == 0,
                    reason="root can read files without read permission")
def test_kdb_without_read_permission_exit_1(saved, capsys):
    kdb = saved / "good.kdb"
    kdb.chmod(0)
    try:
        assert_fails(run_cli(capsys, "-t", "dumpref", "-r", kdb), 1,
                     "good.kdb: Permission denied")
    finally:
        kdb.chmod(0o644)


@pytest.mark.parametrize("make, expected", [
    (lambda d: (d / "x.aln").write_bytes((d / "good.kdb").read_bytes()),
     "does not contain an alignment (it holds a KmerCollection)"),
    (lambda d: (d / "x.aln").write_bytes(b"junk"),
     "is not a file saved by genomic-pseudo-aligner"),
])
def test_unreadable_aln_exit_1(saved, capsys, make, expected):
    make(saved)
    assert_fails(run_cli(capsys, "-t", "dumpalign", "-a", saved / "x.aln"),
                 1, "x.aln", expected)


def test_valid_saved_files_still_work(saved, capsys):
    status, out, err = run_cli(capsys, "-t", "dumpalign",
                               "-a", saved / "good.aln")
    assert (status, err) == (0, "")
    assert json.loads(out)["Statistics"]["unique_mapped_reads"] == 1


####################   BINARY / BAD GZIP INPUT   ####################

@pytest.mark.parametrize("name, content, expected", [
    ("fake.fa.gz", f">g1\n{GENOME}\n".encode(),
     "fake.fa.gz is not valid gzip-compressed FASTA"),
    ("cut.fa.gz", gzip.compress(f">g1\n{GENOME}\n".encode())[:25],
     "cut.fa.gz is not valid gzip-compressed FASTA"),
    ("binary.fa", b"\x89PNG\r\n\x1a\n\xff\xfe\x00",
     "binary.fa is not a FASTA text file"),
])
def test_undecodable_fasta_exit_1(files, capsys, name, content, expected):
    (files / name).write_bytes(content)
    assert_fails(_reference(capsys, files / name), 1, expected)


@pytest.mark.parametrize("name, content, expected", [
    ("fake.fq.gz", b"@r1\nACGT\n+\nIIII\n",
     "fake.fq.gz is not valid gzip-compressed FASTQ"),
    ("binary.fq", b"\xff\xfe\x00\x01", "binary.fq is not a FASTQ text file"),
])
def test_undecodable_fastq_exit_1(files, capsys, name, content, expected):
    (files / name).write_bytes(content)
    assert_fails(_dumpalign(capsys, files, files / name), 1, expected)


def test_gzipped_inputs_still_work(files, capsys):
    (files / "ref.fa.gz").write_bytes(
        gzip.compress((files / "ref.fa").read_bytes()))
    (files / "reads.fq.gz").write_bytes(
        gzip.compress((files / "reads.fq").read_bytes()))
    status, out, err = run_cli(capsys, "-t", "dumpalign",
                               "-g", files / "ref.fa.gz", "-k", 5,
                               "--reads", files / "reads.fq.gz")
    assert (status, err) == (0, "")
    assert json.loads(out)["Statistics"]["unique_mapped_reads"] == 1
