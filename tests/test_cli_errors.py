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
    """Checks the exit status, a single error on stderr, and empty stdout."""
    status, out, err = result
    assert status == code, err
    assert out == ""
    # argparse prints a usage block first; our own errors start the output
    assert err.startswith(ERROR) or err.startswith("usage:"), err
    assert err.count(ERROR) == 1, err
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
