####################        Imports         ####################

import sys
import argparse
from .tasks import (reference_task, dumpref_task,
                         align_task, dumpalign_task)
from .helper_functions import AppErr, UsageErr, PROG

TASKS = {'reference': reference_task,
         'dumpref': dumpref_task,
         'align': align_task,
         'dumpalign': dumpalign_task}

####################         Program          ####################


def _number(value, convert, kind, check):
    """Converts an argparse value, raising a clear ArgumentTypeError if it
    isn't a number or fails check."""
    try:
        number = convert(value)
    except ValueError:
        number = None
    if number is None or not check(number):
        raise argparse.ArgumentTypeError(f"must be {kind}, got '{value}'")
    return number

def positive_int(value):
    """argparse type: an integer > 0."""
    return _number(value, int, "a positive integer", lambda n: n > 0)

def non_negative_int(value):
    """argparse type: an integer >= 0."""
    return _number(value, int, "an integer >= 0", lambda n: n >= 0)

def non_negative_float(value):
    """argparse type: a number >= 0."""
    return _number(value, float, "a number >= 0", lambda n: n >= 0)

def fraction(value):
    """argparse type: a float between 0 and 1."""
    return _number(value, float, "a number between 0 and 1",
                   lambda n: 0.0 <= n <= 1.0)

def readargs(args=None):
    parser = argparse.ArgumentParser(
                    prog=PROG,
    )
    parser.add_argument('-t', '--task',
                        help="task to run",
                        choices=list(TASKS),
                        required=True
                        )
    parser.add_argument('-g', '--genomefile',
                        help="Path to the reference genome FASTA file.",
                        )
    parser.add_argument('-r', '--referencefile',
                        help="Path to the reference k-mer database (.kdb).",
                        )
    parser.add_argument('-k', '--kmer-size',
                        help="length of kmers",
                        type=positive_int
                        )
    parser.add_argument('-a', '--alignfile',
                        help="Path to the alignment output file (.aln).",
                        )
    parser.add_argument('--reads',
                        help="fastq reads file",
                        )
    parser.add_argument('-m', '--unique-threshold',
                        help="unique k-mer threshold (default 1)",
                        type=positive_int,
                        )
    parser.add_argument('-p', '--ambiguous-threshold',
                        help="ambiguous k-mer threshold (default 1)",
                        type=positive_int,
                        )
    parser.add_argument('--min-read-quality',
                        type=non_negative_int,
                        help="drop reads with a lower mean Phred quality",
                        )
    parser.add_argument('--min-kmer-quality',
                        type=non_negative_int,
                        help="drop k-mers with a lower mean Phred quality",
                        )
    parser.add_argument('--max-genomes',
                        type=positive_int,
                        help="drop k-mers found in more genomes than this",
                        )
    parser.add_argument('--genomes',
                        help="comma-separated genome ids to limit coverage "
                             "and variant output to",
                        )
    parser.add_argument('--coverage',
                        action='store_true',
                        help="Enable coverage extension")
    parser.add_argument('--min-coverage',
                        type=positive_int,
                        default=1,
                        )
    parser.add_argument('--full-coverage',
                        action='store_true',
                        )
    parser.add_argument('--detect-variants',
                        action='store_true',
                        help="Enable variant detection (EXTVARTRACK)")
    parser.add_argument('--min-variant-quality',
                        type=non_negative_float,
                        help="Minimal quality score for variants")
    parser.add_argument('--min-variant-coverage',
                        type=non_negative_int,
                        help="Minimal depth (uniquely mapped reads covering "
                             "the position) to report a variant (default 10)")
    parser.add_argument('--min-variant-af',
                        type=fraction,
                        help="Minimal allele frequency (alt_count / depth, "
                             "0-1) to report a variant (default 0)")
    return parser.parse_args(args)

def main(argv=None):
    """Runs the CLI. Results go to stdout; errors and warnings go to stderr.
    Exit status: 0 on success, 1 when the task fails (missing or invalid
    input file, ...), 2 for invalid command-line arguments."""
    try:
        args = readargs(argv) # argparse exits with status 2 on bad arguments
        TASKS[args.task](args)
    except UsageErr as e:
        print(f"{PROG}: error: {e}\nRun '{PROG} --help' for usage.",
              file=sys.stderr)
        sys.exit(2)
    except AppErr as e: # catch the errors raised
        print(f"{PROG}: error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"{PROG}: error: unexpected {type(e).__name__}: {e}",
              file=sys.stderr)
        sys.exit(1)

if __name__=="__main__":
    main()
