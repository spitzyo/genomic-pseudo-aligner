####################        Imports         ####################

import sys
import argparse
from basic_tasks import (reference_task, dumpref_task,
                         align_task, dumpalign_task)
from helper_functions import AppErr

####################         Program          ####################


def readargs(args=None):
    parser = argparse.ArgumentParser(
                    prog='genomic-pseudo-aligner',
    )
    parser.add_argument('-t', '--task',
                        help="task",
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
                        type=int
                        )
    parser.add_argument('-a', '--alignfile',
                        help="Path to the alignment output file (.aln).",
                        )
    parser.add_argument('--reads',
                        help="fastq reads file",
                        )
    parser.add_argument('-m', '--unique-threshold',
                        help="unique k-mer threshold",
                        )
    parser.add_argument('-p', '--ambiguous-threshold',
                        help="ambiguous k-mer threshold",
                        )
    parser.add_argument('--min-read-quality',
                        type=int,
                        )
    parser.add_argument('--min-kmer-quality',
                        type=int,
                        )
    parser.add_argument('--max-genomes',
                        type=int,
                        )
    parser.add_argument('--genomes',
                        )
    parser.add_argument('--coverage',
                        action='store_true',
                        help="Enable coverage extension")
    parser.add_argument('--min-coverage',
                        type=int,
                        default=1,
                        )
    parser.add_argument('--full-coverage',
                        action='store_true',
                        )
    parser.add_argument('--detect-variants',
                        action='store_true',
                        help="Enable variant detection (EXTVARTRACK)")
    parser.add_argument('--min-variant-quality',
                        type=float,
                        help="Minimal quality score for variants")
    parser.add_argument('--min-variant-coverage',
                        type=int,
                        help="Minimal coverage to call changes- variants")
    return parser.parse_args(args)

def main():
    try:
        args = readargs()

        if args.task == "reference":
            reference_task(args)
        elif args.task == "dumpref":
            dumpref_task(args)
        elif args.task == "align":
            align_task(args)
        elif args.task == "dumpalign":
            dumpalign_task(args)
        else:
            print("Error: unknown task")
            sys.exit(1)

    except AppErr as e: # catch the errors raised
        sys.stderr.write(f"{e}\n")
        sys.exit(1)
    except Exception as e:
        sys.stderr.write(f"Unexpected System Error: {e}\n")
        sys.exit(1)

if __name__=="__main__":
    main()