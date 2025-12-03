#!/usr/bin/env python3

####################        Imports         ####################

import argparse
from basic_tasks import (reference_task, dumpref_task,
                         align_task, dumpalign_task)

####################        Constants         ####################

nucleo_letters = {'A', 'C', 'G', 'T', 'N'} # a set of possible letters

####################         Program          ####################


def readargs(args=None):
    parser = argparse.ArgumentParser(
                    prog='Biosequence project',
    )
    #General arguments
    parser.add_argument('-t', '--task',
                        help="task",
                        required=True
                        )
    parser.add_argument('-g', '--genomefile',
                        help="Genome fasta file (multiple records)",
                        )
    parser.add_argument('-r', '--referencefile',
                        help="kdb file. Can be either input or"
                             "name for output file",
                        )
    parser.add_argument('-k', '--kmer-size',
                        help="length of kmers",
                        type=int #ensuring k is an int
                        )
    #Task specific arguments
    #align
    parser.add_argument('-a', '--alignfile',
                        help="aln file. Can be either input "
                             "or name for output file",
                        )
    parser.add_argument('--reads',
                        help="fastq reads file",
                        )
    parser.add_argument('-m', '--unique-threshold',
                        help="unique k-mer threshold",
                        )
    parser.add_argument('-p', '--ambiguous-threhold',
                        help="ambiguous k-mer threshold",
                        )
    #align+quality
    parser.add_argument('--min-read-quality',
                        type=int,
                        )
    parser.add_argument('--min-kmer-quality',
                        type=int,
                        )
    parser.add_argument('--max-genomes',
                        type=int,
                        )
    #coverage
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
    # variant detection (new EXTVARTRACK)
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
    args=readargs()

    # the following if statements execute the relevant functions for the tasks:
    if args.task == "reference":
        reference_task(args)
    elif args.task == "dumpref":
        dumpref_task(args)
    elif args.task == "align":
        align_task(args)
    elif args.task == "dumpalign":
        dumpalign_task(args)
    else: # if an invalid task has been called:
        print("Error: unknown task")

if __name__=="__main__":
    main()