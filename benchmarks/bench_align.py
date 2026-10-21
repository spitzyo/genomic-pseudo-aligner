"""Benchmark the aligner on a synthetic genome.

Builds a random reference genome (1 Mb by default), simulates reads with
substitution errors, then times index building and read alignment with
variant tracking enabled. Only the standard library is used.

Usage (from the repository root):
    python benchmarks/bench_align.py
    python benchmarks/bench_align.py --genome-len 1000000 --reads 1000 --minus-fraction 0.5
"""
import argparse
import os
import random
import sys
import time

# Make the repository root importable when run as a script.
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from aligner import Aligner, Read  # noqa: E402
from kmers import KmerCollection, Reference  # noqa: E402

_COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def simulate(genome_len, n_reads, read_len, error_rate, minus_fraction, seed):
    """Return (genome sequence, list of Read) for a seeded random workload."""
    rng = random.Random(seed)
    genome = "".join(rng.choice("ACGT") for _ in range(genome_len))
    reads = []
    for n in range(n_reads):
        start = rng.randrange(0, genome_len - read_len + 1)
        bases = list(genome[start:start + read_len])
        for i in range(read_len):
            if rng.random() < error_rate:
                bases[i] = rng.choice([b for b in "ACGT" if b != bases[i]])
        seq = "".join(bases)
        quality = [rng.randint(25, 40) for _ in range(read_len)]
        if rng.random() < minus_fraction:
            seq = seq.translate(_COMPLEMENT)[::-1]
            quality = quality[::-1]
        reads.append(Read(f"read{n}", seq, quality))
    return genome, reads


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--genome-len", type=int, default=1_000_000)
    parser.add_argument("--reads", type=int, default=1000)
    parser.add_argument("--read-len", type=int, default=150)
    parser.add_argument("-k", type=int, default=21)
    parser.add_argument("--error-rate", type=float, default=0.01)
    parser.add_argument("--minus-fraction", type=float, default=0.0,
                        help="fraction of reads simulated on the minus strand")
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    genome_seq, reads = simulate(args.genome_len, args.reads, args.read_len,
                                 args.error_rate, args.minus_fraction,
                                 args.seed)

    t0 = time.perf_counter()
    collection = KmerCollection()
    Reference("synthetic", genome_seq).add_ref_kmers(args.k, collection)
    build_s = time.perf_counter() - t0

    aligner = Aligner(collection)
    aligner.enable_vartracking(min_quality=20, min_coverage=1)
    t0 = time.perf_counter()
    for read in reads:
        aligner.align_read(read, args.k)
    align_s = time.perf_counter() - t0

    statuses = {}
    for read in reads:
        statuses[read.status] = statuses.get(read.status, 0) + 1

    print(f"genome length:   {args.genome_len:,} bp (k={args.k})")
    print(f"reads:           {len(reads):,} x {args.read_len} bp, "
          f"{args.minus_fraction:.0%} minus strand")
    print(f"index build:     {build_s:.2f} s")
    print(f"alignment:       {align_s:.2f} s "
          f"({align_s / len(reads) * 1000:.2f} ms/read, "
          f"{len(reads) / align_s:,.0f} reads/s)")
    print(f"read statuses:   {dict(sorted(statuses.items()))}")


if __name__ == "__main__":
    main()
