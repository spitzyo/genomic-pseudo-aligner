"""k-mer based pseudo-aligner with coverage and variant tracking.

Modules:
    kmers            Reference genomes and the k-mer index (KmerCollection)
    aligner          Read and Aligner
    extcoverage      Coverage extension
    extvartrack      Variant tracking extension
    helper_functions File input/output and error handling
    tasks            The reference / dumpref / align / dumpalign tasks
    cli              Command-line interface (entry point: cli.main)
"""
