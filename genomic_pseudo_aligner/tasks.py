####################           Imports            ####################
import gzip
import pickle
import json
from typing import Dict, List, Optional, Tuple
from .aligner import Aligner, Read
from .kmers import KmerCollection
from .helper_functions import (load_kdb_file, load_saved_file, import_fasta,
                              save_kmer_collection, is_pos_int,
                              import_fastq, get_kmer_size_from_collection,
                              handle_f_write, handle_f_read, AppErr,
                              UsageErr, warn)

####################      REFERENCE + DUMPREF     ####################

def build_reference(filename: str, k: int,
                    kmer_collection: KmerCollection) -> None:
    """Builds the reference kmer collection from a FASTA file. Warns about
    genomes too short for k, and raises AppErr if no k-mer was built."""
    genomes = import_fasta(filename)
    for ref_genome in genomes:
        if ref_genome.total_bases < k:
            warn(f"genome '{ref_genome.identifier}' ({ref_genome.total_bases}"
                 f" bases) is shorter than k={k}, so it has no k-mers")
        ref_genome.add_ref_kmers(k, kmer_collection)
    if len(kmer_collection.get_all_kmers()) == 0:
        raise AppErr(f"no k-mers could be built from {filename} with k={k} "
                     f"(genomes shorter than k, or only N bases)")

# Dump helper functions:
def create_kmer_details(kmer_collection: KmerCollection) -> Dict:
    """Generates kmer details for each from collection."""
    kmer_details = {}
    all_kmers = kmer_collection.get_all_kmers()
    for kmer_obj in all_kmers:
        kmer_seq = kmer_obj.sequence
        kmer_details[kmer_seq] = {}
        for genome in kmer_obj.get_genomes():
            genome_id = genome.identifier
            positions = kmer_obj.get_positions(genome)
            if positions:
                if kmer_seq not in kmer_details:
                    kmer_details[kmer_seq] = {}
                kmer_details[kmer_seq][genome_id] = sorted(set(positions))
    return kmer_details

def create_genome_summary(kmer_collection: KmerCollection) -> Dict:
    """This helper function generates genome summary and statistics."""
    genome_summary = {}
    for genome in kmer_collection.get_ordered_genomes():
        genome_id = genome.identifier
        base_length = genome.total_bases

        genome_kmers = kmer_collection.get_kmers_for_genome(genome)
        unique_kmers = 0
        multi_map_kmers = 0

        for kmer_seq in genome_kmers:
            kmer_obj = kmer_collection.get_kmer(kmer_seq)
            if kmer_obj.is_specific():
                unique_kmers += 1
            else:
                multi_map_kmers += 1
        genome_summary[genome_id] = {"total_bases": base_length,
                                     "unique_kmers": unique_kmers,
                                     "multi_mapping_kmers": multi_map_kmers}
    return genome_summary

@handle_f_write("writing reference dump")
def dump_reference(collection, is_file=True, output_file=None) -> None:
    """Dumps the reference k-mer collection to a JSON file or console."""

    if is_file:
        kmer_collection = load_kdb_file(collection)
    else:
        kmer_collection = collection

    kmer_details = create_kmer_details(kmer_collection)
    genome_summary = create_genome_summary(kmer_collection)
    output = {"Kmers": kmer_details, "Summary": genome_summary}

    if output_file:
        with open(output_file, 'w') as outfile:
            json.dump(output, outfile, indent=4) # type: ignore
    else:
        print(json.dumps(output, indent=4))

####################      DUMPALIGN HELPERS      ####################

@handle_f_read("loading alignment (.aln) file")
@handle_f_write("writing alignment results")
def dump_alignment(alignfile=None, reads=None, aligner=None,
                   output_file=None, quality_filtering=None) -> None:
    """This function performs the dumping task for alignment."""

    reads_stats = {'unique_mapped_reads': 0,
                   'ambiguous_mapped_reads': 0,
                   'unmapped_reads': 0}
    genome_mapping_summary = {}

    alignment_summary = None
    if alignfile:
        data = load_saved_file(alignfile)
        if isinstance(data, tuple) and len(data) == 2:
            reads, alignment_summary = data
        else:
            reads = data
        if not (isinstance(reads, list)
                and all(isinstance(read, Read) for read in reads)
                and isinstance(alignment_summary, (dict, type(None)))):
            raise AppErr(f"{alignfile} does not contain an alignment (it "
                         f"holds a {type(data).__name__}); create .aln "
                         f"files with the align task")
    if not alignfile and not reads:
        raise AppErr("no alignment file and no reads to dump")

    if quality_filtering and aligner is not None:
        quality_stats = aligner.get_quality_stats()
        reads_stats.update({
            'filtered_quality_reads': quality_stats.get(
                'filtered_quality_reads', 0),
            'filtered_quality_kmers': quality_stats.get(
                'filtered_quality_kmers', 0),
            'filtered_hr_kmers': quality_stats.get('filtered_hr_kmers', 0)
        })

    for read in reads:
        status = read.status
        mapped_genomes = read.mapped_genomes

        if status == 'unique':
            reads_stats['unique_mapped_reads'] += 1
        elif status == 'ambiguous':
            reads_stats['ambiguous_mapped_reads'] += 1
        elif status == 'unmapped':
            reads_stats['unmapped_reads'] += 1

        for genome in mapped_genomes:
            if genome not in genome_mapping_summary:
                genome_mapping_summary[genome] = {'unique_reads': 0,
                                                  'ambiguous_reads': 0}
            if status == 'unique':
                genome_mapping_summary[genome]['unique_reads'] += 1
            elif status == 'ambiguous':
                genome_mapping_summary[genome]['ambiguous_reads'] += 1

    # ensuring all genomes, even unmapped ones, while preserving the order
    if alignment_summary:
        for genome_id, stats in alignment_summary.items():
            if genome_id not in genome_mapping_summary:
                genome_mapping_summary[genome_id] = {
                    'unique_reads': 0,
                    'ambiguous_reads': 0
                }
    elif aligner is not None:
        alignment_summary = aligner.alignment_summary
        for genome_id, stats in alignment_summary.items():
            if genome_id not in genome_mapping_summary:
                genome_mapping_summary[genome_id] = {'unique_reads': 0,
                                                     'ambiguous_reads': 0}

    result = {'Statistics': reads_stats, 'Summary': genome_mapping_summary}

    if output_file:
        with open(output_file, "w") as outfile:
            json.dump(result, outfile, indent=4) # type: ignore
    else:
        print(json.dumps(result, indent=4))

@handle_f_read("loading or building reference db")
def load_or_build_reference(args) -> Tuple[KmerCollection, int]:
    """Loads (-r) or builds (-g with -k) the kmer collection, once.
    Returns (kmer collection, k). If both are given, the pre-built .kdb
    file (-r) wins and k comes from it; -g and -k are then ignored.
    Raises UsageErr when neither is given."""
    if args.referencefile: # pre-built reference
        kmer_collection = load_kdb_file(args.referencefile)
        k = get_kmer_size_from_collection(kmer_collection)
        if not k:
            raise AppErr(f"{args.referencefile} contains no k-mers; "
                         f"rebuild it with the reference task")
    elif args.genomefile and is_pos_int(args.kmer_size): # building it
        kmer_collection = KmerCollection()
        build_reference(args.genomefile, args.kmer_size, kmer_collection)
        k = args.kmer_size
    else:
        raise UsageErr("a reference is needed: either -r (a .kdb file), "
                       "or -g (a FASTA file) with -k")
    return kmer_collection, k

def run_alignment(args, kmer_collection: KmerCollection,
                  k: int) -> Tuple[List, Aligner]:
    """Aligns the reads in args.reads against an already loaded kmer
    collection, with the extensions enabled in args.
    Returns (reads, aligner)."""
    # thresholds are validated by the CLI; default 1 when not given
    unique_diff = args.unique_threshold \
        if is_pos_int(args.unique_threshold) else 1
    total_diff = args.ambiguous_threshold \
        if is_pos_int(args.ambiguous_threshold) else 1

    # FASTA order, so the output lists genomes in the same order every run
    genomes = kmer_collection.get_ordered_genomes()
    aligner = Aligner(kmer_collection)
    for genome in genomes:
        aligner.add_genome_to_summary(genome.identifier)

    # EXTCOVERAGE
    if args.coverage:
        genome_lens = {}
        for genome in genomes:
            genome_lens[genome.identifier] = genome.total_bases
        aligner.enable_coverage(genome_lens, args.min_coverage)

    # EXTVARTRACK
    if args.detect_variants:
        aligner.enable_vartracking(
            min_quality=args.min_variant_quality,
            min_coverage=args.min_variant_coverage,
            min_af=args.min_variant_af)

    reads = import_fastq(args.reads) # loading reads from file
    for read in reads:
        aligner.align_read(read, k, unique_diff, total_diff,
                           min_read_quality=args.min_read_quality,
                           min_kmer_quality=args.min_kmer_quality,
                           max_genomes=args.max_genomes)
        # EXTQUALITY: quality filtering would be executed if required
    return reads, aligner

def _selected_genomes(args) -> Optional[List[str]]:
    """The genome ids given with --genomes, or None for all genomes."""
    return args.genomes.split(',') if args.genomes else None

def _check_selected_genomes(args, kmer_collection: KmerCollection) -> None:
    """Fails before any output if --genomes names a genome that isn't in
    the reference, and warns if --genomes would have no effect."""
    selected = _selected_genomes(args)
    if selected is None:
        return
    if not (args.coverage or args.detect_variants):
        warn("--genomes has no effect without --coverage or "
             "--detect-variants")
    known = [g.identifier for g in kmer_collection.get_ordered_genomes()]
    unknown = [g for g in selected if g not in known]
    if unknown:
        shown = ", ".join(known[:10]) + (", ..." if len(known) > 10 else "")
        raise UsageErr(f"--genomes: not in the reference: "
                       f"{', '.join(unknown)} (the reference has: {shown})")


####################     Final Tasks      ####################

@handle_f_write("saving reference database")
def reference_task(args) -> None:
    """Executes the REFERENCE task."""
    if not args.genomefile or not args.referencefile or not args.kmer_size:
        raise UsageErr("the reference task needs -g (input FASTA), -k and "
                       "-r (output .kdb)")
    if not is_pos_int(args.kmer_size):
        raise UsageErr(f"-k must be a positive integer, got "
                       f"{args.kmer_size}")
    if not args.referencefile.endswith('.kdb'):
        raise UsageErr(f"-r {args.referencefile}: the output file name "
                       f"must end in .kdb")

    kmer_collection = KmerCollection()
    build_reference(args.genomefile, args.kmer_size, kmer_collection)
    save_kmer_collection(kmer_collection, args.referencefile)

@handle_f_read("loading reference for dumping")
def dumpref_task(args) -> None:
    """Performs the DUMPREF task: prints the reference loaded from -r, or
    built from -g and -k (-r wins if both are given)."""
    kmer_collection, _ = load_or_build_reference(args)
    dump_reference(kmer_collection, False)

@handle_f_read("loading alignment (.aln) file")
@handle_f_write("saving alignment results")
def align_task(args) -> Tuple[List, Aligner]:
    """Executes the ALIGN task. Returns (reads, aligner)."""
    if not args.reads:
        raise UsageErr("the align task needs --reads (a FASTQ file)")

    reference = load_or_build_reference(args)
    if not args.alignfile:
        warn("no -a given, so the alignment won't be saved; use dumpalign "
             "to print the results")
    reads, aligner = run_alignment(args, *reference)

    if args.alignfile:
        with gzip.open(args.alignfile, "wb") as f:
            pickle.dump((
                reads, aligner.alignment_summary), f) # type: ignore
    return reads, aligner

@handle_f_read("loading data for alignment dump")
def dumpalign_task(args) -> None:
    """Executes the DUMPALIGN task: either dumps a saved alignment (-a), or
    aligns the reads (loading the reference once) and dumps the result."""

    if args.alignfile:
        dump_alignment(alignfile=args.alignfile)
        return

    if not args.reads:
        raise UsageErr("the dumpalign task needs either -a (a saved .aln "
                       "file), or --reads with a reference")

    # EXTQUALITY
    quality_filtering = args.min_read_quality is not None or \
                             args.min_kmer_quality is not None or \
                             args.max_genomes is not None

    # kmer collection is retrieved from genomefile+kmer_s, or .kdb file:
    reference = load_or_build_reference(args) # sections 3.7 / 3.8
    _check_selected_genomes(args, reference[0])
    reads, aligner = run_alignment(args, *reference)
    dump_alignment(reads=reads, aligner=aligner,
                   quality_filtering=quality_filtering)

    if args.coverage:
        coverage_output = aligner.dump_coverage(_selected_genomes(args),
                                                args.full_coverage)
        print(json.dumps(coverage_output, indent=4)) # dump coverage

    if args.detect_variants:
        variants_output = aligner.dump_variants(_selected_genomes(args))
        if variants_output:
            print(json.dumps(variants_output, indent=4))
