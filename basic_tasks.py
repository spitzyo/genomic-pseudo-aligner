####################           Imports            ####################
import gzip
import pickle
import json
from typing import Dict
from aligner import Aligner
from kmers import KmerCollection
from helper_functions import (load_kdb_file, import_fasta,
                              save_kmer_collection, is_pos_int,
                              import_fastq, get_kmer_size_from_collection,
                              handle_f_write, handle_f_read)

####################      REFERENCE + DUMPREF     ####################

def build_reference(filename: str, k: int,
                    kmer_collection: KmerCollection) -> None:
    """
    This function builds the reference kmer collection from a FASTA file.
    :param filename: The name of the FASTA file.
    :param k: The kmer size.
    :param kmer_collection: The KmerCollection instance to be added.
    """
    genomes = import_fasta(filename)
    for ref_genome in genomes:
        ref_genome.add_ref_kmers(k, kmer_collection)
        # extracts kmers from every reference genome

# Dump helper functions:
def create_kmer_details(kmer_collection: KmerCollection) -> Dict:
    """This helper function generates kmer details for each from collection."""
    kmer_details = {}
    all_kmers = kmer_collection.get_all_kmers()
    for kmer_obj in all_kmers:
        kmer_seq = kmer_obj.sequence  # getting the str from the instance
        kmer_details[kmer_seq] = {} # inits for every kmer
        for genome in kmer_obj.get_genomes():
            genome_id = genome.identifier
            positions = kmer_obj.get_positions(genome)
            if positions: # if positions exist for the genomes
                if kmer_seq not in kmer_details: # avoiding duplicates
                    kmer_details[kmer_seq] = {}
                kmer_details[kmer_seq][genome_id] = sorted(set(positions))
           # sort numerically and use a set for uniqueness of pos.
    return kmer_details

def create_genome_summary(kmer_collection: KmerCollection) -> Dict:
    """This helper function generates genome summary and statistics."""
    genome_summary = {}  # initializing the final genome details
    for genome in kmer_collection.get_ordered_genomes():
        genome_id = genome.identifier
        base_length = genome.total_bases

        genome_kmers = kmer_collection.get_kmers_for_genome(genome)
        unique_kmers = 0
        multi_map_kmers = 0

        # iterating and counting unique and multi-mapping kmers:
        for kmer_seq in genome_kmers:
            kmer_obj = kmer_collection.get_kmer(kmer_seq)
            # getting the instance from sequence (which is str)
            if kmer_obj.is_specific():  # using kmer's specificity method
                unique_kmers += 1
            else:
                multi_map_kmers += 1
        genome_summary[genome_id] = {"total_bases": base_length,
                                     "unique_kmers": unique_kmers,
                                     "multi_mapping_kmers": multi_map_kmers}
    return genome_summary

@handle_f_write("writing reference dump")
def dump_reference(collection, is_file=True, output_file=None) -> None:
    """
    This function dumps the reference kmer collection to a JSON file.
    :param collection: reference file name (.kdb file) or KmerCollection ins.
    :param is_file: whether the collection arg is a file or not.
    :param output_file: output file name (JSON file).
    """

    if is_file: # if KmerCollection is on a pickle file
        with gzip.open(collection, 'rb'):
            kmer_collection = load_kdb_file(collection)
    else:  # if it was passed directly as an argument
        kmer_collection = collection

    kmer_details = create_kmer_details(kmer_collection) # helper func
    genome_summary = create_genome_summary(kmer_collection) # helper func
    output = {"Kmers": kmer_details, "Summary": genome_summary}

    if output_file: # in case the output should be written to a file
        with open(output_file, 'w') as outfile:
            json.dump(output, outfile, indent=4) # type: ignore
    else: # dumping the output to console, like instructed
        print(json.dumps(output, indent=4))  # Dump the output to the console

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
    if alignfile:  # dumping from an .aln file
        with gzip.open(alignfile, "rb") as f:
            data = pickle.load(f)
            if isinstance(data, tuple) and len(data) == 2:
                reads, alignment_summary = data
            else:
                reads = data
    if not alignfile and not reads:
        print("No alignment file and no reads provided.")
        return # exiting if no reads (to avoid iterator problems)

    # EXTQUALITY (dump filtering statistics) - only included if exist
    if (quality_filtering and aligner and
            hasattr(aligner, '_Aligner__quality_stats')):
        quality_stats = aligner.get_quality_stats()
        reads_stats.update({
            'filtered_quality_reads': quality_stats.get(
                'filtered_quality_reads', 0),
            'filtered_quality_kmers': quality_stats.get(
                'filtered_quality_kmers', 0),
            'filtered_hr_kmers': quality_stats.get('filtered_hr_kmers', 0)
        })

    for read in reads: # calcs stats for every read
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
                # inits a genome_mapping summary
                genome_mapping_summary[genome] = {'unique_reads': 0,
                                                  'ambiguous_reads': 0}
            if status == 'unique':
                genome_mapping_summary[genome]['unique_reads'] += 1
            elif status == 'ambiguous':
                genome_mapping_summary[genome]['ambiguous_reads'] += 1

    # ensuring all genomes, even unmapped ones, while preserving the order
    if alignment_summary:
        # Use the saved alignment summary from the file
        for genome_id, stats in alignment_summary.items():
            if genome_id not in genome_mapping_summary:
                genome_mapping_summary[genome_id] = {
                    'unique_reads': 0,
                    'ambiguous_reads': 0
                }
    elif aligner and hasattr(aligner, 'get_alignment_summary'):
        # Use the aligner's alignment summary
        alignment_summary = aligner.alignment_summary
        for genome_id, stats in alignment_summary.items():
            if genome_id not in genome_mapping_summary:
                genome_mapping_summary[genome_id] = {'unique_reads': 0,
                                                     'ambiguous_reads': 0}

    result = {'Statistics': reads_stats, 'Summary': genome_mapping_summary}

    if output_file: # saving to .aln file if filename is provided
        with open(output_file, "w") as outfile:
            json.dump(result, outfile, indent=4) # type: ignore
    else: # dumping to console otherwise
        print(json.dumps(result, indent=4))

@handle_f_read("loading or building reference db")
def load_or_build_reference(args) -> [KmerCollection, None]:
    """This function loads or builds a kmer collection from ref genomes."""
    # building reference (no file):
    if args.genomefile and args.kmer_size and is_pos_int(args.kmer_size):
        kmer_collection = KmerCollection()  # inits an empty KmerCollection ins
        build_reference(args.genomefile, args.kmer_size, kmer_collection)
    elif args.referencefile: # loading reference
        kmer_collection = load_kdb_file(args.referencefile)  # import
        if not kmer_collection:
            print(f"Error: failed to load reference kmers from "
                  f"kdb file {args.referencefile}.")
            return  # exit
    else:
        print(f"Error: missing arguments to build the kmer reference.")
        return  # exit
    return kmer_collection


####################     Final Tasks      ####################

@handle_f_write("saving reference database")
def reference_task(args) -> None:
    """This function handles the REFERENCE task."""
    # checking all arguments were given to this function:
    if not args.genomefile or not args.referencefile or not args.kmer_size:
        print("Error: missing or invalid arguments for 'reference' task.")
        return
    if not is_pos_int(args.kmer_size):
        print("Error: kmer_size must be a positive integer.")

    kmer_collection = KmerCollection() # initialize kmers obj empty collection
    build_reference(args.genomefile, args.kmer_size, kmer_collection)
    save_kmer_collection(kmer_collection, args.referencefile) # saving the ref

@handle_f_read("loading reference for dumping")
def dumpref_task(args) -> None:
    """This function performs the DUMPREF task."""
    if args.referencefile and args.genomefile is None:
        dump_reference(args.referencefile) # only dumping
        return

    if args.genomefile and args.kmer_size: # building and dumping
        kmer_collection = KmerCollection()
        build_reference(args.genomefile, args.kmer_size, kmer_collection)
        dump_reference(kmer_collection, False)

@handle_f_read("loading alignment (.aln) file")
@handle_f_write("saving alignment results")
def align_task(args):
    """This function performs the ALIGN task."""

    # in any case, the reads should be provided
    if not args.reads:
        print("Error: missing reads file for 'align' task.")
        return

    # Section 3.4: pre-built reference
    if args.referencefile:
        kmer_collection = load_kdb_file(args.referencefile)  # loading kmer ref
        k = get_kmer_size_from_collection(kmer_collection)  # get kmer size

    # Section 3.5: building the reference
    elif args.genomefile and args.kmer_size:
        if not is_pos_int(args.kmer_size):
            return # if k size is not a positive int
        kmer_collection = KmerCollection() # inits a KmerCollection instance
        build_reference(args.genomefile, args.kmer_size, kmer_collection)
        k = args.kmer_size
    else: # in case unable to build or import an existing kdb reference
        print("Error: Unable to build or import a kmer reference.")
        return

    m = int(args.unique_threshold) \
        if (is_pos_int(args.unique_threshold)) else 1
    p = int(args.ambiguous_threhold) \
        if (is_pos_int(args.ambiguous_threhold)) else 1
    # m, p set to default values if not ints / smaller than 0

    aligner = Aligner(kmer_collection) # inits an aligner instance with kmers
    for genome in kmer_collection.get_all_genomes():
        aligner.add_genome_to_summary(genome.identifier)
        # this step ensures unmapped genomes appear in the summary output

    # EXTCOVERAGE
    if args.coverage:
        genome_lens = {}
        for genome in kmer_collection.get_all_genomes():
            genome_lens[genome.identifier] = genome.total_bases
        aligner.enable_coverage(genome_lens, args.min_coverage)

    # EXTVARTRACK
    if args.detect_variants:
        # if no other arguments are given -> default values will be used:
        aligner.enable_vartracking(
            min_quality=args.min_variant_quality if
            args.min_variant_quality else None,
            min_coverage=args.min_variant_coverage
            if args.min_variant_coverage else None)

    reads = import_fastq(args.reads) # loading reads from file
    for read in reads: # iterating and aligning each read
        aligner.align_read(read, k, m, p,
                           MRQ=args.min_read_quality, MKQ=args.min_kmer_quality
                           , MG=args.max_genomes) # execute alignment of read
        # EXTQUALITY: quality filtering would be executed if required

    if args.alignfile:
        with gzip.open(args.alignfile, "wb") as f:
            pickle.dump((
                reads, aligner.alignment_summary), f) # type: ignore

    if args.coverage: # output the coverage
        selected_genomes = args.genomes.split(',') if args.genomes else None
        coverage = aligner.dump_coverage(selected_genomes, args.full_coverage)
        return reads, coverage, aligner

    return reads, aligner # in case needed not in an alignfile

@handle_f_read("loading data for alignment dump")
def dumpalign_task(args) -> None:
    """This function performs the DUMPALIGN task (three diff ways)."""

    # EXTQUALITY
    quality_filtering = args.min_read_quality is not None or \
                             args.min_kmer_quality is not None or \
                             args.max_genomes is not None

    if args.alignfile: # section 3.6
        dump_alignment(alignfile=args.alignfile) # dumping .aln file to console
        return

    # else - kmer collection is retrieved from genomefile+kmer_s, or .kdb file:
    kmer_collection = load_or_build_reference(args) # sections 3.7 / 3.8
    if not kmer_collection:
        return # if kmer collection is invalid

    alignment_results = align_task(args) # perform basic alignment
    reads = alignment_results[0] # as defined in align_task
    aligner = alignment_results[-1] # last item in any case is aligner
    dump_alignment(reads=reads, aligner=aligner,
                   quality_filtering=quality_filtering) # dump

    # Extension output:
    if args.coverage:
        coverage_output = alignment_results[1]
        print(json.dumps(coverage_output, indent=4)) # dump coverage

    if args.detect_variants:
        selected_genomes = args.genomes.split(',') if args.genomes else None
        variants_output = aligner.dump_variants(selected_genomes)
        if variants_output:
            print(json.dumps(variants_output, indent=4))