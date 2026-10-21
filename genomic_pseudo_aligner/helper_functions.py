####################        Imports         ####################
import functools
import gzip
import pickle
from typing import List, Callable, Optional
from .aligner import Read
from .kmers import Reference, KmerCollection
import sys


####################        Constants         ####################
NUCLEO_LETTERS = {'A', 'C', 'G', 'T', 'N'} # a set of possible letters
PROG = "genomic-pseudo-aligner" # name used in error and warning messages

####################     Error Handling       ####################

class AppErr(Exception):
    """An error to report to the user. The CLI prints the message to stderr
    and exits with status 1."""
    pass

class UsageErr(AppErr):
    """A problem with the command-line arguments (exit status 2)."""
    pass

def warn(message: str) -> None:
    """Prints a warning to stderr, keeping stdout for results only."""
    print(f"{PROG}: warning: {message}", file=sys.stderr)

def handle_file_errors(operation_desc: str) -> Callable:
    """
    This decorator turns file errors into an AppErr with a clear message,
    e.g. "ref.fa: No such file or directory (while loading FASTA file)".
    An AppErr raised inside is passed through unchanged, so nested
    decorated calls don't stack their descriptions.
    :param operation_desc: for example, "loading reference databases", etc.
    """
    def decorator(func):
        @functools.wraps(func) # maintain function attrs after wrapping
        def wrapper(*args, **kwargs):
            try:
                return func(*args, **kwargs)
            except AppErr:
                raise
            except OSError as e:
                reason = e.strerror or str(e)
                target = f"{e.filename}: " if e.filename else ""
                raise AppErr(f"{target}{reason} (while {operation_desc})") \
                    from e
            except Exception as e:
                raise AppErr(f"unexpected {type(e).__name__} while "
                             f"{operation_desc}: {e}") from e
        return wrapper
    return decorator

# Older names, kept for the decorated functions and any external callers.
handle_f_read = handle_file_errors
handle_f_write = handle_file_errors

####################        File Reading         ####################

def open_gz_file(filename: str, mode='r'):
    """This helper function is used to open a gzip or regular file."""
    if filename.endswith('.gz'):
        return gzip.open(filename, mode + 't') # mode in gzip format
    else:
        return open(filename, mode)

@handle_f_read("loading FASTA file")
def import_fasta(filename) -> List[Reference]:
    """This function imports sequences from a FASTA-formatted file.
    This format is usually used for bacterial genomes.
    Records with an empty header, no sequence, or characters other than
    A, C, G, T, N are skipped with a warning on stderr. Raises AppErr if
    the file isn't FASTA or has no usable record."""
    if not filename.endswith(('.fa', '.fasta', '.fa.gz', '.fasta.gz')):
        raise UsageErr(f"{filename} is not a FASTA file (expected a name "
                       f"ending in .fa, .fasta, .fa.gz or .fasta.gz)")

    # each record: [identifier, sequence lines, first bad line number or
    # None, header line number]
    records = []
    with open_gz_file(filename, "r") as f:
        for line_no, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                records.append([line[1:].strip(), [], None, line_no])
            elif not records:
                raise AppErr(f"{filename} is not in FASTA format: line "
                             f"{line_no} comes before the first '>' header")
            else:
                record = records[-1]
                if record[2] is None and not set(line) <= NUCLEO_LETTERS:
                    record[2] = line_no
                record[1].append(line)

    genomes = []
    for identifier, lines, bad_line, header_line in records:
        if not identifier:
            warn(f"{filename}: skipped the record at line {header_line}: "
                 f"empty '>' header")
        elif bad_line is not None:
            warn(f"{filename}: skipped '{identifier}': line {bad_line} has "
                 f"characters other than A, C, G, T, N")
        elif not lines:
            warn(f"{filename}: skipped '{identifier}': no sequence")
        else:
            genomes.append(Reference(identifier, ''.join(lines)))
    if not genomes:
        reason = "see the warnings above" if records else "the file is empty"
        raise AppErr(f"{filename} has no usable FASTA records ({reason})")
    return genomes

@handle_file_errors("reading a saved .kdb/.aln file")
def load_saved_file(filename: str):
    """Unpickles a gzipped .kdb / .aln file. Files saved before version
    0.3.0 refer to the old top-level modules (kmers, aligner) and can't be
    read any more; they raise an AppErr saying to rebuild them."""
    with gzip.open(filename, "rb") as f:
        try:
            return pickle.load(f)
        except (ModuleNotFoundError, AttributeError) as e:
            raise AppErr(f"Error: {filename} was saved by an older version "
                         f"of genomic-pseudo-aligner and can't be read "
                         f"({e}). Rebuild it with this version.")

@handle_f_read("loading reference databases")
def load_kdb_file(filename: str) -> KmerCollection:
    """
    This function loads a previously saved kmer collection from a pickle file,
    and returns it as a KmerCollection instance if not corrupt.
    """
    if not filename.endswith('.kdb'):
        raise AppErr(f"Error: {filename} is not a KDB file.")

    return load_saved_file(filename)

@handle_f_read("loading FASTQ file")
def import_fastq(filename) -> List[Read]:
    """This function imports sequences from a FASTQ-formatted file.
    This format is usually used for NGS-sequenced genomes.
    Malformed records are skipped, with one warning on stderr giving the
    count and the first problem. Raises AppErr if no read is usable."""
    if not filename.endswith(('.fq', '.fastq', '.fq.gz', '.fastq.gz')):
        raise UsageErr(f"{filename} is not a FASTQ file (expected a name "
                       f"ending in .fq, .fastq, .fq.gz or .fastq.gz)")

    reads = []
    skipped = 0
    first_problem = None
    with open_gz_file(filename, "r") as f:
        line_no = 0
        while True:
            header = f.readline()
            if not header: # end of file
                break
            line_no += 1
            header = header.strip()
            if not header: # blank line between records
                continue
            record_line = line_no
            sequence, plus, quality = (f.readline().strip()
                                       for _ in range(3))
            line_no += 3

            problem = None
            if not header.startswith("@"):
                problem = "the header doesn't start with '@'"
            elif not plus.startswith("+"):
                problem = "the third line doesn't start with '+'"
            elif not sequence:
                problem = "the sequence is empty"
            elif not set(sequence) <= NUCLEO_LETTERS:
                problem = ("the sequence has characters other than "
                           "A, C, G, T, N")
            elif len(quality) != len(sequence):
                problem = "the quality and sequence lengths differ"
            if problem:
                skipped += 1
                if first_problem is None:
                    first_problem = f"line {record_line}: {problem}"
                continue
            # converting using a Phred33 format
            quality_str = [ord(char) - 33 for char in quality]
            reads.append(Read(header[1:], sequence, quality_str))

    if skipped:
        warn(f"{filename}: skipped {skipped} malformed record(s); first at "
             f"{first_problem}")
    if not reads:
        reason = "see the warning above" if skipped else "the file is empty"
        raise AppErr(f"{filename} has no valid FASTQ reads ({reason})")
    return reads

####################        File Writing         ####################

@handle_f_write("saving kmer collection")
def save_kmer_collection(kmer_collection: KmerCollection,
                         output_file: str) -> None:
    """This function saves the kmer collection in a pickle file.
    It uses gzip compression to open the file in binary mode."""
    if not output_file.endswith('.kdb'):
        raise UsageErr(f"{output_file}: the output file name must end in "
                       f".kdb")
    with gzip.open(output_file, "wb") as f:
        pickle.dump(kmer_collection, f) # type: ignore


####################            Other              ####################

def is_pos_int(var) -> bool:
    """This function checks if a variable is a positive integer."""
    return isinstance(var, int) and var > 0

def get_kmer_size_from_collection(
        kmer_collection: KmerCollection) -> Optional[int]:
    """This function extracts the kmer size from the collection, or returns
    None if the collection has no k-mers."""
    if len(kmer_collection.get_all_kmers()) > 0:
        return len(next(iter(kmer_collection.get_all_kmers())).sequence)
    return None