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
    This format is usually used for bacterial genomes."""
    if not filename.endswith(('.fa', '.fasta', '.fa.gz', '.fasta.gz')):
        raise AppErr(f"Error: {filename} is not a FASTA-formatted file.")

    def read_fasta_helper():
        current_id = ''
        sequence = []
        with open_gz_file(filename, "r") as f:
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if current_id and sequence: # if valid
                        sequence = ''.join(sequence)
                        if sequence:
                            yield Reference(current_id, sequence)
                            # only if something was added -> create Ref inst
                    current_id = line[1:].strip()
                    sequence = [] # after creating instance -> inits a new
                else:
                    if all(char in NUCLEO_LETTERS for char in line):
                        sequence.append(line) # only if valid DNA letters
            if current_id and sequence:
                sequence = ''.join(sequence) #last
                if sequence:
                    yield Reference(current_id, sequence)
    return list(read_fasta_helper())

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
    This format is usually used for NGS-sequenced genomes."""
    if not filename.endswith(('.fq', '.fastq', '.fq.gz', '.fastq.gz')):
        raise AppErr(f"Error: {filename} is not a FASTQ-formatted file.")

    def read_fastq_chunks():
        """This helper function basically iterates over 4 line chunks
        from the FASTQ file in an efficient manner."""
        with open_gz_file(filename, "r") as f:
            while True:
                header = f.readline().strip()
                if not header:
                    break
                sequence = f.readline().strip()
                plus = f.readline().strip()
                quality = f.readline().strip()

                if not header.startswith("@") or not plus.startswith("+"):
                    continue  # invalid format
                if not all(char in NUCLEO_LETTERS for char in sequence):
                    continue  # invalid sequence
                quality_str = [ord(char) - 33 for char in quality.strip()]
                # converting using a Phred33 format
                if len(sequence) != len(quality_str):
                    continue  # quality str and seq don't match in length
                yield Read(header[1:], sequence, quality_str)
                # generator: returns and will continue to next in the next call

    return list(read_fastq_chunks())

####################        File Writing         ####################

@handle_f_write("saving kmer collection")
def save_kmer_collection(kmer_collection: KmerCollection,
                         output_file: str) -> None:
    """This function saves the kmer collection in a pickle file.
    It uses gzip compression to open the file in binary mode."""
    if not output_file.endswith('.kdb'):
        raise AppErr(f"Error: {output_file} is not a KDB file.")
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