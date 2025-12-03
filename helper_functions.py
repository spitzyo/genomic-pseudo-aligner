####################        Imports         ####################
import functools
import gzip
import pickle
from typing import List, Callable
from aligner import Read
from kmers import Reference, KmerCollection
import sys


####################        Constants         ####################
nucleo_letters = {'A', 'C', 'G', 'T', 'N'} # a set of possible letters

####################     Error Handling       ####################

def handle_f_read(operation_desc: str) -> Callable:
    """
    This decorator is used to handle errors of file read operations.
    It is used throughout the program with different files that need opening.
    It treats two common errors, as well as a general exception, and exits
    the program "gracefully" with exit code 1 (that indicates a failure).
    :param operation_desc: for example, "loading reference databases", etc.
    :return: the decorated function with proper file reading errors' handling.
    """
    def decorator(func):
        @functools.wraps(func) # maintain function attrs after wrapping
        def wrapper(*args, **kwargs):
            try:
                return func(*args, **kwargs) # if no error
            except FileNotFoundError:
                print(f"Error: File not found while {operation_desc}.")
                sys.exit(1)
            except PermissionError:
                print(f"Error: Permission denied to "
                      f"file while {operation_desc}.")
                sys.exit(1)
            except Exception as e:
                print(f"Error while {operation_desc}: {str(e)}")
                sys.exit(1)
        return wrapper
    return decorator

def handle_f_write(operation_desc: str) -> Callable:
    """
    This decorator is used to handle errors of file write operations.
    It treats two common errors, as well as a general exception, and exits
    the program gracefully with exit code 1.
    :param operation_desc: for example, "writing alignment results", etc.
    :return: the decorated function with proper file writing errors' handling.
    """
    def decorator(func):
        @functools.wraps(func) # maintain function attrs after wrapping
        def wrapper(*args, **kwargs):
            try:
                return func(*args, **kwargs) # no error
            except PermissionError:
                print(f"Error: Permission "
                      f"denied to file while {operation_desc}.")
                sys.exit(1)
            except OSError: # a broader exception than PermissionError
                print(f"Error: OS error occurred while {operation_desc}.")
                sys.exit(1)
            except Exception as e:
                print(f"Error while {operation_desc}: {str(e)}")
                sys.exit(1)
        return wrapper
    return decorator

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
    # check if it is indeed a FASTA file:
    if not filename.endswith(('.fa', '.fasta', '.fa.gz', '.fasta.gz')):
        print(f"Error: {filename} is not a FASTA-formatted file.")
        sys.exit(1) # exit the program

    def read_fasta_helper():
        # initializers:
        current_id = ''
        sequence = []
        with open_gz_file(filename, "r") as f:
            for line in f:
                line = line.strip() # removes whitespaces and so on
                if line.startswith(">"):
                    if current_id and sequence: # if valid
                        sequence = ''.join(sequence)
                        if sequence:
                            yield Reference(current_id, sequence)
                            # only if something was added -> create Ref inst
                    current_id = line[1:].strip()
                    sequence = [] # after creating instance -> inits a new
                else:
                    if all(char in nucleo_letters for char in line):
                        sequence.append(line) # only if valid DNA letters
            if current_id and sequence:
                sequence = ''.join(sequence) #last
                if sequence:
                    yield Reference(current_id, sequence)
    return list(read_fasta_helper())

@handle_f_read("loading reference databases")
def load_kdb_file(filename: str) -> KmerCollection:
    """
    This function loads a previously saved kmer collection from a pickle file,
    and returns it as a KmerCollection instance if not corrupt.
    """
    # check if it is indeed a FASTA file:
    if not filename.endswith('.kdb'):
        print(f"Error: {filename} is not a KDB file.")
        sys.exit(1) # exit the program

    with gzip.open(filename, "rb") as f:
        kmer_collection = pickle.load(f)
        return kmer_collection

@handle_f_read("loading FASTQ file")
def import_fastq(filename) -> List[Read]:
    """This function imports sequences from a FASTQ-formatted file.
    This format is usually used for NGS-sequenced genomes."""

    # check if it is indeed a FASTA file:
    if not filename.endswith(('.fq', '.fastq', '.fq.gz', '.fastq.gz')):
        print(f"Error: {filename} is not a FASTQ-formatted file.")
        sys.exit(1) # exit the program

    def read_fastq_chunks():
        """This helper function basically iterates over 4 line chunks
        from the FASTQ file in an efficient manner."""
        with open_gz_file(filename, "r") as f:
            while True:
                header = f.readline().strip()
                if not header:
                    break # end of file (EOF)
                sequence = f.readline().strip()
                plus = f.readline().strip()
                quality = f.readline().strip()

                if not header.startswith("@") or not plus.startswith("+"):
                    continue  # invalid format
                if not all(char in nucleo_letters for char in sequence):
                    continue  # invalid sequence
                quality_str = [ord(char) - 33 for char in quality.strip()]
                # converting using a Phred33 format
                if len(sequence) != len(quality_str):
                    continue  # quality str and seq don't match in length
                yield Read(header[1:], sequence, quality_str)
                # generator: returns and will continue to next in the next call

    return list(read_fastq_chunks()) # returning a list of Read-class items

####################        File Writing         ####################

@handle_f_write("saving kmer collection")
def save_kmer_collection(kmer_collection: KmerCollection,
                         output_file: str) -> None:
    """This function saves the kmer collection in a pickle file.
    It uses gzip compression to open the file in binary mode."""
    if not output_file.endswith('.kdb'):
        print(f"Error: {output_file} is not a KDB file.")
        sys.exit(1) # exit the program
    with gzip.open(output_file, "wb") as f:
        pickle.dump(kmer_collection, f) # type: ignore


####################            Other              ####################

def is_pos_int(var) -> bool:
    """This function checks if a variable is a positive integer."""
    return isinstance(var, int) and var > 0

def get_kmer_size_from_collection(kmer_collection: KmerCollection) -> int:
    """This function extracts the kmer size from the collection."""
    if len(kmer_collection.get_all_kmers()) > 0:
        return len(next(iter(kmer_collection.get_all_kmers())).get_sequence())
    else:
        print("Collection does not contain any kmers.")