import os
import sys
from Bio import SeqIO

def count_first_sequence_length(folder_path):
    """
    Counts and prints the length of the first sequence in each FASTA file within a given directory.

    Args:
        folder_path (str): The path to the directory containing FASTA files.

    If no FASTA files or sequences are found, the function prints an appropriate message. It also
    handles cases where the directory path is incorrect or empty.
    """
    # Validate the directory
    if not os.path.isdir(folder_path):
        print("Error: The provided path is not a directory.")
        return

    # Find all FASTA files in the directory
    fasta_files = [file for file in os.listdir(folder_path) if file.endswith(".fasta")]

    if not fasta_files:
        print("No FASTA files found in the directory.")
        return

    # Process each FASTA file
    for file in fasta_files:
        file_path = os.path.join(folder_path, file)
        with open(file_path, "r") as fasta_file:
            records = list(SeqIO.parse(fasta_file, "fasta"))

            if not records:
                print(f"No sequences found in {file}")
            else:
                first_sequence_length = len(records[0].seq)
                print(f"Length of the first sequence in {file}: {first_sequence_length}")

if __name__ == "__main__":
    # Check command-line arguments
    if len(sys.argv) != 2:
        print("Usage: python script.py /path/to/folder")
        sys.exit(1)

    folder_path = sys.argv[1]
    count_first_sequence_length(folder_path)
