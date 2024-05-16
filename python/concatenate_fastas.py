import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def preprocess_sequence(sequence):
    """
    Replaces all '-' characters with 'n' at the edges of the sequence until
        a nucleotide is encountered.

    Args:
        sequence (str): The original sequence from the FASTA file.

    Returns:
        str: The preprocessed sequence with '-' at the edges replaced by 'n'.
    """
    import re

    # Regex to find leading and trailing '-' characters
    leading_pattern = re.compile(r'^-+')
    trailing_pattern = re.compile(r'-+$')

    # Replace leading '-' with 'n'
    sequence = leading_pattern.sub(lambda x: 'n' * len(x.group()), sequence)
    sequence = trailing_pattern.sub(lambda x: 'n' * len(x.group()), sequence)

    return sequence


def concatenate_fastas(folder_path, output_file):
    """
    Concatenates sequences from multiple FASTA files within a specified folder
        and writes a new FASTA file, after replacing all '-' characters with
        'n' on the edges of the sequences.

    Args:
        folder_path (str): The path to the folder containing the FASTA files.
        output_file (str): The path to the output file to save the
            concatenated sequences.
    """
    all_sequences = {}
    sequence_order = []

    # List and sort FASTA files in the folder
    fasta_files = [f for f in os.listdir(folder_path)
                   if f.endswith((".fasta", ".fa"))]
    fasta_files.sort()

    # Read sequences from each file
    for file_name in fasta_files:
        file_path = os.path.join(folder_path, file_name)
        max_length = 0
        with open(file_path, 'r') as f:
            for record in SeqIO.parse(f, "fasta"):
                seq_id = record.id
                sequence = preprocess_sequence(str(record.seq))
                if seq_id not in all_sequences:
                    all_sequences[seq_id] = {}
                    sequence_order.append(seq_id)
                all_sequences[seq_id][file_name] = sequence
                max_length = max(max_length, len(sequence))
            all_sequences['_max_lengths_' + file_name] = max_length

    # Construct concatenated sequences
    concatenated_sequences = {}
    for seq_id in sequence_order:
        concatenated_sequence = ''
        for file_name in fasta_files:
            sequence = all_sequences[seq_id].get(
                file_name, 'n' * all_sequences['_max_lengths_' + file_name])
            concatenated_sequence += sequence
        concatenated_sequences[seq_id] = concatenated_sequence

    # Write concatenated sequences to the output file
    with open(output_file, 'w') as output:
        for seq_id in sequence_order:
            record = SeqRecord(Seq(concatenated_sequences[seq_id]),
                               id=seq_id, description="")
            SeqIO.write(record, output, "fasta")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <folder_path> <output_file>")
        sys.exit(1)

    folder_path = sys.argv[1]
    output_file = sys.argv[2]

    concatenate_fastas(folder_path, output_file)
