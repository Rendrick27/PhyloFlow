import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def concatenate_fastas(folder_path, output_file):
    """
    Concatenates sequences from multiple FASTA files within a specified folder and writes a new FASTA file.

    Args:
        folder_path (str): The path to the folder containing the FASTA files.
        output_file (str): The path to the output file to save the concatenated sequences.

    This function reads all FASTA files in the given directory, concatenates sequences based on their IDs,
    and writes the concatenated sequences to an output file. If a sequence ID is not found in a file, it
    pads the sequence with 'n' corresponding to the maximum length of sequences in that file.
    """
    all_sequences = {}
    sequence_order = []

    # List and sort FASTA files in the folder
    fasta_files = [f for f in os.listdir(folder_path) if f.endswith((".fasta", ".fa"))]
    fasta_files.sort()

    # Read sequences from each file
    for file_name in fasta_files:
        file_path = os.path.join(folder_path, file_name)
        max_length = 0
        with open(file_path, 'r') as f:
            for record in SeqIO.parse(f, "fasta"):
                seq_id = record.id
                sequence = str(record.seq)
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
            sequence = all_sequences[seq_id].get(file_name, 'n' * all_sequences['_max_lengths_' + file_name])
            concatenated_sequence += sequence
        concatenated_sequences[seq_id] = concatenated_sequence

    # Write concatenated sequences to the output file
    with open(output_file, 'w') as output:
        for seq_id in sequence_order:
            record = SeqRecord(Seq(concatenated_sequences[seq_id]), id=seq_id, description="")
            SeqIO.write(record, output, "fasta")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python concatenate_fastas.py <folder_path> <output_file>")
        sys.exit(1)

    folder_path = sys.argv[1]
    output_file = sys.argv[2]

    concatenate_fastas(folder_path, output_file)
