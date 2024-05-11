import sys
import os
import glob

def get_sequence_length(fasta_file):
    """
    Calculate the length of the sequence in a FASTA file.

    Args:
        fasta_file (str): Path to the FASTA file.

    Returns:
        int: The total length of the sequence in the file, excluding header lines and whitespace.
    """
    seq_length = 0
    with open(fasta_file, 'r') as file:
        for line in file:
            if line.startswith('>'):
                if seq_length > 0:
                    break  # Stop if we reach a new sequence or repeat header
                continue
            seq_length += len(line.strip())
    return seq_length

def process_sequences(directory, output_file):
    """
    Processes each FASTA file in the directory to calculate sequence lengths and generate output
    based on corresponding model data.

    Args:
        directory (str): Directory containing FASTA files and corresponding model files.
        output_file (str): File to which the formatted output will be written.
    """
    fasta_files = sorted(glob.glob(os.path.join(directory, '*.fasta')))
    outputs = []
    start = 1

    for fasta_file in fasta_files:
        base_name = os.path.splitext(os.path.basename(fasta_file))[0]
        model_file = os.path.join(directory, f"{base_name}Model.txt")

        if os.path.exists(model_file):
            seq_length = get_sequence_length(fasta_file)

            with open(model_file, 'r') as file:
                model_data = file.read().strip()

            end = start + seq_length - 1
            formatted_output = f"{model_data}, {base_name}={start}-{end}"
            outputs.append(formatted_output)
            start = end + 1  # Update start for the next sequence
        else:
            print(f"Warning: Model file not found for {base_name}")

    with open(output_file, 'w') as file:
        for output in outputs:
            file.write(output + "\n")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python script.py <directory> <output_file>")
    else:
        process_sequences(sys.argv[1], sys.argv[2])
