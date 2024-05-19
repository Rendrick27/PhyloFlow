import sys


def clean_fasta(input_file, output_file):
    """
    Processes a FASTA file to remove all occurrences of the letter 'n'
        (uppercase) from the DNA sequences.

    Args:
    input_file (str): The path to the input FASTA file.
    output_file (str): The path to the output FASTA file
        where the cleaned sequences will be saved.

    Each sequence in the file is expected to start
        with a header line that begins with '>'.
    This script retains the header lines and writes the cleaned sequences
        to the output file.
    """
    with open(input_file, 'r') as f, open(output_file, 'w') as fout:
        for line in f:
            if line.startswith('>'):
                fout.write(line)
            else:
                cleaned_sequence = line.strip().replace('N', '')
                fout.write(cleaned_sequence + '\n')


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 script_name.py <input_fasta> <output_fasta>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    clean_fasta(input_file, output_file)
