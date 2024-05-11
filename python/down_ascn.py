from Bio import Entrez, SeqIO
import os
import sys

def download_sequences_from_asn_file(input_file, output_file):
    """
    Downloads DNA sequences from the NCBI database using ASN numbers provided in an input file and writes them to an output file.

    Args:
        input_file (str): Path to a file containing species names and ASN numbers, separated by a semicolon.
        output_file (str): Path to the output FASTA file where downloaded sequences will be saved.

    Each line in the input file should contain a species name followed by its ASN number, separated by a semicolon.
    The function modifies the sequence ID and description to use the species name and writes each sequence to the output file.
    """
    Entrez.email = "fakemail@gmail.com"  # Set your email address for Entrez

    with open(input_file, "r") as f, open(output_file, "w") as out_f:
        for line in f:
            line = line.strip()
            if ";" not in line:
                continue

            species_name, asn_number = line.split(";")

            try:
                handle = Entrez.efetch(db="nucleotide", id=asn_number, rettype="fasta", retmode="text")
                record = SeqIO.read(handle, "fasta")
                handle.close()

                record.id = record.description = species_name
                SeqIO.write(record, out_f, "fasta")

                print(f"Downloaded sequence for {species_name}")
            except Exception as e:
                print(f"Error downloading sequence for {species_name}: {e}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_file> <output_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    download_sequences_from_asn_file(input_file, output_file)
