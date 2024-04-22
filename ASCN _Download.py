#!/usr/bin/env python3

from Bio import Entrez
import sys

def fetch_sequences(database, accession_numbers):
    Entrez.email = "fakeemail@exemplo.com"
    Entrez.tool = "fetch_sequences.py"
    
    # Join the accession numbers into a comma-separated string
    accession_list = ",".join(accession_numbers)
    
    fetch_results = Entrez.efetch(db=database, id=accession_list, rettype="fasta", retmode="text")
    sys.stdout.write(fetch_results.read())
    fetch_results.close()

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python fetch_sequences.py [database] [file_path_to_accession_numbers]")
        sys.exit(1)

    database = sys.argv[1]
    file_path = sys.argv[2]

    try:
        with open(file_path) as file:
            accession_numbers = [line.strip() for line in file if line.strip()]
    except FileNotFoundError:
        print(f"Error: File '{file_path}' not found.")
        sys.exit(1)

    if not accession_numbers:
        print("Error: No accession numbers found in the file.")
        sys.exit(1)

    fetch_sequences(database, accession_numbers)

