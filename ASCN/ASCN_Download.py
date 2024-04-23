from Bio import Entrez
import sys


def fetch_sequences(database, accession_numbers):
    """
    Fetch sequences from NCBI database using accession numbers.
    
    Args:
        database (str): Database name to fetch sequences from.
        accession_numbers (list): List of accession numbers.
        
    Returns:
        None
    """
    Entrez.email = "fakeemail@exemplo.com"
    Entrez.tool = "fetch_sequences.py"

    # Join the accession numbers into a comma-separated string
    accession_list = ",".join(accession_numbers)

    try:
        # Fetch the sequences
        fetch_results = Entrez.efetch(db=database, id=accession_list, rettype="fasta", retmode="text")
        sys.stdout.write(fetch_results.read())
        fetch_results.close()
    except Exception as e:
        print(f"Error fetching sequences: {e}")
        sys.exit(1)


def read_accession_numbers(file_path):
    """
    Read accession numbers from a file.
    
    Args:
        file_path (str): Path to the file containing accession numbers.
        
    Returns:
        list: List of accession numbers.
    """
    try:
        with open(file_path) as file:
            # Filter out empty lines and strip whitespace from each line
            accession_numbers = [line.strip() for line in file if line.strip()]
        return accession_numbers
    except FileNotFoundError:
        print(f"Error: File '{file_path}' not found.")
        sys.exit(1)


def main():
    """
    Main function to execute the script.
    """
    if len(sys.argv) != 3:
        print("Usage: python fetch_sequences.py [database] "
        "[file_path_to_accession_numbers]")
        sys.exit(1)

    database = sys.argv[1]
    file_path = sys.argv[2]

    accession_numbers = read_accession_numbers(file_path)

    if not accession_numbers:
        print("Error: No accession numbers found in the file.")
        sys.exit(1)

    fetch_sequences(database, accession_numbers)


if __name__ == '__main__':
    main()
