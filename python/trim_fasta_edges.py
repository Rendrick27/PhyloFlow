import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def trim_edges(sequences, threshold):
    """
    Trims the edges of aligned sequences based on the threshold
        of nucleotide presence.

    Args:
        sequences (list of SeqRecord): The list of sequence records
            to be trimmed.
        threshold (float): The percentage of sequences that must have
            a nucleotide (not a gap '-') at a position to include that
            column in the output.

    Returns:
        list of SeqRecord: The trimmed sequence records.
    """
    num_seqs = len(sequences)
    seq_length = len(sequences[0].seq)
    min_presence = (threshold / 100.0) * num_seqs

    # Find the first and last indices that meet the threshold
    start_index = next((idx for idx in range(seq_length)
                        if sum(1 for seq in sequences if seq.seq[idx] != '-')
                        >= min_presence), seq_length)
    end_index = next((idx for idx in reversed(range(seq_length))
                      if sum(1 for seq in sequences if seq.seq[idx] != '-')
                      >= min_presence), -1)

    # Trim sequences and create new SeqRecords
    trimmed_sequences = [SeqRecord(seq.seq[start_index:end_index+1],
                                   id=seq.id, description="")
                         for seq in sequences if end_index >= start_index]

    return trimmed_sequences


def main():
    if len(sys.argv) != 4:
        print("Usage: python script.py <input_fasta> <output_fasta>"
              "<threshold>")
        sys.exit(1)

    input_fasta = sys.argv[1]
    output_fasta = sys.argv[2]
    threshold = float(sys.argv[3])

    sequences = list(SeqIO.parse(input_fasta, 'fasta'))
    trimmed_sequences = trim_edges(sequences, threshold)
    SeqIO.write(trimmed_sequences, output_fasta, 'fasta')

    print(f"Trimmed sequences have been saved to {output_fasta}")


if __name__ == "__main__":
    main()
