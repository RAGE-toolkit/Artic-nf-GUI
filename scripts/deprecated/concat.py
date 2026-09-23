from Bio import SeqIO
import argparse

def main():
    parser = argparse.ArgumentParser(
        description="Concatenate multiple FASTA files into one non-redundant file."
    )
    parser.add_argument(
        "-i", "--input",
        nargs="+",
        required=True,
        help="Input FASTA files (space separated)"
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output FASTA file"
    )
    args = parser.parse_args()

    unique_records = {}

    for fasta in args.input:
        for record in SeqIO.parse(fasta, "fasta"):
            seq_str = str(record.seq)
            if seq_str not in unique_records:
                unique_records[seq_str] = record

    with open(args.output, "w") as out_f:
        SeqIO.write(unique_records.values(), out_f, "fasta")

    print(f"✅ Non-redundant FASTA written to {args.output}")

if __name__ == "__main__":
    main()

