"""
fetch_sequences.py

Fetches sequences from NCBI GenBank by accession number and saves them
as a single FASTA file, ready for build_tree.py.

Usage:
    python fetch_sequences.py

Edit the ACCESSIONS list and EMAIL below before running.
NCBI requires a real email address for API access (it's just used to
contact you if your usage causes problems — this is standard practice
for the Entrez API, not optional).
"""

from Bio import Entrez, SeqIO

# --- Configuration ---
Entrez.email = "email@example.com"  
OUTPUT_FILE = "data/sample_sequences.fasta"


ACCESSIONS = [
    "NC_012920.1",
    "NC_001643.1",
    "NC_005089.1",
      
]


def fetch_sequences(accessions, output_file):
    records = []
    for acc in accessions:
        print(f"Fetching {acc}...")
        handle = Entrez.efetch(db="nucleotide", id=acc, rettype="fasta", retmode="text")
        record = SeqIO.read(handle, "fasta")
        handle.close()
        records.append(record)

    SeqIO.write(records, output_file, "fasta")
    print(f"\nSaved {len(records)} sequences to {output_file}")


if __name__ == "__main__":
    if Entrez.email == "email@example.com":
        print("Please set your real email address in Entrez.email before running.")
    else:
        fetch_sequences(ACCESSIONS, OUTPUT_FILE)
