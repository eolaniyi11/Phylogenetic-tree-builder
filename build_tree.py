"""
build_tree.py

Builds a phylogenetic tree from a set of DNA sequences using a
distance-based method:

  1. Load sequences from a FASTA file
  2. Compute pairwise genetic distances (via global alignment identity)
  3. Build a distance matrix
  4. Construct a tree using Neighbor-Joining
  5. Visualize the tree and save it as an image + Newick file

Usage:
    python build_tree.py
"""

from Bio import SeqIO, Align
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
from Bio import Phylo
import matplotlib.pyplot as plt
import itertools
import os

INPUT_FILE = "data/sample_sequences.fasta"
OUTPUT_DIR = "output"


def load_sequences(fasta_path):
    records = list(SeqIO.parse(fasta_path, "fasta"))
    if len(records) < 3:
        raise ValueError("Need at least 3 sequences to build a meaningful tree.")
    return records


def compute_distance_matrix(records):
    """
    Aligns every pair of sequences and converts alignment identity into
    a genetic distance (0 = identical, 1 = completely different).
    """
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 1
    aligner.mismatch_score = 0
    aligner.open_gap_score = -1
    aligner.extend_gap_score = -0.5

    names = [rec.id for rec in records]
    seqs = [str(rec.seq) for rec in records]
    n = len(seqs)

    # DistanceMatrix expects a lower-triangular matrix (list of lists)
    matrix = [[0.0] * (i + 1) for i in range(n)]

    print("Computing pairwise distances...")
    for i, j in itertools.combinations(range(n), 2):
        alignment = aligner.align(seqs[i], seqs[j])[0]
        score = alignment.score
        max_possible = max(len(seqs[i]), len(seqs[j]))
        similarity = score / max_possible
        distance = 1 - similarity
        matrix[j][i] = round(distance, 5)
        print(f"  {names[i]} vs {names[j]}: distance = {distance:.4f}")

    return DistanceMatrix(names, matrix)


def build_tree(distance_matrix):
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(distance_matrix)  # Neighbor-Joining
    tree.root_at_midpoint()
    return tree


def save_outputs(tree, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    # Save tree in Newick format (standard for phylogenetics — can be
    # opened in tools like iTOL, FigTree, or shared with other software)
    newick_path = os.path.join(output_dir, "tree.newick")
    Phylo.write(tree, newick_path, "newick")
    print(f"\nSaved Newick tree to {newick_path}")

    # Draw and save a visual figure
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(1, 1, 1)
    Phylo.draw(tree, axes=ax, do_show=False)
    ax.set_title("Phylogenetic Tree (Neighbor-Joining)")
    plt.tight_layout()
    image_path = os.path.join(output_dir, "tree.png")
    plt.savefig(image_path, dpi=150)
    print(f"Saved tree image to {image_path}")


def main():
    records = load_sequences(INPUT_FILE)
    print(f"Loaded {len(records)} sequences: {[r.id for r in records]}\n")

    dm = compute_distance_matrix(records)
    tree = build_tree(dm)

    print("\nTree structure (ASCII):")
    Phylo.draw_ascii(tree)

    save_outputs(tree, OUTPUT_DIR)


if __name__ == "__main__":
    main()
