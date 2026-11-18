import csv
from pathlib import Path

import pandas as pd

csv.field_size_limit(10_000_000)

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"

SEQUENCE_FILE = DATA / "simple_model_sequences.csv"
COMPARISON_FILE = DATA / "frameshift_vs_wildtype.csv"

OUTPUT_CSV = DATA / "frameshift_biology_candidates.csv"
PLUS1_FASTA = DATA / "plus1_frameshift_candidates.fasta"
PLUS2_FASTA = DATA / "plus2_frameshift_candidates.fasta"

STOP_CODONS = {"TAA", "TAG", "TGA"}

CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TGT": "C", "TGC": "C",
    "TGG": "W",

    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",

    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",

    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G"
}


def translate_until_stop(sequence, offset):
    """
    Translate +1 or +2 frame until the first stop codon.

    offset = 1 for +1
    offset = 2 for +2
    """

    shifted = sequence[offset:]

    protein = []
    stop_position = None
    stop_codon = None

    for i in range(0, len(shifted) - 2, 3):

        codon = shifted[i:i + 3]

        if codon in STOP_CODONS:
            stop_position = (i // 3) + 1
            stop_codon = codon
            break

        amino_acid = CODON_TABLE.get(
            codon,
            "X"
        )

        protein.append(amino_acid)

    return (
        "".join(protein),
        stop_position,
        stop_codon
    )


# --------------------------------------------------
# Load candidate classifications
# --------------------------------------------------

comparison = pd.read_csv(
    COMPARISON_FILE
)

interesting_categories = {
    "Longer than wild type",
    "No stop in available transcript"
}


# --------------------------------------------------
# Load transcript sequences
# --------------------------------------------------

sequences = {}

with SEQUENCE_FILE.open(newline="") as f:

    reader = csv.DictReader(f)

    for row in reader:

        transcript_id = row[
            "canonical_transcript_id"
        ]

        sequences[transcript_id] = row


# --------------------------------------------------
# Generate candidate proteins
# --------------------------------------------------

results = []

for _, row in comparison.iterrows():

    transcript_id = row["transcript_id"]

    if transcript_id not in sequences:
        continue

    sequence_row = sequences[
        transcript_id
    ]

    sequence = (
        sequence_row[
            "simple_model_sequence"
        ]
        .strip()
        .upper()
        .replace("U", "T")
    )

    reference_protein = (
        sequence_row[
            "protein_sequence"
        ]
        .strip()
        .rstrip("*")
    )

    for frame, offset, category_column in [
        (
            "+1",
            1,
            "plus1_vs_wildtype"
        ),
        (
            "+2",
            2,
            "plus2_vs_wildtype"
        )
    ]:

        category = row[
            category_column
        ]

        if category not in interesting_categories:
            continue

        (
            protein,
            stop_position,
            stop_codon
        ) = translate_until_stop(
            sequence,
            offset
        )

        results.append({
            "gene_id":
                row["gene_id"],

            "gene_name":
                row["gene_name"],

            "transcript_id":
                transcript_id,

            "frame":
                frame,

            "classification":
                category,

            "wildtype_protein_length":
                len(reference_protein),

            "frameshift_protein_length":
                len(protein),

            "first_stop_position":
                stop_position,

            "first_stop_codon":
                stop_codon,

            "length_difference_vs_wildtype":
                (
                    len(protein)
                    - len(reference_protein)
                ),

            "frameshift_protein_sequence":
                protein
        })


# --------------------------------------------------
# Save candidate table
# --------------------------------------------------

candidate_df = pd.DataFrame(
    results
)

candidate_df.to_csv(
    OUTPUT_CSV,
    index=False
)


# --------------------------------------------------
# Write FASTA files
# --------------------------------------------------

with PLUS1_FASTA.open("w") as plus1_file, \
     PLUS2_FASTA.open("w") as plus2_file:

    for _, row in candidate_df.iterrows():

        header = (
            f">{row['gene_name']}|"
            f"{row['transcript_id']}|"
            f"{row['frame']}|"
            f"{row['classification'].replace(' ', '_')}|"
            f"WT={row['wildtype_protein_length']}|"
            f"FS={row['frameshift_protein_length']}"
        )

        sequence = row[
            "frameshift_protein_sequence"
        ]

        # Split FASTA sequence into readable lines
        wrapped = "\n".join(
            sequence[i:i + 60]
            for i in range(
                0,
                len(sequence),
                60
            )
        )

        if row["frame"] == "+1":

            plus1_file.write(
                header
                + "\n"
                + wrapped
                + "\n"
            )

        else:

            plus2_file.write(
                header
                + "\n"
                + wrapped
                + "\n"
            )


# --------------------------------------------------
# Summary
# --------------------------------------------------

print(
    "Total protein candidates:",
    len(candidate_df)
)

print("\n+1:")

print(
    candidate_df[
        candidate_df["frame"] == "+1"
    ]["classification"].value_counts()
)

print("\n+2:")

print(
    candidate_df[
        candidate_df["frame"] == "+2"
    ]["classification"].value_counts()
)


print(
    "\nLongest +1 products:"
)

print(
    candidate_df[
        candidate_df["frame"] == "+1"
    ]
    .sort_values(
        "frameshift_protein_length",
        ascending=False
    )
    [
        [
            "gene_name",
            "transcript_id",
            "wildtype_protein_length",
            "frameshift_protein_length",
            "classification"
        ]
    ]
    .head(10)
    .to_string(index=False)
)


print(
    "\nLongest +2 products:"
)

print(
    candidate_df[
        candidate_df["frame"] == "+2"
    ]
    .sort_values(
        "frameshift_protein_length",
        ascending=False
    )
    [
        [
            "gene_name",
            "transcript_id",
            "wildtype_protein_length",
            "frameshift_protein_length",
            "classification"
        ]
    ]
    .head(10)
    .to_string(index=False)
)


print("\nSaved:")
print(OUTPUT_CSV)
print(PLUS1_FASTA)
print(PLUS2_FASTA)