import pandas as pd
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"

INPUT = DATA / "frameshift_biology_candidates.csv"

OUTPUT_CSV = DATA / "interproscan_candidates_over100aa.csv"
OUTPUT_FASTA = DATA / "interproscan_candidates_over100aa.fasta"


# --------------------------------------------------
# Load biology candidates
# --------------------------------------------------

df = pd.read_csv(INPUT)


# --------------------------------------------------
# Filter for proteins > 100 amino acids
# --------------------------------------------------

filtered = df[
    df["frameshift_protein_length"] > 100
].copy()


# --------------------------------------------------
# Clean frame labels
# --------------------------------------------------

filtered["frame_label"] = filtered["frame"].apply(
    lambda x: "+1" if int(x) == 1 else "+2"
)


# --------------------------------------------------
# Save CSV
# --------------------------------------------------

filtered.to_csv(
    OUTPUT_CSV,
    index=False
)


# --------------------------------------------------
# Save FASTA
# --------------------------------------------------

with OUTPUT_FASTA.open("w") as f:

    for _, row in filtered.iterrows():

        gene_name = row["gene_name"]

        if pd.isna(gene_name) or str(gene_name).strip() == "":
            display_name = row["transcript_id"]
        else:
            display_name = str(gene_name)

        header = (
            f">{display_name}|"
            f"{row['transcript_id']}|"
            f"{row['frame_label']}|"
            f"WT={int(row['wildtype_protein_length'])}|"
            f"FS={int(row['frameshift_protein_length'])}"
        )

        sequence = str(
            row["frameshift_protein_sequence"]
        )

        wrapped = "\n".join(
            sequence[i:i + 60]
            for i in range(
                0,
                len(sequence),
                60
            )
        )

        f.write(
            header
            + "\n"
            + wrapped
            + "\n"
        )


# --------------------------------------------------
# Summary
# --------------------------------------------------

print(
    "Total biology candidates:",
    len(df)
)

print(
    "Candidates >100 aa:",
    len(filtered)
)

print("\nBy frame:")

print(
    filtered[
        "frame_label"
    ].value_counts()
)

print(
    "\nProtein length summary:"
)

print(
    filtered[
        "frameshift_protein_length"
    ].describe()
)

print("\nSaved:")
print(OUTPUT_CSV)
print(OUTPUT_FASTA)