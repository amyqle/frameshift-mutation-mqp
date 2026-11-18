import pandas as pd
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"

INPUT = DATA / "top_frameshift_blast_candidates_named.csv"

OUTPUT_FASTA = DATA / "top_frameshift_blast_candidates_clean.fasta"


df = pd.read_csv(INPUT)


with OUTPUT_FASTA.open("w") as f:

    for _, row in df.iterrows():

        gene_name = row["gene_name"]
        transcript_id = row["transcript_id"]

        if pd.isna(gene_name) or str(gene_name).strip() == "":
            display_name = transcript_id
        else:
            display_name = str(gene_name)

        frame_label = (
            "+1"
            if int(row["frame"]) == 1
            else "+2"
        )

        header = (
            f">{display_name}|"
            f"{transcript_id}|"
            f"{frame_label}|"
            f"WT={int(row['wildtype_protein_length'])}|"
            f"FS={int(row['frameshift_protein_length'])}|"
            f"EXT={int(row['extension_aa'])}"
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


print("Saved clean FASTA:")
print(OUTPUT_FASTA)

print("\nHeaders:")

for _, row in df.iterrows():

    gene_name = row["gene_name"]
    transcript_id = row["transcript_id"]

    if pd.isna(gene_name) or str(gene_name).strip() == "":
        display_name = transcript_id
    else:
        display_name = str(gene_name)

    frame_label = (
        "+1"
        if int(row["frame"]) == 1
        else "+2"
    )

    print(
        f"{display_name} | "
        f"{frame_label} | "
        f"WT={int(row['wildtype_protein_length'])} | "
        f"FS={int(row['frameshift_protein_length'])} | "
        f"EXT={int(row['extension_aa'])}"
    )