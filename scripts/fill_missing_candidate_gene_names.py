import csv
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"

INPUT = DATA / "top_frameshift_blast_candidates.csv"
GTF_FILE = DATA / "genes.gtf"

OUTPUT = DATA / "top_frameshift_blast_candidates_named.csv"


def parse_attributes(attribute_text):
    attributes = {}

    for item in attribute_text.strip().split(";"):
        item = item.strip()

        if not item:
            continue

        if " " in item:
            key, value = item.split(" ", 1)

            attributes[key] = (
                value.strip()
                .strip('"')
            )

    return attributes


def normalize_id(value):
    return str(value).split(".")[0]


# --------------------------------------------------
# Load candidate transcripts
# --------------------------------------------------

df = pd.read_csv(INPUT)

candidate_ids = {
    normalize_id(x)
    for x in df["transcript_id"]
}

print(
    "Candidate transcripts:",
    len(candidate_ids)
)


# --------------------------------------------------
# Read transcript -> gene mappings from GTF
# --------------------------------------------------

transcript_map = {}

with open(
    GTF_FILE,
    "r",
    encoding="utf-8"
) as handle:

    for line in handle:

        if line.startswith("#"):
            continue

        fields = line.rstrip("\n").split("\t")

        if len(fields) != 9:
            continue

        attributes = parse_attributes(
            fields[8]
        )

        transcript_id = attributes.get(
            "transcript_id"
        )

        if not transcript_id:
            continue

        normalized = normalize_id(
            transcript_id
        )

        if normalized not in candidate_ids:
            continue

        transcript_map[normalized] = {
            "gene_name":
                attributes.get(
                    "gene_name",
                    ""
                ),

            "gene_id":
                attributes.get(
                    "gene_id",
                    ""
                )
        }


print(
    "Transcript mappings found:",
    len(transcript_map)
)


# --------------------------------------------------
# Fill missing gene names
# --------------------------------------------------

for index, row in df.iterrows():

    transcript_id = normalize_id(
        row["transcript_id"]
    )

    current_name = row["gene_name"]

    if (
        pd.isna(current_name)
        or str(current_name).strip() == ""
    ):

        info = transcript_map.get(
            transcript_id
        )

        if info:

            df.at[
                index,
                "gene_name"
            ] = info["gene_name"]

            if (
                "gene_id" in df.columns
                and (
                    pd.isna(
                        row.get(
                            "gene_id"
                        )
                    )
                    or str(
                        row.get(
                            "gene_id",
                            ""
                        )
                    ).strip() == ""
                )
            ):
                df.at[
                    index,
                    "gene_id"
                ] = info["gene_id"]


# --------------------------------------------------
# Save updated file
# --------------------------------------------------

df.to_csv(
    OUTPUT,
    index=False
)


# --------------------------------------------------
# Print results
# --------------------------------------------------

print("\nUpdated candidates:")

print(
    df[
        [
            "gene_name",
            "transcript_id",
            "frame",
            "wildtype_protein_length",
            "frameshift_protein_length",
            "extension_aa"
        ]
    ].to_string(index=False)
)


# --------------------------------------------------
# Show any names still missing
# --------------------------------------------------

missing = df[
    df["gene_name"].isna()
    | (
        df["gene_name"]
        .astype(str)
        .str.strip()
        == ""
    )
]

print(
    "\nCandidates still missing gene names:",
    len(missing)
)

if len(missing) > 0:

    print(
        missing[
            [
                "transcript_id",
                "frame"
            ]
        ].to_string(index=False)
    )


print("\nSaved:")
print(OUTPUT)