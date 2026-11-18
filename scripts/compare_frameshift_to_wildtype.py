import csv
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt

csv.field_size_limit(10_000_000)

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
FIGURES = ROOT / "figures"

INPUT = DATA / "simple_model_sequences.csv"
VALIDATION = DATA / "cds_translation_validation.csv"

OUTPUT = DATA / "frameshift_vs_wildtype.csv"
FIGURE = FIGURES / "frameshift_vs_wildtype_survival.png"

STOP_CODONS = {"TAA", "TAG", "TGA"}


def find_first_stop(sequence, offset):
    """
    Read sequence in the requested shifted frame.

    offset = 1 -> +1 frame
    offset = 2 -> +2 frame

    Returns the first stop position in amino-acid/codon units
    measured from the CDS start.

    Returns None if no stop is found in the available sequence.
    """

    shifted = sequence[offset:]

    for i in range(0, len(shifted) - 2, 3):

        codon = shifted[i:i + 3]

        if codon in STOP_CODONS:

            # Number of translated codons before the stop,
            # plus 1 for the stop codon itself
            return (i // 3) + 1

    return None


# --------------------------------------------------
# Load validated transcript IDs
# --------------------------------------------------

validated_ids = set()

with VALIDATION.open(newline="") as f:

    reader = csv.DictReader(f)

    for row in reader:

        if row["match"].strip().lower() == "true":

            validated_ids.add(
                row["canonical_transcript_id"]
            )


print(
    "Validated transcripts:",
    len(validated_ids)
)


# --------------------------------------------------
# Compare +1 and +2 against wild type
# --------------------------------------------------

results = []

with INPUT.open(newline="") as f:

    reader = csv.DictReader(f)

    for row in reader:

        transcript_id = row[
            "canonical_transcript_id"
        ]

        if transcript_id not in validated_ids:
            continue

        sequence = (
            row["simple_model_sequence"]
            .strip()
            .upper()
            .replace("U", "T")
        )

        cds_length = int(
            row["cds_length"]
        )

        reference_protein = (
            row["protein_sequence"]
            .strip()
            .rstrip("*")
        )

        # ------------------------------------------
        # Wild-type stop
        #
        # Reference protein length tells us how many
        # amino acids are produced before the normal
        # stop codon.
        # ------------------------------------------

        wildtype_stop_position = (
            len(reference_protein) + 1
        )

        plus1_stop = find_first_stop(
            sequence,
            offset=1
        )

        plus2_stop = find_first_stop(
            sequence,
            offset=2
        )

        # ------------------------------------------
        # Classification
        # ------------------------------------------

        def classify(stop_position):

            if stop_position is None:
                return "No stop in available transcript"

            if stop_position > wildtype_stop_position:
                return "Longer than wild type"

            if stop_position == wildtype_stop_position:
                return "Same as wild type"

            return "Shorter than wild type"

        plus1_class = classify(
            plus1_stop
        )

        plus2_class = classify(
            plus2_stop
        )

        results.append({
            "gene_id":
                row.get("gene_id", ""),

            "gene_name":
                row.get("gene_name", ""),

            "transcript_id":
                transcript_id,

            "cds_length_nt":
                cds_length,

            "reference_protein_length_aa":
                len(reference_protein),

            "wildtype_stop_position":
                wildtype_stop_position,

            "plus1_first_stop_position":
                plus1_stop,

            "plus1_vs_wildtype":
                plus1_class,

            "plus2_first_stop_position":
                plus2_stop,

            "plus2_vs_wildtype":
                plus2_class
        })


# --------------------------------------------------
# Save detailed results
# --------------------------------------------------

df = pd.DataFrame(results)

df.to_csv(
    OUTPUT,
    index=False
)


print(
    "\nTranscripts analyzed:",
    len(df)
)


# --------------------------------------------------
# Summary
# --------------------------------------------------

categories = [
    "Shorter than wild type",
    "Same as wild type",
    "Longer than wild type",
    "No stop in available transcript"
]

plus1_counts = (
    df["plus1_vs_wildtype"]
    .value_counts()
    .reindex(categories, fill_value=0)
)

plus2_counts = (
    df["plus2_vs_wildtype"]
    .value_counts()
    .reindex(categories, fill_value=0)
)


print("\n+1 vs wild type:")

for category, count in plus1_counts.items():

    print(
        f"{category}: {count}"
    )


print("\n+2 vs wild type:")

for category, count in plus2_counts.items():

    print(
        f"{category}: {count}"
    )


# --------------------------------------------------
# Specifically print survival beyond wild type
# --------------------------------------------------

plus1_longer = (
    plus1_counts[
        "Longer than wild type"
    ]
)

plus2_longer = (
    plus2_counts[
        "Longer than wild type"
    ]
)


print(
    "\n+1 survived beyond wild type:",
    plus1_longer
)

print(
    "+2 survived beyond wild type:",
    plus2_longer
)


# --------------------------------------------------
# Plot
# --------------------------------------------------

FIGURES.mkdir(
    parents=True,
    exist_ok=True
)

plot_df = pd.DataFrame({
    "Category": categories,
    "+1": plus1_counts.values,
    "+2": plus2_counts.values
})

ax = plot_df.set_index(
    "Category"
).plot(
    kind="bar"
)

ax.set_title(
    "+1 and +2 Survival Compared with Wild Type"
)

ax.set_xlabel("")

ax.set_ylabel(
    "Number of transcripts"
)

plt.xticks(
    rotation=20,
    ha="right"
)

plt.tight_layout()

plt.savefig(
    FIGURE,
    dpi=300,
    bbox_inches="tight"
)

plt.close()


print(
    "\nSaved:",
    OUTPUT
)

print(
    "Saved:",
    FIGURE
)


# --------------------------------------------------
# Show examples that survive beyond WT
# --------------------------------------------------

print(
    "\nExamples where +1 survives beyond wild type:"
)

examples_plus1 = df[
    df["plus1_vs_wildtype"]
    == "Longer than wild type"
].sort_values(
    "plus1_first_stop_position",
    ascending=False
).head(10)

print(
    examples_plus1[
        [
            "gene_name",
            "transcript_id",
            "wildtype_stop_position",
            "plus1_first_stop_position"
        ]
    ].to_string(index=False)
)


print(
    "\nExamples where +2 survives beyond wild type:"
)

examples_plus2 = df[
    df["plus2_vs_wildtype"]
    == "Longer than wild type"
].sort_values(
    "plus2_first_stop_position",
    ascending=False
).head(10)

print(
    examples_plus2[
        [
            "gene_name",
            "transcript_id",
            "wildtype_stop_position",
            "plus2_first_stop_position"
        ]
    ].to_string(index=False)
)