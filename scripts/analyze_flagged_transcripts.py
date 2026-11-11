import csv
from collections import Counter

import matplotlib.pyplot as plt

csv.field_size_limit(10_000_000)

validation_file = "data/cds_translation_validation.csv"
genes_file = "data/genes_with_cds.csv"

output_file = "data/flagged_transcript_analysis.csv"
figure_file = "figures/flagged_transcript_reasons.png"


# --------------------------------------------------
# Load all gene/CDS information
# --------------------------------------------------

genes = {}

with open(genes_file, "r") as f:
    reader = csv.DictReader(f)

    for row in reader:
        genes[row["canonical_transcript_id"]] = row


# --------------------------------------------------
# Find transcripts that failed validation
# --------------------------------------------------

flagged = []

with open(validation_file, "r") as f:
    reader = csv.DictReader(f)

    for row in reader:

        if row["match"].lower() == "true":
            continue

        transcript_id = row["canonical_transcript_id"]

        if transcript_id not in genes:
            continue

        gene = genes[transcript_id]

        cds = gene["cds_sequence"].upper()
        reference_protein = gene["protein_sequence"]

        translated_length = int(
            row["translated_protein_length"]
        )

        reference_length = len(reference_protein)

        start_codon = cds[:3] if len(cds) >= 3 else ""

        reasons = []

        # ------------------------------------------
        # Classify possible reasons
        # ------------------------------------------

        if translated_length < reference_length:
            reasons.append("translated_protein_shorter")

        if translated_length == reference_length:
            reasons.append("same_length_sequence_differs")

        if start_codon != "ATG":
            reasons.append("non_ATG_start")

        if len(cds) % 3 != 0:
            reasons.append("CDS_not_divisible_by_3")

        if not reasons:
            reasons.append("other")

        flagged.append({
            "gene_id": gene["gene_id"],
            "gene_name": gene["gene_name"],
            "transcript_id": transcript_id,

            "start_codon": start_codon,

            "cds_length": len(cds),

            "reference_protein_length":
                reference_length,

            "translated_protein_length":
                translated_length,

            "length_difference":
                reference_length - translated_length,

            "cds_divisible_by_3":
                len(cds) % 3 == 0,

            "reason_flagged":
                "; ".join(reasons)
        })


# --------------------------------------------------
# Save detailed table
# --------------------------------------------------

with open(output_file, "w", newline="") as f:

    fieldnames = flagged[0].keys()

    writer = csv.DictWriter(
        f,
        fieldnames=fieldnames
    )

    writer.writeheader()
    writer.writerows(flagged)


# --------------------------------------------------
# Count reasons
# --------------------------------------------------

reason_counts = Counter()

for row in flagged:

    for reason in row["reason_flagged"].split("; "):
        reason_counts[reason] += 1


print(f"Total flagged transcripts: {len(flagged)}")

print("\nReasons:")

for reason, count in reason_counts.most_common():
    print(f"{reason}: {count}")


# --------------------------------------------------
# Print first flagged genes
# --------------------------------------------------

print("\nFirst 20 flagged genes:")

for row in flagged[:20]:
    print(
        row["gene_name"],
        row["transcript_id"],
        row["reason_flagged"]
    )


# --------------------------------------------------
# Make visualization
# --------------------------------------------------

labels = list(reason_counts.keys())
values = list(reason_counts.values())

plt.figure(figsize=(9, 5))

bars = plt.bar(labels, values)

plt.ylabel("Number of Flagged Transcripts")
plt.title("Reasons for CDS-to-Protein Validation Mismatches")

plt.xticks(
    rotation=25,
    ha="right"
)

for bar, value in zip(bars, values):
    plt.text(
        bar.get_x() + bar.get_width() / 2,
        value + 3,
        str(value),
        ha="center"
    )

plt.tight_layout()

plt.savefig(
    figure_file,
    dpi=300
)

plt.close()


print(f"\nSaved detailed results to {output_file}")
print(f"Saved figure to {figure_file}")
