import csv
from collections import Counter

import matplotlib.pyplot as plt

csv.field_size_limit(10_000_000)

genes_file = "data/genes_with_cds.csv"
validation_file = "data/cds_translation_validation.csv"

output_file = "data/frameshift_stop_positions.csv"

figure_plus1 = "figures/plus1_stop_position_distribution.png"
figure_plus2 = "figures/plus2_stop_position_distribution.png"
figure_compare = "figures/frameshift_relative_stop_comparison.png"
figure_codons = "figures/frameshift_stop_codon_distribution.png"


STOP_CODONS = {"TAA", "TAG", "TGA"}


# --------------------------------------------------
# Find first stop codon in a shifted reading frame
# --------------------------------------------------

def find_first_stop(sequence, shift):

    shifted_sequence = sequence[shift:]

    for i in range(0, len(shifted_sequence) - 2, 3):

        codon = shifted_sequence[i:i + 3]

        if codon in STOP_CODONS:

            # Position within shifted reading frame
            codon_number = (i // 3) + 1

            # Convert back to position in original CDS
            stop_base_position = shift + i + 1

            return {
                "stop_codon": codon,
                "stop_codon_number": codon_number,
                "stop_base_position": stop_base_position
            }

    return {
        "stop_codon": "",
        "stop_codon_number": "",
        "stop_base_position": ""
    }


# --------------------------------------------------
# Load validated transcript IDs
# --------------------------------------------------

validated_ids = set()

with open(validation_file, "r") as f:
    reader = csv.DictReader(f)

    for row in reader:

        if row["match"].lower() == "true":
            validated_ids.add(
                row["canonical_transcript_id"]
            )


print(
    f"Validated transcripts: "
    f"{len(validated_ids)}"
)


# --------------------------------------------------
# Analyze +1 and +2 frames
# --------------------------------------------------

results = []

with open(genes_file, "r") as f:

    reader = csv.DictReader(f)

    for row in reader:

        transcript_id = row[
            "canonical_transcript_id"
        ]

        if transcript_id not in validated_ids:
            continue

        cds = row["cds_sequence"].upper()

        if len(cds) < 3:
            continue

        for shift in [1, 2]:

            stop = find_first_stop(
                cds,
                shift
            )

            if stop["stop_base_position"] != "":

                relative_stop_position = (
                    stop["stop_base_position"]
                    / len(cds)
                )

                distance_from_start_codons = (
                    stop["stop_codon_number"]
                )

            else:

                relative_stop_position = ""
                distance_from_start_codons = ""

            results.append({
                "gene_id":
                    row["gene_id"],

                "gene_name":
                    row["gene_name"],

                "transcript_id":
                    transcript_id,

                "frameshift":
                    f"+{shift}",

                "cds_length":
                    len(cds),

                "stop_codon":
                    stop["stop_codon"],

                "stop_codon_number":
                    stop["stop_codon_number"],

                "stop_base_position":
                    stop["stop_base_position"],

                "relative_stop_position":
                    relative_stop_position,

                "no_stop_found":
                    stop["stop_codon"] == ""
            })


# --------------------------------------------------
# Save detailed results
# --------------------------------------------------

with open(output_file, "w", newline="") as f:

    fieldnames = results[0].keys()

    writer = csv.DictWriter(
        f,
        fieldnames=fieldnames
    )

    writer.writeheader()
    writer.writerows(results)


print(
    f"Total frame analyses: "
    f"{len(results)}"
)


# --------------------------------------------------
# Summaries
# --------------------------------------------------

for shift in ["+1", "+2"]:

    subset = [
        row
        for row in results
        if row["frameshift"] == shift
    ]

    with_stop = [
        row
        for row in subset
        if not row["no_stop_found"]
    ]

    no_stop = [
        row
        for row in subset
        if row["no_stop_found"]
    ]

    stop_positions = [
        int(row["stop_codon_number"])
        for row in with_stop
    ]

    relative_positions = [
        float(row["relative_stop_position"])
        for row in with_stop
    ]

    codon_counts = Counter(
        row["stop_codon"]
        for row in with_stop
    )

    print(f"\nFrameshift {shift}")

    print(
        f"Transcripts analyzed: "
        f"{len(subset)}"
    )

    print(
        f"Stop codon found: "
        f"{len(with_stop)}"
    )

    print(
        f"No stop codon found: "
        f"{len(no_stop)}"
    )

    if stop_positions:

        avg_stop = (
            sum(stop_positions)
            / len(stop_positions)
        )

        median_stop = sorted(
            stop_positions
        )[len(stop_positions) // 2]

        avg_relative = (
            sum(relative_positions)
            / len(relative_positions)
        )

        print(
            f"Average first stop: "
            f"{avg_stop:.1f} codons after shift"
        )

        print(
            f"Median first stop: "
            f"{median_stop} codons after shift"
        )

        print(
            f"Average relative stop position: "
            f"{avg_relative:.3f}"
        )

    print("Stop codons:")

    for codon in ["TAA", "TAG", "TGA"]:

        print(
            f"  {codon}: "
            f"{codon_counts[codon]}"
        )


# --------------------------------------------------
# Helper data for plots
# --------------------------------------------------

plus1 = [
    row
    for row in results
    if (
        row["frameshift"] == "+1"
        and not row["no_stop_found"]
    )
]

plus2 = [
    row
    for row in results
    if (
        row["frameshift"] == "+2"
        and not row["no_stop_found"]
    )
]


plus1_codons = [
    int(row["stop_codon_number"])
    for row in plus1
]

plus2_codons = [
    int(row["stop_codon_number"])
    for row in plus2
]

plus1_relative = [
    float(row["relative_stop_position"])
    for row in plus1
]

plus2_relative = [
    float(row["relative_stop_position"])
    for row in plus2
]


# --------------------------------------------------
# Figure 1: +1 first stop distribution
# --------------------------------------------------

plt.figure(figsize=(8, 5))

plt.hist(
    plus1_codons,
    bins=50
)

plt.xlabel(
    "Codons Until First Stop"
)

plt.ylabel(
    "Number of Transcripts"
)

plt.title(
    "+1 Frameshift: First Stop Codon Distribution"
)

plt.tight_layout()

plt.savefig(
    figure_plus1,
    dpi=300
)

plt.close()


# --------------------------------------------------
# Figure 2: +2 first stop distribution
# --------------------------------------------------

plt.figure(figsize=(8, 5))

plt.hist(
    plus2_codons,
    bins=50
)

plt.xlabel(
    "Codons Until First Stop"
)

plt.ylabel(
    "Number of Transcripts"
)

plt.title(
    "+2 Frameshift: First Stop Codon Distribution"
)

plt.tight_layout()

plt.savefig(
    figure_plus2,
    dpi=300
)

plt.close()


# --------------------------------------------------
# Figure 3: Compare normalized stop positions
# --------------------------------------------------

plt.figure(figsize=(8, 5))

plt.hist(
    plus1_relative,
    bins=40,
    alpha=0.6,
    label="+1"
)

plt.hist(
    plus2_relative,
    bins=40,
    alpha=0.6,
    label="+2"
)

plt.xlabel(
    "Relative Position of First Stop Codon"
)

plt.ylabel(
    "Number of Transcripts"
)

plt.title(
    "Relative Stop Position: +1 vs +2 Frameshift"
)

plt.legend()

plt.tight_layout()

plt.savefig(
    figure_compare,
    dpi=300
)

plt.close()


# --------------------------------------------------
# Figure 4: Stop codon type distribution
# --------------------------------------------------

plus1_counts = Counter(
    row["stop_codon"]
    for row in plus1
)

plus2_counts = Counter(
    row["stop_codon"]
    for row in plus2
)

codons = ["TAA", "TAG", "TGA"]

x = range(len(codons))
width = 0.35

plt.figure(figsize=(8, 5))

plt.bar(
    [i - width / 2 for i in x],
    [plus1_counts[c] for c in codons],
    width,
    label="+1"
)

plt.bar(
    [i + width / 2 for i in x],
    [plus2_counts[c] for c in codons],
    width,
    label="+2"
)

plt.xticks(
    list(x),
    codons
)

plt.xlabel(
    "Stop Codon"
)

plt.ylabel(
    "Number of Transcripts"
)

plt.title(
    "Stop Codon Distribution After Frameshift"
)

plt.legend()

plt.tight_layout()

plt.savefig(
    figure_codons,
    dpi=300
)

plt.close()


print(
    f"\nSaved results to {output_file}"
)

print("\nFigures created:")

print(figure_plus1)
print(figure_plus2)
print(figure_compare)
print(figure_codons)