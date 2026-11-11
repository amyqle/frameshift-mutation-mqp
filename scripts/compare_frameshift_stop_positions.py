import csv
import math

import matplotlib.pyplot as plt

csv.field_size_limit(10_000_000)

input_file = "data/frameshift_stop_positions.csv"

output_file = "data/frameshift_stop_comparison.csv"

figure_difference = "figures/frameshift_stop_position_difference.png"
figure_scatter = "figures/plus1_vs_plus2_stop_position.png"


# --------------------------------------------------
# Load +1 and +2 results by transcript
# --------------------------------------------------

transcripts = {}

with open(input_file, "r") as f:
    reader = csv.DictReader(f)

    for row in reader:

        transcript_id = row["transcript_id"]

        if transcript_id not in transcripts:
            transcripts[transcript_id] = {
                "gene_id": row["gene_id"],
                "gene_name": row["gene_name"],
                "transcript_id": transcript_id,
                "cds_length": int(row["cds_length"])
            }

        transcripts[transcript_id][row["frameshift"]] = row


# --------------------------------------------------
# Build gene-by-gene comparison
# --------------------------------------------------

results = []

for transcript_id, data in transcripts.items():

    if "+1" not in data or "+2" not in data:
        continue

    plus1 = data["+1"]
    plus2 = data["+2"]

    plus1_no_stop = (
        plus1["no_stop_found"].lower() == "true"
    )

    plus2_no_stop = (
        plus2["no_stop_found"].lower() == "true"
    )

    plus1_stop = (
        int(plus1["stop_codon_number"])
        if not plus1_no_stop
        else None
    )

    plus2_stop = (
        int(plus2["stop_codon_number"])
        if not plus2_no_stop
        else None
    )

    plus1_relative = (
        float(plus1["relative_stop_position"])
        if not plus1_no_stop
        else None
    )

    plus2_relative = (
        float(plus2["relative_stop_position"])
        if not plus2_no_stop
        else None
    )


    # --------------------------------------------------
    # Determine which frame survives longer
    # --------------------------------------------------

    if plus1_no_stop and plus2_no_stop:
        longer_frame = "both_no_stop"
        stop_difference = ""

    elif plus1_no_stop:
        longer_frame = "+1"
        stop_difference = ""

    elif plus2_no_stop:
        longer_frame = "+2"
        stop_difference = ""

    else:

        stop_difference = (
            plus2_stop - plus1_stop
        )

        if plus1_stop > plus2_stop:
            longer_frame = "+1"

        elif plus2_stop > plus1_stop:
            longer_frame = "+2"

        else:
            longer_frame = "same"


    results.append({
        "gene_id":
            data["gene_id"],

        "gene_name":
            data["gene_name"],

        "transcript_id":
            transcript_id,

        "cds_length":
            data["cds_length"],

        "plus1_first_stop_codon":
            plus1["stop_codon"],

        "plus1_stop_codons_after_shift":
            "" if plus1_stop is None else plus1_stop,

        "plus1_relative_stop_position":
            "" if plus1_relative is None else plus1_relative,

        "plus1_no_stop":
            plus1_no_stop,

        "plus2_first_stop_codon":
            plus2["stop_codon"],

        "plus2_stop_codons_after_shift":
            "" if plus2_stop is None else plus2_stop,

        "plus2_relative_stop_position":
            "" if plus2_relative is None else plus2_relative,

        "plus2_no_stop":
            plus2_no_stop,

        "stop_difference_plus2_minus_plus1":
            stop_difference,

        "longer_surviving_frame":
            longer_frame
    })


# --------------------------------------------------
# Save comparison table
# --------------------------------------------------

with open(output_file, "w", newline="") as f:

    fieldnames = results[0].keys()

    writer = csv.DictWriter(
        f,
        fieldnames=fieldnames
    )

    writer.writeheader()
    writer.writerows(results)


# --------------------------------------------------
# Summary counts
# --------------------------------------------------

plus1_longer = sum(
    row["longer_surviving_frame"] == "+1"
    for row in results
)

plus2_longer = sum(
    row["longer_surviving_frame"] == "+2"
    for row in results
)

same = sum(
    row["longer_surviving_frame"] == "same"
    for row in results
)

plus1_no_stop = sum(
    row["plus1_no_stop"]
    for row in results
)

plus2_no_stop = sum(
    row["plus2_no_stop"]
    for row in results
)

both_no_stop = sum(
    row["longer_surviving_frame"] == "both_no_stop"
    for row in results
)


print(
    f"Transcripts compared: "
    f"{len(results)}"
)

print("\nWhich frame survives longer?")

print(
    f"+1 survives longer: "
    f"{plus1_longer}"
)

print(
    f"+2 survives longer: "
    f"{plus2_longer}"
)

print(
    f"Same stop position: "
    f"{same}"
)

print("\nNo-stop cases")

print(
    f"+1 no stop: "
    f"{plus1_no_stop}"
)

print(
    f"+2 no stop: "
    f"{plus2_no_stop}"
)

print(
    f"No stop in either frame: "
    f"{both_no_stop}"
)


# --------------------------------------------------
# Largest +1 vs +2 differences
# --------------------------------------------------

comparable = [
    row
    for row in results
    if row["stop_difference_plus2_minus_plus1"] != ""
]


sorted_positive = sorted(
    comparable,
    key=lambda row:
        row["stop_difference_plus2_minus_plus1"],
    reverse=True
)

sorted_negative = sorted(
    comparable,
    key=lambda row:
        row["stop_difference_plus2_minus_plus1"]
)


print(
    "\nTop 10 genes where +2 survives "
    "much longer than +1:"
)

for row in sorted_positive[:10]:

    print(
        row["gene_name"],
        " +1:",
        row["plus1_stop_codons_after_shift"],
        " +2:",
        row["plus2_stop_codons_after_shift"],
        " difference:",
        row["stop_difference_plus2_minus_plus1"]
    )


print(
    "\nTop 10 genes where +1 survives "
    "much longer than +2:"
)

for row in sorted_negative[:10]:

    print(
        row["gene_name"],
        " +1:",
        row["plus1_stop_codons_after_shift"],
        " +2:",
        row["plus2_stop_codons_after_shift"],
        " difference:",
        row["stop_difference_plus2_minus_plus1"]
    )


# --------------------------------------------------
# Figure 1:
# Distribution of +2 minus +1 stop difference
# --------------------------------------------------

differences = [
    row["stop_difference_plus2_minus_plus1"]
    for row in comparable
]

plt.figure(figsize=(8, 5))

plt.hist(
    differences,
    bins=60
)

plt.axvline(
    0,
    linestyle="--"
)

plt.xlabel(
    "+2 Stop Position − +1 Stop Position (codons)"
)

plt.ylabel(
    "Number of Transcripts"
)

plt.title(
    "Difference Between +1 and +2 First Stop Positions"
)

plt.tight_layout()

plt.savefig(
    figure_difference,
    dpi=300
)

plt.close()


# --------------------------------------------------
# Figure 2:
# +1 vs +2 stop positions
# --------------------------------------------------

x = [
    row["plus1_stop_codons_after_shift"]
    for row in comparable
]

y = [
    row["plus2_stop_codons_after_shift"]
    for row in comparable
]

plt.figure(figsize=(7, 7))

plt.scatter(
    x,
    y,
    alpha=0.3,
    s=10
)

max_value = max(
    max(x),
    max(y)
)

plt.plot(
    [0, max_value],
    [0, max_value],
    linestyle="--"
)

plt.xlabel(
    "+1 First Stop Position (codons)"
)

plt.ylabel(
    "+2 First Stop Position (codons)"
)

plt.title(
    "+1 vs +2 First Stop Position by Transcript"
)

plt.tight_layout()

plt.savefig(
    figure_scatter,
    dpi=300
)

plt.close()


print(
    f"\nSaved comparison results to "
    f"{output_file}"
)

print("\nFigures created:")

print(figure_difference)
print(figure_scatter)