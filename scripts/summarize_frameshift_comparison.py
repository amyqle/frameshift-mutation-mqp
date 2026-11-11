import csv
import matplotlib.pyplot as plt

input_file = "data/frameshift_stop_comparison.csv"

bar_chart_file = "figures/frameshift_survival_comparison.png"
extreme_file = "data/extreme_frameshift_examples.csv"


# --------------------------------------------------
# Load comparison results
# --------------------------------------------------

rows = []

with open(input_file, "r") as f:
    reader = csv.DictReader(f)

    for row in reader:
        rows.append(row)


# --------------------------------------------------
# Count which frame survives longer
# --------------------------------------------------

plus1_longer = sum(
    row["longer_surviving_frame"] == "+1"
    for row in rows
)

plus2_longer = sum(
    row["longer_surviving_frame"] == "+2"
    for row in rows
)

both_no_stop = sum(
    row["longer_surviving_frame"] == "both_no_stop"
    for row in rows
)


# --------------------------------------------------
# Bar chart
# --------------------------------------------------

labels = [
    "+1 survives longer",
    "+2 survives longer",
    "No stop in either"
]

values = [
    plus1_longer,
    plus2_longer,
    both_no_stop
]

plt.figure(figsize=(8, 5))

bars = plt.bar(labels, values)

plt.ylabel("Number of Transcripts")
plt.title("Which Frameshift Reading Frame Survives Longer?")

for bar, value in zip(bars, values):
    plt.text(
        bar.get_x() + bar.get_width() / 2,
        value + 150,
        str(value),
        ha="center"
    )

plt.tight_layout()

plt.savefig(
    bar_chart_file,
    dpi=300
)

plt.close()


# --------------------------------------------------
# Find extreme examples
# --------------------------------------------------

comparable = [
    row for row in rows
    if row["stop_difference_plus2_minus_plus1"] != ""
]

for row in comparable:
    row["stop_difference_plus2_minus_plus1"] = int(
        row["stop_difference_plus2_minus_plus1"]
    )

plus2_extreme = sorted(
    comparable,
    key=lambda row:
        row["stop_difference_plus2_minus_plus1"],
    reverse=True
)[:5]

plus1_extreme = sorted(
    comparable,
    key=lambda row:
        row["stop_difference_plus2_minus_plus1"]
)[:5]


extreme_rows = []

for row in plus2_extreme:
    extreme_rows.append({
        "gene_name": row["gene_name"],
        "plus1_stop": row[
            "plus1_stop_codons_after_shift"
        ],
        "plus2_stop": row[
            "plus2_stop_codons_after_shift"
        ],
        "difference": row[
            "stop_difference_plus2_minus_plus1"
        ],
        "longer_frame": "+2"
    })

for row in plus1_extreme:
    extreme_rows.append({
        "gene_name": row["gene_name"],
        "plus1_stop": row[
            "plus1_stop_codons_after_shift"
        ],
        "plus2_stop": row[
            "plus2_stop_codons_after_shift"
        ],
        "difference": row[
            "stop_difference_plus2_minus_plus1"
        ],
        "longer_frame": "+1"
    })


# --------------------------------------------------
# Save extreme examples
# --------------------------------------------------

with open(extreme_file, "w", newline="") as f:

    fieldnames = [
        "gene_name",
        "plus1_stop",
        "plus2_stop",
        "difference",
        "longer_frame"
    ]

    writer = csv.DictWriter(
        f,
        fieldnames=fieldnames
    )

    writer.writeheader()
    writer.writerows(extreme_rows)


# --------------------------------------------------
# Print summary
# --------------------------------------------------

print("Frameshift survival summary")

print(
    f"+1 survives longer: {plus1_longer}"
)

print(
    f"+2 survives longer: {plus2_longer}"
)

print(
    f"No stop in either frame: {both_no_stop}"
)


print(
    "\nTop 5 examples where +2 survives longer:"
)

for row in plus2_extreme:
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
    "\nTop 5 examples where +1 survives longer:"
)

for row in plus1_extreme:
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
    f"\nSaved chart to {bar_chart_file}"
)

print(
    f"Saved examples to {extreme_file}"
)