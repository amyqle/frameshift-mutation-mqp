import csv
from collections import Counter

csv.field_size_limit(10_000_000)

input_file = "data/cds_mismatches.csv"
output_file = "data/cds_mismatch_classification.csv"

classified = []
counts = Counter()

with open(input_file, "r") as f:
    reader = csv.DictReader(f)

    for row in reader:
        cds_length = int(row["cds_length"])
        start_codon = row["start_codon"]
        reference_length = int(row["reference_protein_length"])
        translated_length = int(row["translated_protein_length"])

        categories = []

        if cds_length % 3 != 0:
            categories.append("CDS_not_divisible_by_3")

        if start_codon != "ATG":
            categories.append("non_ATG_start")

        if translated_length < reference_length:
            categories.append("translated_protein_shorter")

        if translated_length > reference_length:
            categories.append("translated_protein_longer")

        if translated_length == reference_length:
            categories.append("same_length_but_sequence_differs")

        if not categories:
            categories.append("other")

        for category in categories:
            counts[category] += 1

        row["classification"] = ";".join(categories)
        classified.append(row)


with open(output_file, "w", newline="") as f:
    fieldnames = list(classified[0].keys())

    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(classified)


print("Mismatch classifications:\n")

for category, count in counts.most_common():
    print(f"{category}: {count}")

print(f"\nSaved to {output_file}")