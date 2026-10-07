import csv

# Allow very large sequence fields
csv.field_size_limit(10_000_000)

gtf_file = "data/genes.gtf"
input_csv = "data/final_gene_dataset.csv"
output_csv = "data/genes_with_cds.csv"

canonical_transcripts = set()
rows = []

with open(input_csv, "r") as f:
    reader = csv.DictReader(f)

    for row in reader:
        rows.append(row)
        canonical_transcripts.add(row["canonical_transcript_id"])


# Collect CDS coordinates
cds_parts = {}

with open(gtf_file, "r") as f:
    for line in f:
        if line.startswith("#"):
            continue

        parts = line.strip().split("\t")

        if len(parts) < 9 or parts[2] != "CDS":
            continue

        start = int(parts[3])
        end = int(parts[4])
        attributes = parts[8]

        transcript_id = ""

        for item in attributes.split(";"):
            item = item.strip()

            if item.startswith("transcript_id"):
                transcript_id = item.split('"')[1]

        if transcript_id in canonical_transcripts:
            cds_parts.setdefault(transcript_id, []).append((start, end))


# Build CDS sequences
matched = 0

for row in rows:

    transcript_id = row["canonical_transcript_id"]

    gene_start = int(row["start"])
    gene_end = int(row["end"])

    strand = row["strand"]
    genomic_sequence = row["genomic_sequence"]

    pieces = cds_parts.get(transcript_id, [])

    if not pieces:
        row["cds_sequence"] = ""
        continue

    # Put CDS pieces in transcript order
    if strand == "+":
        pieces.sort()
    else:
        pieces.sort(reverse=True)

    cds_sequence = ""

    for start, end in pieces:

        if strand == "+":
            relative_start = start - gene_start
            relative_end = end - gene_start + 1

        else:
            relative_start = gene_end - end
            relative_end = gene_end - start + 1

        cds_sequence += genomic_sequence[
            relative_start:relative_end
        ]

    row["cds_sequence"] = cds_sequence
    matched += 1


fieldnames = list(rows[0].keys())

with open(output_csv, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)


print(f"Matched {matched} CDS sequences.")
print(f"Saved to {output_csv}")