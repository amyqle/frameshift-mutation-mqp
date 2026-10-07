import csv
from collections import Counter

csv.field_size_limit(10_000_000)

gtf_file = "data/genes.gtf"
csv_file = "data/final_gene_dataset.csv"

# Load gene information
genes = {}

with open(csv_file, "r") as f:
    reader = csv.DictReader(f)

    for row in reader:
        genes[row["canonical_transcript_id"]] = {
            "start": int(row["start"]),
            "end": int(row["end"]),
            "strand": row["strand"],
            "sequence": row["genomic_sequence"]
        }

# Collect start/stop codon coordinates
codon_parts = {
    "start_codon": {},
    "stop_codon": {}
}

with open(gtf_file, "r") as f:
    for line in f:
        if line.startswith("#"):
            continue

        parts = line.strip().split("\t")

        if len(parts) < 9:
            continue

        feature = parts[2]

        if feature not in ("start_codon", "stop_codon"):
            continue

        start = int(parts[3])
        end = int(parts[4])
        attributes = parts[8]

        transcript_id = ""

        for item in attributes.split(";"):
            item = item.strip()

            if item.startswith("transcript_id"):
                transcript_id = item.split('"')[1]

        if transcript_id in genes:
            codon_parts[feature].setdefault(
                transcript_id, []
            ).append((start, end))


def extract_sequence(transcript_id, pieces):

    gene = genes[transcript_id]

    gene_start = gene["start"]
    gene_end = gene["end"]
    strand = gene["strand"]
    sequence = gene["sequence"]

    if strand == "+":
        pieces.sort()
    else:
        pieces.sort(reverse=True)

    result = ""

    for start, end in pieces:

        if strand == "+":
            relative_start = start - gene_start
            relative_end = end - gene_start + 1

        else:
            relative_start = gene_end - end
            relative_end = gene_end - start + 1

        result += sequence[relative_start:relative_end]

    return result.upper()


start_counts = Counter()
stop_counts = Counter()

for transcript_id, pieces in codon_parts["start_codon"].items():
    codon = extract_sequence(transcript_id, pieces)
    start_counts[codon] += 1

for transcript_id, pieces in codon_parts["stop_codon"].items():
    codon = extract_sequence(transcript_id, pieces)
    stop_counts[codon] += 1


print("Start codons:")
for codon, count in start_counts.most_common():
    print(f"{codon}: {count}")

print("\nStop codons:")
for codon, count in stop_counts.most_common():
    print(f"{codon}: {count}")