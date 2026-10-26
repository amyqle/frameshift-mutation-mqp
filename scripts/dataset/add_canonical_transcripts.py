import csv

gtf_file = "data/genes.gtf"
gene_csv = "data/protein_coding_genes.csv"
output_file = "data/genes_with_canonical_transcripts.csv"

# Get protein-coding gene IDs
protein_coding_genes = set()

with open(gene_csv, "r") as f:
    reader = csv.DictReader(f)
    for row in reader:
        protein_coding_genes.add(row["gene_id"])

canonical = {}

with open(gtf_file, "r") as f:
    for line in f:
        if line.startswith("#"):
            continue

        parts = line.strip().split("\t")

        if len(parts) < 9 or parts[2] != "transcript":
            continue

        attributes = parts[8]

        if 'tag "Ensembl_canonical"' not in attributes:
            continue

        gene_id = ""
        transcript_id = ""

        for item in attributes.split(";"):
            item = item.strip()

            if item.startswith("gene_id"):
                gene_id = item.split('"')[1]

            elif item.startswith("transcript_id"):
                transcript_id = item.split('"')[1]

        if gene_id in protein_coding_genes:
            canonical[gene_id] = transcript_id


with open(gene_csv, "r") as infile, open(output_file, "w", newline="") as outfile:
    reader = csv.DictReader(infile)

    fieldnames = reader.fieldnames + ["canonical_transcript_id"]
    writer = csv.DictWriter(outfile, fieldnames=fieldnames)

    writer.writeheader()

    for row in reader:
        row["canonical_transcript_id"] = canonical.get(row["gene_id"], "")
        writer.writerow(row)


print(f"Found {len(canonical)} canonical protein-coding transcripts.")
print(f"Saved to {output_file}")