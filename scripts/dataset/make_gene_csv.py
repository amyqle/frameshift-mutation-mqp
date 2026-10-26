import csv

input_file = "data/genes.gtf"
output_file = "data/protein_coding_genes.csv"

rows = []

with open(input_file, "r") as f:
    for line in f:
        if line.startswith("#"):
            continue

        parts = line.strip().split("\t")

        if len(parts) < 9:
            continue

        chromosome = parts[0]
        feature = parts[2]
        start = parts[3]
        end = parts[4]
        strand = parts[6]
        attributes = parts[8]

        if feature == "gene" and 'gene_biotype "protein_coding"' in attributes:
            gene_id = ""
            gene_name = ""

            for item in attributes.split(";"):
                item = item.strip()

                if item.startswith("gene_id"):
                    gene_id = item.split('"')[1]

                elif item.startswith("gene_name"):
                    gene_name = item.split('"')[1]

            rows.append([
                gene_id,
                gene_name,
                chromosome,
                start,
                end,
                strand
            ])

with open(output_file, "w", newline="") as f:
    writer = csv.writer(f)

    writer.writerow([
        "gene_id",
        "gene_name",
        "chromosome",
        "start",
        "end",
        "strand"
    ])

    writer.writerows(rows)

print(f"Saved {len(rows)} genes to {output_file}")
