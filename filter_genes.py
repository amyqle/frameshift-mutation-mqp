input_file = "genes.gtf"
output_file = "protein_coding_genes.txt"

genes = set()

with open(input_file, "r") as f:
    for line in f:
        if line.startswith("#"):
            continue

        if '\tgene\t' in line and 'gene_biotype "protein_coding"' in line:
            parts = line.strip().split("\t")
            attributes = parts[8]

            for item in attributes.split(";"):
                item = item.strip()
                if item.startswith("gene_name"):
                    gene_name = item.split('"')[1]
                    genes.add(gene_name)

with open(output_file, "w") as f:
    for gene in sorted(genes):
        f.write(gene + "\n")

print(f"Found {len(genes)} protein-coding genes.")
