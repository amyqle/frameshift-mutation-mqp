import csv
import gzip

csv_file = "genes_with_cdna_and_protein.csv"
fasta_file = "unmasked.fa.bgz"
output_file = "final_gene_dataset.csv"

# Load chromosome sequences
chromosomes = {}

with gzip.open(fasta_file, "rt") as f:
    chromosome = None
    sequence_parts = []

    for line in f:
        line = line.strip()

        if line.startswith(">"):
            if chromosome is not None:
                chromosomes[chromosome] = "".join(sequence_parts)

            chromosome = line.split()[0][1:]
            sequence_parts = []

        else:
            sequence_parts.append(line)

    if chromosome is not None:
        chromosomes[chromosome] = "".join(sequence_parts)

# Add genomic sequence for each gene
with open(csv_file, "r") as infile, open(output_file, "w", newline="") as outfile:
    reader = csv.DictReader(infile)

    fieldnames = reader.fieldnames + ["genomic_sequence"]
    writer = csv.DictWriter(outfile, fieldnames=fieldnames)

    writer.writeheader()

    matched = 0

    for row in reader:
        chromosome = row["chromosome"]
        start = int(row["start"])
        end = int(row["end"])
        strand = row["strand"]

        if chromosome in chromosomes:
            # GTF coordinates are 1-based
            sequence = chromosomes[chromosome][start - 1:end]

            # Reverse complement genes on the negative strand
            if strand == "-":
                complement = str.maketrans("ACGTNacgtn", "TGCANtgcan")
                sequence = sequence.translate(complement)[::-1]

            row["genomic_sequence"] = sequence
            matched += 1
        else:
            row["genomic_sequence"] = ""

        writer.writerow(row)

print(f"Matched {matched} genomic sequences.")
print(f"Saved to {output_file}")