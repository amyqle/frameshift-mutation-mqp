import csv

csv.field_size_limit(10_000_000)

input_file = "data/genes_with_cds.csv"
output_file = "data/cds_translation_validation.csv"

# Standard genetic code
codon_table = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L",
    "TCT":"S","TCC":"S","TCA":"S","TCG":"S",
    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*",
    "TGT":"C","TGC":"C","TGA":"*","TGG":"W",

    "CTT":"L","CTC":"L","CTA":"L","CTG":"L",
    "CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
    "CGT":"R","CGC":"R","CGA":"R","CGG":"R",

    "ATT":"I","ATC":"I","ATA":"I","ATG":"M",
    "ACT":"T","ACC":"T","ACA":"T","ACG":"T",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K",
    "AGT":"S","AGC":"S","AGA":"R","AGG":"R",

    "GTT":"V","GTC":"V","GTA":"V","GTG":"V",
    "GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "GAT":"D","GAC":"D","GAA":"E","GAG":"E",
    "GGT":"G","GGC":"G","GGA":"G","GGG":"G"
}


def translate(sequence):
    protein = ""

    for i in range(0, len(sequence) - 2, 3):
        codon = sequence[i:i+3].upper()
        amino_acid = codon_table.get(codon, "X")

        if amino_acid == "*":
            break

        protein += amino_acid

    return protein


rows = []

with open(input_file, "r") as f:
    reader = csv.DictReader(f)

    for row in reader:
        cds = row["cds_sequence"]
        reference_protein = row["protein_sequence"]

        translated = translate(cds)

        match = translated == reference_protein

        rows.append({
            "gene_id": row["gene_id"],
            "gene_name": row["gene_name"],
            "canonical_transcript_id": row["canonical_transcript_id"],
            "cds_length": len(cds),
            "reference_protein_length": len(reference_protein),
            "translated_protein_length": len(translated),
            "match": match
        })


with open(output_file, "w", newline="") as f:
    fieldnames = rows[0].keys()

    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)


matches = sum(row["match"] for row in rows)

print(f"Total genes checked: {len(rows)}")
print(f"Exact matches: {matches}")
print(f"Non-matches: {len(rows) - matches}")

print("\nFirst 10 genes:")
for row in rows[:10]:
    print(row["gene_name"], row["match"])

print("\nLast 10 genes:")
for row in rows[-10:]:
    print(row["gene_name"], row["match"])

print(f"\nSaved to {output_file}")