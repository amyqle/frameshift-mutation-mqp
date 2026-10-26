import csv

csv.field_size_limit(10_000_000)

validation_file = "data/cds_translation_validation.csv"
genes_file = "data/genes_with_cds.csv"
output_file = "data/test_frameshift_deletion_results.csv"

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


# 1. Get only validated transcript IDs
validated_ids = []

with open(validation_file, "r") as f:
    reader = csv.DictReader(f)

    for row in reader:
        if row["match"].lower() == "true":
            validated_ids.append(row["canonical_transcript_id"])


print(f"Validated genes available: {len(validated_ids)}")


# 2. Select first 10 + last 10 validated genes
selected_ids = set(
    validated_ids[:10] +
    validated_ids[-10:]
)


# 3. Load those genes
selected_genes = []

with open(genes_file, "r") as f:
    reader = csv.DictReader(f)

    for row in reader:
        if row["canonical_transcript_id"] in selected_ids:
            selected_genes.append(row)


# 4. Apply a 1-base deletion
results = []

for row in selected_genes:

    cds = row["cds_sequence"]
    reference_protein = row["protein_sequence"]

    # Delete one base near the beginning,
    # but AFTER the start codon so ATG stays intact
    deletion_position = 3

    mutated_cds = (
        cds[:deletion_position] +
        cds[deletion_position + 1:]
    )

    mutated_protein = translate(mutated_cds)

    # Find first amino-acid difference
    first_difference = None

    for i, (ref_aa, mut_aa) in enumerate(
        zip(reference_protein, mutated_protein)
    ):
        if ref_aa != mut_aa:
            first_difference = i + 1
            break

    if first_difference is None:
        if len(reference_protein) != len(mutated_protein):
            first_difference = min(
                len(reference_protein),
                len(mutated_protein)
            ) + 1

    results.append({
        "gene_id": row["gene_id"],
        "gene_name": row["gene_name"],
        "transcript_id": row["canonical_transcript_id"],
        "deletion_position": deletion_position + 1,
        "reference_cds_length": len(cds),
        "mutated_cds_length": len(mutated_cds),
        "reference_protein_length": len(reference_protein),
        "mutated_protein_length": len(mutated_protein),
        "first_difference_position": first_difference,
        "reference_protein_first_20": reference_protein[:20],
        "mutated_protein_first_20": mutated_protein[:20]
    })


# 5. Save results
with open(output_file, "w", newline="") as f:
    fieldnames = results[0].keys()

    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(results)


print(f"Test genes used: {len(results)}")
print("\nResults:")

for row in results:
    print(
        row["gene_name"],
        "| ref protein:", row["reference_protein_length"],
        "| mutated:", row["mutated_protein_length"],
        "| first difference:", row["first_difference_position"]
    )

print(f"\nSaved to {output_file}")