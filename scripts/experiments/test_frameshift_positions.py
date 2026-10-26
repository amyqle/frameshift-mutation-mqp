import csv

csv.field_size_limit(10_000_000)

validation_file = "data/cds_translation_validation.csv"
genes_file = "data/genes_with_cds.csv"
output_file = "data/test_frameshift_positions_results.csv"

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


# Get validated transcript IDs
validated_ids = []

with open(validation_file, "r") as f:
    reader = csv.DictReader(f)

    for row in reader:
        if row["match"].lower() == "true":
            validated_ids.append(row["canonical_transcript_id"])


# First 10 + last 10 validated genes
selected_ids = set(
    validated_ids[:10] +
    validated_ids[-10:]
)


# Load selected genes
selected_genes = []

with open(genes_file, "r") as f:
    reader = csv.DictReader(f)

    for row in reader:
        if row["canonical_transcript_id"] in selected_ids:
            selected_genes.append(row)


positions = {
    "25%": 0.25,
    "50%": 0.50,
    "75%": 0.75
}

results = []


for row in selected_genes:

    cds = row["cds_sequence"]
    reference_protein = row["protein_sequence"]

    for position_name, fraction in positions.items():

        # Calculate deletion position
        deletion_position = int(len(cds) * fraction)

        # Avoid deleting inside the start codon
        if deletion_position < 3:
            deletion_position = 3

        # Delete 1 nucleotide
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

        premature_stop = len(mutated_protein) < len(reference_protein)

        results.append({
            "gene_id": row["gene_id"],
            "gene_name": row["gene_name"],
            "transcript_id": row["canonical_transcript_id"],
            "mutation_position": position_name,
            "deleted_base_position": deletion_position + 1,
            "reference_cds_length": len(cds),
            "mutated_cds_length": len(mutated_cds),
            "reference_protein_length": len(reference_protein),
            "mutated_protein_length": len(mutated_protein),
            "protein_length_change":
                len(mutated_protein) - len(reference_protein),
            "first_difference_position": first_difference,
            "premature_stop": premature_stop
        })


with open(output_file, "w", newline="") as f:
    fieldnames = results[0].keys()

    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(results)


print(f"Validated genes available: {len(validated_ids)}")
print(f"Test genes used: {len(selected_genes)}")
print(f"Total mutations tested: {len(results)}")


for position in positions:

    subset = [
        row for row in results
        if row["mutation_position"] == position
    ]

    avg_length = sum(
        row["mutated_protein_length"] for row in subset
    ) / len(subset)

    premature_count = sum(
        row["premature_stop"] for row in subset
    )

    print(f"\nMutation at {position}")
    print(f"Average mutated protein length: {avg_length:.1f}")
    print(
        f"Premature stops: "
        f"{premature_count}/{len(subset)}"
    )


print(f"\nSaved to {output_file}")