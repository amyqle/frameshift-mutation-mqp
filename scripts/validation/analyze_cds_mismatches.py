import csv

csv.field_size_limit(10_000_000)

input_file = "data/genes_with_cds.csv"
output_file = "data/cds_mismatches.csv"

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


mismatches = []

with open(input_file, "r") as f:
    reader = csv.DictReader(f)

    for row in reader:
        cds = row["cds_sequence"]
        reference = row["protein_sequence"]
        translated = translate(cds)

        if translated != reference:

            first_difference = None

            for i, (a, b) in enumerate(zip(translated, reference)):
                if a != b:
                    first_difference = i + 1
                    break

            if first_difference is None and len(translated) != len(reference):
                first_difference = min(len(translated), len(reference)) + 1

            mismatches.append({
                "gene_id": row["gene_id"],
                "gene_name": row["gene_name"],
                "transcript_id": row["canonical_transcript_id"],
                "cds_length": len(cds),
                "cds_divisible_by_3": len(cds) % 3 == 0,
                "start_codon": cds[:3],
                "end_codon": cds[-3:],
                "reference_protein_length": len(reference),
                "translated_protein_length": len(translated),
                "first_difference_position": first_difference
            })


with open(output_file, "w", newline="") as f:
    fieldnames = mismatches[0].keys()

    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(mismatches)


print(f"Total mismatches: {len(mismatches)}")

print("\nCDS length divisible by 3:")
divisible = sum(row["cds_divisible_by_3"] for row in mismatches)
print(f"Yes: {divisible}")
print(f"No: {len(mismatches) - divisible}")

print("\nFirst 10 mismatches:")
for row in mismatches[:10]:
    print(
        row["gene_name"],
        "start:", row["start_codon"],
        "CDS:", row["cds_length"],
        "reference:", row["reference_protein_length"],
        "translated:", row["translated_protein_length"]
    )

print(f"\nSaved to {output_file}")