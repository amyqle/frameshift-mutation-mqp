import csv

csv.field_size_limit(10_000_000)

validation_file = "data/cds_translation_validation.csv"
genes_file = "data/genes_with_cds.csv"
output_file = "data/simple_frameshift_model_results.csv"

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
    stop_position = None

    for i in range(0, len(sequence) - 2, 3):
        codon = sequence[i:i+3].upper()
        amino_acid = codon_table.get(codon, "X")

        if amino_acid == "*":
            stop_position = (i // 3) + 1
            break

        protein += amino_acid

    return protein, stop_position


def find_first_difference(reference, mutated):
    for i, (a, b) in enumerate(zip(reference, mutated)):
        if a != b:
            return i + 1

    if len(reference) != len(mutated):
        return min(len(reference), len(mutated)) + 1

    return None


def find_downstream_atg(sequence, start_position):
    pos = sequence.find("ATG", start_position)

    if pos == -1:
        return None

    return pos


# Load validated transcript IDs
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


selected_genes = []

with open(genes_file, "r") as f:
    reader = csv.DictReader(f)

    for row in reader:
        if row["canonical_transcript_id"] in selected_ids:
            selected_genes.append(row)


results = []

for row in selected_genes:

    cds = row["cds_sequence"].upper()
    reference_protein = row["protein_sequence"]

    # Midpoint mutation for pilot model
    mutation_position = len(cds) // 2

    if mutation_position < 3:
        mutation_position = 3

    for shift_size in [1, 2]:

        mutated_cds = (
            cds[:mutation_position] +
            cds[mutation_position + shift_size:]
        )

        mutated_protein, stop_position = translate(mutated_cds)

        first_difference = find_first_difference(
            reference_protein,
            mutated_protein
        )

        # Find first downstream ATG
        downstream_atg = find_downstream_atg(
            mutated_cds,
            mutation_position
        )

        downstream_atg_position = None
        distance_to_downstream_atg = None
        candidate_protein_length = None
        candidate_stop_position = None

        if downstream_atg is not None:

            downstream_atg_position = downstream_atg + 1

            distance_to_downstream_atg = (
                downstream_atg - mutation_position
            )

            candidate_sequence = mutated_cds[downstream_atg:]

            candidate_protein, candidate_stop_position = translate(
                candidate_sequence
            )

            candidate_protein_length = len(candidate_protein)

        results.append({
            "gene_id": row["gene_id"],
            "gene_name": row["gene_name"],
            "transcript_id": row["canonical_transcript_id"],
            "frameshift": f"+{shift_size}",
            "mutation_base_position": mutation_position + 1,
            "reference_cds_length": len(cds),
            "mutated_cds_length": len(mutated_cds),
            "reference_protein_length": len(reference_protein),
            "mutated_protein_length": len(mutated_protein),
            "first_difference_position": first_difference,
            "premature_stop":
                len(mutated_protein) < len(reference_protein),
            "stop_position_in_mutated_protein": stop_position,

            "downstream_ATG_found":
                downstream_atg is not None,

            "downstream_ATG_position":
                downstream_atg_position,

            "distance_from_mutation_to_ATG":
                distance_to_downstream_atg,

            "candidate_protein_length":
                candidate_protein_length,

            "candidate_stop_position":
                candidate_stop_position
        })


with open(output_file, "w", newline="") as f:
    fieldnames = results[0].keys()

    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(results)


print(f"Validated genes available: {len(validated_ids)}")
print(f"Test genes used: {len(selected_genes)}")
print(f"Total model outputs: {len(results)}")


for shift in ["+1", "+2"]:

    subset = [
        row for row in results
        if row["frameshift"] == shift
    ]

    premature = sum(
        row["premature_stop"] for row in subset
    )

    with_atg = [
        row for row in subset
        if row["downstream_ATG_found"]
    ]

    avg_mutated_length = sum(
        row["mutated_protein_length"]
        for row in subset
    ) / len(subset)

    print(f"\nFrameshift {shift}")
    print(
        f"Average mutated protein length: "
        f"{avg_mutated_length:.1f}"
    )
    print(
        f"Premature stops: "
        f"{premature}/{len(subset)}"
    )
    print(
        f"Downstream ATG found: "
        f"{len(with_atg)}/{len(subset)}"
    )

    if with_atg:

        avg_distance = sum(
            row["distance_from_mutation_to_ATG"]
            for row in with_atg
        ) / len(with_atg)

        avg_candidate_length = sum(
            row["candidate_protein_length"]
            for row in with_atg
        ) / len(with_atg)

        print(
            f"Average distance to downstream ATG: "
            f"{avg_distance:.1f} bases"
        )

        print(
            f"Average candidate protein length: "
            f"{avg_candidate_length:.1f} amino acids"
        )


print(f"\nSaved to {output_file}")