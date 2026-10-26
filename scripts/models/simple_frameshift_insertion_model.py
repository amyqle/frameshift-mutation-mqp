import csv
import itertools

csv.field_size_limit(10_000_000)

input_file = "data/simple_model_sequences.csv"
validation_file = "data/cds_translation_validation.csv"
output_file = "data/simple_frameshift_insertion_results.csv"


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
    stop_aa_position = None
    stop_nt_position = None

    for i in range(0, len(sequence) - 2, 3):

        codon = sequence[i:i+3].upper()
        amino_acid = codon_table.get(codon, "X")

        if amino_acid == "*":
            stop_aa_position = (i // 3) + 1
            stop_nt_position = i + 1
            break

        protein += amino_acid

    return protein, stop_aa_position, stop_nt_position


def find_first_difference(reference, mutated):

    for i, (ref_aa, mut_aa) in enumerate(
        zip(reference, mutated)
    ):
        if ref_aa != mut_aa:
            return i + 1

    if len(reference) != len(mutated):
        return min(
            len(reference),
            len(mutated)
        ) + 1

    return None


# --------------------------------------------------
# Load validated transcript IDs
# --------------------------------------------------

validated_ids = set()

with open(validation_file, "r") as f:
    reader = csv.DictReader(f)

    for row in reader:
        if row["match"].lower() == "true":
            validated_ids.add(
                row["canonical_transcript_id"]
            )


print(
    f"Validated transcripts available: "
    f"{len(validated_ids)}"
)


# --------------------------------------------------
# Load simple-model sequences
# --------------------------------------------------

genes = []

with open(input_file, "r") as f:
    reader = csv.DictReader(f)

    for row in reader:
        if (
            row["canonical_transcript_id"]
            in validated_ids
        ):
            genes.append(row)


# --------------------------------------------------
# Full dataset: all validated genes
# --------------------------------------------------

selected_genes = genes


# --------------------------------------------------
# All possible insertions
# --------------------------------------------------

bases = ["A", "T", "C", "G"]

plus1_insertions = bases

plus2_insertions = [
    "".join(combo)
    for combo in itertools.product(
        bases,
        repeat=2
    )
]


# --------------------------------------------------
# Simulate frameshift insertions
# --------------------------------------------------

results = []

for row in selected_genes:

    transcript_id = (
        row["canonical_transcript_id"]
    )

    sequence = (
        row["simple_model_sequence"].upper()
    )

    cds_length = int(
        row["cds_length"]
    )

    reference_protein = (
        row["protein_sequence"]
    )

    # Midpoint of CDS for pilot
    mutation_position = cds_length // 2

    insertion_sets = {
        "+1": plus1_insertions,
        "+2": plus2_insertions
    }

    for frameshift_type, insertions \
            in insertion_sets.items():

        for inserted_sequence in insertions:

            insertion_length = len(
                inserted_sequence
            )

            mutated_sequence = (
                sequence[:mutation_position]
                + inserted_sequence
                + sequence[mutation_position:]
            )

            (
                mutated_protein,
                stop_aa_position,
                stop_nt_position
            ) = translate(mutated_sequence)

            difference_position = \
                find_first_difference(
                    reference_protein,
                    mutated_protein
                )

            # The original CDS boundary moves
            # because bases were inserted
            mutated_cds_boundary = (
                cds_length + insertion_length
            )

            # True only if the new stop codon
            # begins after the original coding region
            if stop_nt_position is not None:
                entered_downstream_region = (
                    stop_nt_position
                    > mutated_cds_boundary
                )
            else:
                # No stop codon was found before
                # the available transcript sequence ended
                entered_downstream_region = (
                    len(mutated_sequence)
                    > mutated_cds_boundary
                )

            premature_stop = (
                stop_aa_position is not None
                and stop_aa_position
                <= len(reference_protein)
            )

            results.append({
                "gene_id":
                    row["gene_id"],

                "gene_name":
                    row["gene_name"],

                "transcript_id":
                    transcript_id,

                "frameshift":
                    frameshift_type,

                "inserted_sequence":
                    inserted_sequence,

                "mutation_base_position":
                    mutation_position + 1,

                "reference_protein_length":
                    len(reference_protein),

                "mutated_protein_length":
                    len(mutated_protein),

                "first_difference_position":
                    difference_position,

                "stop_aa_position":
                    stop_aa_position,

                "stop_nt_position":
                    stop_nt_position,

                "premature_stop":
                    premature_stop,

                "entered_downstream_region":
                    entered_downstream_region,

                "mutated_protein_sequence":
                    mutated_protein
            })


# --------------------------------------------------
# Save results
# --------------------------------------------------

with open(output_file, "w", newline="") as f:

    fieldnames = results[0].keys()

    writer = csv.DictWriter(
        f,
        fieldnames=fieldnames
    )

    writer.writeheader()
    writer.writerows(results)


# --------------------------------------------------
# Summary
# --------------------------------------------------

print(
    f"Genes modeled: "
    f"{len(selected_genes)}"
)

print(
    f"Total model outputs: "
    f"{len(results)}"
)


for shift in ["+1", "+2"]:

    subset = [
        row
        for row in results
        if row["frameshift"] == shift
    ]

    premature = sum(
        row["premature_stop"]
        for row in subset
    )

    entered_downstream = sum(
        row["entered_downstream_region"]
        for row in subset
    )

    avg_length = sum(
        row["mutated_protein_length"]
        for row in subset
    ) / len(subset)

    no_stop = sum(
        row["stop_nt_position"] is None
        for row in subset
    )

    print(f"\nFrameshift {shift}")

    print(
        f"Outputs: "
        f"{len(subset)}"
    )

    print(
        f"Average mutated protein length: "
        f"{avg_length:.1f}"
    )

    print(
        f"Premature stops: "
        f"{premature}/{len(subset)}"
    )

    print(
        f"Entered downstream region: "
        f"{entered_downstream}/{len(subset)}"
    )

    print(
        f"No stop codon found in available sequence: "
        f"{no_stop}/{len(subset)}"
    )


print(
    f"\nSaved to {output_file}"
)