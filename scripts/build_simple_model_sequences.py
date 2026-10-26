import csv

csv.field_size_limit(10_000_000)

input_file = "data/genes_with_cds.csv"
validation_file = "data/cds_translation_validation.csv"
output_file = "data/simple_model_sequences.csv"


# --------------------------------------------------
# Load validated transcript IDs
# --------------------------------------------------

validated_ids = set()

with open(validation_file, "r") as f:
    reader = csv.DictReader(f)

    for row in reader:
        if row["match"].lower() == "true":
            validated_ids.add(row["canonical_transcript_id"])


print(f"Validated transcripts: {len(validated_ids)}")


# --------------------------------------------------
# Build CDS + downstream transcript sequence
# --------------------------------------------------

results = []
not_found = []

with open(input_file, "r") as f:
    reader = csv.DictReader(f)

    for row in reader:

        transcript_id = row["canonical_transcript_id"]

        if transcript_id not in validated_ids:
            continue

        cdna = row["cdna_sequence"].upper()
        cds = row["cds_sequence"].upper()

        # Find the CDS within the spliced cDNA
        cds_start = cdna.find(cds)

        if cds_start == -1:
            not_found.append(
                (row["gene_name"], transcript_id)
            )
            continue

        cds_end = cds_start + len(cds)

        # Sequence after the CDS
        downstream_sequence = cdna[cds_end:]

        # Simple model input:
        # CDS + downstream transcript sequence
        simple_model_sequence = cdna[cds_start:]

        results.append({
            "gene_id": row["gene_id"],
            "gene_name": row["gene_name"],
            "canonical_transcript_id": transcript_id,

            "cdna_length": len(cdna),
            "cds_length": len(cds),

            "cds_start_in_cdna": cds_start + 1,
            "cds_end_in_cdna": cds_end,

            "downstream_sequence_length":
                len(downstream_sequence),

            "simple_model_sequence_length":
                len(simple_model_sequence),

            "protein_sequence":
                row["protein_sequence"],

            "cds_sequence":
                cds,

            "downstream_sequence":
                downstream_sequence,

            "simple_model_sequence":
                simple_model_sequence
        })


# --------------------------------------------------
# Save output
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

print(f"Simple-model sequences created: {len(results)}")
print(f"CDS not found inside cDNA: {len(not_found)}")

if results:

    avg_downstream = sum(
        row["downstream_sequence_length"]
        for row in results
    ) / len(results)

    print(
        f"Average downstream sequence length: "
        f"{avg_downstream:.1f} bp"
    )


if not_found:

    print("\nFirst CDS/cDNA mismatches:")

    for gene, transcript in not_found[:10]:
        print(gene, transcript)


print(f"\nSaved to {output_file}")