import csv
import gzip

csv_file = "genes_with_cdna.csv"
fasta_file = "pep.fa.bgz"
output_file = "genes_with_cdna_and_protein.csv"

protein_sequences = {}

with gzip.open(fasta_file, "rt") as f:
    transcript_id = None
    sequence_parts = []

    for line in f:
        line = line.strip()

        if line.startswith(">"):
            if transcript_id is not None:
                protein_sequences[transcript_id] = "".join(sequence_parts)

            header = line

            transcript_id = None

            for item in header.split():
                if item.startswith("transcript:"):
                    transcript_id = item.split(":")[1].split(".")[0]

            sequence_parts = []

        else:
            sequence_parts.append(line)

    if transcript_id is not None:
        protein_sequences[transcript_id] = "".join(sequence_parts)


with open(csv_file, "r") as infile, open(output_file, "w", newline="") as outfile:
    reader = csv.DictReader(infile)

    fieldnames = reader.fieldnames + ["protein_sequence"]
    writer = csv.DictWriter(outfile, fieldnames=fieldnames)

    writer.writeheader()

    matched = 0

    for row in reader:
        transcript_id = row["canonical_transcript_id"]

        sequence = protein_sequences.get(transcript_id, "")

        if sequence:
            matched += 1

        row["protein_sequence"] = sequence
        writer.writerow(row)


print(f"Matched {matched} protein sequences.")
print(f"Saved to {output_file}")