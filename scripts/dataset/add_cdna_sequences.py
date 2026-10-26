import csv
import gzip

csv_file = "data/genes_with_canonical_transcripts.csv"
fasta_file = "data/cdna.fa.bgz"
output_file = "data/genes_with_cdna.csv"

cdna_sequences = {}

with gzip.open(fasta_file, "rt") as f:
    transcript_id = None
    sequence_parts = []

    for line in f:
        line = line.strip()

        if line.startswith(">"):
            if transcript_id is not None:
                cdna_sequences[transcript_id] = "".join(sequence_parts)

            transcript_id = line.split()[0][1:].split(".")[0]
            sequence_parts = []

        else:
            sequence_parts.append(line)

    if transcript_id is not None:
        cdna_sequences[transcript_id] = "".join(sequence_parts)


with open(csv_file, "r") as infile, open(output_file, "w", newline="") as outfile:
    reader = csv.DictReader(infile)

    fieldnames = reader.fieldnames + ["cdna_sequence"]
    writer = csv.DictWriter(outfile, fieldnames=fieldnames)

    writer.writeheader()

    matched = 0

    for row in reader:
        transcript_id = row["canonical_transcript_id"]

        sequence = cdna_sequences.get(transcript_id, "")

        if sequence:
            matched += 1

        row["cdna_sequence"] = sequence
        writer.writerow(row)


print(f"Matched {matched} cDNA sequences.")
print(f"Saved to {output_file}")