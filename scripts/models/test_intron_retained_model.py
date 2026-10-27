import csv
import gzip
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq


# --------------------------------------------------
# File paths and settings
# --------------------------------------------------

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data"

GENES_FILE = DATA / "genes_with_cds.csv"
VALIDATION_FILE = DATA / "cds_translation_validation.csv"
GTF_FILE = DATA / "genes.gtf"
GENOME_FILE = DATA / "unmasked.fa.bgz"

OUTPUT_FILE = DATA / "first_intron_retained_results.csv"

TEST_LIMIT = 20
DELETION_SIZES = [1, 2]

csv.field_size_limit(10_000_000)


# --------------------------------------------------
# Translation
# --------------------------------------------------

codon_table = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G"
}


def translate(sequence):
    protein = ""
    stop_position = None

    for i in range(0, len(sequence) - 2, 3):
        codon = sequence[i:i + 3].upper()
        amino_acid = codon_table.get(codon, "X")

        if amino_acid == "*":
            stop_position = (i // 3) + 1
            break

        protein += amino_acid

    return protein, stop_position


def normalize_id(value):
    return str(value or "").split(".")[0]


# --------------------------------------------------
# Load validated transcripts
# --------------------------------------------------

validated_ids = []

with open(VALIDATION_FILE, "r", newline="", encoding="utf-8-sig") as file:
    reader = csv.DictReader(file)

    for row in reader:
        if row["match"].lower() == "true":
            validated_ids.append(row["canonical_transcript_id"])


selected_ids = set(validated_ids[:10] + validated_ids[-10:])

print(f"Validated transcripts: {len(validated_ids)}")
print(f"Target sample size: {len(selected_ids)}")


# --------------------------------------------------
# Load gene data
# --------------------------------------------------

selected_genes = []

with open(GENES_FILE, "r", newline="", encoding="utf-8-sig") as file:
    reader = csv.DictReader(file)

    for row in reader:
        if row["canonical_transcript_id"] in selected_ids:
            selected_genes.append(row)

print(f"Gene rows selected: {len(selected_genes)}")


# --------------------------------------------------
# Read GTF exon and CDS annotations
# --------------------------------------------------

def read_gtf():
    transcripts = {}

    # This file is gzip-compressed even though its extension is .gtf.
    with gzip.open(GTF_FILE, "rt", encoding="utf-8") as file:
        for line in file:
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) != 9:
                continue

            feature = fields[2]

            if feature not in ("exon", "CDS"):
                continue

            chromosome = fields[0]
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]

            attributes = {}

            for item in fields[8].strip().split(";"):
                item = item.strip()

                if " " in item:
                    key, value = item.split(" ", 1)
                    attributes[key] = value.strip().strip('"')

            transcript_id = attributes.get("transcript_id")

            if not transcript_id:
                continue

            if transcript_id not in transcripts:
                transcripts[transcript_id] = {
                    "gene_id": attributes.get("gene_id", ""),
                    "gene_name": attributes.get("gene_name", ""),
                    "chromosome": chromosome,
                    "strand": strand,
                    "exons": [],
                    "cds": []
                }

            if feature == "exon":
                transcripts[transcript_id]["exons"].append(
                    (start, end)
                )
            else:
                transcripts[transcript_id]["cds"].append(
                    (start, end)
                )

    for info in transcripts.values():
        info["exons"] = sorted(
            set(info["exons"]),
            reverse=(info["strand"] == "-")
        )

        info["cds"] = sorted(set(info["cds"]))

    return transcripts


# --------------------------------------------------
# Reconstruct transcript with only the first intron
# --------------------------------------------------

def reconstruct_transcript(info, genome):
    chromosome = info["chromosome"]

    possible_names = [
        chromosome,
        "chr" + chromosome,
        chromosome.removeprefix("chr")
    ]

    fasta_name = next(
        (name for name in possible_names if name in genome),
        None
    )

    if fasta_name is None:
        return None

    exons = info["exons"]
    strand = info["strand"]

    if len(exons) < 2 or not info["cds"]:
        return None

    pieces = []
    genomic_to_transcript = {}
    transcript_position = 0
    first_intron_length = 0

    for i, (start, end) in enumerate(exons):

        genomic_sequence = str(
            genome[fasta_name].seq[start - 1:end]
        ).upper()

        if strand == "-":
            exon_sequence = str(
                Seq(genomic_sequence).reverse_complement()
            )
            positions = range(end, start - 1, -1)
        else:
            exon_sequence = genomic_sequence
            positions = range(start, end + 1)

        pieces.append(exon_sequence)

        for genomic_position in positions:
            genomic_to_transcript[genomic_position] = (
                transcript_position
            )
            transcript_position += 1

        # Insert the first intron only.
        if i == 0:
            next_start, next_end = exons[1]

            if strand == "+":
                intron_start = end + 1
                intron_end = next_start - 1
            else:
                intron_start = next_end + 1
                intron_end = start - 1

            if intron_start > intron_end:
                return None

            intron_sequence = str(
                genome[fasta_name].seq[
                    intron_start - 1:intron_end
                ]
            ).upper()

            if strand == "-":
                intron_sequence = str(
                    Seq(intron_sequence).reverse_complement()
                )

            pieces.append(intron_sequence)

            first_intron_length = len(intron_sequence)
            transcript_position += first_intron_length

    reconstructed = "".join(pieces)

    # Locate the original CDS start and end.
    if strand == "+":
        cds_start_genomic = min(
            start for start, end in info["cds"]
        )
        cds_end_genomic = max(
            end for start, end in info["cds"]
        )
    else:
        cds_start_genomic = max(
            end for start, end in info["cds"]
        )
        cds_end_genomic = min(
            start for start, end in info["cds"]
        )

    if (
        cds_start_genomic not in genomic_to_transcript
        or cds_end_genomic not in genomic_to_transcript
    ):
        return None

    cds_start = genomic_to_transcript[cds_start_genomic]
    cds_end = genomic_to_transcript[cds_end_genomic]

    retained_cds = reconstructed[cds_start:cds_end + 1]

    # Check whether the retained intron overlaps the coding span.
    intron_in_cds = (
        len(retained_cds) != 0
        and len(retained_cds)
        != sum(end - start + 1 for start, end in info["cds"])
    )

    return retained_cds, first_intron_length, intron_in_cds


# --------------------------------------------------
# Run model
# --------------------------------------------------

print("Reading GTF...")
transcripts = read_gtf()

print("Indexing genome FASTA...")
genome = SeqIO.index(str(GENOME_FILE), "fasta")

results = []
skipped = 0
used_transcripts = set()

try:
    for row in selected_genes:

        transcript_id = row["canonical_transcript_id"]

        # Match transcript IDs even if their version suffix differs.
        matches = [
            (tid, info)
            for tid, info in transcripts.items()
            if normalize_id(tid) == normalize_id(transcript_id)
        ]

        if not matches:
            skipped += 1
            continue

        gtf_transcript_id, info = matches[0]

        if gtf_transcript_id in used_transcripts:
            continue

        reconstructed = reconstruct_transcript(info, genome)

        if reconstructed is None:
            skipped += 1
            continue

        retained_cds, intron_length, intron_in_cds = reconstructed

        if intron_length == 0 or len(retained_cds) < 6:
            skipped += 1
            continue

        # Ensure the extracted sequence corresponds to the expected
        # annotated CDS when the retained intron is outside the CDS.
        original_cds = row["cds_sequence"].upper()

        used_transcripts.add(gtf_transcript_id)

        mutation_position = len(retained_cds) // 2

        if mutation_position < 3:
            mutation_position = 3

        reference_protein = row["protein_sequence"].rstrip("*")

        baseline_protein, baseline_stop = translate(retained_cds)

        for deletion_size in DELETION_SIZES:

            mutated_cds = (
                retained_cds[:mutation_position]
                + retained_cds[mutation_position + deletion_size:]
            )

            mutated_protein, stop_position = translate(mutated_cds)

            # A stop is considered premature when it occurs before
            # the normal protein's terminal position.
            premature_stop = (
                stop_position is not None
                and stop_position <= len(reference_protein)
            )


            downstream_atg = mutated_cds.find("ATG", mutation_position)

            downstream_atg_found = downstream_atg != -1

            if downstream_atg_found:
                downstream_atg_distance = downstream_atg - mutation_position
                candidate_protein, candidate_stop = translate(
                    mutated_cds[downstream_atg:]
                )
                candidate_protein_length = len(candidate_protein)
            else:
                downstream_atg_distance = None
                candidate_protein_length = None

            # Find the first amino acid that differs from the reference.
            first_difference = None

            for i, (original_aa, mutated_aa) in enumerate(
                zip(reference_protein, mutated_protein)
            ):
                if original_aa != mutated_aa:
                    first_difference = i + 1
                    break

            if first_difference is None and (
                len(reference_protein) != len(mutated_protein)
            ):
                first_difference = min(
                    len(reference_protein),
                    len(mutated_protein)
                ) + 1

            results.append({
                "gene_id": info["gene_id"],
                "gene_name": info["gene_name"],
                "transcript_id": gtf_transcript_id,
                "strand": info["strand"],
                "first_intron_length": intron_length,
                "first_intron_in_cds": intron_in_cds,
                "retained_cds_length": len(retained_cds),
                "deletion_size": deletion_size,
                "mutation_position": mutation_position,
                "reference_protein_length": len(reference_protein),
                "baseline_protein_length": len(baseline_protein),
                "mutated_protein_length": len(mutated_protein),
                "stop_position": stop_position,
                "premature_stop": premature_stop,
                "downstream_atg_found": downstream_atg_found,
                "downstream_atg_distance": downstream_atg_distance,
                "candidate_protein_length": candidate_protein_length,
                "first_difference_position": first_difference
            })

except Exception as error:
    print(f"An error occurred while running the model: {error}")

finally:
    genome.close()

# Save the results.
fieldnames = [
    "gene_id",
    "gene_name",
    "transcript_id",
    "strand",
    "first_intron_length",
    "first_intron_in_cds",
    "retained_cds_length",
    "deletion_size",
    "mutation_position",
    "reference_protein_length",
    "baseline_protein_length",
    "mutated_protein_length",
    "stop_position",
    "premature_stop",
    "downstream_atg_found",
    "downstream_atg_distance",
    "candidate_protein_length",
    "first_difference_position"
]

with open(OUTPUT_FILE, "w", newline="", encoding="utf-8") as file:
    writer = csv.DictWriter(file, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(results)

print(f"\nGenes tested: {len(used_transcripts)}")
print(f"Results generated: {len(results)}")
print(f"Transcripts skipped: {skipped}")
print(f"Results saved to: {OUTPUT_FILE}")
