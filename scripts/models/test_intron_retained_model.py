
import csv
import gzip
import sys
from pathlib import Path

limit = sys.maxsize

while True:
    try:
        csv.field_size_limit(limit)
        break
    except OverflowError:
        limit //= 10

from Bio import SeqIO
from Bio.Seq import Seq


# Project folders
ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data"

GENE_CSV = DATA / "genes_with_cds.csv"
GTF_FILE = DATA / "genes.gtf"
GENOME_FILE = DATA / "unmasked.fa.bgz"

OUTPUT_FILE = DATA / "intron_retained_model_results.csv"

TEST_LIMIT = 20
DELETION_SIZES = [1, 2]


def find_column(columns, possible_names):
    """Find a column even if its name differs slightly."""
    lookup = {c.lower().strip(): c for c in columns}

    for name in possible_names:
        if name.lower() in lookup:
            return lookup[name.lower()]

    return None


def normalize_id(value):
    """Allow IDs with or without their version suffix."""
    return str(value or "").split(".")[0]


def read_gtf():
    """Collect CDS coordinates for each transcript."""
    transcripts = {}

    # The GTF file is gzip-compressed, despite its .gtf extension.
    with gzip.open(GTF_FILE, "rt", encoding="utf-8") as file:
        for line in file:
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) != 9 or fields[2] != "CDS":
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
                    "start": start,
                    "end": end,
                }
            else:
                transcripts[transcript_id]["start"] = min(
                    transcripts[transcript_id]["start"], start
                )
                transcripts[transcript_id]["end"] = max(
                    transcripts[transcript_id]["end"], end
                )

    return transcripts


def find_transcript(row, transcripts, columns):
    """Match a CSV gene to its canonical transcript when available."""
    gene_col = find_column(
        columns, ["gene_name", "gene_symbol", "gene", "external_gene_name"]
    )
    gene_id_col = find_column(columns, ["gene_id", "ensembl_gene_id"])
    transcript_col = find_column(
        columns,
        [
            "canonical_transcript_id",
            "canonical_transcript",
            "transcript_id",
        ],
    )

    transcript_id = row.get(transcript_col, "") if transcript_col else ""
    gene_name = row.get(gene_col, "") if gene_col else ""
    gene_id = row.get(gene_id_col, "") if gene_id_col else ""

    # Prefer the canonical transcript if the CSV contains it.
    if transcript_id:
        for tid, info in transcripts.items():
            if normalize_id(tid) == normalize_id(transcript_id):
                return tid, info

    # Otherwise match by gene name or gene ID.
    matches = []
    for tid, info in transcripts.items():
        name_match = (
            gene_name
            and info["gene_name"].upper() == gene_name.upper()
        )
        id_match = (
            gene_id
            and normalize_id(info["gene_id"]) == normalize_id(gene_id)
        )

        if name_match or id_match:
            matches.append((tid, info))

    if not matches:
        return None, None

    # If there are multiple transcripts, use the longest genomic span.
    return max(
        matches,
        key=lambda item: item[1]["end"] - item[1]["start"]
    )


def get_intron_retained_sequence(info, genome):
    """
    Extract the full genomic span of the coding region.

    Unlike the normal CDS, this includes the introns between
    the coding exons.
    """
    chromosome = info["chromosome"]

    # Handle FASTA headers that use either "1" or "chr1".
    candidates = [
        chromosome,
        "chr" + chromosome,
        chromosome.removeprefix("chr"),
    ]

    fasta_name = next(
        (name for name in candidates if name in genome),
        None
    )

    if fasta_name is None:
        raise KeyError(
            f"Chromosome {chromosome} was not found in the genome FASTA."
        )

    start = info["start"]
    end = info["end"]

    # GTF coordinates start at 1; Python sequence slices start at 0.
    sequence = str(genome[fasta_name].seq[start - 1:end]).upper()

    # Put the sequence in the direction of transcription.
    if info["strand"] == "-":
        sequence = str(Seq(sequence).reverse_complement())

    return sequence


def translate_sequence(sequence):
    """Translate complete codons and preserve stop codons."""
    complete_length = len(sequence) - (len(sequence) % 3)

    if complete_length == 0:
        return ""

    return str(
        Seq(sequence[:complete_length]).translate(to_stop=False)
    )


def analyze_mutation(sequence, deletion_size):
    """Delete bases, translate, and measure the predicted outcome."""
    mutation_position = len(sequence) // 2

    # Delete one or two bases near the middle of the sequence.
    mutated_sequence = (
        sequence[:mutation_position]
        + sequence[mutation_position + deletion_size:]
    )

    protein = translate_sequence(mutated_sequence)

    first_stop = protein.find("*")
    premature_stop = (
        first_stop != -1 and first_stop < len(protein) - 1
    )

    if first_stop != -1:
        mutated_protein_length = first_stop
        stop_base_position = first_stop * 3
        search_start = stop_base_position + 3
    else:
        mutated_protein_length = len(protein)
        search_start = mutation_position

    # Find the next ATG downstream of the premature stop.
    atg_position = mutated_sequence.find("ATG", search_start)

    if atg_position != -1:
        candidate_protein = translate_sequence(
            mutated_sequence[atg_position:]
        )
        candidate_stop = candidate_protein.find("*")

        if candidate_stop != -1:
            candidate_length = candidate_stop
        else:
            candidate_length = len(candidate_protein)

        distance_to_atg = atg_position - search_start
    else:
        candidate_length = 0
        distance_to_atg = ""

    return {
        "mutation_position": mutation_position,
        "mutated_sequence_length": len(mutated_sequence),
        "mutated_protein_length": mutated_protein_length,
        "premature_stop": premature_stop,
        "downstream_atg_found": atg_position != -1,
        "distance_to_downstream_atg": distance_to_atg,
        "candidate_protein_length": candidate_length,
    }


def main():
    if not GENE_CSV.exists():
        raise FileNotFoundError(f"Missing file: {GENE_CSV}")

    if not GTF_FILE.exists():
        raise FileNotFoundError(f"Missing file: {GTF_FILE}")

    if not GENOME_FILE.exists():
        raise FileNotFoundError(f"Missing file: {GENOME_FILE}")

    transcripts = read_gtf()
    print(f"Transcripts with CDS annotations: {len(transcripts)}")

    with open(GENE_CSV, "r", newline="", encoding="utf-8-sig") as file:
        reader = csv.DictReader(file)
        columns = reader.fieldnames
        rows = list(reader)

    if not columns:
        raise ValueError("The gene CSV has no column headers.")

    # Index the BGZF-compressed genome without loading the whole
    # human genome into memory.
    genome = SeqIO.index(str(GENOME_FILE), "fasta")

    results = []
    used_genes = set()

    try:
        for row in rows:
            transcript_id, info = find_transcript(
                row, transcripts, columns
            )

            if info is None:
                continue

            gene_name = info["gene_name"] or info["gene_id"] or transcript_id

            if gene_name in used_genes:
                continue

            sequence = get_intron_retained_sequence(info, genome)

            if len(sequence) < 6:
                continue

            used_genes.add(gene_name)

            for deletion_size in DELETION_SIZES:
                result = analyze_mutation(sequence, deletion_size)

                results.append({
                    "gene_name": gene_name,
                    "transcript_id": transcript_id,
                    "model": "intron_retained",
                    "frameshift": f"+{deletion_size}",
                    "retained_sequence_length": len(sequence),
                    **result,
                })

            if len(used_genes) >= TEST_LIMIT:
                break

    finally:
        genome.close()

    if not results:
        raise ValueError(
            "No genes were matched. Check the CSV column names, "
            "GTF transcript IDs, and genome FASTA chromosome names."
        )

    OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)

    with open(OUTPUT_FILE, "w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=results[0].keys())
        writer.writeheader()
        writer.writerows(results)

    print(f"\nGenes tested: {len(used_genes)}")
    print(f"Total model outputs: {len(results)}")

    for deletion_size in DELETION_SIZES:
        subset = [
            row for row in results
            if row["frameshift"] == f"+{deletion_size}"
        ]

        if not subset:
            continue

        print(f"\nFrameshift +{deletion_size}")
        print(
            "Average mutated protein length:",
            round(
                sum(r["mutated_protein_length"] for r in subset)
                / len(subset), 1
            ),
        )
        print(
            "Premature stops:",
            f"{sum(r['premature_stop'] for r in subset)}/{len(subset)}",
        )
        print(
            "Downstream ATG found:",
            f"{sum(r['downstream_atg_found'] for r in subset)}/{len(subset)}",
        )

        distances = [
            r["distance_to_downstream_atg"]
            for r in subset
            if r["distance_to_downstream_atg"] != ""
        ]

        if distances:
            print(
                "Average distance to downstream ATG:",
                round(sum(distances) / len(distances), 1),
                "bases",
            )

        print(
            "Average candidate protein length:",
            round(
                sum(r["candidate_protein_length"] for r in subset)
                / len(subset), 1
            ),
            "amino acids",
        )

    print(f"\nSaved to {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
