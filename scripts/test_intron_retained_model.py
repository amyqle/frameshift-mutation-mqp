
import csv
import gzip
import sys
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq


# --------------------------------------------------
# File paths and settings
# --------------------------------------------------

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"

GENE_CSV = DATA / "genes_with_cds.csv"
GTF_FILE = DATA / "genes.gtf"
GENOME_FILE = DATA / "unmasked.fa.bgz"

OUTPUT_FILE = DATA / "first_intron_retained_results.csv"

TEST_LIMIT = 20
DELETION_SIZES = [1, 2]


# Allow CSV fields larger than the default 131072 characters.
limit = sys.maxsize

while True:
    try:
        csv.field_size_limit(limit)
        break
    except OverflowError:
        limit //= 10


# --------------------------------------------------
# Helper functions
# --------------------------------------------------

def find_column(columns, possible_names):
    lookup = {c.lower().strip(): c for c in columns if c}

    for name in possible_names:
        if name.lower() in lookup:
            return lookup[name.lower()]

    return None


def normalize_id(value):
    return str(value or "").split(".")[0]


# --------------------------------------------------
# Read GTF annotations
# --------------------------------------------------

def read_gtf():
    transcripts = {}

    # genes.gtf is gzip-compressed even though its extension is .gtf
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
                    "cds": [],
                }

            info = transcripts[transcript_id]

            if feature == "exon":
                info["exons"].append((start, end))
            else:
                info["cds"].append((start, end))

    # Remove duplicate exon/CDS intervals and sort exons in
    # transcript order (5' to 3').
    for info in transcripts.values():
        info["exons"] = sorted(
            set(info["exons"]),
            reverse=(info["strand"] == "-"),
        )

        info["cds"] = sorted(set(info["cds"]))

    return transcripts


# --------------------------------------------------
# Match CSV genes to GTF transcripts
# --------------------------------------------------

def find_transcript(row, columns, transcripts):
    gene_col = find_column(
        columns,
        ["gene_name", "gene_symbol", "gene", "external_gene_name"],
    )

    gene_id_col = find_column(
        columns,
        ["gene_id", "ensembl_gene_id"],
    )

    transcript_col = find_column(
        columns,
        [
            "canonical_transcript_id",
            "canonical_transcript",
            "transcript_id",
            "ensembl_transcript_id",
        ],
    )

    transcript_id = row.get(transcript_col, "") if transcript_col else ""
    gene_name = row.get(gene_col, "") if gene_col else ""
    gene_id = row.get(gene_id_col, "") if gene_id_col else ""

    # Prefer an exact transcript ID match.
    if transcript_id:
        for tid, info in transcripts.items():
            if normalize_id(tid) == normalize_id(transcript_id):
                if info["exons"] and info["cds"]:
                    return tid, info

    # Otherwise, match by gene ID or gene name.
    matches = []

    for tid, info in transcripts.items():
        if not info["exons"] or not info["cds"]:
            continue

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

    # Choose the transcript with the greatest annotated CDS span.
    return max(
        matches,
        key=lambda item: (
            max(end for start, end in item[1]["cds"])
            - min(start for start, end in item[1]["cds"])
        ),
    )


# --------------------------------------------------
# Construct a transcript retaining only its first intron
# --------------------------------------------------

def get_first_intron_retained_cds(info, genome):
    chromosome = info["chromosome"]

    candidates = [
        chromosome,
        "chr" + chromosome,
        chromosome.removeprefix("chr"),
    ]

    fasta_name = next(
        (name for name in candidates if name in genome),
        None,
    )

    if fasta_name is None:
        raise KeyError(
            f"Chromosome {chromosome} was not found in the genome FASTA."
        )

    exons = info["exons"]
    cds_intervals = info["cds"]
    strand = info["strand"]

    if len(exons) < 2 or not cds_intervals:
        return None

    # Store every genomic base's position in the reconstructed
    # transcript. This allows the annotated CDS boundaries to
    # be located after inserting the first intron.
    pieces = []
    genomic_to_transcript = {}
    transcript_position = 0
    retained_intron_length = 0

    for i, (start, end) in enumerate(exons):
        genomic_sequence = str(
            genome[fasta_name].seq[start - 1:end]
        ).upper()

        if strand == "-":
            exon_sequence = str(
                Seq(genomic_sequence).reverse_complement()
            )
            genomic_positions = range(end, start - 1, -1)
        else:
            exon_sequence = genomic_sequence
            genomic_positions = range(start, end + 1)

        pieces.append(exon_sequence)

        for genomic_position in genomic_positions:
            genomic_to_transcript[genomic_position] = transcript_position
            transcript_position += 1

        # Insert only the first intron, between the first and
        # second exons in transcript (5' to 3') order.
        if i == 0:
            next_start, next_end = exons[1]

            if strand == "+":
                intron_start = end + 1
                intron_end = next_start - 1
            else:
                intron_start = next_end + 1
                intron_end = start - 1

            if intron_start <= intron_end:
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
                retained_intron_length = len(intron_sequence)
                transcript_position += retained_intron_length

    retained_transcript = "".join(pieces)

    # Find the annotated translation start and end using the CDS
    # coordinates, accounting for the transcript strand.
    if strand == "+":
        cds_start_genomic = min(start for start, end in cds_intervals)
        cds_end_genomic = max(end for start, end in cds_intervals)
    else:
        cds_start_genomic = max(end for start, end in cds_intervals)
        cds_end_genomic = min(start for start, end in cds_intervals)

    if (
        cds_start_genomic not in genomic_to_transcript
        or cds_end_genomic not in genomic_to_transcript
    ):
        return None

    cds_start_position = genomic_to_transcript[cds_start_genomic]
    cds_end_position = genomic_to_transcript[cds_end_genomic]

    # Extract the annotated coding span. If the retained first
    # intron lies inside the CDS, it remains in this sequence.
    coding_sequence = retained_transcript[
        cds_start_position:cds_end_position + 1
    ]

    if len(coding_sequence) < 6:
        return None

    return coding_sequence, retained_intron_length


# --------------------------------------------------
# Translation and frameshift analysis
# --------------------------------------------------

def translate_sequence(sequence):
    complete_length = len(sequence) - (len(sequence) % 3)

    if complete_length == 0:
        return ""

    return str(
        Seq(sequence[:complete_length]).translate(to_stop=False)
    )


def analyze_mutation(sequence, deletion_size):
    # Keep the same midpoint deletion assumption as the previous model.
    mutation_position = len(sequence) // 2

    mutated_sequence = (
        sequence[:mutation_position]
        + sequence[mutation_position + deletion_size:]
    )

    protein = translate_sequence(mutated_sequence)

    first_stop = protein.find("*")

    premature_stop = (
        first_stop != -1
        and first_stop < len(protein) - 1
    )

    if first_stop != -1:
        mutated_protein_length = first_stop
        stop_base_position = first_stop * 3
        search_start = stop_base_position + 3
    else:
        mutated_protein_length = len(protein)
        search_start = mutation_position

    # Search for a potential downstream start codon.
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
        "retained_cds_length": len(sequence),
        "mutated_sequence_length": len(mutated_sequence),
        "mutated_protein_length": mutated_protein_length,
        "premature_stop": premature_stop,
        "downstream_atg_found": atg_position != -1,
        "distance_to_downstream_atg": distance_to_atg,
        "candidate_protein_length": candidate_length,
    }


# --------------------------------------------------
# Main program
# --------------------------------------------------

def main():
    for path in (GENE_CSV, GTF_FILE, GENOME_FILE):
        if not path.exists():
            raise FileNotFoundError(f"Required file not found: {path}")

    print("Reading GTF annotations...", flush=True)
    transcripts = read_gtf()

    print(
        f"Transcripts with exon and CDS annotations: "
        f"{sum(bool(t['exons'] and t['cds']) for t in transcripts.values())}",
        flush=True,
    )

    print("Reading gene CSV...", flush=True)

    with open(
        GENE_CSV, "r", newline="", encoding="utf-8-sig"
    ) as file:
        reader = csv.DictReader(file)
        columns = reader.fieldnames

        if not columns:
            raise ValueError("The gene CSV has no column headers.")

        rows = list(reader)

    print(f"Loaded {len(rows)} CSV rows.", flush=True)
    print("Indexing genome FASTA...", flush=True)

    genome = SeqIO.index(str(GENOME_FILE), "fasta")

    results = []
    used_transcripts = set()
    skipped = 0

    try:
        for row in rows:
            transcript_id, info = find_transcript(
                row, columns, transcripts
            )

            if info is None or transcript_id in used_transcripts:
                continue

            extracted = get_first_intron_retained_cds(info, genome)

            if extracted is None:
                skipped += 1
                continue

            sequence, intron_length = extracted

            # Do not include transcripts with no first intron.
            if intron_length == 0:
                skipped += 1
                continue

            used_transcripts.add(transcript_id)

            gene_name = (
                info["gene_name"]
                or info["gene_id"]
                or transcript_id
            )

            for deletion_size in DELETION_SIZES:
                result = analyze_mutation(sequence, deletion_size)

                results.append({
                    "gene_name": gene_name,
                    "gene_id": info["gene_id"],
                    "transcript_id": transcript_id,
                    "model": "first_intron_retained",
                    "frameshift": f"+{deletion_size}",
                    "first_intron_length": intron_length,
                    **result,
                })

            if len(used_transcripts) >= TEST_LIMIT:
                break

    finally:
        genome.close()

    if not results:
        raise ValueError(
            "No results were generated. Check the CSV column names "
            "and transcript matching."
        )

    OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)

    with open(
        OUTPUT_FILE, "w", newline="", encoding="utf-8"
    ) as file:
        writer = csv.DictWriter(
            file, fieldnames=list(results[0].keys())
        )
        writer.writeheader()
        writer.writerows(results)

    result_df = None

    print("\nFirst-intron-retained model results", flush=True)
    print(f"Genes tested: {len(used_transcripts)}")
    print(f"Total model outputs: {len(results)}")
    print(f"Transcripts skipped: {skipped}")

    for frameshift in ["+1", "+2"]:
        group = [
            r for r in results if r["frameshift"] == frameshift
        ]

        if not group:
            continue

        premature_count = sum(
            r["premature_stop"] for r in group
        )

        atg_count = sum(
            r["downstream_atg_found"] for r in group
        )

        distances = [
            r["distance_to_downstream_atg"]
            for r in group
            if r["distance_to_downstream_atg"] != ""
        ]

        print(f"\nFrameshift {frameshift}")
        print(
            "Average mutated protein length: "
            f"{sum(r['mutated_protein_length'] for r in group) / len(group):.1f}"
        )
        print(
            f"Premature stops: {premature_count}/{len(group)}"
        )
        print(
            f"Downstream ATG found: {atg_count}/{len(group)}"
        )

        if distances:
            print(
                "Average distance to downstream ATG: "
                f"{sum(distances) / len(distances):.1f} bases"
            )

        print(
            "Average candidate protein length: "
            f"{sum(r['candidate_protein_length'] for r in group) / len(group):.1f} amino acids"
        )

    print(f"\nSaved to: {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
