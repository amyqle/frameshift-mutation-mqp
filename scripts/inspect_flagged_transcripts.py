import csv
from pathlib import Path
from Bio.Seq import Seq

csv.field_size_limit(10_000_000)

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"

GENES_FILE = DATA / "genes_with_cds.csv"
VALIDATION_FILE = DATA / "cds_translation_validation.csv"
OUTPUT_FILE = DATA / "flagged_transcript_details.csv"


def normalize_protein(sequence):
    return str(sequence or "").strip().upper()


def remove_terminal_stop(sequence):
    """
    Remove only a stop codon at the very end.
    Internal * characters are preserved.
    """
    sequence = normalize_protein(sequence)

    if sequence.endswith("*"):
        return sequence[:-1]

    return sequence


def translate_cds(cds):
    """
    Translate the CDS using the standard genetic code.
    Keep * for stop codons so internal stops remain visible.
    """
    cds = str(cds or "").strip().upper().replace("U", "T")

    complete_length = len(cds) - (len(cds) % 3)

    if complete_length == 0:
        return ""

    return str(
        Seq(cds[:complete_length]).translate(to_stop=False)
    )


def first_difference(seq1, seq2):
    """
    Return the 1-based amino-acid position of the first mismatch.
    """
    for i, (a, b) in enumerate(zip(seq1, seq2)):
        if a != b:
            return i + 1

    if len(seq1) != len(seq2):
        return min(len(seq1), len(seq2)) + 1

    return None


# --------------------------------------------------
# Load flagged transcript IDs
# --------------------------------------------------

flagged_ids = set()

with VALIDATION_FILE.open(newline="") as f:
    reader = csv.DictReader(f)

    for row in reader:
        if row["match"].strip().lower() == "false":
            flagged_ids.add(
                row["canonical_transcript_id"]
            )

print(f"Flagged transcripts found: {len(flagged_ids)}")


# --------------------------------------------------
# Analyze flagged transcripts
# --------------------------------------------------

results = []

with GENES_FILE.open(newline="") as f:
    reader = csv.DictReader(f)

    for row in reader:

        transcript_id = row["canonical_transcript_id"]

        if transcript_id not in flagged_ids:
            continue

        cds = str(
            row["cds_sequence"]
        ).strip().upper().replace("U", "T")

        reference_raw = normalize_protein(
            row["protein_sequence"]
        )

        translated_raw = translate_cds(cds)

        # Remove only normal terminal stop
        reference = remove_terminal_stop(
            reference_raw
        )

        translated = remove_terminal_stop(
            translated_raw
        )

        # --------------------------------------------------
        # Detect internal stop codons
        # --------------------------------------------------

        first_internal_stop = None
        first_internal_stop_codon = ""

        if "*" in translated:

            first_internal_stop = (
                translated.index("*") + 1
            )

            nt_start = (
                first_internal_stop - 1
            ) * 3

            first_internal_stop_codon = cds[
                nt_start:nt_start + 3
            ]

        # --------------------------------------------------
        # Find first actual protein mismatch
        # --------------------------------------------------

        diff = first_difference(
            translated,
            reference
        )

        our_aa = ""
        reference_aa = ""
        codon = ""

        protein_context_ours = ""
        protein_context_reference = ""

        if diff is not None:

            index = diff - 1

            if index < len(translated):
                our_aa = translated[index]

            if index < len(reference):
                reference_aa = reference[index]

            nt_start = index * 3

            if nt_start + 3 <= len(cds):
                codon = cds[
                    nt_start:nt_start + 3
                ]

            context_start = max(
                0,
                index - 10
            )

            context_end = (
                index + 11
            )

            protein_context_ours = translated[
                context_start:context_end
            ]

            protein_context_reference = reference[
                context_start:context_end
            ]

        # --------------------------------------------------
        # Classification
        # --------------------------------------------------

        reasons = []

        if len(cds) % 3 != 0:
            reasons.append(
                "CDS length not divisible by 3"
            )

        if len(cds) >= 3 and not cds.startswith("ATG"):
            reasons.append(
                "Non-ATG start"
            )

        if first_internal_stop is not None:
            reasons.append(
                "Internal stop codon"
            )

        if len(translated) < len(reference):
            reasons.append(
                "Translated protein shorter"
            )

        elif len(translated) > len(reference):
            reasons.append(
                "Translated protein longer"
            )

        if (
            len(translated) == len(reference)
            and translated != reference
        ):
            reasons.append(
                "Same length, sequence differs"
            )

        if not reasons and translated != reference:
            reasons.append(
                "Other sequence mismatch"
            )

        results.append({
            "gene_id":
                row.get("gene_id", ""),

            "gene_name":
                row.get("gene_name", ""),

            "canonical_transcript_id":
                transcript_id,

            "cds_length_nt":
                len(cds),

            "cds_mod_3":
                len(cds) % 3,

            "start_codon":
                cds[:3],

            "our_translation_length":
                len(translated),

            "ensembl_protein_length":
                len(reference),

            "first_difference_aa":
                diff,

            "our_amino_acid":
                our_aa,

            "ensembl_amino_acid":
                reference_aa,

            "codon_at_difference":
                codon,

            "first_internal_stop_aa":
                first_internal_stop,

            "first_internal_stop_codon":
                first_internal_stop_codon,

            "reason":
                "; ".join(reasons),

            "our_translation_context":
                protein_context_ours,

            "ensembl_context":
                protein_context_reference,

            "our_full_translation":
                translated,

            "ensembl_reference_protein":
                reference,

            "cds_sequence":
                cds
        })


# --------------------------------------------------
# Save detailed file
# --------------------------------------------------

fieldnames = [
    "gene_id",
    "gene_name",
    "canonical_transcript_id",
    "cds_length_nt",
    "cds_mod_3",
    "start_codon",
    "our_translation_length",
    "ensembl_protein_length",
    "first_difference_aa",
    "our_amino_acid",
    "ensembl_amino_acid",
    "codon_at_difference",
    "first_internal_stop_aa",
    "first_internal_stop_codon",
    "reason",
    "our_translation_context",
    "ensembl_context",
    "our_full_translation",
    "ensembl_reference_protein",
    "cds_sequence"
]

with OUTPUT_FILE.open(
    "w",
    newline=""
) as f:

    writer = csv.DictWriter(
        f,
        fieldnames=fieldnames
    )

    writer.writeheader()
    writer.writerows(results)


print(
    f"Detailed flagged transcripts: "
    f"{len(results)}"
)

print(
    f"Saved: {OUTPUT_FILE}"
)


# --------------------------------------------------
# Print examples
# --------------------------------------------------

print("\nExamples:\n")

for row in results[:10]:

    print(
        row["gene_name"],
        row["canonical_transcript_id"]
    )

    print(
        "Reason:",
        row["reason"]
    )

    print(
        "Lengths:",
        row["our_translation_length"],
        "vs",
        row["ensembl_protein_length"]
    )

    print(
        "First difference:",
        row["first_difference_aa"]
    )

    print(
        "Codon:",
        row["codon_at_difference"]
    )

    print(
        "Ours:   ",
        row["our_translation_context"]
    )

    print(
        "Ensembl:",
        row["ensembl_context"]
    )

    print(
        "Internal stop:",
        row["first_internal_stop_aa"],
        row["first_internal_stop_codon"]
    )

    print("-" * 70)