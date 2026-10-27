import argparse
import csv
import gzip
import itertools
import re
import tempfile
from collections import Counter, defaultdict
from contextlib import ExitStack
from pathlib import Path

from Bio.Data import CodonTable
from Bio.Seq import Seq


ROOT = Path(__file__).resolve().parents[2]
TABLE = CodonTable.unambiguous_dna_by_id[1]
INSERTIONS = ["A", "T", "C", "G"] + [
    "".join(bases) for bases in itertools.product("ATCG", repeat=2)
]
csv.field_size_limit(10_000_000)


def normalize_id(value):
    return str(value or "").split(".")[0]


def translate(sequence):
    codons = []
    stop_nt = None
    for position in range(0, len(sequence) - 2, 3):
        codon = sequence[position:position + 3].upper()
        if codon in TABLE.stop_codons:
            stop_nt = position + 1
            break
        codons.append(codon if codon in TABLE.forward_table else "NNN")
    protein = str(Seq("".join(codons)).translate(table=1))
    stop_aa = None if stop_nt is None else (stop_nt - 1) // 3 + 1
    return protein, stop_aa, stop_nt


def first_difference(reference, changed):
    for position, (left, right) in enumerate(zip(reference, changed), 1):
        if left != right:
            return position
    return min(len(reference), len(changed)) + 1 if len(reference) != len(changed) else None


def measure(sequence, reference, cds_boundary):
    protein, stop_aa, stop_nt = translate(sequence)
    return {
        "mutated_protein_sequence": protein,
        "mutated_protein_length": len(protein),
        "first_difference_position": first_difference(reference, protein),
        "stop_aa_position": stop_aa,
        "stop_nt_position": stop_nt,
        "premature_stop": stop_aa is not None and stop_aa <= len(reference),
        "entered_downstream_region": (
            stop_nt > cds_boundary if stop_nt is not None else len(sequence) > cds_boundary
        ),
    }


def load_inputs(data, limit):
    with (data / "cds_translation_validation.csv").open(newline="", encoding="utf-8-sig") as file:
        validated = {
            normalize_id(row["canonical_transcript_id"])
            for row in csv.DictReader(file) if row["match"].lower() == "true"
        }
    inputs = {}
    with (data / "simple_model_sequences.csv").open(newline="", encoding="utf-8-sig") as file:
        for row in csv.DictReader(file):
            transcript = normalize_id(row["canonical_transcript_id"])
            if transcript not in validated:
                continue
            if transcript in inputs:
                raise ValueError(f"Duplicate model transcript: {transcript}")
            inputs[transcript] = row
            if limit is not None and len(inputs) >= limit:
                break
    if not inputs:
        raise ValueError("No validated simple-model inputs were found.")
    return inputs


def read_gtf(path, selected):
    with path.open("rb") as file:
        opener = gzip.open if file.read(2) == b"\x1f\x8b" else open
    transcripts = {}
    with opener(path, "rt", encoding="utf-8") as file:
        for line in file:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] not in {"exon", "CDS"}:
                continue
            match = re.search(r'transcript_id "([^"]+)"', fields[8])
            if match is None:
                continue
            transcript = normalize_id(match.group(1))
            if transcript not in selected:
                continue
            info = transcripts.setdefault(transcript, {
                "chromosome": fields[0], "strand": fields[6], "exons": set(), "cds": set(),
            })
            if (fields[0], fields[6]) != (info["chromosome"], info["strand"]):
                raise ValueError(f"Inconsistent annotation for {transcript}")
            info["exons" if fields[2] == "exon" else "cds"].add((int(fields[3]), int(fields[4])))
    for info in transcripts.values():
        for name in ("exons", "cds"):
            info[name] = sorted(info[name], reverse=info["strand"] == "-")
    return transcripts


def reconstruct_transcript(row, model, info):
    exons, cds = info["exons"], info["cds"]
    if not exons or not cds:
        raise ValueError("missing_exon_or_cds_annotation")
    if row["strand"] != info["strand"] or row["chromosome"].removeprefix("chr") != info["chromosome"].removeprefix("chr"):
        raise ValueError("annotation_location_mismatch")
    gene_start, gene_end = int(row["start"]), int(row["end"])
    genomic = row["genomic_sequence"].upper()
    if len(genomic) != gene_end - gene_start + 1:
        raise ValueError("genomic_length_mismatch")

    def extract(start, end):
        if not gene_start <= start <= end <= gene_end:
            raise ValueError("annotation_outside_gene")
        if info["strand"] == "+":
            return genomic[start - gene_start:end - gene_start + 1]
        return genomic[gene_end - end:gene_end - start + 1]

    def transcript_offset(coordinate):
        offset = 0
        for start, end in exons:
            if start <= coordinate <= end:
                return offset + (coordinate - start if info["strand"] == "+" else end - coordinate)
            offset += end - start + 1
        raise ValueError("cds_boundary_outside_exons")

    cdna = "".join(extract(start, end) for start, end in exons)
    original_cds = "".join(extract(start, end) for start, end in cds)
    if cdna != row["cdna_sequence"].upper():
        raise ValueError("reconstructed_cdna_mismatch")
    if original_cds != model["cds_sequence"].upper() or original_cds != row["cds_sequence"].upper():
        raise ValueError("reconstructed_cds_mismatch")
    start_coordinate = cds[0][0] if info["strand"] == "+" else cds[0][1]
    end_coordinate = cds[-1][1] if info["strand"] == "+" else cds[-1][0]
    cds_start = transcript_offset(start_coordinate)
    cds_end = transcript_offset(end_coordinate) + 1
    if cds_start != int(model["cds_start_in_cdna"]) - 1 or cds_end - cds_start != len(original_cds):
        raise ValueError("cds_boundary_mismatch")
    simple = cdna[cds_start:]
    if simple != model["simple_model_sequence"].upper() or len(original_cds) != int(model["cds_length"]):
        raise ValueError("simple_model_input_mismatch")
    reference = model["protein_sequence"]
    if reference != row["protein_sequence"] or translate(simple)[0] != reference:
        raise ValueError("reference_translation_mismatch")

    first_intron = ""
    location = "none"
    point = None
    if len(exons) > 1:
        left, right = exons[0], exons[1]
        start, end = ((left[1] + 1, right[0] - 1) if info["strand"] == "+"
                      else (right[1] + 1, left[0] - 1))
        if start > end:
            raise ValueError("overlapping_or_adjacent_first_exons")
        first_intron = extract(start, end)
        point = left[1] - left[0] + 1 - cds_start
        location = "5prime_UTR" if point <= 0 else ("CDS" if point < len(original_cds) else "3prime_UTR")

    included = len(first_intron) if point is not None and point > 0 else 0
    retained = simple[:point] + first_intron + simple[point:] if included else simple
    mutation_position = len(original_cds) // 2
    # Keep the insertion at the same coding position.
    retained_position = mutation_position + (included if included and point <= mutation_position else 0)
    retained_boundary = len(original_cds) + (included if location == "CDS" else 0)
    return {
        "simple": simple, "retained": retained, "reference": reference,
        "cds_length": len(original_cds), "retained_cds_length": retained_boundary,
        "mutation_position": mutation_position, "retained_position": retained_position,
        "first_intron_length": len(first_intron), "retained_bases_in_model": included,
        "first_intron_location": location,
        "cds_has_partial_final_codon": len(original_cds) % 3 != 0,
        "retained_input_has_ambiguous_bases": any(base not in "ACGT" for base in retained),
    }


def simulate(info, inserted, baseline=None):
    baseline = baseline or measure(info["retained"], info["reference"], info["retained_cds_length"])
    position, mapped = info["mutation_position"], info["retained_position"]
    size = len(inserted)
    simple = measure(info["simple"][:position] + inserted + info["simple"][position:],
                     info["reference"], info["cds_length"] + size)
    retained = measure(info["retained"][:mapped] + inserted + info["retained"][mapped:],
                       info["reference"], info["retained_cds_length"] + size)
    stop = baseline["stop_nt_position"]
    before_mutation = stop is not None and stop - 1 + 3 <= mapped
    return {
        "frameshift": f"+{size}", "inserted_sequence": inserted,
        "mutation_base_position": position + 1, "retained_mutation_base_position": mapped + 1,
        **{f"simple_{name}": value for name, value in simple.items()}, **retained,
        "baseline_retained_protein_length": baseline["mutated_protein_length"],
        "baseline_retained_stop_nt_position": stop,
        "baseline_retained_premature_stop": baseline["premature_stop"],
        "baseline_stop_before_mutation": before_mutation,
        "protein_changed_from_retained_baseline": retained["mutated_protein_sequence"] != baseline["mutated_protein_sequence"],
        "protein_length_change_from_retained_baseline": retained["mutated_protein_length"] - baseline["mutated_protein_length"],
    }


META_FIELDS = ["gene_id", "gene_name", "transcript_id", "strand", "first_intron_length",
               "retained_bases_in_model", "first_intron_location", "reference_protein_length",
               "cds_length", "retained_cds_length", "cds_has_partial_final_codon",
               "retained_input_has_ambiguous_bases"]
METRIC_FIELDS = ["mutated_protein_sequence", "mutated_protein_length", "first_difference_position",
                 "stop_aa_position", "stop_nt_position", "premature_stop", "entered_downstream_region"]
RESULT_FIELDS = META_FIELDS + ["frameshift", "inserted_sequence", "mutation_base_position",
    "retained_mutation_base_position"] + [f"simple_{name}" for name in METRIC_FIELDS] + METRIC_FIELDS + [
    "baseline_retained_protein_length", "baseline_retained_stop_nt_position", "baseline_retained_premature_stop",
    "baseline_stop_before_mutation", "protein_changed_from_retained_baseline", "protein_length_change_from_retained_baseline"]
CONTROL_FIELDS = META_FIELDS + ["simple_control_protein_sequence", "simple_control_protein_length",
    "simple_control_stop_nt_position", "retained_control_protein_sequence", "retained_control_protein_length",
    "retained_control_stop_nt_position", "retained_control_premature_stop"]
SUMMARY_FIELDS = ["scope", "frameshift", "cases", "simple_premature_stops", "retained_premature_stops",
    "simple_downstream", "retained_downstream", "simple_no_stop", "retained_no_stop",
    "baseline_stop_before_mutation", "protein_changed_from_retained_baseline", "simple_length_sum", "retained_length_sum"]


def main(argv=None):
    parser = argparse.ArgumentParser(description="Compare first-intron retention with the simple insertion model.")
    parser.add_argument("--data-dir", type=Path, default=ROOT / "data")
    parser.add_argument("--limit", type=int)
    parser.add_argument("--output-prefix", default="first_intron_retained")
    args = parser.parse_args(argv)
    if args.limit is not None and args.limit < 1:
        parser.error("--limit must be positive")
    if not re.fullmatch(r"[A-Za-z0-9_-]+", args.output_prefix):
        parser.error("--output-prefix must be a filename prefix")
    data = args.data_dir.resolve()
    for name in ["genes_with_cds.csv", "cds_translation_validation.csv", "simple_model_sequences.csv", "genes.gtf"]:
        if not (data / name).is_file():
            raise FileNotFoundError(f"Missing {data / name}. Run the reference, validation and model-input steps first.")
    inputs = load_inputs(data, args.limit)
    print(f"Selected validated inputs: {len(inputs)}", flush=True)
    print("Reading transcript annotations...", flush=True)
    annotations = read_gtf(data / "genes.gtf", inputs)
    groups = defaultdict(Counter)
    processed, skipped, cases = 0, 0, 0
    seen = set()
    staged = []
    try:
        with ExitStack() as stack:
            writers = {}
            schemas = {"results": RESULT_FIELDS, "controls": CONTROL_FIELDS,
                       "skipped": ["transcript_id", "reason"], "summary": SUMMARY_FIELDS}
            for label, fields in schemas.items():
                file = stack.enter_context(tempfile.NamedTemporaryFile(
                    mode="w", newline="", encoding="utf-8", dir=data, suffix=".tmp", delete=False))
                staged.append((Path(file.name), data / f"{args.output_prefix}_{label}.csv"))
                writers[label] = csv.DictWriter(file, fieldnames=fields)
                writers[label].writeheader()
            with (data / "genes_with_cds.csv").open(newline="", encoding="utf-8-sig") as file:
                for row in csv.DictReader(file):
                    transcript = normalize_id(row["canonical_transcript_id"])
                    if transcript not in inputs:
                        continue
                    if transcript in seen:
                        raise ValueError(f"Duplicate gene transcript: {transcript}")
                    seen.add(transcript)
                    try:
                        if transcript not in annotations:
                            raise ValueError("missing_canonical_annotation")
                        info = reconstruct_transcript(row, inputs[transcript], annotations[transcript])
                    except ValueError as error:
                        writers["skipped"].writerow({"transcript_id": transcript, "reason": str(error)})
                        skipped += 1
                        continue
                    meta = {"gene_id": row["gene_id"], "gene_name": row["gene_name"],
                            "transcript_id": transcript, "strand": row["strand"],
                            "reference_protein_length": len(info["reference"]),
                            **{name: info[name] for name in META_FIELDS if name in info}}
                    simple_control = measure(info["simple"], info["reference"], info["cds_length"])
                    retained_control = measure(info["retained"], info["reference"], info["retained_cds_length"])
                    writers["controls"].writerow({**meta,
                        "simple_control_protein_sequence": simple_control["mutated_protein_sequence"],
                        "simple_control_protein_length": simple_control["mutated_protein_length"],
                        "simple_control_stop_nt_position": simple_control["stop_nt_position"],
                        "retained_control_protein_sequence": retained_control["mutated_protein_sequence"],
                        "retained_control_protein_length": retained_control["mutated_protein_length"],
                        "retained_control_stop_nt_position": retained_control["stop_nt_position"],
                        "retained_control_premature_stop": retained_control["premature_stop"]})
                    for inserted in INSERTIONS:
                        outcome = simulate(info, inserted, retained_control)
                        writers["results"].writerow({**meta, **outcome})
                        cases += 1
                        for scope in ["all", info["first_intron_location"]]:
                            group = groups[scope, outcome["frameshift"]]
                            group["cases"] += 1
                            for prefix in ["simple_", ""]:
                                label = "simple" if prefix else "retained"
                                group[f"{label}_premature_stops"] += outcome[prefix + "premature_stop"]
                                group[f"{label}_downstream"] += outcome[prefix + "entered_downstream_region"]
                                group[f"{label}_no_stop"] += outcome[prefix + "stop_nt_position"] is None
                                group[f"{label}_length_sum"] += outcome[prefix + "mutated_protein_length"]
                            for name in ["baseline_stop_before_mutation", "protein_changed_from_retained_baseline"]:
                                group[name] += outcome[name]
                    processed += 1
                    if processed % 1000 == 0:
                        print(f"Compared {processed} transcripts", flush=True)
            for transcript in inputs.keys() - seen:
                writers["skipped"].writerow({"transcript_id": transcript, "reason": "missing_gene_row"})
                skipped += 1
            if not processed:
                raise ValueError("No comparable transcripts; no result files replaced.")
            for (scope, frameshift), counts in sorted(groups.items()):
                writers["summary"].writerow({"scope": scope, "frameshift": frameshift,
                    **{field: counts[field] for field in SUMMARY_FIELDS if field not in {"scope", "frameshift"}}})
        for temporary, target in staged:
            temporary.replace(target)
    finally:
        for temporary, target in staged:
            if temporary.exists():
                temporary.unlink()
    print(f"Compared transcripts: {processed}; skipped: {skipped}; paired insertion cases: {cases}", flush=True)
    print(f"Saved {args.output_prefix}_results/controls/skipped/summary.csv in {data}", flush=True)
    return {"transcripts": processed, "skipped": skipped, "cases": cases}


if __name__ == "__main__":
    main()
