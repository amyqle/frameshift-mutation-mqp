"""Build the reference CSV files using the documented pipeline."""
import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
STEPS = [
    "dataset/make_gene_csv.py",
    "dataset/add_canonical_transcripts.py",
    "dataset/add_cdna_sequences.py",
    "dataset/add_protein_sequences.py",
    "dataset/add_genomic_sequences.py",
    "dataset/extract_cds_sequences.py",
]
SIMPLE_MODEL_STEPS = [
    "validation/validate_cds_translation.py",
    "models/build_simple_model_sequences.py",
    "models/simple_frameshift_insertion_model.py",
    "analysis/analyze_simple_frameshift_results.py",
]


def main(argv=None):
    parser = argparse.ArgumentParser(description="Run the MQP dataset or simple-model steps.")
    parser.add_argument("--simple-model", action="store_true",
                        help="Run validation, build model inputs, run the insertion model and summarize results.")
    args = parser.parse_args(argv)
    steps = SIMPLE_MODEL_STEPS if args.simple_model else STEPS
    required = (["genes_with_cds.csv"] if args.simple_model else
                ["genes.gtf", "cdna.fa.bgz", "pep.fa.bgz", "unmasked.fa.bgz"])
    missing = [name for name in required if not (ROOT / "data" / name).is_file()]
    if missing:
        raise SystemExit("Missing files in data/: " + ", ".join(missing))
    for index, step in enumerate(steps, 1):
        print(f"\n[{index}/{len(steps)}] {step}", flush=True)
        subprocess.run([sys.executable, "-u", str(ROOT / "scripts" / step)],
                       cwd=ROOT, check=True)
    if args.simple_model:
        print("\nComplete: data/simple_frameshift_insertion_results.csv. Summary printed above.", flush=True)
    else:
        print("\nComplete: data/final_gene_dataset.csv and data/genes_with_cds.csv", flush=True)


if __name__ == "__main__":
    main()
