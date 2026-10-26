"""Build the reference CSV files using the documented pipeline."""
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
STEPS = [
    "make_gene_csv.py",
    "add_canonical_transcripts.py",
    "add_cdna_sequences.py",
    "add_protein_sequences.py",
    "add_genomic_sequences.py",
    "extract_cds_sequences.py",
]

def main():
    required = ["genes.gtf", "cdna.fa.bgz", "pep.fa.bgz", "unmasked.fa.bgz"]
    missing = [name for name in required if not (ROOT / "data" / name).is_file()]
    if missing:
        raise SystemExit("Missing files in data/: " + ", ".join(missing))
    for index, step in enumerate(STEPS, 1):
        print(f"\n[{index}/{len(STEPS)}] {step}", flush=True)
        subprocess.run([sys.executable, "-u", str(ROOT / "scripts" / step)],
                       cwd=ROOT, check=True)
    print("\nComplete: data/final_gene_dataset.csv and data/genes_with_cds.csv", flush=True)

if __name__ == "__main__":
    main()
