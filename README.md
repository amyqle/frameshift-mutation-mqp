# Frameshift Mutation MQP

We are building a reference dataset of human protein-coding genes for frameshift mutation analysis. The current dataset has 20,131 genes with canonical transcripts, cDNA, protein, genomic and CDS sequences.

## Data source

We used Ensembl GRCh38.p14, assembly accession **GCA_000001405.29**, from the **2026_04** dataset.

- [Gene annotations, cDNA and protein downloads](https://ftp.ebi.ac.uk/pub/ensemblorganisms/GCA/000/001/405/29/ensembl/2026_04/geneset/)
- [Genome download](https://ftp.ebi.ac.uk/pub/ensemblorganisms/GCA/000/001/405/29/ensembl/2026_04/genome/)

## Setup

1. Download or clone this repository and open `frameshift-mutation-mqp` in VS Code. We used Python 3.12.
2. Create a `data` folder inside the project.
3. Download the files from the sources above and put them in `data`:

```text
data/
  genes.gtf
  cdna.fa.bgz
  pep.fa.bgz
  unmasked.fa.bgz
```

The data files are stored locally and are not included in this repository.

## Project folders

```text
scripts/
  run_pipeline.py  # Build the reference CSVs
  dataset/         # Prepare genes and sequences
  validation/      # Check CDS translations and mismatches
  models/          # Build model inputs and simulate frameshifts
  analysis/        # Statistics and graphs
  experiments/     # Initial filtering and deletion/position experiments
data/              # Local input files and generated CSVs
figures/           # Saved graphs
```

Run commands from the project folder, not from inside a script subfolder.

## Build the CSV files

Run this from the project folder:

```powershell
python scripts/run_pipeline.py
```

The runner runs the six dataset scripts in order. Running it again replaces the generated CSV files.

Main outputs:

- `data/final_gene_dataset.csv`: reference genes with cDNA, protein and genomic sequences.
- `data/genes_with_cds.csv`: the reference dataset with coding sequences added.

## Statistics and graphs

Install the analysis libraries, then run:

```powershell
python -m pip install -r requirements.txt
python scripts/analysis/basic_statistics.py
python scripts/analysis/make_histograms.py
python scripts/analysis/analyze_codons.py
python scripts/analysis/make_codon_charts.py
```

## Validation and simple model

After building the reference CSVs, run all four steps with:

```powershell
python scripts/run_pipeline.py --simple-model
```

This validates CDS translations, builds the model inputs, runs the insertion model and prints the results summary. It stops if a step fails. Results are saved in `data/`; rerunning replaces these generated files.

See [CODING_BOOK.md](CODING_BOOK.md) for the script notes.
