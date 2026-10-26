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
python scripts/basic_statistics.py
python scripts/make_histograms.py
python scripts/analyze_codons.py
python scripts/make_codon_charts.py
```
