# Frameshift Mutation MQP

This project builds a reference dataset of human protein-coding genes from Ensembl GRCh38.p14 for frameshift mutation analysis.

## Project Goal

The current pipeline:

1. Identifies protein-coding genes
2. Matches each gene to its canonical transcript
3. Adds cDNA sequences
4. Adds protein sequences
5. Adds genomic sequences
6. Extracts CDS sequences
7. Calculates sequence-length statistics
8. Analyzes start and stop codons

## Data Source

Human genome data were downloaded from Ensembl using:

- Assembly: GRCh38.p14
- Genome accession: GCA_000001405.29

Required files:

- `genes.gtf`
- `cdna.fa.bgz`
- `pep.fa.bgz`
- `unmasked.fa.bgz`

These files are too large to store in this GitHub repository.

## Folder Structure

```text
frameshift-mutation-mqp/
├── scripts/
├── figures/
├── data/
├── CODING_BOOK.md
├── requirements.txt
└── .gitignore

The data/ folder is stored locally and is not uploaded to GitHub.
Setup
1. Clone the repository
git clone https://github.com/amyqle/frameshift-mutation-mqp.git
cd frameshift-mutation-mqp

2. Install required Python libraries
python3 -m pip install -r requirements.txt

3. Create a data folder
mkdir data

Download the following Ensembl files and place them inside the data/ folder:
- genes.gtf
- cdna.fa.bgz
- pep.fa.bgz
- unmasked.fa.bgz

Run Order
Run the scripts from the main project folder in this order:
python3 scripts/make_gene_csv.py
python3 scripts/add_canonical_transcripts.py
python3 scripts/add_cdna_sequences.py
python3 scripts/add_protein_sequences.py
python3 scripts/add_genomic_sequences.py
python3 scripts/extract_cds_sequences.py

For analysis:
python3 scripts/basic_statistics.py
python3 scripts/make_histograms.py
python3 scripts/analyze_codons.py
python3 scripts/make_codon_charts.py

filter_genes.py was used as an initial filtering test and is not required for the main pipeline.

Current Results

The pipeline produced:
- 20,131 protein-coding gene records
- 20,131 canonical transcript matches
- 20,131 cDNA sequence matches
- 20,131 protein sequence matches
- 20,131 genomic sequence matches
- 20,131 CDS sequence matches

The current analysis includes:
- Genomic sequence length distribution
- cDNA sequence length distribution
- Protein sequence length distribution
- Basic sequence-length statistics
- Start codon distribution
- Stop codon distribution

Documentation
See CODING_BOOK.md for an explanation of what each Python script does.

Next Steps
- Validate CDS-to-protein translation
- Simulate frameshift mutations
- Translate mutated CDS sequences
- Compare reference and mutated proteins
