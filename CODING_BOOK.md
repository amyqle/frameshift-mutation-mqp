# Frameshift Mutation Project Coding Book

## filter_genes.py

**Purpose:**  
Finds protein-coding genes in the Ensembl GTF file.

**Input:**  
`genes.gtf`

**Output:**  
`protein_coding_genes.txt`

**What it does:**  
- Reads the GTF file
- Looks for rows where `gene_biotype` is `protein_coding`
- Collects the gene names
- Saves one gene name per line

**Result:**  
19,452 unique protein-coding gene names

**Why we used it:**  
This was our first test to confirm we could successfully filter the Ensembl annotation file.

---

## 1. make_gene_csv.py

**Purpose:**  
Creates a CSV of protein-coding genes from the Ensembl GTF file.

**Input:**  
`genes.gtf`

**Output:**  
`protein_coding_genes.csv`

**What it does:**  
- Filters for protein-coding genes
- Extracts gene ID, gene name, chromosome, start, end, and strand

**Result:**  
20,131 protein-coding gene records

---

## 2. add_canonical_transcripts.py

**Purpose:**  
Adds the canonical transcript for each protein-coding gene.

**Inputs:**  
- `genes.gtf`
- `protein_coding_genes.csv`

**Output:**  
`genes_with_canonical_transcripts.csv`

**Result:**  
20,131 canonical transcripts matched

---

## 3. add_cdna_sequences.py

**Purpose:**  
Matches each canonical transcript to its cDNA sequence.

**Inputs:**  
- `genes_with_canonical_transcripts.csv`
- `cdna.fa.bgz`

**Output:**  
`genes_with_cdna.csv`

**Result:**  
20,131 cDNA sequences matched

---

## 4. add_protein_sequences.py

**Purpose:**  
Matches each canonical transcript to its protein sequence.

**Inputs:**  
- `genes_with_cdna.csv`
- `pep.fa.bgz`

**Output:**  
`genes_with_cdna_and_protein.csv`

**Result:**  
20,131 protein sequences matched

---

## 5. add_genomic_sequences.py

**Purpose:**  
Extracts the genomic DNA sequence for each gene.

**Inputs:**  
- `genes_with_cdna_and_protein.csv`
- `unmasked.fa.bgz`

**Output:**  
`final_gene_dataset.csv`

**Result:**  
20,131 genomic sequences matched

---

## 6. make_histograms.py

**Purpose:**  
Creates histograms of:
- genomic sequence length
- cDNA sequence length
- protein sequence length

**Input:**  
`final_gene_dataset.csv`

**Outputs:**  
- `genomic_length_histogram_clean.png`
- `cdna_length_histogram_clean.png`
- `protein_length_histogram_clean.png`

---

## 7. basic_statistics.py

**Purpose:**  
Calculates basic statistics for sequence lengths.

**Statistics:**  
- Count
- Mean
- Median
- Standard deviation
- Minimum
- Maximum

**Input:**  
`final_gene_dataset.csv`

---

## 8. extract_cds_sequences.py

**Purpose:**  
Extracts the coding sequence, or CDS, for each canonical transcript.

**Inputs:**  
- `genes.gtf`
- `final_gene_dataset.csv`

**Output:**  
`genes_with_cds.csv`

**Result:**  
20,131 CDS sequences matched

---

## 9. analyze_codons.py

**Purpose:**  
Analyzes annotated start and stop codons from the GTF file.

**Inputs:**  
- `genes.gtf`
- `final_gene_dataset.csv`

**Start codon results:**  
- ATG: 19,894
- Other: 36

**Stop codon results:**  
- TGA: 9,900
- TAA: 5,633
- TAG: 4,465
- Other: 2

---

## 10. make_codon_charts.py

**Purpose:**  
Creates charts showing the start and stop codon distributions.

**Outputs:**  
- `start_codon_distribution_labeled.png`
- `stop_codon_distribution_labeled.png`

---

# Pipeline Order

1. `make_gene_csv.py`
2. `add_canonical_transcripts.py`
3. `add_cdna_sequences.py`
4. `add_protein_sequences.py`
5. `add_genomic_sequences.py`
6. `make_histograms.py`
7. `basic_statistics.py`
8. `extract_cds_sequences.py`
9. `analyze_codons.py`
10. `make_codon_charts.py`

# Current Status

Completed:
- Reference gene dataset
- Canonical transcript matching
- cDNA matching
- Protein matching
- Genomic sequence extraction
- Sequence-length analysis
- CDS extraction
- Start/stop codon analysis

Next:
- Validate CDS translation
- Simulate frameshift mutations
- Compare reference and mutated proteins