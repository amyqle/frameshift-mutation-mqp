import pandas as pd

df = pd.read_csv("data/final_gene_dataset.csv")

df["genomic_length"] = df["genomic_sequence"].str.len()
df["cdna_length"] = df["cdna_sequence"].str.len()
df["protein_length"] = df["protein_sequence"].str.len()

for column, label in [
    ("genomic_length", "Genomic Sequence"),
    ("cdna_length", "cDNA Sequence"),
    ("protein_length", "Protein Sequence")
]:
    print(f"\n{label}")
    print(f"Count: {df[column].count()}")
    print(f"Mean: {df[column].mean():,.2f}")
    print(f"Median: {df[column].median():,.2f}")
    print(f"Standard Deviation: {df[column].std():,.2f}")
    print(f"Minimum: {df[column].min():,.0f}")
    print(f"Maximum: {df[column].max():,.0f}")