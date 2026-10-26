import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("data/final_gene_dataset.csv")

# Calculate lengths
df["genomic_length"] = df["genomic_sequence"].str.len()
df["cdna_length"] = df["cdna_sequence"].str.len()
df["protein_length"] = df["protein_sequence"].str.len()


def make_histogram(data, title, xlabel, filename):
    cutoff = data.quantile(0.99)
    median = data.median()

    filtered = data[data <= cutoff]

    plt.figure(figsize=(8, 5))
    plt.hist(filtered, bins=50)

    plt.axvline(
        median,
        linestyle="--",
        linewidth=2,
        label=f"Median = {median:,.0f}"
    )

    plt.xlabel(xlabel)
    plt.ylabel("Number of Genes")
    plt.title(title)
    plt.legend()
    plt.tight_layout()

    plt.savefig(filename, dpi=300)
    plt.close()


make_histogram(
    df["genomic_length"],
    "Distribution of Genomic Sequence Lengths",
    "Genomic Sequence Length (bp)",
    "figures/genomic_length_histogram_clean.png"
)

make_histogram(
    df["cdna_length"],
    "Distribution of cDNA Sequence Lengths",
    "cDNA Sequence Length (bp)",
    "figures/cdna_length_histogram_clean.png"
)

make_histogram(
    df["protein_length"],
    "Distribution of Protein Sequence Lengths",
    "Protein Sequence Length (amino acids)",
    "figures/protein_length_histogram_clean.png"
)

print("Saved 3 histograms to the figures folder.")