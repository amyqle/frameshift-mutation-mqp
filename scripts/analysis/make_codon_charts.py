import matplotlib.pyplot as plt

# Start codon counts
start_counts = {
    "ATG": 19894,
    "Other": 36
}

# Stop codon counts
stop_counts = {
    "TGA": 9900,
    "TAA": 5633,
    "TAG": 4465,
    "Other": 2
}


def add_labels(bars):
    for bar in bars:
        height = bar.get_height()

        plt.text(
            bar.get_x() + bar.get_width() / 2,
            height,
            f"{int(height):,}",
            ha="center",
            va="bottom"
        )


# Start codon chart
plt.figure(figsize=(7, 5))

bars = plt.bar(
    start_counts.keys(),
    start_counts.values()
)

add_labels(bars)

plt.xlabel("Start Codon")
plt.ylabel("Number of Transcripts")
plt.title("Start Codon Distribution")
plt.tight_layout()

plt.savefig(
    "figures/start_codon_distribution_labeled.png",
    dpi=300
)

plt.close()


# Stop codon chart
plt.figure(figsize=(7, 5))

bars = plt.bar(
    stop_counts.keys(),
    stop_counts.values()
)

add_labels(bars)

plt.xlabel("Stop Codon")
plt.ylabel("Number of Transcripts")
plt.title("Stop Codon Distribution")
plt.tight_layout()

plt.savefig(
    "figures/stop_codon_distribution_labeled.png",
    dpi=300
)

plt.close()

print("Saved codon charts to the figures folder.")