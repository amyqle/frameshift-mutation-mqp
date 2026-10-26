import matplotlib.pyplot as plt

positions = ["25%", "50%", "75%"]
avg_lengths = [164.4, 297.4, 401.5]

plt.figure(figsize=(7, 5))

bars = plt.bar(positions, avg_lengths)

for bar, value in zip(bars, avg_lengths):
    plt.text(
        bar.get_x() + bar.get_width() / 2,
        value,
        f"{value:.1f}",
        ha="center",
        va="bottom"
    )

plt.xlabel("Frameshift Position in CDS")
plt.ylabel("Average Mutated Protein Length (amino acids)")
plt.title("Effect of Frameshift Position on Protein Length")
plt.tight_layout()

plt.savefig(
    "figures/frameshift_position_protein_length.png",
    dpi=300
)

plt.close()

print("Saved frameshift position chart.")