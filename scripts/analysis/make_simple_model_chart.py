import matplotlib.pyplot as plt


labels = ["+1 Frameshift", "+2 Frameshift"]

premature_stop = [96.79, 97.18]
downstream = [3.16, 2.79]


# --------------------------------------------------
# Chart 1: Premature Stop Percentage
# --------------------------------------------------

plt.figure(figsize=(7, 5))

bars = plt.bar(
    labels,
    premature_stop
)

plt.ylabel("Percentage of Simulations (%)")
plt.title("Premature Stop Rate by Frameshift Type")
plt.ylim(90, 100)

for bar, value in zip(bars, premature_stop):
    plt.text(
        bar.get_x() + bar.get_width() / 2,
        value + 0.15,
        f"{value:.2f}%",
        ha="center"
    )

plt.tight_layout()

output1 = "figures/simple_model_premature_stop.png"
plt.savefig(output1, dpi=300)

print(f"Saved to {output1}")

plt.close()


# --------------------------------------------------
# Chart 2: Downstream Region Percentage
# --------------------------------------------------

plt.figure(figsize=(7, 5))

bars = plt.bar(
    labels,
    downstream
)

plt.ylabel("Percentage of Simulations (%)")
plt.title("Frameshifts Extending Into Downstream Region")
plt.ylim(0, 4)

for bar, value in zip(bars, downstream):
    plt.text(
        bar.get_x() + bar.get_width() / 2,
        value + 0.08,
        f"{value:.2f}%",
        ha="center"
    )

plt.tight_layout()

output2 = "figures/simple_model_downstream_region.png"
plt.savefig(output2, dpi=300)

print(f"Saved to {output2}")

plt.close()