import pandas as pd
from scipy.stats import wilcoxon

INPUT = "data/frameshift_stop_comparison.csv"

df = pd.read_csv(INPUT)

# Keep only transcripts where BOTH +1 and +2 found a stop codon
paired = df.dropna(
    subset=[
        "plus1_stop_codons_after_shift",
        "plus2_stop_codons_after_shift"
    ]
).copy()

plus1 = paired["plus1_stop_codons_after_shift"]
plus2 = paired["plus2_stop_codons_after_shift"]

print("Transcripts included:", len(paired))

print("\n+1:")
print("Mean:", plus1.mean())
print("Median:", plus1.median())

print("\n+2:")
print("Mean:", plus2.mean())
print("Median:", plus2.median())

# Wilcoxon signed-rank test
statistic, p_value = wilcoxon(
    plus1,
    plus2,
    alternative="two-sided"
)

print("\nWilcoxon signed-rank test")
print("Statistic:", statistic)
print("P-value:", p_value)

if p_value < 0.05:
    print(
        "\nResult: The difference between +1 and +2 "
        "first-stop positions is statistically significant."
    )
else:
    print(
        "\nResult: The difference between +1 and +2 "
        "first-stop positions is not statistically significant."
    )

# Direction
plus1_longer = (plus1 > plus2).sum()
plus2_longer = (plus2 > plus1).sum()
same = (plus1 == plus2).sum()

print("\nDirection:")
print("+1 survives longer:", plus1_longer)
print("+2 survives longer:", plus2_longer)
print("Same stop position:", same)

difference = plus2 - plus1

print("\nDifference (+2 - +1):")
print("Mean difference:", difference.mean())
print("Median difference:", difference.median())