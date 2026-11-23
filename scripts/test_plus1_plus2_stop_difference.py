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

# Each wild-type comparison uses the matched positions in its own CSV.
wildtype_results = pd.read_csv("data/frameshift_vs_wildtype.csv")
wildtype_column = "wildtype_stop_position"
p_values = {"+1 vs. +2": p_value}

for frame, stop_column in (
    ("+1", "plus1_first_stop_position"),
    ("+2", "plus2_first_stop_position"),
):
    label = f"{frame} vs. wild type"
    matched = wildtype_results.dropna(subset=[wildtype_column, stop_column])
    differences = matched[stop_column] - matched[wildtype_column]

    print(f"\n{label}")
    print("Transcripts included:", len(matched))
    print("Earlier stops:", (differences < 0).sum())
    print("Same stop position:", (differences == 0).sum())
    print("Later stops:", (differences > 0).sum())
    extensions = differences[differences > 0]
    if not extensions.empty:
        print(f"Average extension beyond wild type: {extensions.mean():.2f} aa")

    if differences.empty or not differences.ne(0).any():
        p_values[label] = None
        continue

    result = wilcoxon(differences, alternative="two-sided")
    p_values[label] = result.pvalue
    print("Wilcoxon statistic:", result.statistic)

print("\nP-value summary")
for label, comparison_p in p_values.items():
    if comparison_p is None:
        print(f"{label}: no nonzero paired differences to test")
    elif comparison_p < 0.05:
        print(f"{label}: p < 0.05 (statistically significant)")
    elif comparison_p > 0.05:
        print(f"{label}: p > 0.05 (not statistically significant)")
    else:
        print(f"{label}: p = 0.05 (does not meet p < 0.05)")