import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
csv.field_size_limit(10_000_000)


def read_rows(path):
    with path.open(newline="", encoding="utf-8-sig") as file:
        return list(csv.DictReader(file))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, default=ROOT / "data")
    parser.add_argument("--figure-dir", type=Path, default=ROOT / "figures")
    args = parser.parse_args()
    args.figure_dir.mkdir(parents=True, exist_ok=True)
    controls = read_rows(args.data_dir / "first_intron_retained_controls.csv")
    summaries = {(r["scope"], r["frameshift"]): r for r in read_rows(args.data_dir / "first_intron_retained_summary.csv")}
    palette = ["#366C8D", "#B26D24"]
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11,
                         "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(1, 2, figsize=(11.8, 5.7), sharey=True)
    largest = 0
    for ax, scope, title in zip(axes, ["all", "CDS"], ["All matched transcripts", "First intron inside the CDS"]):
        cohort = controls if scope == "all" else [r for r in controls if r["first_intron_location"] == scope]
        if not cohort:
            raise ValueError(f"No transcripts in {scope}")
        simple = [sum(int(r["simple_control_protein_length"]) for r in cohort) / len(cohort)]
        retained = [sum(int(r["retained_control_protein_length"]) for r in cohort) / len(cohort)]
        for shift in ["+1", "+2"]:
            summary = summaries[scope, shift]
            simple.append(int(summary["simple_length_sum"]) / int(summary["cases"]))
            retained.append(int(summary["retained_length_sum"]) / int(summary["cases"]))
        x = np.arange(3)
        for delta, values, color, label in [(-0.18, simple, palette[0], "Simple sequence"),
                                             (0.18, retained, palette[1], "First intron retained")]:
            bars = ax.bar(x + delta, values, width=0.34, color=color, label=label)
            ax.bar_label(bars, labels=[f"{v:.0f}" for v in values], padding=5)
        ax.set_xticks(x, ["No insertion", "+1 base", "+2 bases"])
        ax.set_title(f"{title}\n{len(cohort):,} transcripts", fontsize=12, pad=15)
        largest = max(largest, *simple, *retained)
        ax.set_axisbelow(True)
        ax.yaxis.grid(True, color="#E3E7EA", linewidth=0.7)
    axes[0].set_ylim(0, max(largest, 1) * 1.22)
    axes[0].set_ylabel("Mean translated length (amino acids)")
    axes[1].legend(frameon=False, loc="upper right")
    fig.suptitle("First-intron retention with matched insertions", fontsize=16, y=0.98)
    fig.text(0.5, 0.91, "Same transcripts, inserted bases and original coding positions", ha="center")
    summary = summaries["CDS", "+1"]
    rate = int(summary["baseline_stop_before_mutation"]) / int(summary["cases"]) * 100
    fig.text(0.08, 0.08, f"In the CDS-retention group, {rate:.1f}% already stop before the insertion site.", fontsize=10)
    fig.text(0.08, 0.035, "Standard-code sequence predictions; no-stop cases end at the available sequence boundary. Protein function is not tested.", fontsize=9)
    fig.subplots_adjust(left=0.08, right=0.98, top=0.77, bottom=0.22, wspace=0.12)
    fig.savefig(args.figure_dir / "first_intron_comparison.png", dpi=180)
    plt.close(fig)

    by_id = {row["transcript_id"]: row for row in controls}
    example = None
    with (args.data_dir / "first_intron_retained_results.csv").open(newline="", encoding="utf-8-sig") as file:
        for row in csv.DictReader(file):
            if (row["inserted_sequence"] == "A" and row["first_intron_location"] == "CDS"
                    and row["baseline_stop_before_mutation"] == "True"
                    and 100 <= int(row["reference_protein_length"]) <= 800
                    and row["simple_stop_nt_position"]
                    and by_id[row["transcript_id"]]["simple_control_stop_nt_position"]):
                example = row
                break
    if example is not None:
        control = by_id[example["transcript_id"]]
        values = [int(control["simple_control_protein_length"]), int(example["simple_mutated_protein_length"]),
                  int(control["retained_control_protein_length"]), int(example["mutated_protein_length"])]
        labels = ["Simple control", "Simple +A", "Retained control", "Retained +A"]
        fig, ax = plt.subplots(figsize=(9.6, 5.0))
        bars = ax.barh(labels, values, color=[palette[0], palette[0], palette[1], palette[1]], height=0.6)
        ax.bar_label(bars, padding=5)
        ax.invert_yaxis()
        ax.set_xlim(0, max(values) * 1.17)
        ax.set_xlabel("Translated length (amino acids)")
        ax.set_title(f"One matched example: {example['gene_name']}\n{example['transcript_id']}", fontsize=14, pad=18)
        ax.set_axisbelow(True)
        ax.xaxis.grid(True, color="#E3E7EA", linewidth=0.7)
        fig.text(0.23, 0.09, f"Same A insertion: original CDS base {int(example['mutation_base_position']):,}; retained-input base {int(example['retained_mutation_base_position']):,}.", fontsize=10)
        fig.text(0.23, 0.035, "The retained control already stops before this site; adding A leaves that predicted protein unchanged.", fontsize=9)
        fig.subplots_adjust(left=0.23, right=0.95, top=0.78, bottom=0.24)
        fig.savefig(args.figure_dir / "first_intron_example.png", dpi=180)
        plt.close(fig)
    print(f"Saved comparison figures in {args.figure_dir}")


if __name__ == "__main__":
    main()
