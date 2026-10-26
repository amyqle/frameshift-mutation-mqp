import csv

input_file = "data/simple_frameshift_insertion_results.csv"

rows = []

with open(input_file, "r") as f:
    reader = csv.DictReader(f)

    for row in reader:

        row["mutated_protein_length"] = int(
            row["mutated_protein_length"]
        )

        row["premature_stop"] = (
            row["premature_stop"].lower() == "true"
        )

        row["entered_downstream_region"] = (
            row["entered_downstream_region"].lower() == "true"
        )

        row["no_stop"] = (
            row["stop_nt_position"] == ""
        )

        rows.append(row)


for shift in ["+1", "+2"]:

    subset = [
        row for row in rows
        if row["frameshift"] == shift
    ]

    total = len(subset)

    premature = sum(
        row["premature_stop"]
        for row in subset
    )

    downstream = sum(
        row["entered_downstream_region"]
        for row in subset
    )

    no_stop = sum(
        row["no_stop"]
        for row in subset
    )

    avg_length = sum(
        row["mutated_protein_length"]
        for row in subset
    ) / total

    print(f"\nFrameshift {shift}")
    print(f"Total simulations: {total}")

    print(
        f"Average mutated protein length: "
        f"{avg_length:.1f} aa"
    )

    print(
        f"Premature stops: "
        f"{premature} "
        f"({premature / total * 100:.2f}%)"
    )

    print(
        f"Entered downstream region: "
        f"{downstream} "
        f"({downstream / total * 100:.2f}%)"
    )

    print(
        f"No stop found: "
        f"{no_stop} "
        f"({no_stop / total * 100:.2f}%)"
    )