from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"

INPUT = DATA / "interproscan_candidates_over100aa.fasta"

BATCH_SIZE = 100


def read_fasta(path):
    records = []
    header = None
    sequence = []

    with path.open() as f:
        for line in f:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    records.append(
                        (header, "".join(sequence))
                    )

                header = line
                sequence = []

            else:
                sequence.append(line)

        if header is not None:
            records.append(
                (header, "".join(sequence))
            )

    return records


records = read_fasta(INPUT)

print("Total sequences:", len(records))

for start in range(
    0,
    len(records),
    BATCH_SIZE
):

    batch_number = (
        start // BATCH_SIZE
    ) + 1

    batch = records[
        start:start + BATCH_SIZE
    ]

    output = DATA / (
        f"interproscan_batch_{batch_number}.fasta"
    )

    with output.open("w") as f:

        for header, sequence in batch:

            f.write(header + "\n")

            for i in range(
                0,
                len(sequence),
                60
            ):
                f.write(
                    sequence[i:i + 60]
                    + "\n"
                )

    print(
        f"Batch {batch_number}: "
        f"{len(batch)} sequences"
    )

    print(
        "Saved:",
        output
    )