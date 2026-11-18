import csv
from collections import defaultdict, Counter
from pathlib import Path
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
FIGURES = ROOT / "figures"

PATHS = [
    DATA / "interproscan_batch_1.tsv",
    DATA / "interproscan_batch_2.tsv",
    DATA / "interproscan_batch_3.tsv"
]

OUTPUT = DATA / "interproscan_combined_summary.csv"


# --------------------------------------------------
# Read all InterProScan TSV files
# --------------------------------------------------

rows = []

for path in PATHS:

    with path.open(
        newline="",
        encoding="utf-8"
    ) as f:

        reader = csv.reader(
            f,
            delimiter="\t"
        )

        rows.extend(reader)


# --------------------------------------------------
# Group matches by protein
# --------------------------------------------------

by_protein = defaultdict(list)

database_counts = Counter()

for row in rows:

    protein = row[0]

    by_protein[
        protein
    ].append(row)

    database_counts[
        row[3]
    ] += 1


# --------------------------------------------------
# Summarize each protein
# --------------------------------------------------

summary = []

for protein, hits in by_protein.items():

    parts = protein.split("|")

    display_name = (
        parts[0]
        if len(parts) > 0
        else protein
    )

    transcript_id = (
        parts[1]
        if len(parts) > 1
        else ""
    )

    frame = (
        parts[2]
        if len(parts) > 2
        else ""
    )

    wt_length = None
    fs_length = None

    if len(parts) > 3:

        wt_length = int(
            parts[3].replace(
                "WT=",
                ""
            )
        )

    if len(parts) > 4:

        fs_length = int(
            parts[4].replace(
                "FS=",
                ""
            )
        )

    # ----------------------------------------------
    # Unique InterPro IDs
    #
    # TSV column 11 = InterPro accession
    # ----------------------------------------------

    unique_interpro_ids = sorted({
        hit[11]
        for hit in hits
        if len(hit) > 11
        and hit[11] not in (
            "",
            "-"
        )
    })


    # ----------------------------------------------
    # Count member database matches
    # ----------------------------------------------

    pfam_matches = sum(
        1
        for hit in hits
        if hit[3] == "Pfam"
    )

    cath_matches = sum(
        1
        for hit in hits
        if hit[3] == "CATH-Gene3D"
    )

    superfamily_matches = sum(
        1
        for hit in hits
        if hit[3] == "SUPERFAMILY"
    )


    summary.append({

        "display_name":
            display_name,

        "transcript_id":
            transcript_id,

        "frame":
            frame,

        "wildtype_length_aa":
            wt_length,

        "frameshift_length_aa":
            fs_length,

        "extension_aa":
            (
                fs_length - wt_length
                if wt_length is not None
                and fs_length is not None
                else ""
            ),

        "total_domain_matches":
            len(hits),

        "unique_interpro_ids":
            len(
                unique_interpro_ids
            ),

        "interpro_ids":
            ";".join(
                unique_interpro_ids
            ),

        "pfam_matches":
            pfam_matches,

        "cath_gene3d_matches":
            cath_matches,

        "superfamily_matches":
            superfamily_matches
    })


# --------------------------------------------------
# Sort strongest candidates first
# --------------------------------------------------

summary.sort(
    key=lambda row: (
        row[
            "total_domain_matches"
        ],
        row[
            "unique_interpro_ids"
        ]
    ),
    reverse=True
)


# --------------------------------------------------
# Save combined summary CSV
# --------------------------------------------------

with OUTPUT.open(
    "w",
    newline="",
    encoding="utf-8"
) as f:

    writer = csv.DictWriter(
        f,
        fieldnames=summary[0].keys()
    )

    writer.writeheader()
    writer.writerows(
        summary
    )


# --------------------------------------------------
# Overall results
# --------------------------------------------------

print(
    "Total InterProScan matches:",
    len(rows)
)

print(
    "Proteins with at least one match:",
    len(summary)
)

print(
    "\nMatches by database:"
)

print(
    "Pfam:",
    database_counts["Pfam"]
)

print(
    "CATH-Gene3D:",
    database_counts[
        "CATH-Gene3D"
    ]
)

print(
    "SUPERFAMILY:",
    database_counts[
        "SUPERFAMILY"
    ]
)


# --------------------------------------------------
# +1 vs +2 summary
# --------------------------------------------------

frame_proteins = Counter()
frame_matches = Counter()
frame_unique_interpro = Counter()

for row in summary:

    frame = row["frame"]

    frame_proteins[
        frame
    ] += 1

    frame_matches[
        frame
    ] += row[
        "total_domain_matches"
    ]

    frame_unique_interpro[
        frame
    ] += row[
        "unique_interpro_ids"
    ]


print(
    "\n+1 vs +2:"
)

for frame in [
    "+1",
    "+2"
]:

    print(
        frame
    )

    print(
        "Proteins with hits:",
        frame_proteins[
            frame
        ]
    )

    print(
        "Total matches:",
        frame_matches[
            frame
        ]
    )

    print(
        "Summed unique InterPro IDs:",
        frame_unique_interpro[
            frame
        ]
    )


# --------------------------------------------------
# Print strongest candidates
# --------------------------------------------------

print(
    "\nTop proteins:"
)

for row in summary[:15]:

    print(
        row[
            "display_name"
        ],
        row[
            "transcript_id"
        ],
        row[
            "frame"
        ],
        "| matches =",
        row[
            "total_domain_matches"
        ],
        "| unique InterPro IDs =",
        row[
            "unique_interpro_ids"
        ]
    )


# --------------------------------------------------
# Create figures folder
# --------------------------------------------------

FIGURES.mkdir(
    parents=True,
    exist_ok=True
)


# --------------------------------------------------
# Labels for figures
# --------------------------------------------------

labels = []

for row in summary:

    if (
        row["display_name"]
        == row["transcript_id"]
    ):

        name = row[
            "transcript_id"
        ]

    else:

        name = row[
            "display_name"
        ]

    labels.append(
        f"{name} {row['frame']}"
    )


# --------------------------------------------------
# Figure 1:
# Total domain matches per protein
# --------------------------------------------------

values = [
    row[
        "total_domain_matches"
    ]
    for row in summary
]

fig, ax = plt.subplots(
    figsize=(11, 7)
)

ax.barh(
    labels[::-1],
    values[::-1]
)

ax.set_xlabel(
    "Total domain/signature matches"
)

ax.set_title(
    "InterProScan Domain Matches per Frameshift Protein"
)

plt.tight_layout()

plt.savefig(
    FIGURES
    / "interproscan_total_domain_matches.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close()


# --------------------------------------------------
# Figure 2:
# Unique InterPro IDs per protein
# --------------------------------------------------

unique_values = [
    row[
        "unique_interpro_ids"
    ]
    for row in summary
]

fig, ax = plt.subplots(
    figsize=(11, 7)
)

ax.barh(
    labels[::-1],
    unique_values[::-1]
)

ax.set_xlabel(
    "Number of unique InterPro IDs"
)

ax.set_title(
    "Unique InterPro Domains per Frameshift Protein"
)

plt.tight_layout()

plt.savefig(
    FIGURES
    / "interproscan_unique_interpro_ids.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close()


# --------------------------------------------------
# Figure 3:
# +1 vs +2 comparison
# --------------------------------------------------

frame_labels = [
    "+1",
    "+2"
]

total_matches = [
    frame_matches[
        frame
    ]
    for frame in frame_labels
]

unique_counts = [
    frame_unique_interpro[
        frame
    ]
    for frame in frame_labels
]

x = range(
    len(frame_labels)
)

width = 0.35

fig, ax = plt.subplots(
    figsize=(7, 5)
)

ax.bar(
    [
        i - width / 2
        for i in x
    ],
    total_matches,
    width,
    label="Total matches"
)

ax.bar(
    [
        i + width / 2
        for i in x
    ],
    unique_counts,
    width,
    label="Summed unique InterPro IDs"
)

ax.set_xticks(
    list(x)
)

ax.set_xticklabels(
    frame_labels
)

ax.set_ylabel(
    "Count"
)

ax.set_title(
    "+1 vs +2 InterProScan Results"
)

ax.legend()

plt.tight_layout()

plt.savefig(
    FIGURES
    / "interproscan_plus1_vs_plus2.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close()


# --------------------------------------------------
# Finished
# --------------------------------------------------

print(
    "\nSaved:"
)

print(
    OUTPUT
)

print(
    FIGURES
    / "interproscan_total_domain_matches.png"
)

print(
    FIGURES
    / "interproscan_unique_interpro_ids.png"
)

print(
    FIGURES
    / "interproscan_plus1_vs_plus2.png"
)