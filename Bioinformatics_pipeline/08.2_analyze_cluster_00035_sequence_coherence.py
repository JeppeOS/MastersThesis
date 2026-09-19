#!/usr/bin/env python3

from pathlib import Path
import shutil
import subprocess
import sys

import pandas as pd


## ================================================================== ##
## Configuration
## ================================================================== ##

TARGET_CLUSTER = "Cluster_00035"

HERE = Path(__file__).resolve().parent
WORKFLOW = HERE.parent

CLUSTER_MEMBERSHIP = (
    WORKFLOW
    / "06_clustering"
    / "cluster_membership.tsv"
)

OUTDIR = (
    HERE
    / "cluster_00035_sequence_analysis"
)

OUTDIR.mkdir(
    parents=True,
    exist_ok=True
)

OUT_FASTA = OUTDIR / "Cluster_00035.faa"
OUT_METADATA = OUTDIR / "Cluster_00035_members.tsv"
OUT_LENGTH_SUMMARY = OUTDIR / "Cluster_00035_length_summary.tsv"

OUT_MMSEQS_RAW = OUTDIR / "Cluster_00035_mmseqs_all_vs_all.tsv"
MMSEQS_TMP = OUTDIR / "mmseqs_tmp"

OUT_PAIRWISE = OUTDIR / "Cluster_00035_pairwise_identities.tsv"
OUT_PAIRWISE_SUMMARY = OUTDIR / "Cluster_00035_pairwise_identity_summary.tsv"

OUT_NON_MTOA = (
    OUTDIR
    / "Cluster_00035_non_MtoA_vs_MtoA.tsv"
)


## ================================================================== ##
## Utilities
## ================================================================== ##

def fail(message):
    print(
        f"\nERROR: {message}",
        file=sys.stderr
    )
    sys.exit(1)


def clean(value):
    if pd.isna(value):
        return ""

    value = str(value).strip()

    if value.lower() in {
        "",
        "na",
        "nan",
        "none",
        "null",
    }:
        return ""

    return value


def classify_fegenie(row):
    """
    Assign each protein to its FeGenie annotation class.

    For Cluster_00035 we expect:
        MtoA
        MtrA
        FeGenie_negative
    """

    positive = int(row["fegenie_positive"])

    if positive == 0:
        return "FeGenie_negative"

    hmm = clean(
        row["fegenie_HMMs"]
    )

    if hmm == "":
        fail(
            f"{row['protein_id']} is FeGenie-positive "
            "but has no HMM."
        )

    return hmm


def resolve_source_proteome(path_string):
    """
    Resolve source_proteome.

    First try the path stored in cluster_membership.tsv.

    If that /home/... path is unavailable, try the corresponding
    /faststorage/project/... path.
    """

    path = Path(path_string)

    if path.exists():
        return path

    old_prefix = (
        "/home/jeppeos/methanotrophs/"
    )

    new_prefix = (
        "/faststorage/project/methanotrophs/"
    )

    text = str(path)

    if text.startswith(old_prefix):
        alternative = Path(
            text.replace(
                old_prefix,
                new_prefix,
                1
            )
        )

        if alternative.exists():
            return alternative

    fail(
        "Cannot locate source proteome:\n"
        f"{path}"
    )


def read_fasta(path):
    """
    Simple FASTA reader.

    Protein ID is taken as the first whitespace-delimited token
    after '>'.
    """

    sequences = {}

    current_id = None
    current_sequence = []

    with open(path) as handle:

        for line in handle:

            line = line.rstrip()

            if not line:
                continue

            if line.startswith(">"):

                if current_id is not None:
                    sequences[current_id] = "".join(
                        current_sequence
                    )

                current_id = (
                    line[1:]
                    .split()[0]
                )

                current_sequence = []

            else:
                current_sequence.append(
                    line.strip()
                )

    if current_id is not None:
        sequences[current_id] = "".join(
            current_sequence
        )

    return sequences


def pair_type(a, b):
    """
    Deterministic category ordering.
    """

    order = {
        "MtoA": 1,
        "MtrA": 2,
        "FeGenie_negative": 3,
    }

    pair = sorted(
        [a, b],
        key=lambda x: (
            order.get(x, 99),
            x
        )
    )

    return f"{pair[0]}__vs__{pair[1]}"


## ================================================================== ##
## 1. Read target-cluster membership
## ================================================================== ##

print("=" * 80)
print("CLUSTER_00035 SEQUENCE COHERENCE ANALYSIS")
print("=" * 80)

if not CLUSTER_MEMBERSHIP.exists():
    fail(
        f"Cannot find:\n{CLUSTER_MEMBERSHIP}"
    )

members = pd.read_csv(
    CLUSTER_MEMBERSHIP,
    sep="\t",
    dtype=str
)

required = {
    "cluster",
    "genome",
    "protein_id",
    "length",
    "clustering_sequence_length",
    "fegenie_positive",
    "fegenie_HMMs",
    "number_of_hemes",
    "localization_class",
    "source_proteome",
}

missing = required - set(
    members.columns
)

if missing:
    fail(
        "Missing column(s): "
        + ", ".join(
            sorted(missing)
        )
    )

x = members[
    members["cluster"] == TARGET_CLUSTER
].copy()

if len(x) == 0:
    fail(
        f"No proteins found for {TARGET_CLUSTER}."
    )

x["length"] = pd.to_numeric(
    x["length"],
    errors="raise"
)

x["clustering_sequence_length"] = pd.to_numeric(
    x["clustering_sequence_length"],
    errors="raise"
)

x["fegenie_positive"] = pd.to_numeric(
    x["fegenie_positive"],
    errors="raise"
).astype(int)

x["fegenie_class"] = x.apply(
    classify_fegenie,
    axis=1
)

print()
print(
    f"Proteins: {len(x):,}"
)

print(
    f"Genomes:  {x['genome'].nunique():,}"
)

print()
print("FeGenie classes:")

print(
    x["fegenie_class"]
    .value_counts()
    .to_string()
)


## ================================================================== ##
## 2. Length summary
## ================================================================== ##

length_summary = (
    x.groupby(
        "fegenie_class",
        as_index=False
    )
    .agg(
        n_proteins=(
            "protein_id",
            "size"
        ),
        n_genomes=(
            "genome",
            "nunique"
        ),
        mean_length=(
            "clustering_sequence_length",
            "mean"
        ),
        sd_length=(
            "clustering_sequence_length",
            "std"
        ),
        median_length=(
            "clustering_sequence_length",
            "median"
        ),
        min_length=(
            "clustering_sequence_length",
            "min"
        ),
        max_length=(
            "clustering_sequence_length",
            "max"
        ),
    )
)

length_summary.to_csv(
    OUT_LENGTH_SUMMARY,
    sep="\t",
    index=False,
    float_format="%.2f"
)

print()
print("Length summary")
print(
    length_summary.to_string(
        index=False
    )
)


## ================================================================== ##
## 3. Extract exact protein sequences
## ================================================================== ##

proteome_cache = {}

sequence_records = []

for row in x.itertuples(
    index=False
):

    source = resolve_source_proteome(
        row.source_proteome
    )

    source_key = str(source)

    if source_key not in proteome_cache:

        proteome_cache[source_key] = (
            read_fasta(source)
        )

    proteome_sequences = (
        proteome_cache[source_key]
    )

    if row.protein_id not in proteome_sequences:
        fail(
            f"Protein {row.protein_id} not found in:\n"
            f"{source}"
        )

    raw_sequence = (
        proteome_sequences[
            row.protein_id
        ]
    )

    ## Prodigal proteins can contain a terminal '*'.
    sequence = raw_sequence.rstrip("*")

    if (
        len(sequence)
        != row.clustering_sequence_length
    ):
        fail(
            f"Sequence length mismatch for "
            f"{row.protein_id}:\n"
            f"FASTA = {len(sequence)}\n"
            f"cluster_membership = "
            f"{row.clustering_sequence_length}"
        )

    sequence_records.append(
        (
            row.protein_id,
            sequence
        )
    )


## Write FASTA ##

with open(
    OUT_FASTA,
    "w"
) as handle:

    for protein_id, sequence in sequence_records:

        handle.write(
            f">{protein_id}\n"
        )

        for start in range(
            0,
            len(sequence),
            80
        ):
            handle.write(
                sequence[
                    start:start + 80
                ]
                + "\n"
            )


## Save member metadata ##

metadata_cols = [
    "genome",
    "protein_id",
    "fegenie_class",
    "fegenie_positive",
    "fegenie_HMMs",
    "length",
    "clustering_sequence_length",
    "number_of_hemes",
    "localization_class",
]

x[
    metadata_cols
].sort_values(
    [
        "fegenie_class",
        "genome"
    ]
).to_csv(
    OUT_METADATA,
    sep="\t",
    index=False
)

print()
print(
    f"Extracted {len(sequence_records):,} "
    "protein sequences."
)


## ================================================================== ##
## 4. Locate MMseqs2
## ================================================================== ##

MMSEQS_EXECUTABLE = Path(
    "~/miniforge3/envs/MMseqs2/bin/mmseqs"
).expanduser()

if not MMSEQS_EXECUTABLE.exists():
    fail(
        f"Cannot find MMseqs2 executable:\n"
        f"{MMSEQS_EXECUTABLE}"
    )

mmseqs = str(MMSEQS_EXECUTABLE)

print()
print(
    f"MMseqs2 found: {mmseqs}"
)


## ================================================================== ##
## 5. All-vs-all MMseqs2 search
##
## Do NOT impose the old 40% / 80% threshold here.
## We want to observe the actual within-cluster similarity first.
## ================================================================== ##

MMSEQS_TMP.mkdir(
    parents=True,
    exist_ok=True
)

cmd = [
    mmseqs,
    "easy-search",
    str(OUT_FASTA),
    str(OUT_FASTA),
    str(OUT_MMSEQS_RAW),
    str(MMSEQS_TMP),

    "--max-seqs",
    "1000",

    "-s",
    "7.5",

    "--threads",
    "4",

    "--format-output",
    (
        "query,target,pident,alnlen,"
        "qcov,tcov,evalue,bits"
    ),
]

print()
print("Running all-vs-all MMseqs2...")

subprocess.run(
    cmd,
    check=True
)


## ================================================================== ##
## 6. Parse MMseqs2 pairwise results
## ================================================================== ##

columns = [
    "query",
    "target",
    "pident",
    "alnlen",
    "qcov",
    "tcov",
    "evalue",
    "bits",
]

pairs = pd.read_csv(
    OUT_MMSEQS_RAW,
    sep="\t",
    names=columns
)

for col in [
    "pident",
    "alnlen",
    "qcov",
    "tcov",
    "evalue",
    "bits",
]:
    pairs[col] = pd.to_numeric(
        pairs[col],
        errors="raise"
    )


## Remove self-hits ##

pairs = pairs[
    pairs["query"]
    != pairs["target"]
].copy()


## ================================================================== ##
## 7. Collapse reciprocal hits into unordered protein pairs
## ================================================================== ##

pairs["protein_a"] = pairs[
    ["query", "target"]
].min(
    axis=1
)

pairs["protein_b"] = pairs[
    ["query", "target"]
].max(
    axis=1
)

## If both A→B and B→A exist, retain the stronger alignment. ##

pairs = pairs.sort_values(
    "bits",
    ascending=False
)

pairs = pairs.drop_duplicates(
    [
        "protein_a",
        "protein_b"
    ],
    keep="first"
).copy()


## ================================================================== ##
## 8. Add annotations
## ================================================================== ##

annotation_lookup = dict(
    zip(
        x["protein_id"],
        x["fegenie_class"]
    )
)

genome_lookup = dict(
    zip(
        x["protein_id"],
        x["genome"]
    )
)

pairs["class_a"] = (
    pairs["protein_a"]
    .map(annotation_lookup)
)

pairs["class_b"] = (
    pairs["protein_b"]
    .map(annotation_lookup)
)

pairs["genome_a"] = (
    pairs["protein_a"]
    .map(genome_lookup)
)

pairs["genome_b"] = (
    pairs["protein_b"]
    .map(genome_lookup)
)

if (
    pairs["class_a"].isna().any()
    or
    pairs["class_b"].isna().any()
):
    fail(
        "Failed to annotate one or more MMseqs pairs."
    )

pairs["comparison"] = pairs.apply(
    lambda r: pair_type(
        r["class_a"],
        r["class_b"]
    ),
    axis=1
)

pairs["min_coverage"] = pairs[
    ["qcov", "tcov"]
].min(
    axis=1
)


## Original MMseqs clustering criterion ##

pairs[
    "meets_40id_80cov"
] = (
    (pairs["pident"] >= 40)
    &
    (pairs["min_coverage"] >= 0.80)
)


## ================================================================== ##
## 9. Check how complete the all-vs-all comparison is
## ================================================================== ##

n = len(x)

expected_pairs = (
    n * (n - 1) // 2
)

observed_pairs = len(pairs)

print()
print("Pairwise search coverage")
print(
    f"  Expected unordered pairs: {expected_pairs:,}"
)
print(
    f"  Observed MMseqs pairs:     {observed_pairs:,}"
)

if observed_pairs < expected_pairs:

    print(
        f"  No MMseqs alignment for:   "
        f"{expected_pairs - observed_pairs:,} pair(s)"
    )


## ================================================================== ##
## 10. Save full pairwise table
## ================================================================== ##

pair_columns = [
    "protein_a",
    "genome_a",
    "class_a",
    "protein_b",
    "genome_b",
    "class_b",
    "comparison",
    "pident",
    "alnlen",
    "qcov",
    "tcov",
    "min_coverage",
    "evalue",
    "bits",
    "meets_40id_80cov",
]

pairs[
    pair_columns
].sort_values(
    [
        "comparison",
        "pident"
    ],
    ascending=[
        True,
        False
    ]
).to_csv(
    OUT_PAIRWISE,
    sep="\t",
    index=False,
    float_format="%.4f"
)


## ================================================================== ##
## 11. Summarize pairwise identity by annotation combination
## ================================================================== ##

def summarize_pair_group(group):

    meets = group[
        group["meets_40id_80cov"]
    ]

    return pd.Series(
        {
            "n_pairs": len(group),

            "mean_pident":
                group["pident"].mean(),

            "median_pident":
                group["pident"].median(),

            "min_pident":
                group["pident"].min(),

            "max_pident":
                group["pident"].max(),

            "mean_min_coverage":
                group["min_coverage"].mean(),

            "n_pairs_meeting_40id_80cov":
                len(meets),

            "pct_pairs_meeting_40id_80cov":
                (
                    100
                    * len(meets)
                    / len(group)
                ),

            "mean_pident_among_40id_80cov":
                (
                    meets["pident"].mean()
                    if len(meets) > 0
                    else float("nan")
                ),
        }
    )


pair_summary = (
    pairs
    .groupby(
        "comparison"
    )
    .apply(
        summarize_pair_group,
        include_groups=False
    )
    .reset_index()
)

pair_summary.to_csv(
    OUT_PAIRWISE_SUMMARY,
    sep="\t",
    index=False,
    float_format="%.3f"
)

print()
print("Pairwise identity by FeGenie class")
print(
    pair_summary.to_string(
        index=False
    )
)


## ================================================================== ##
## 12. Specifically measure every non-MtoA protein against MtoA
##
## This directly asks:
##
##   "Does each MtrA / FeGenie-negative member have close MtoA
##    homologues within this sequence family?"
## ================================================================== ##

non_mtoa_records = []

non_mtoa = x[
    x["fegenie_class"] != "MtoA"
]

for row in non_mtoa.itertuples(
    index=False
):

    protein = row.protein_id

    relevant = pairs[
        (
            (pairs["protein_a"] == protein)
            &
            (pairs["class_b"] == "MtoA")
        )
        |
        (
            (pairs["protein_b"] == protein)
            &
            (pairs["class_a"] == "MtoA")
        )
    ].copy()

    if len(relevant) == 0:

        non_mtoa_records.append(
            {
                "protein_id": protein,
                "genome": row.genome,
                "fegenie_class": row.fegenie_class,
                "length": row.clustering_sequence_length,

                "n_MtoA_alignments": 0,

                "max_identity_to_MtoA":
                    float("nan"),

                "median_identity_to_MtoA":
                    float("nan"),

                "mean_identity_to_MtoA":
                    float("nan"),

                "max_min_coverage_to_MtoA":
                    float("nan"),

                "n_MtoA_pairs_meeting_40id_80cov":
                    0,
            }
        )

        continue

    qualifying = relevant[
        relevant["meets_40id_80cov"]
    ]

    non_mtoa_records.append(
        {
            "protein_id": protein,
            "genome": row.genome,
            "fegenie_class": row.fegenie_class,
            "length": row.clustering_sequence_length,

            "n_MtoA_alignments":
                len(relevant),

            "max_identity_to_MtoA":
                relevant["pident"].max(),

            "median_identity_to_MtoA":
                relevant["pident"].median(),

            "mean_identity_to_MtoA":
                relevant["pident"].mean(),

            "max_min_coverage_to_MtoA":
                relevant["min_coverage"].max(),

            "n_MtoA_pairs_meeting_40id_80cov":
                len(qualifying),
        }
    )

non_mtoa_summary = pd.DataFrame(
    non_mtoa_records
)

non_mtoa_summary = (
    non_mtoa_summary
    .sort_values(
        [
            "fegenie_class",
            "max_identity_to_MtoA"
        ],
        ascending=[
            True,
            False
        ]
    )
)

non_mtoa_summary.to_csv(
    OUT_NON_MTOA,
    sep="\t",
    index=False,
    float_format="%.3f"
)

print()
print("Non-MtoA members compared with MtoA")
print(
    non_mtoa_summary.to_string(
        index=False
    )
)


## ================================================================== ##
## Final report
## ================================================================== ##

print()
print("=" * 80)
print("SUCCESS")
print("=" * 80)

print(f"FASTA:                    {OUT_FASTA}")
print(f"Member metadata:          {OUT_METADATA}")
print(f"Length summary:           {OUT_LENGTH_SUMMARY}")
print(f"Raw MMseqs results:       {OUT_MMSEQS_RAW}")
print(f"Pairwise identities:      {OUT_PAIRWISE}")
print(f"Pairwise class summary:   {OUT_PAIRWISE_SUMMARY}")
print(f"Non-MtoA vs MtoA:         {OUT_NON_MTOA}")
