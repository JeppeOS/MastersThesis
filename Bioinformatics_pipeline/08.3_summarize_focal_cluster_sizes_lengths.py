#!/usr/bin/env python3

from pathlib import Path
import pandas as pd


###############################################################################
## Paths
###############################################################################

## This script is expected to be located in:
## genome_analysis_workflow/08_fegenie_composition/
SCRIPT_DIR = Path(__file__).resolve().parent
WORKFLOW_DIR = SCRIPT_DIR.parent

MEMBERSHIP_FILE = (
    WORKFLOW_DIR
    / "06_clustering"
    / "cluster_membership.tsv"
)

OUTPUT_FILE = (
    SCRIPT_DIR
    / "focal_cluster_size_length_summary.tsv"
)


###############################################################################
## Focal clusters
###############################################################################

CLUSTERS = [
    "Cluster_00035",
    "Cluster_00050",
    "Cluster_00206",
]


###############################################################################
## Read cluster membership table
###############################################################################

print(f"Reading: {MEMBERSHIP_FILE}")

df = pd.read_csv(
    MEMBERSHIP_FILE,
    sep="\t"
)

required_columns = {
    "cluster",
    "protein_id",
    "genome",
    "length",
}

missing = required_columns - set(df.columns)

if missing:
    raise ValueError(
        "cluster_membership.tsv is missing the following required columns: "
        + ", ".join(sorted(missing))
    )


###############################################################################
## Retain focal clusters
###############################################################################

focal = df[
    df["cluster"].isin(CLUSTERS)
].copy()

if focal.empty:
    raise ValueError(
        "None of the focal clusters were found in cluster_membership.tsv"
    )


###############################################################################
## Remove duplicated protein assignments
###############################################################################

## Each protein should contribute only once to a cluster summary.
focal = focal.drop_duplicates(
    subset=["cluster", "protein_id"]
).copy()

## Ensure protein length is numeric.
focal["length"] = pd.to_numeric(
    focal["length"],
    errors="coerce"
)

if focal["length"].isna().any():
    bad = focal.loc[
        focal["length"].isna(),
        ["cluster", "protein_id"]
    ]

    raise ValueError(
        "Non-numeric or missing protein lengths were found:\n"
        + bad.to_string(index=False)
    )


###############################################################################
## Summarize cluster size and protein lengths
###############################################################################

summary = (
    focal
    .groupby("cluster", as_index=False)
    .agg(
        n_proteins=("protein_id", "nunique"),
        n_genomes=("genome", "nunique"),
        min_length=("length", "min"),
        max_length=("length", "max"),
        mean_length=("length", "mean"),
        median_length=("length", "median"),
    )
)

## Keep the clusters in the requested biological order rather than sorting
## alphabetically.
summary["cluster"] = pd.Categorical(
    summary["cluster"],
    categories=CLUSTERS,
    ordered=True
)

summary = (
    summary
    .sort_values("cluster")
    .reset_index(drop=True)
)

## Make the output easier to read.
summary["mean_length"] = summary["mean_length"].round(1)
summary["median_length"] = summary["median_length"].round(1)

summary["length_interval"] = (
    summary["min_length"].astype(int).astype(str)
    + "-"
    + summary["max_length"].astype(int).astype(str)
    + " aa"
)


###############################################################################
## Write output
###############################################################################

summary.to_csv(
    OUTPUT_FILE,
    sep="\t",
    index=False
)

print()
print("Focal cluster summary:")
print(summary.to_string(index=False))

print()
print(f"Written to: {OUTPUT_FILE}")
