#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import numpy as np


###############################################################################
## Configuration
###############################################################################

WF = (
    Path.home()
    / "methanotrophs"
    / "methanotroph_project"
    / "jeppe"
    / "genome_analysis_workflow"
)

CLUSTER = "Cluster_00063"

OUTDIR = (
    WF
    / "08_fegenie_composition"
    / "Cluster_00063_MCA0421"
)

OUTDIR.mkdir(
    parents=True,
    exist_ok=True
)


###############################################################################
## Input files
###############################################################################

CLUSTER_MEMBERSHIP = (
    WF
    / "06_clustering"
    / "cluster_membership.tsv"
)

PROTEIN_MODULE_MEMBERSHIP = (
    WF
    / "09_module_occurrence"
    / "protein_module_membership.tsv"
)

TAXONOMY = (
    WF
    / "16_visualization"
    / "globdb_r226_taxonomy.tsv"
)


###############################################################################
## Output files
###############################################################################

SUMMARY_OUT = (
    OUTDIR
    / "cluster_00063_summary.tsv"
)

MODULE_OUT = (
    OUTDIR
    / "cluster_00063_module_assignment.tsv"
)

PROTEINS_OUT = (
    OUTDIR
    / "cluster_00063_proteins.tsv"
)

GENOMES_OUT = (
    OUTDIR
    / "cluster_00063_genomes.tsv"
)

COPY_NUMBER_OUT = (
    OUTDIR
    / "cluster_00063_copy_number_summary.tsv"
)


###############################################################################
## Helpers
###############################################################################

def require_columns(df, required, name):

    missing = set(required) - set(df.columns)

    if missing:
        raise RuntimeError(
            f"{name} is missing required columns: "
            f"{sorted(missing)}"
        )


def taxonomy_rank(taxonomy, prefix):

    if pd.isna(taxonomy):
        return np.nan

    for item in str(taxonomy).split(";"):

        if item.startswith(prefix):
            return item[len(prefix):]

    return np.nan


def format_value(value):

    if pd.isna(value):
        return ""

    text = str(value)

    try:
        number = float(text)

        if number.is_integer():
            return str(int(number))

    except ValueError:
        pass

    return text


def join_values(series):

    return ";".join(
        format_value(x)
        for x in series
        if not pd.isna(x)
        and str(x) != ""
    )


###############################################################################
## Check inputs
###############################################################################

for path in [
    CLUSTER_MEMBERSHIP,
    PROTEIN_MODULE_MEMBERSHIP,
    TAXONOMY,
]:

    if not path.is_file():
        raise FileNotFoundError(
            f"Required input not found: {path}"
        )


###############################################################################
## Read Cluster_00063 membership
###############################################################################

membership = pd.read_csv(
    CLUSTER_MEMBERSHIP,
    sep="\t",
    low_memory=False
)

require_columns(
    membership,
    [
        "cluster",
        "genome",
        "protein_id",
        "length",
        "number_of_hemes",
    ],
    "cluster_membership.tsv"
)

cluster = membership[
    membership["cluster"] == CLUSTER
].copy()

print("=" * 100)
print(
    "Cluster_00063 — MCA0421 / c553O-LIKE FAMILY"
)
print("=" * 100)

print(
    f"Raw rows: {len(cluster)}"
)


###############################################################################
## Ensure one row per protein
###############################################################################

cluster = (
    cluster
    .drop_duplicates(
        subset=["protein_id"]
    )
    .copy()
)

print(
    f"Unique proteins: "
    f"{cluster['protein_id'].nunique()}"
)

print(
    f"Unique genomes: "
    f"{cluster['genome'].nunique()}"
)

if cluster.empty:

    raise RuntimeError(
        f"{CLUSTER} was not found in "
        f"{CLUSTER_MEMBERSHIP}"
    )


###############################################################################
## Numeric versions of length and heme count
###############################################################################

cluster["_length_num"] = pd.to_numeric(
    cluster["length"],
    errors="coerce"
)

cluster["_heme_num"] = pd.to_numeric(
    cluster["number_of_hemes"],
    errors="coerce"
)

if cluster["_length_num"].isna().any():

    bad = cluster.loc[
        cluster["_length_num"].isna(),
        ["protein_id", "length"]
    ]

    raise RuntimeError(
        "Non-numeric protein lengths found:\n"
        + bad.to_string(index=False)
    )


###############################################################################
## Basic Cluster_00063 summary
###############################################################################

min_length = int(
    cluster["_length_num"].min()
)

max_length = int(
    cluster["_length_num"].max()
)

summary = pd.DataFrame(
    [
        {
            "cluster":
                CLUSTER,

            "n_proteins":
                cluster["protein_id"].nunique(),

            "n_genomes":
                cluster["genome"].nunique(),

            "length_interval":
                f"{min_length}\u2013{max_length} aa",

            "mean_length":
                round(
                    cluster["_length_num"].mean(),
                    1
                ),

            "median_length":
                round(
                    cluster["_length_num"].median(),
                    1
                ),

            "min_hemes":
                int(
                    cluster["_heme_num"].min()
                ),

            "max_hemes":
                int(
                    cluster["_heme_num"].max()
                ),

            "median_hemes":
                round(
                    cluster["_heme_num"].median(),
                    1
                ),
        }
    ]
)

summary.to_csv(
    SUMMARY_OUT,
    sep="\t",
    index=False
)

print()
print("=" * 100)
print("BASIC SUMMARY")
print("=" * 100)
print(
    summary.to_string(index=False)
)


###############################################################################
## Determine MCL module assignment
###############################################################################

modules = pd.read_csv(
    PROTEIN_MODULE_MEMBERSHIP,
    sep="\t",
    low_memory=False
)

require_columns(
    modules,
    [
        "protein_id",
        "module",
    ],
    "protein_module_membership.tsv"
)

print()
print("=" * 100)
print("MODULE ASSIGNMENT")
print("=" * 100)

print(
    "protein_module_membership.tsv columns:"
)

print(
    modules.columns.tolist()
)


## Prefer the direct cluster assignment when available.
if "cluster" in modules.columns:

    cluster_modules = (
        modules[
            modules["cluster"] == CLUSTER
        ]
        [["cluster", "module"]]
        .drop_duplicates()
        .copy()
    )

else:

    cluster_modules = pd.DataFrame()


## Fallback: recover module assignment from the protein IDs.
if cluster_modules.empty:

    target_proteins = set(
        cluster["protein_id"]
    )

    cluster_modules = (
        modules[
            modules["protein_id"].isin(
                target_proteins
            )
        ]
        [["module"]]
        .drop_duplicates()
        .copy()
    )

    cluster_modules.insert(
        0,
        "cluster",
        CLUSTER
    )


module_values = sorted(
    cluster_modules[
        "module"
    ]
    .dropna()
    .astype(str)
    .unique()
)

if len(module_values) != 1:

    raise RuntimeError(
        f"Expected exactly one MCL module for "
        f"{CLUSTER}, found: {module_values}"
    )

TARGET_MODULE = module_values[0]

cluster_modules = pd.DataFrame(
    {
        "cluster": [CLUSTER],
        "module": [TARGET_MODULE],
    }
)

cluster_modules.to_csv(
    MODULE_OUT,
    sep="\t",
    index=False
)

print()
print(
    cluster_modules.to_string(index=False)
)

print()
print(
    f"Unique module assignments: "
    f"{module_values}"
)

print()
print(
    f">>> {CLUSTER} belongs to "
    f"{TARGET_MODULE} <<<"
)


###############################################################################
## Add module to individual proteins
###############################################################################

protein_modules = (
    modules[
        modules["protein_id"].isin(
            set(cluster["protein_id"])
        )
    ]
    [["protein_id", "module"]]
    .drop_duplicates(
        subset=["protein_id"]
    )
)

cluster = cluster.merge(
    protein_modules,
    on="protein_id",
    how="left"
)

cluster["module"] = (
    cluster["module"]
    .fillna(TARGET_MODULE)
)


###############################################################################
## Read GTDB r226 taxonomy
###############################################################################

taxonomy = pd.read_csv(
    TAXONOMY,
    sep="\t",
    header=None,
    dtype=str
)

taxonomy = (
    taxonomy
    .iloc[:, :2]
    .copy()
)

taxonomy.columns = [
    "genome",
    "taxonomy",
]

taxonomy["genus"] = (
    taxonomy["taxonomy"]
    .map(
        lambda x:
        taxonomy_rank(
            x,
            "g__"
        )
    )
)

taxonomy["species"] = (
    taxonomy["taxonomy"]
    .map(
        lambda x:
        taxonomy_rank(
            x,
            "s__"
        )
    )
)

taxonomy = (
    taxonomy[
        [
            "genome",
            "genus",
            "species",
        ]
    ]
    .drop_duplicates(
        subset=["genome"]
    )
)


###############################################################################
## Add taxonomy
###############################################################################

cluster = cluster.merge(
    taxonomy,
    on="genome",
    how="left"
)


###############################################################################
## Protein-level output
###############################################################################

preferred_columns = [
    "genome",
    "genus",
    "species",
    "protein_id",
    "cluster",
    "module",
    "representative_protein",
    "length",
    "candidate_source",
    "fegenie_positive",
    "fegenie_HMMs",
    "fegenie_categories",
    "findmehemes_positive",
    "number_of_hemes",
    "signalp_prediction",
    "deeptmhmm_class",
    "deeptmhmm_n_tm_helices",
    "export_evidence",
    "localization_class",
]

protein_columns = [
    column
    for column in preferred_columns
    if column in cluster.columns
]

proteins = (
    cluster[
        protein_columns
    ]
    .sort_values(
        [
            "genome",
            "protein_id",
        ]
    )
    .reset_index(drop=True)
)

proteins.to_csv(
    PROTEINS_OUT,
    sep="\t",
    index=False
)


###############################################################################
## Genome-level summary
###############################################################################

genomes = (
    cluster
    .sort_values(
        [
            "genome",
            "protein_id",
        ]
    )
    .groupby(
        "genome",
        as_index=False,
        dropna=False
    )
    .agg(
        genus=(
            "genus",
            "first"
        ),

        species=(
            "species",
            "first"
        ),

        n_cluster_proteins=(
            "protein_id",
            "nunique"
        ),

        lengths=(
            "_length_num",
            join_values
        ),

        heme_counts=(
            "_heme_num",
            join_values
        ),
    )
)

genomes = (
    genomes
    .sort_values(
        [
            "genus",
            "species",
            "genome",
        ],
        na_position="last"
    )
    .reset_index(drop=True)
)

genomes.to_csv(
    GENOMES_OUT,
    sep="\t",
    index=False
)

print()
print("=" * 120)
print(
    "GENOMES CARRYING CLUSTER_00063"
)
print("=" * 120)

print(
    genomes.to_string(
        index=False
    )
)


###############################################################################
## Copy-number distribution
###############################################################################

copy_number = (
    genomes
    .groupby(
        "n_cluster_proteins"
    )
    .size()
    .reset_index(
        name="n_genomes"
    )
    .rename(
        columns={
            "n_cluster_proteins":
                "copies_per_genome"
        }
    )
    .sort_values(
        "copies_per_genome"
    )
    .reset_index(drop=True)
)

copy_number.to_csv(
    COPY_NUMBER_OUT,
    sep="\t",
    index=False
)

print()
print("=" * 100)
print(
    "COPY NUMBER PER GENOME"
)
print("=" * 100)

print(
    copy_number.to_string(
        index=False
    )
)


###############################################################################
## Final checks
###############################################################################

expected_n_proteins = 46
expected_n_genomes = 20

observed_n_proteins = (
    cluster["protein_id"]
    .nunique()
)

observed_n_genomes = (
    cluster["genome"]
    .nunique()
)

if (
    observed_n_proteins
    != expected_n_proteins
):

    print(
        "WARNING: expected 46 Cluster_00063 "
        f"proteins but found "
        f"{observed_n_proteins}."
    )

if (
    observed_n_genomes
    != expected_n_genomes
):

    print(
        "WARNING: expected 20 genomes but "
        f"found {observed_n_genomes}."
    )


###############################################################################
## Report outputs
###############################################################################

print()
print("=" * 100)
print("OUTPUT FILES")
print("=" * 100)

for path in sorted(
    OUTDIR.glob("*.tsv")
):

    print(path)
