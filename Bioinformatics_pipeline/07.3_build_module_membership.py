#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import re
import sys

## Paths ##
HERE = Path(__file__).resolve().parent

CYTOSCAPE_CSV = HERE / "cooccurence_network.csv"

CLUSTER_SUMMARY = (
    HERE.parent
    / "06_clustering"
    / "cluster_summary.tsv"
)

GENOME_CLUSTER_PAIRS = (
    HERE.parent
    / "06_clustering"
    / "genome_cluster_pairs.tsv"
)

OUT_MEMBERSHIP = HERE / "module_membership.tsv"
OUT_METADATA = HERE / "module_membership_metadata.tsv"
OUT_ID_MAP = HERE / "module_id_map.tsv"
OUT_SUMMARY = HERE / "module_summary.tsv"

## Expected values from the finalized network/MCL analysis ##
EXPECTED_NETWORK_NODES = 156
EXPECTED_MCL_MODULES = 35


def fail(message):
    print(f"\nERROR: {message}", file=sys.stderr)
    sys.exit(1)


def normalize_mcl(value):
    """
    Normalize Cytoscape MCL labels.

    Examples:
        "6"   -> "6"
        6     -> "6"
        6.0   -> "6"
    """
    if pd.isna(value):
        return None

    text = str(value).strip()

    if text == "":
        return None

    try:
        number = float(text)
        if number.is_integer():
            return str(int(number))
    except ValueError:
        pass

    return text


print("=" * 72)
print("BUILD MODULE MEMBERSHIP")
print("=" * 72)

## ------------------------------------------------------------------ ##
## 1. Read Cytoscape node export
## ------------------------------------------------------------------ ##

if not CYTOSCAPE_CSV.exists():
    fail(f"Cannot find Cytoscape export:\n{CYTOSCAPE_CSV}")

df = pd.read_csv(CYTOSCAPE_CSV, dtype=str)

required = {
    "__mclCluster",
    "component",
    "component_size",
    "name",
}

missing = required - set(df.columns)

if missing:
    fail(
        "Cytoscape CSV is missing required column(s): "
        + ", ".join(sorted(missing))
    )

print(f"Cytoscape rows read: {len(df):,}")


## ------------------------------------------------------------------ ##
## 2. Clean the essential identifiers
## ------------------------------------------------------------------ ##

df["cluster"] = df["name"].astype(str).str.strip()
df["mcl_cluster_original"] = df["__mclCluster"].apply(normalize_mcl)
df["component"] = df["component"].astype(str).str.strip()

## Validate cluster names ##
bad_cluster_names = df.loc[
    ~df["cluster"].str.match(r"^Cluster_\d+$", na=False),
    "cluster"
].tolist()

if bad_cluster_names:
    fail(
        "Unexpected cluster name(s), for example: "
        + ", ".join(bad_cluster_names[:10])
    )

## Each protein-family cluster must occur exactly once ##
duplicates = df.loc[
    df["cluster"].duplicated(keep=False),
    "cluster"
].unique()

if len(duplicates) > 0:
    fail(
        "Duplicate cluster rows found: "
        + ", ".join(duplicates[:20])
    )

## Every connected node must have an MCL assignment ##
missing_mcl = df.loc[
    df["mcl_cluster_original"].isna(),
    "cluster"
].tolist()

if missing_mcl:
    fail(
        f"{len(missing_mcl)} cluster(s) have no MCL assignment:\n"
        + "\n".join(missing_mcl[:20])
    )


## ------------------------------------------------------------------ ##
## 3. Network-level QC
## ------------------------------------------------------------------ ##

n_nodes = df["cluster"].nunique()
n_mcl = df["mcl_cluster_original"].nunique()
n_components = df["component"].nunique()

print()
print("Network QC")
print(f"  Unique network nodes:      {n_nodes:,}")
print(f"  Connected components:      {n_components:,}")
print(f"  Original MCL groups:       {n_mcl:,}")

if n_nodes != EXPECTED_NETWORK_NODES:
    fail(
        f"Expected {EXPECTED_NETWORK_NODES} network nodes, "
        f"but found {n_nodes}."
    )

if n_mcl != EXPECTED_MCL_MODULES:
    fail(
        f"Expected {EXPECTED_MCL_MODULES} MCL modules, "
        f"but found {n_mcl}."
    )

## If this column exists, all exported rows should be connected ##
if "is_connected_at_threshold" in df.columns:

    connected = pd.to_numeric(
        df["is_connected_at_threshold"],
        errors="coerce"
    )

    bad_connected = df.loc[connected != 1, "cluster"].tolist()

    if bad_connected:
        fail(
            "The Cytoscape export contains nodes that are not marked "
            "as connected at the final threshold."
        )


## ------------------------------------------------------------------ ##
## 4. Validate that MCL modules do not cross connected components
## ------------------------------------------------------------------ ##

mcl_component_check = (
    df.groupby("mcl_cluster_original")["component"]
      .nunique()
)

cross_component = mcl_component_check[
    mcl_component_check != 1
]

if len(cross_component) > 0:
    fail(
        "At least one MCL group spans multiple connected components:\n"
        + cross_component.to_string()
    )

print("  MCL modules crossing components: 0")


## ------------------------------------------------------------------ ##
## 5. Cross-check against authoritative MMseqs cluster summary
## ------------------------------------------------------------------ ##

if not CLUSTER_SUMMARY.exists():
    fail(f"Cannot find cluster summary:\n{CLUSTER_SUMMARY}")

cluster_summary = pd.read_csv(
    CLUSTER_SUMMARY,
    sep="\t",
    dtype=str
)

if "cluster" not in cluster_summary.columns:
    fail(
        f"{CLUSTER_SUMMARY.name} does not contain a 'cluster' column."
    )

known_clusters = set(
    cluster_summary["cluster"].dropna().astype(str).str.strip()
)

missing_from_clustering = sorted(
    set(df["cluster"]) - known_clusters
)

if missing_from_clustering:
    fail(
        "Network clusters missing from authoritative cluster_summary.tsv:\n"
        + "\n".join(missing_from_clustering[:20])
    )

print("  All network nodes found in cluster_summary.tsv: YES")


## ------------------------------------------------------------------ ##
## 6. Determine deterministic module IDs
##
## Rule:
##   1. Larger MCL modules receive lower Module numbers.
##   2. Ties are broken by the lexicographically lowest Cluster_XXXXX.
##
## This avoids treating Cytoscape's arbitrary MCL numbering as the
## permanent biological module identifier.
## ------------------------------------------------------------------ ##

mcl_stats = (
    df.groupby("mcl_cluster_original", as_index=False)
      .agg(
          module_size_clusters=("cluster", "nunique"),
          component=("component", "first"),
          minimum_cluster=("cluster", "min"),
      )
)

mcl_stats = mcl_stats.sort_values(
    ["module_size_clusters", "minimum_cluster"],
    ascending=[False, True],
    kind="stable"
).reset_index(drop=True)

width = max(2, len(str(len(mcl_stats))))

mcl_stats["module"] = [
    f"Module_{i:0{width}d}"
    for i in range(1, len(mcl_stats) + 1)
]

mcl_to_module = dict(
    zip(
        mcl_stats["mcl_cluster_original"],
        mcl_stats["module"]
    )
)

mcl_to_size = dict(
    zip(
        mcl_stats["mcl_cluster_original"],
        mcl_stats["module_size_clusters"]
    )
)

df["module"] = df["mcl_cluster_original"].map(mcl_to_module)

df["module_size_clusters"] = (
    df["mcl_cluster_original"]
    .map(mcl_to_size)
    .astype(int)
)


## ------------------------------------------------------------------ ##
## 7. Validate module assignment
## ------------------------------------------------------------------ ##

if df["module"].isna().any():
    fail("At least one network cluster failed to receive a module ID.")

assignment_counts = (
    df.groupby("cluster")["module"]
      .nunique()
)

bad_assignments = assignment_counts[
    assignment_counts != 1
]

if len(bad_assignments) > 0:
    fail("Some clusters received more than one module assignment.")

print()
print("Module assignment QC")
print(f"  Clusters assigned:         {df['cluster'].nunique():,}")
print(f"  Final modules:             {df['module'].nunique():,}")
print("  Multiple assignments:      0")


## ------------------------------------------------------------------ ##
## 8. Authoritative minimal membership table
##
## This is deliberately kept simple for downstream joins:
##
##       cluster -> module
## ------------------------------------------------------------------ ##

membership = (
    df[
        [
            "cluster",
            "module",
        ]
    ]
    .sort_values(
        ["module", "cluster"],
        kind="stable"
    )
    .reset_index(drop=True)
)

membership.to_csv(
    OUT_MEMBERSHIP,
    sep="\t",
    index=False
)


## ------------------------------------------------------------------ ##
## 9. Detailed module-membership metadata
## ------------------------------------------------------------------ ##

metadata_columns = [
    "cluster",
    "module",
    "mcl_cluster_original",
    "module_size_clusters",
    "component",
    "component_size",
]

## Preserve useful cluster-level attributes from the Cytoscape export ##
optional_columns = [
    "n_genomes",
    "n_proteins",
    "degree",
    "weighted_degree",
    "dominant_fegenie_HMM",
    "dominant_signalp_prediction",
    "dominant_deeptmhmm_class",
    "dominant_localization_class",
    "mean_heme_count",
    "max_heme_count",
    "n_fegenie_positive",
    "n_fegenie_negative",
    "pct_fegenie_negative",
    "n_strong_unknown_soluble_periplasmic_like",
    "pct_strong_unknown_soluble_periplasmic_like",
    "representative_protein",
]

metadata_columns += [
    col
    for col in optional_columns
    if col in df.columns
]

metadata = (
    df[metadata_columns]
    .sort_values(
        ["module", "cluster"],
        kind="stable"
    )
    .reset_index(drop=True)
)

metadata.to_csv(
    OUT_METADATA,
    sep="\t",
    index=False
)


## ------------------------------------------------------------------ ##
## 10. MCL -> deterministic module ID map
## ------------------------------------------------------------------ ##

module_map = (
    mcl_stats[
        [
            "module",
            "mcl_cluster_original",
            "module_size_clusters",
            "component",
            "minimum_cluster",
        ]
    ]
    .sort_values("module")
    .reset_index(drop=True)
)

module_map.to_csv(
    OUT_ID_MAP,
    sep="\t",
    index=False
)


## ------------------------------------------------------------------ ##
## 11. Build module-level summary
##
## n_genomes_with_module is calculated from the authoritative
## genome_cluster_pairs.tsv rather than summing cluster prevalence.
## ------------------------------------------------------------------ ##

if not GENOME_CLUSTER_PAIRS.exists():
    fail(
        f"Cannot find genome-cluster pairs:\n{GENOME_CLUSTER_PAIRS}"
    )

pairs = pd.read_csv(
    GENOME_CLUSTER_PAIRS,
    sep="\t",
    dtype=str
)

required_pair_cols = {"genome", "cluster"}

missing_pair_cols = required_pair_cols - set(pairs.columns)

if missing_pair_cols:
    fail(
        "genome_cluster_pairs.tsv is missing required column(s): "
        + ", ".join(sorted(missing_pair_cols))
    )

pairs = pairs[
    ["genome", "cluster"]
].drop_duplicates()

pairs_with_module = pairs.merge(
    membership,
    on="cluster",
    how="inner",
    validate="many_to_one"
)

module_genome_counts = (
    pairs_with_module
    .groupby("module")["genome"]
    .nunique()
    .rename("n_genomes_with_module")
    .reset_index()
)

summary = module_map.merge(
    module_genome_counts,
    on="module",
    how="left",
    validate="one_to_one"
)

summary["n_genomes_with_module"] = (
    summary["n_genomes_with_module"]
    .fillna(0)
    .astype(int)
)

## Optional total protein count across the member sequence clusters ##
if "n_proteins" in df.columns:

    tmp = df[
        ["module", "n_proteins"]
    ].copy()

    tmp["n_proteins"] = pd.to_numeric(
        tmp["n_proteins"],
        errors="coerce"
    )

    protein_totals = (
        tmp.groupby("module")["n_proteins"]
           .sum(min_count=1)
           .rename("n_cluster_member_proteins")
           .reset_index()
    )

    summary = summary.merge(
        protein_totals,
        on="module",
        how="left",
        validate="one_to_one"
    )

summary = summary.sort_values(
    "module"
).reset_index(drop=True)

summary.to_csv(
    OUT_SUMMARY,
    sep="\t",
    index=False
)


## ------------------------------------------------------------------ ##
## 12. Final report
## ------------------------------------------------------------------ ##

print()
print("Final module sizes")
print(
    summary[
        [
            "module",
            "mcl_cluster_original",
            "module_size_clusters",
            "component",
            "n_genomes_with_module",
        ]
    ].to_string(index=False)
)

print()
print("=" * 72)
print("SUCCESS")
print("=" * 72)
print(f"Authoritative membership : {OUT_MEMBERSHIP}")
print(f"Membership metadata      : {OUT_METADATA}")
print(f"Module ID map            : {OUT_ID_MAP}")
print(f"Module summary           : {OUT_SUMMARY}")
print()
print(
    f"{len(membership):,} protein-family clusters assigned to "
    f"{membership['module'].nunique():,} modules."
)
