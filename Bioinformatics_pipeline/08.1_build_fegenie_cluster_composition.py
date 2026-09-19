#!/usr/bin/env python3

from pathlib import Path
from collections import Counter
import math
import re
import sys

import pandas as pd


## ================================================================== ##
## Paths
## ================================================================== ##

HERE = Path(__file__).resolve().parent
WORKFLOW = HERE.parent

CLUSTER_MEMBERSHIP = (
    WORKFLOW
    / "06_clustering"
    / "cluster_membership.tsv"
)

CLUSTER_SUMMARY = (
    WORKFLOW
    / "06_clustering"
    / "cluster_summary.tsv"
)

MODULE_MEMBERSHIP = (
    WORKFLOW
    / "07_cooccurrence_network"
    / "module_membership.tsv"
)


## ================================================================== ##
## Outputs
## ================================================================== ##

OUT_LONG = HERE / "cluster_fegenie_hmm_long.tsv"

OUT_SIGNATURES = (
    HERE
    / "cluster_fegenie_signature_composition.tsv"
)

OUT_SUMMARY = (
    HERE
    / "cluster_fegenie_composition_summary.tsv"
)

OUT_MATRIX_COUNTS = (
    HERE
    / "cluster_fegenie_hmm_matrix_counts.tsv"
)

OUT_MATRIX_PCT_PROTEINS = (
    HERE
    / "cluster_fegenie_hmm_matrix_pct_proteins.tsv"
)

OUT_MATRIX_PCT_GENOMES = (
    HERE
    / "cluster_fegenie_hmm_matrix_pct_genomes.tsv"
)

OUT_MODULE_NETWORK_SUMMARY = (
    HERE
    / "network_module_cluster_fegenie_composition.tsv"
)

OUT_QC = (
    HERE
    / "fegenie_composition_qc.tsv"
)


## ================================================================== ##
## Expected values from the locked upstream workflow
## ================================================================== ##

EXPECTED_CLUSTERED_PROTEINS = 14308
EXPECTED_CLUSTERS = 1390

EXPECTED_MODULE_CLUSTERS = 156
EXPECTED_MODULES = 35


## ================================================================== ##
## Constants
## ================================================================== ##

FEGENIE_NEGATIVE_LABEL = "FeGenie_negative"


## ================================================================== ##
## Utility functions
## ================================================================== ##

def fail(message):
    print(f"\nERROR: {message}", file=sys.stderr)
    sys.exit(1)


def clean_string(value):
    """
    Convert missing/string-like values into a clean string.
    """
    if pd.isna(value):
        return ""

    value = str(value).strip()

    if value.lower() in {
        "",
        "na",
        "nan",
        "none",
        "null"
    }:
        return ""

    return value


def parse_hmms(value):
    """
    Parse the integrated FeGenie HMM field.

    The parser tolerates common list delimiters:
        ;
        |
        ,

    It also strips brackets/quotes in case values were stored
    in a list-like representation.

    Returns a sorted tuple of unique HMM names.
    """
    value = clean_string(value)

    if value == "":
        return tuple()

    value = value.strip("[]")
    value = value.replace('"', "")
    value = value.replace("'", "")

    parts = re.split(r"[;,|]+", value)

    hmms = sorted({
        x.strip()
        for x in parts
        if x.strip() != ""
    })

    return tuple(hmms)


def hmm_signature(hmms):
    """
    Produce a mutually exclusive annotation signature.

    Examples:
        no FeGenie hit
            -> FeGenie_negative

        one HMM
            -> MtoA

        multiple HMMs
            -> MtoA|MtrA
    """
    if len(hmms) == 0:
        return FEGENIE_NEGATIVE_LABEL

    return "|".join(sorted(hmms))


def shannon_entropy(counts):
    """
    Shannon entropy using natural logarithms.

    H = -sum(p_i * ln(p_i))

    Returns 0 for a single-category distribution.
    """
    counts = [
        x
        for x in counts
        if x > 0
    ]

    total = sum(counts)

    if total == 0:
        return float("nan")

    entropy = 0.0

    for count in counts:
        p = count / total
        entropy -= p * math.log(p)

    return entropy


def normalized_entropy(counts):
    """
    Shannon entropy normalized to 0-1.

    0 = all observations have the same signature.
    1 = observations are evenly distributed among
        all observed signatures.
    """
    counts = [
        x
        for x in counts
        if x > 0
    ]

    k = len(counts)

    if k == 0:
        return float("nan")

    if k == 1:
        return 0.0

    h = shannon_entropy(counts)

    return h / math.log(k)


def dominant_labels(counter):
    """
    Return all tied dominant labels, sorted.

    This avoids silently breaking ties.
    """
    if len(counter) == 0:
        return tuple(), 0

    maximum = max(counter.values())

    labels = tuple(sorted(
        key
        for key, value in counter.items()
        if value == maximum
    ))

    return labels, maximum


## ================================================================== ##
## 1. Read authoritative clustered-protein membership
## ================================================================== ##

print("=" * 78)
print("BUILD FEGENIE CLUSTER COMPOSITION")
print("=" * 78)

if not CLUSTER_MEMBERSHIP.exists():
    fail(
        f"Cannot find:\n{CLUSTER_MEMBERSHIP}"
    )

members = pd.read_csv(
    CLUSTER_MEMBERSHIP,
    sep="\t",
    dtype=str
)

required_columns = {
    "cluster",
    "genome",
    "protein_id",
    "fegenie_positive",
    "fegenie_HMMs",
}

missing = required_columns - set(members.columns)

if missing:
    fail(
        "cluster_membership.tsv is missing required columns: "
        + ", ".join(sorted(missing))
    )

print()
print(f"Clustered protein rows read: {len(members):,}")


## ================================================================== ##
## 2. Basic identifiers and uniqueness QC
## ================================================================== ##

members["cluster"] = (
    members["cluster"]
    .astype(str)
    .str.strip()
)

members["genome"] = (
    members["genome"]
    .astype(str)
    .str.strip()
)

members["protein_id"] = (
    members["protein_id"]
    .astype(str)
    .str.strip()
)

duplicate_keys = members.duplicated(
    ["genome", "protein_id"],
    keep=False
)

if duplicate_keys.any():
    examples = (
        members.loc[
            duplicate_keys,
            ["genome", "protein_id", "cluster"]
        ]
        .head(20)
    )

    fail(
        "Duplicate genome + protein_id keys detected:\n"
        + examples.to_string(index=False)
    )

n_clustered_proteins = len(members)
n_clusters = members["cluster"].nunique()

print(f"Unique protein keys:       {n_clustered_proteins:,}")
print(f"Unique MMseqs clusters:    {n_clusters:,}")

if n_clustered_proteins != EXPECTED_CLUSTERED_PROTEINS:
    fail(
        f"Expected {EXPECTED_CLUSTERED_PROTEINS:,} clustered proteins "
        f"but found {n_clustered_proteins:,}."
    )

if n_clusters != EXPECTED_CLUSTERS:
    fail(
        f"Expected {EXPECTED_CLUSTERS:,} MMseqs clusters "
        f"but found {n_clusters:,}."
    )


## ================================================================== ##
## 3. Parse FeGenie status and HMM sets
## ================================================================== ##

members["fegenie_positive_num"] = pd.to_numeric(
    members["fegenie_positive"],
    errors="coerce"
)

bad_status = members[
    ~members["fegenie_positive_num"].isin([0, 1])
]

if len(bad_status) > 0:
    fail(
        "Unexpected values in fegenie_positive."
    )

members["fegenie_positive_num"] = (
    members["fegenie_positive_num"]
    .astype(int)
)

members["hmm_tuple"] = (
    members["fegenie_HMMs"]
    .apply(parse_hmms)
)

members["n_hmms"] = (
    members["hmm_tuple"]
    .apply(len)
)

members["signature"] = (
    members["hmm_tuple"]
    .apply(hmm_signature)
)


## ================================================================== ##
## 4. Strict internal FeGenie QC
## ================================================================== ##

positive_without_hmm = members[
    (members["fegenie_positive_num"] == 1)
    &
    (members["n_hmms"] == 0)
]

negative_with_hmm = members[
    (members["fegenie_positive_num"] == 0)
    &
    (members["n_hmms"] > 0)
]

if len(positive_without_hmm) > 0:
    fail(
        f"{len(positive_without_hmm):,} FeGenie-positive proteins "
        "have no FeGenie HMM listed."
    )

if len(negative_with_hmm) > 0:
    fail(
        f"{len(negative_with_hmm):,} FeGenie-negative proteins "
        "nevertheless have a FeGenie HMM listed."
    )

n_positive = int(
    members["fegenie_positive_num"].sum()
)

n_negative = (
    len(members) - n_positive
)

multi_hmm_proteins = int(
    (members["n_hmms"] > 1).sum()
)

unique_hmms = sorted({
    hmm
    for hmms in members["hmm_tuple"]
    for hmm in hmms
})

print()
print("FeGenie protein-level QC")
print(f"  FeGenie-positive proteins: {n_positive:,}")
print(f"  FeGenie-negative proteins: {n_negative:,}")
print(f"  Multi-HMM proteins:         {multi_hmm_proteins:,}")
print(f"  Distinct HMM names:         {len(unique_hmms):,}")


## ================================================================== ##
## 5. Attach final MCL module membership
##
## Most clusters do not belong to the co-occurrence network.
## Those clusters intentionally retain module = NA.
## ================================================================== ##

if not MODULE_MEMBERSHIP.exists():
    fail(
        f"Cannot find:\n{MODULE_MEMBERSHIP}"
    )

modules = pd.read_csv(
    MODULE_MEMBERSHIP,
    sep="\t",
    dtype=str
)

required_module_columns = {
    "cluster",
    "module"
}

missing_module = (
    required_module_columns
    - set(modules.columns)
)

if missing_module:
    fail(
        "module_membership.tsv is missing: "
        + ", ".join(sorted(missing_module))
    )

if modules["cluster"].duplicated().any():
    fail(
        "Duplicate cluster IDs detected in module_membership.tsv."
    )

if modules["module"].nunique() != EXPECTED_MODULES:
    fail(
        f"Expected {EXPECTED_MODULES} modules but found "
        f"{modules['module'].nunique()}."
    )

if len(modules) != EXPECTED_MODULE_CLUSTERS:
    fail(
        f"Expected {EXPECTED_MODULE_CLUSTERS} module-assigned clusters "
        f"but found {len(modules)}."
    )

members = members.merge(
    modules,
    on="cluster",
    how="left",
    validate="many_to_one"
)

print()
print("Module linkage QC")
print(
    f"  Module-assigned clusters: "
    f"{members.loc[members['module'].notna(), 'cluster'].nunique():,}"
)
print(
    f"  Final modules represented: "
    f"{members['module'].nunique():,}"
)


## ================================================================== ##
## 6. Cross-check totals against existing cluster_summary.tsv
## ================================================================== ##

if not CLUSTER_SUMMARY.exists():
    fail(
        f"Cannot find:\n{CLUSTER_SUMMARY}"
    )

old_summary = pd.read_csv(
    CLUSTER_SUMMARY,
    sep="\t",
    dtype=str
)

summary_required = {
    "cluster",
    "n_proteins",
    "n_genomes",
    "n_fegenie_positive",
    "n_fegenie_negative",
    "dominant_fegenie_HMM",
}

missing_summary = (
    summary_required
    - set(old_summary.columns)
)

if missing_summary:
    fail(
        "cluster_summary.tsv is missing: "
        + ", ".join(sorted(missing_summary))
    )

check = (
    members
    .groupby("cluster", as_index=False)
    .agg(
        n_proteins_new=("protein_id", "size"),
        n_genomes_new=("genome", "nunique"),
        n_fegenie_positive_new=(
            "fegenie_positive_num",
            "sum"
        ),
    )
)

check["n_fegenie_negative_new"] = (
    check["n_proteins_new"]
    - check["n_fegenie_positive_new"]
)

old_check = old_summary[
    [
        "cluster",
        "n_proteins",
        "n_genomes",
        "n_fegenie_positive",
        "n_fegenie_negative",
    ]
].copy()

for col in [
    "n_proteins",
    "n_genomes",
    "n_fegenie_positive",
    "n_fegenie_negative",
]:
    old_check[col] = pd.to_numeric(
        old_check[col],
        errors="coerce"
    )

check = check.merge(
    old_check,
    on="cluster",
    how="left",
    validate="one_to_one"
)

count_mismatch = check[
    (check["n_proteins_new"] != check["n_proteins"])
    |
    (check["n_genomes_new"] != check["n_genomes"])
    |
    (
        check["n_fegenie_positive_new"]
        != check["n_fegenie_positive"]
    )
    |
    (
        check["n_fegenie_negative_new"]
        != check["n_fegenie_negative"]
    )
]

if len(count_mismatch) > 0:
    fail(
        "New composition counts disagree with cluster_summary.tsv "
        f"for {len(count_mismatch)} cluster(s)."
    )

print(
    "  Existing cluster_summary count cross-check: PASS"
)


## ================================================================== ##
## 7. Build per-cluster HMM presence table
##
## IMPORTANT:
## A protein with multiple HMMs contributes once to EACH HMM.
##
## Therefore pct_cluster_proteins across different HMM rows can sum
## to >100% for a cluster.
##
## FeGenie_negative is mutually exclusive with HMM hits.
## ================================================================== ##

cluster_totals = (
    members
    .groupby("cluster", as_index=False)
    .agg(
        n_cluster_proteins=("protein_id", "size"),
        n_cluster_genomes=("genome", "nunique"),
    )
)

module_lookup = (
    members[
        ["cluster", "module"]
    ]
    .drop_duplicates()
)

long_records = []

for row in members.itertuples(index=False):

    if row.n_hmms == 0:
        labels = [FEGENIE_NEGATIVE_LABEL]
    else:
        labels = list(row.hmm_tuple)

    for hmm in labels:
        long_records.append(
            {
                "cluster": row.cluster,
                "module": row.module,
                "genome": row.genome,
                "protein_id": row.protein_id,
                "HMM": hmm,
            }
        )

hmm_membership = pd.DataFrame(
    long_records
)

hmm_long = (
    hmm_membership
    .groupby(
        ["cluster", "module", "HMM"],
        dropna=False,
        as_index=False
    )
    .agg(
        n_proteins_with_HMM=("protein_id", "nunique"),
        n_genomes_with_HMM=("genome", "nunique"),
    )
)

hmm_long = hmm_long.merge(
    cluster_totals,
    on="cluster",
    how="left",
    validate="many_to_one"
)

hmm_long["pct_cluster_proteins"] = (
    100
    * hmm_long["n_proteins_with_HMM"]
    / hmm_long["n_cluster_proteins"]
)

hmm_long["pct_cluster_genomes"] = (
    100
    * hmm_long["n_genomes_with_HMM"]
    / hmm_long["n_cluster_genomes"]
)

hmm_long = hmm_long[
    [
        "module",
        "cluster",
        "HMM",
        "n_proteins_with_HMM",
        "pct_cluster_proteins",
        "n_genomes_with_HMM",
        "pct_cluster_genomes",
        "n_cluster_proteins",
        "n_cluster_genomes",
    ]
]

hmm_long = hmm_long.sort_values(
    [
        "module",
        "cluster",
        "n_proteins_with_HMM",
        "HMM",
    ],
    ascending=[
        True,
        True,
        False,
        True,
    ],
    na_position="last",
    kind="stable"
)

hmm_long.to_csv(
    OUT_LONG,
    sep="\t",
    index=False,
    float_format="%.3f"
)


## ================================================================== ##
## 8. Build mutually exclusive HMM-signature composition
##
## Examples:
##
##     FeGenie_negative
##     MtoA
##     MtoA|MtrA
##
## Unlike the HMM-presence table, every protein contributes to
## exactly ONE signature.
## ================================================================== ##

signature_composition = (
    members
    .groupby(
        ["cluster", "module", "signature"],
        dropna=False,
        as_index=False
    )
    .agg(
        n_proteins=("protein_id", "size"),
        n_genomes=("genome", "nunique"),
    )
)

signature_composition = signature_composition.merge(
    cluster_totals,
    on="cluster",
    how="left",
    validate="many_to_one"
)

signature_composition["pct_cluster_proteins"] = (
    100
    * signature_composition["n_proteins"]
    / signature_composition["n_cluster_proteins"]
)

signature_composition["pct_cluster_genomes"] = (
    100
    * signature_composition["n_genomes"]
    / signature_composition["n_cluster_genomes"]
)

signature_composition = signature_composition[
    [
        "module",
        "cluster",
        "signature",
        "n_proteins",
        "pct_cluster_proteins",
        "n_genomes",
        "pct_cluster_genomes",
        "n_cluster_proteins",
        "n_cluster_genomes",
    ]
]

signature_composition = signature_composition.sort_values(
    [
        "module",
        "cluster",
        "n_proteins",
        "signature",
    ],
    ascending=[
        True,
        True,
        False,
        True,
    ],
    na_position="last",
    kind="stable"
)

signature_composition.to_csv(
    OUT_SIGNATURES,
    sep="\t",
    index=False,
    float_format="%.3f"
)


## ================================================================== ##
## 9. Build cluster-level heterogeneity summary
## ================================================================== ##

summary_records = []

for cluster, group in members.groupby(
    "cluster",
    sort=False
):

    module_values = (
        group["module"]
        .dropna()
        .unique()
        .tolist()
    )

    if len(module_values) > 1:
        fail(
            f"Cluster {cluster} maps to multiple modules."
        )

    module = (
        module_values[0]
        if len(module_values) == 1
        else pd.NA
    )

    n_proteins = len(group)
    n_genomes = group["genome"].nunique()

    positive_group = group[
        group["fegenie_positive_num"] == 1
    ]

    negative_group = group[
        group["fegenie_positive_num"] == 0
    ]

    n_pos = len(positive_group)
    n_neg = len(negative_group)

    ## -------------------------------------------------------------- ##
    ## Individual HMM composition
    ## -------------------------------------------------------------- ##

    hmm_counter = Counter()

    hmm_genomes = {}

    for row in group.itertuples(index=False):

        for hmm in row.hmm_tuple:
            hmm_counter[hmm] += 1

            if hmm not in hmm_genomes:
                hmm_genomes[hmm] = set()

            hmm_genomes[hmm].add(
                row.genome
            )

    dominant_hmms, dominant_hmm_n = (
        dominant_labels(hmm_counter)
    )

    if len(dominant_hmms) == 0:

        dominant_hmm = pd.NA
        dominant_hmm_n_genomes = 0
        dominant_hmm_pct_all = 0.0
        dominant_hmm_pct_positive = float("nan")
        dominant_hmm_pct_genomes = 0.0

    else:

        ## Preserve ties rather than choosing one arbitrarily.
        dominant_hmm = "|".join(
            dominant_hmms
        )

        dominant_hmm_n_genomes = max(
            len(hmm_genomes[hmm])
            for hmm in dominant_hmms
        )

        dominant_hmm_pct_all = (
            100
            * dominant_hmm_n
            / n_proteins
        )

        dominant_hmm_pct_positive = (
            100
            * dominant_hmm_n
            / n_pos
            if n_pos > 0
            else float("nan")
        )

        dominant_hmm_pct_genomes = (
            100
            * dominant_hmm_n_genomes
            / n_genomes
        )

    ## -------------------------------------------------------------- ##
    ## Signature-based heterogeneity
    ## -------------------------------------------------------------- ##

    signature_counter_all = Counter(
        group["signature"]
    )

    signature_counter_positive = Counter(
        positive_group["signature"]
    )

    dominant_signatures, dominant_signature_n = (
        dominant_labels(
            signature_counter_all
        )
    )

    dominant_signature = (
        "|OR|".join(dominant_signatures)
        if len(dominant_signatures) > 0
        else pd.NA
    )

    dominant_signature_pct = (
        100
        * dominant_signature_n
        / n_proteins
        if n_proteins > 0
        else float("nan")
    )

    entropy_all = shannon_entropy(
        signature_counter_all.values()
    )

    norm_entropy_all = normalized_entropy(
        signature_counter_all.values()
    )

    if n_pos > 0:

        entropy_positive = shannon_entropy(
            signature_counter_positive.values()
        )

        norm_entropy_positive = normalized_entropy(
            signature_counter_positive.values()
        )

    else:

        entropy_positive = float("nan")
        norm_entropy_positive = float("nan")

    summary_records.append(
        {
            "module": module,
            "cluster": cluster,
            "n_proteins": n_proteins,
            "n_genomes": n_genomes,

            "n_fegenie_positive": n_pos,
            "n_fegenie_negative": n_neg,

            "pct_fegenie_positive": (
                100 * n_pos / n_proteins
            ),

            "pct_fegenie_negative": (
                100 * n_neg / n_proteins
            ),

            "n_distinct_fegenie_HMMs": len(
                hmm_counter
            ),

            "n_distinct_HMM_signatures_all": len(
                signature_counter_all
            ),

            "n_distinct_HMM_signatures_positive": len(
                signature_counter_positive
            ),

            "dominant_fegenie_HMM": dominant_hmm,

            "dominant_HMM_n_proteins": (
                dominant_hmm_n
            ),

            "dominant_HMM_pct_all_proteins": (
                dominant_hmm_pct_all
            ),

            "dominant_HMM_pct_fegenie_positive_proteins": (
                dominant_hmm_pct_positive
            ),

            "dominant_HMM_n_genomes": (
                dominant_hmm_n_genomes
            ),

            "dominant_HMM_pct_cluster_genomes": (
                dominant_hmm_pct_genomes
            ),

            "dominant_HMM_signature": (
                dominant_signature
            ),

            "dominant_HMM_signature_n_proteins": (
                dominant_signature_n
            ),

            "dominant_HMM_signature_pct_proteins": (
                dominant_signature_pct
            ),

            "HMM_signature_entropy_all": (
                entropy_all
            ),

            "HMM_signature_normalized_entropy_all": (
                norm_entropy_all
            ),

            "HMM_signature_entropy_fegenie_positive": (
                entropy_positive
            ),

            "HMM_signature_normalized_entropy_fegenie_positive": (
                norm_entropy_positive
            ),
        }
    )

composition_summary = pd.DataFrame(
    summary_records
)

composition_summary = composition_summary.sort_values(
    ["module", "cluster"],
    na_position="last",
    kind="stable"
)

composition_summary.to_csv(
    OUT_SUMMARY,
    sep="\t",
    index=False,
    float_format="%.3f"
)


## ================================================================== ##
## 10. Validate dominant HMM against existing cluster_summary
##
## If cluster_summary has a dominant HMM, that HMM must be among
## our tied maximum HMMs.
## ================================================================== ##

old_dominant = old_summary[
    [
        "cluster",
        "dominant_fegenie_HMM"
    ]
].copy()

old_dominant[
    "dominant_fegenie_HMM"
] = old_dominant[
    "dominant_fegenie_HMM"
].apply(clean_string)

new_dominant = composition_summary[
    [
        "cluster",
        "dominant_fegenie_HMM"
    ]
].copy()

new_dominant = new_dominant.rename(
    columns={
        "dominant_fegenie_HMM":
        "dominant_fegenie_HMM_new"
    }
)

dominant_check = old_dominant.merge(
    new_dominant,
    on="cluster",
    how="left",
    validate="one_to_one"
)

dominant_mismatches = []

for row in dominant_check.itertuples(index=False):

    old = clean_string(
        row.dominant_fegenie_HMM
    )

    new = clean_string(
        row.dominant_fegenie_HMM_new
    )

    if old == "" and new == "":
        continue

    if old == "" and new != "":
        dominant_mismatches.append(
            row.cluster
        )
        continue

    if old != "" and new == "":
        dominant_mismatches.append(
            row.cluster
        )
        continue

    new_tied = set(
        new.split("|")
    )

    if old not in new_tied:
        dominant_mismatches.append(
            row.cluster
        )

if len(dominant_mismatches) > 0:

    fail(
        f"{len(dominant_mismatches)} cluster(s) disagree with the "
        "existing dominant_fegenie_HMM assignment.\n"
        "Examples:\n"
        + "\n".join(
            dominant_mismatches[:20]
        )
    )

print(
    "  Existing dominant FeGenie HMM cross-check: PASS"
)


## ================================================================== ##
## 11. Wide HMM matrices
##
## FeGenie-negative proteins are included as their own column.
##
## Because individual positive proteins can have >1 HMM, the
## percentages across HMM columns can exceed 100%.
## ================================================================== ##

matrix_counts = hmm_long.pivot(
    index="cluster",
    columns="HMM",
    values="n_proteins_with_HMM"
).fillna(0)

matrix_counts = (
    matrix_counts
    .astype(int)
    .reset_index()
)

matrix_pct_proteins = hmm_long.pivot(
    index="cluster",
    columns="HMM",
    values="pct_cluster_proteins"
).fillna(0)

matrix_pct_proteins = (
    matrix_pct_proteins
    .reset_index()
)

matrix_pct_genomes = hmm_long.pivot(
    index="cluster",
    columns="HMM",
    values="pct_cluster_genomes"
).fillna(0)

matrix_pct_genomes = (
    matrix_pct_genomes
    .reset_index()
)


## Attach module as first descriptive column ##

matrix_counts = module_lookup.merge(
    matrix_counts,
    on="cluster",
    how="right",
    validate="one_to_one"
)

matrix_pct_proteins = module_lookup.merge(
    matrix_pct_proteins,
    on="cluster",
    how="right",
    validate="one_to_one"
)

matrix_pct_genomes = module_lookup.merge(
    matrix_pct_genomes,
    on="cluster",
    how="right",
    validate="one_to_one"
)


## Sort module-assigned clusters first ##

matrix_counts = matrix_counts.sort_values(
    ["module", "cluster"],
    na_position="last",
    kind="stable"
)

matrix_pct_proteins = (
    matrix_pct_proteins
    .sort_values(
        ["module", "cluster"],
        na_position="last",
        kind="stable"
    )
)

matrix_pct_genomes = (
    matrix_pct_genomes
    .sort_values(
        ["module", "cluster"],
        na_position="last",
        kind="stable"
    )
)

matrix_counts.to_csv(
    OUT_MATRIX_COUNTS,
    sep="\t",
    index=False
)

matrix_pct_proteins.to_csv(
    OUT_MATRIX_PCT_PROTEINS,
    sep="\t",
    index=False,
    float_format="%.3f"
)

matrix_pct_genomes.to_csv(
    OUT_MATRIX_PCT_GENOMES,
    sep="\t",
    index=False,
    float_format="%.3f"
)


## ================================================================== ##
## 12. Network-module-only summary
##
## This is the natural input for a future heatmap grouped by
## Module_01 ... Module_35.
## ================================================================== ##

network_summary = composition_summary[
    composition_summary["module"].notna()
].copy()

if (
    network_summary["cluster"].nunique()
    != EXPECTED_MODULE_CLUSTERS
):
    fail(
        "Unexpected number of module-assigned clusters in "
        "network summary."
    )

network_summary.to_csv(
    OUT_MODULE_NETWORK_SUMMARY,
    sep="\t",
    index=False,
    float_format="%.3f"
)


## ================================================================== ##
## 13. QC summary
## ================================================================== ##

all_negative_clusters = int(
    (
        composition_summary[
            "n_fegenie_positive"
        ] == 0
    ).sum()
)

positive_clusters = int(
    (
        composition_summary[
            "n_fegenie_positive"
        ] > 0
    ).sum()
)

multi_hmm_clusters = int(
    (
        composition_summary[
            "n_distinct_fegenie_HMMs"
        ] > 1
    ).sum()
)

multi_signature_positive_clusters = int(
    (
        composition_summary[
            "n_distinct_HMM_signatures_positive"
        ] > 1
    ).sum()
)

network_positive_clusters = int(
    (
        network_summary[
            "n_fegenie_positive"
        ] > 0
    ).sum()
)

qc = pd.DataFrame(
    [
        ["clustered_proteins", n_clustered_proteins],
        ["mmseqs_clusters", n_clusters],
        ["fegenie_positive_proteins", n_positive],
        ["fegenie_negative_proteins", n_negative],
        ["proteins_with_multiple_HMMs", multi_hmm_proteins],
        ["distinct_Fegenie_HMMs", len(unique_hmms)],
        ["clusters_all_fegenie_negative", all_negative_clusters],
        ["clusters_with_any_fegenie_positive", positive_clusters],
        ["clusters_with_multiple_distinct_HMMs", multi_hmm_clusters],
        [
            "clusters_with_multiple_positive_HMM_signatures",
            multi_signature_positive_clusters
        ],
        ["module_assigned_clusters", EXPECTED_MODULE_CLUSTERS],
        ["modules", EXPECTED_MODULES],
        [
            "module_clusters_with_any_fegenie_positive",
            network_positive_clusters
        ],
    ],
    columns=[
        "metric",
        "value"
    ]
)

qc.to_csv(
    OUT_QC,
    sep="\t",
    index=False
)


## ================================================================== ##
## 14. Final report
## ================================================================== ##

print()
print("Cluster-level FeGenie composition")
print(
    f"  All-FeGenie-negative clusters:          "
    f"{all_negative_clusters:,}"
)
print(
    f"  Clusters with >=1 FeGenie-positive:     "
    f"{positive_clusters:,}"
)
print(
    f"  Clusters with >1 distinct HMM:          "
    f"{multi_hmm_clusters:,}"
)
print(
    f"  Clusters with >1 positive signature:    "
    f"{multi_signature_positive_clusters:,}"
)

print()
print("Final-network subset")
print(
    f"  MCL module clusters:                     "
    f"{len(network_summary):,}"
)
print(
    f"  Module clusters with FeGenie positives:  "
    f"{network_positive_clusters:,}"
)

print()
print("=" * 78)
print("SUCCESS")
print("=" * 78)

print(f"HMM long table              : {OUT_LONG}")
print(f"Signature composition       : {OUT_SIGNATURES}")
print(f"Composition summary         : {OUT_SUMMARY}")
print(f"HMM count matrix            : {OUT_MATRIX_COUNTS}")
print(f"HMM protein-% matrix        : {OUT_MATRIX_PCT_PROTEINS}")
print(f"HMM genome-% matrix         : {OUT_MATRIX_PCT_GENOMES}")
print(f"Network/module subset       : {OUT_MODULE_NETWORK_SUMMARY}")
print(f"QC summary                  : {OUT_QC}")
