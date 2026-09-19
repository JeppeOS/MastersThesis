#!/usr/bin/env python3

from pathlib import Path
import sys

import numpy as np
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

MODULE_MEMBERSHIP = (
    WORKFLOW
    / "07_cooccurrence_network"
    / "module_membership.tsv"
)

MODULE_SUMMARY_07 = (
    WORKFLOW
    / "07_cooccurrence_network"
    / "module_summary.tsv"
)


## ================================================================== ##
## Outputs
## ================================================================== ##

OUT_PROTEIN_MEMBERSHIP = (
    HERE
    / "protein_module_membership.tsv"
)

OUT_CLUSTER_PRESENCE = (
    HERE
    / "genome_module_cluster_presence.tsv"
)

OUT_GENOME_MODULE = (
    HERE
    / "genome_module_summary.tsv"
)

OUT_GENOME_MODULE_ALL = (
    HERE
    / "genome_module_summary_all.tsv"
)

OUT_COMPLETENESS_MATRIX = (
    HERE
    / "genome_module_completeness_matrix.tsv"
)

OUT_CLUSTER_COUNT_MATRIX = (
    HERE
    / "genome_module_cluster_count_matrix.tsv"
)

OUT_MODULE_SUMMARY = (
    HERE
    / "module_completeness_summary.tsv"
)

OUT_QC = (
    HERE
    / "module_occurrence_qc.tsv"
)


## ================================================================== ##
## Locked expectations from the finalized upstream analysis
## ================================================================== ##

EXPECTED_TOTAL_CLUSTERED_PROTEINS = 14308
EXPECTED_TOTAL_MMSEQS_CLUSTERS = 1390

EXPECTED_NETWORK_CLUSTERS = 156
EXPECTED_MODULES = 35
EXPECTED_GENOMES = 631


## ================================================================== ##
## Utilities
## ================================================================== ##

def fail(message):
    print(
        f"\nERROR: {message}",
        file=sys.stderr
    )
    sys.exit(1)


def require_columns(df, required, filename):
    missing = set(required) - set(df.columns)

    if missing:
        fail(
            f"{filename} is missing required column(s): "
            + ", ".join(sorted(missing))
        )


## ================================================================== ##
## 1. Read upstream tables
## ================================================================== ##

print("=" * 80)
print("BUILD MODULE OCCURRENCE AND COMPLETENESS")
print("=" * 80)

for path in [
    CLUSTER_MEMBERSHIP,
    MODULE_MEMBERSHIP,
    MODULE_SUMMARY_07,
]:
    if not path.exists():
        fail(
            f"Cannot find required input:\n{path}"
        )


cluster_members = pd.read_csv(
    CLUSTER_MEMBERSHIP,
    sep="\t",
    dtype=str
)

module_members = pd.read_csv(
    MODULE_MEMBERSHIP,
    sep="\t",
    dtype=str
)

old_module_summary = pd.read_csv(
    MODULE_SUMMARY_07,
    sep="\t",
    dtype=str
)


## ================================================================== ##
## 2. Validate cluster_membership.tsv
## ================================================================== ##

require_columns(
    cluster_members,
    [
        "cluster",
        "genome",
        "protein_id",
    ],
    CLUSTER_MEMBERSHIP.name
)

cluster_members["cluster"] = (
    cluster_members["cluster"]
    .astype(str)
    .str.strip()
)

cluster_members["genome"] = (
    cluster_members["genome"]
    .astype(str)
    .str.strip()
)

cluster_members["protein_id"] = (
    cluster_members["protein_id"]
    .astype(str)
    .str.strip()
)


## genome + protein_id is authoritative ##
duplicate_proteins = cluster_members.duplicated(
    ["genome", "protein_id"],
    keep=False
)

if duplicate_proteins.any():

    examples = (
        cluster_members.loc[
            duplicate_proteins,
            [
                "genome",
                "protein_id",
                "cluster"
            ]
        ]
        .head(20)
    )

    fail(
        "Duplicate genome + protein_id keys found:\n"
        + examples.to_string(index=False)
    )


n_clustered_proteins = len(
    cluster_members
)

n_mmseqs_clusters = (
    cluster_members["cluster"]
    .nunique()
)

n_all_genomes = (
    cluster_members["genome"]
    .nunique()
)

print()
print("Cluster membership input")
print(
    f"  Protein rows:            "
    f"{n_clustered_proteins:,}"
)
print(
    f"  MMseqs clusters:         "
    f"{n_mmseqs_clusters:,}"
)
print(
    f"  Genomes represented:     "
    f"{n_all_genomes:,}"
)


if (
    n_clustered_proteins
    != EXPECTED_TOTAL_CLUSTERED_PROTEINS
):
    fail(
        f"Expected "
        f"{EXPECTED_TOTAL_CLUSTERED_PROTEINS:,} "
        f"clustered proteins but found "
        f"{n_clustered_proteins:,}."
    )


if (
    n_mmseqs_clusters
    != EXPECTED_TOTAL_MMSEQS_CLUSTERS
):
    fail(
        f"Expected "
        f"{EXPECTED_TOTAL_MMSEQS_CLUSTERS:,} "
        f"MMseqs clusters but found "
        f"{n_mmseqs_clusters:,}."
    )


if n_all_genomes != EXPECTED_GENOMES:
    fail(
        f"Expected {EXPECTED_GENOMES} genomes "
        f"but found {n_all_genomes}."
    )


## ================================================================== ##
## 3. Validate module_membership.tsv
## ================================================================== ##

require_columns(
    module_members,
    [
        "cluster",
        "module",
    ],
    MODULE_MEMBERSHIP.name
)

module_members["cluster"] = (
    module_members["cluster"]
    .astype(str)
    .str.strip()
)

module_members["module"] = (
    module_members["module"]
    .astype(str)
    .str.strip()
)


if module_members["cluster"].duplicated().any():

    duplicate_clusters = (
        module_members.loc[
            module_members[
                "cluster"
            ].duplicated(
                keep=False
            ),
            "cluster"
        ]
        .unique()
    )

    fail(
        "A protein-family cluster has more than one "
        "module assignment:\n"
        + "\n".join(
            duplicate_clusters[:20]
        )
    )


n_network_clusters = (
    module_members["cluster"]
    .nunique()
)

n_modules = (
    module_members["module"]
    .nunique()
)


print()
print("Module membership input")
print(
    f"  Network clusters:        "
    f"{n_network_clusters:,}"
)
print(
    f"  MCL modules:             "
    f"{n_modules:,}"
)


if (
    n_network_clusters
    != EXPECTED_NETWORK_CLUSTERS
):
    fail(
        f"Expected {EXPECTED_NETWORK_CLUSTERS} "
        f"module-assigned clusters but found "
        f"{n_network_clusters}."
    )


if n_modules != EXPECTED_MODULES:
    fail(
        f"Expected {EXPECTED_MODULES} modules "
        f"but found {n_modules}."
    )


## Every module cluster must exist in cluster_membership ##

known_clusters = set(
    cluster_members["cluster"]
)

unknown_module_clusters = sorted(
    set(
        module_members["cluster"]
    )
    - known_clusters
)

if unknown_module_clusters:
    fail(
        "module_membership.tsv contains cluster(s) "
        "not present in cluster_membership.tsv:\n"
        + "\n".join(
            unknown_module_clusters[:20]
        )
    )


## ================================================================== ##
## 4. Calculate module sizes
##
## Module size means number of distinct MMseqs protein-family
## clusters belonging to the MCL module.
## ================================================================== ##

module_sizes = (
    module_members
    .groupby(
        "module",
        as_index=False
    )
    .agg(
        module_size_clusters=(
            "cluster",
            "nunique"
        )
    )
)

module_size_lookup = dict(
    zip(
        module_sizes["module"],
        module_sizes[
            "module_size_clusters"
        ]
    )
)


## ================================================================== ##
## 5. Build protein-level module membership
##
## Only proteins whose MMseqs cluster is one of the 156
## MCL-assigned network clusters are retained.
## ================================================================== ##

protein_module = (
    cluster_members
    .merge(
        module_members,
        on="cluster",
        how="inner",
        validate="many_to_one"
    )
)


if len(protein_module) == 0:
    fail(
        "No clustered proteins joined to module membership."
    )


## Keep authoritative identifiers first. Preserve selected useful
## annotation columns if they already exist upstream. ##

protein_output_columns = [
    "genome",
    "protein_id",
    "cluster",
    "module",
]

optional_annotation_columns = [
    "representative_protein",
    "length",
    "fegenie_positive",
    "fegenie_HMMs",
    "findmehemes_positive",
    "number_of_hemes",
    "signalp_prediction",
    "deeptmhmm_class",
    "deeptmhmm_n_tm_helices",
    "localization_class",
]

protein_output_columns += [
    col
    for col in optional_annotation_columns
    if col in protein_module.columns
]


protein_module = (
    protein_module[
        protein_output_columns
    ]
    .sort_values(
        [
            "module",
            "genome",
            "cluster",
            "protein_id",
        ],
        kind="stable"
    )
    .reset_index(drop=True)
)


protein_module.to_csv(
    OUT_PROTEIN_MEMBERSHIP,
    sep="\t",
    index=False
)


n_module_proteins = len(
    protein_module
)

n_module_genomes = (
    protein_module["genome"]
    .nunique()
)


print()
print("Protein-level module linkage")
print(
    f"  Proteins in MCL modules: "
    f"{n_module_proteins:,}"
)
print(
    f"  Genomes represented:     "
    f"{n_module_genomes:,}"
)


## Because Module_01 occurred in all 631 genomes upstream,
## every genome should be represented here. ##

if n_module_genomes != EXPECTED_GENOMES:
    fail(
        f"Expected module proteins in all "
        f"{EXPECTED_GENOMES} genomes, but only "
        f"{n_module_genomes} genomes are represented."
    )


## ================================================================== ##
## 6. Genome × module × cluster presence
##
## One row = one distinct cluster from a module present in a genome.
##
## n_proteins_in_cluster > 1 identifies within-genome paralogues
## for that sequence family.
## ================================================================== ##

cluster_presence = (
    protein_module
    .groupby(
        [
            "genome",
            "module",
            "cluster",
        ],
        as_index=False
    )
    .agg(
        n_proteins_in_cluster=(
            "protein_id",
            "nunique"
        )
    )
)


cluster_presence[
    "has_cluster_paralogue"
] = (
    cluster_presence[
        "n_proteins_in_cluster"
    ] > 1
).astype(int)


cluster_presence[
    "module_size_clusters"
] = (
    cluster_presence["module"]
    .map(
        module_size_lookup
    )
    .astype(int)
)


cluster_presence = (
    cluster_presence[
        [
            "genome",
            "module",
            "cluster",
            "module_size_clusters",
            "n_proteins_in_cluster",
            "has_cluster_paralogue",
        ]
    ]
    .sort_values(
        [
            "module",
            "genome",
            "cluster",
        ],
        kind="stable"
    )
    .reset_index(drop=True)
)


cluster_presence.to_csv(
    OUT_CLUSTER_PRESENCE,
    sep="\t",
    index=False
)


## ================================================================== ##
## 7. Genome × module summary
##
## IMPORTANT:
##
## module_completeness is based on DISTINCT cluster families,
## not protein-copy count.
##
## paralogue excess is:
##
##     n_module_proteins - n_clusters_present
##
## so a module with one copy of each member has paralogue_excess = 0.
## ================================================================== ##

genome_module = (
    cluster_presence
    .groupby(
        [
            "genome",
            "module",
            "module_size_clusters",
        ],
        as_index=False
    )
    .agg(
        n_clusters_present=(
            "cluster",
            "nunique"
        ),
        n_module_proteins=(
            "n_proteins_in_cluster",
            "sum"
        ),
        n_clusters_with_paralogues=(
            "has_cluster_paralogue",
            "sum"
        ),
    )
)


genome_module[
    "module_completeness"
] = (
    genome_module[
        "n_clusters_present"
    ]
    /
    genome_module[
        "module_size_clusters"
    ]
)


genome_module[
    "module_completeness_pct"
] = (
    100
    * genome_module[
        "module_completeness"
    ]
)


genome_module[
    "n_paralog_excess"
] = (
    genome_module[
        "n_module_proteins"
    ]
    -
    genome_module[
        "n_clusters_present"
    ]
)


genome_module[
    "has_module_paralogues"
] = (
    genome_module[
        "n_paralog_excess"
    ] > 0
).astype(int)


genome_module[
    "is_full_module"
] = (
    genome_module[
        "n_clusters_present"
    ]
    ==
    genome_module[
        "module_size_clusters"
    ]
).astype(int)


genome_module = (
    genome_module[
        [
            "genome",
            "module",
            "module_size_clusters",
            "n_clusters_present",
            "module_completeness",
            "module_completeness_pct",
            "n_module_proteins",
            "n_paralog_excess",
            "n_clusters_with_paralogues",
            "has_module_paralogues",
            "is_full_module",
        ]
    ]
    .sort_values(
        [
            "module",
            "genome",
        ],
        kind="stable"
    )
    .reset_index(drop=True)
)


genome_module.to_csv(
    OUT_GENOME_MODULE,
    sep="\t",
    index=False,
    float_format="%.6f"
)


## ================================================================== ##
## 8. Build complete 631 genomes × 35 modules table
##
## genome_module_summary.tsv contains only observed genome-module
## combinations.
##
## This second table explicitly includes absent modules as zeros,
## which is useful for matrices and visualization.
## ================================================================== ##

all_genomes = (
    pd.DataFrame(
        {
            "genome": sorted(
                cluster_members[
                    "genome"
                ].unique()
            )
        }
    )
)

all_modules = (
    module_sizes[
        [
            "module",
            "module_size_clusters",
        ]
    ]
    .sort_values(
        "module"
    )
    .reset_index(drop=True)
)


all_genomes["_key"] = 1
all_modules["_key"] = 1


full_grid = (
    all_genomes
    .merge(
        all_modules,
        on="_key",
        how="inner"
    )
    .drop(
        columns="_key"
    )
)


full_grid = full_grid.merge(
    genome_module.drop(
        columns=[
            "module_size_clusters"
        ]
    ),
    on=[
        "genome",
        "module",
    ],
    how="left",
    validate="one_to_one"
)


zero_integer_columns = [
    "n_clusters_present",
    "n_module_proteins",
    "n_paralog_excess",
    "n_clusters_with_paralogues",
    "has_module_paralogues",
    "is_full_module",
]


for col in zero_integer_columns:

    full_grid[col] = (
        full_grid[col]
        .fillna(0)
        .astype(int)
    )


full_grid[
    "module_completeness"
] = (
    full_grid[
        "module_completeness"
    ]
    .fillna(0.0)
)


full_grid[
    "module_completeness_pct"
] = (
    full_grid[
        "module_completeness_pct"
    ]
    .fillna(0.0)
)


full_grid[
    "module_present"
] = (
    full_grid[
        "n_clusters_present"
    ] > 0
).astype(int)


full_grid = (
    full_grid[
        [
            "genome",
            "module",
            "module_size_clusters",
            "module_present",
            "n_clusters_present",
            "module_completeness",
            "module_completeness_pct",
            "n_module_proteins",
            "n_paralog_excess",
            "n_clusters_with_paralogues",
            "has_module_paralogues",
            "is_full_module",
        ]
    ]
    .sort_values(
        [
            "genome",
            "module",
        ],
        kind="stable"
    )
    .reset_index(drop=True)
)


expected_full_rows = (
    EXPECTED_GENOMES
    * EXPECTED_MODULES
)


if len(full_grid) != expected_full_rows:
    fail(
        f"Expected {expected_full_rows:,} rows in "
        f"complete genome × module grid, but found "
        f"{len(full_grid):,}."
    )


full_grid.to_csv(
    OUT_GENOME_MODULE_ALL,
    sep="\t",
    index=False,
    float_format="%.6f"
)


## ================================================================== ##
## 9. Wide matrices
##
## Useful later for heatmaps / ordinations.
## ================================================================== ##

completeness_matrix = (
    full_grid
    .pivot(
        index="genome",
        columns="module",
        values="module_completeness"
    )
    .reset_index()
)


completeness_matrix.to_csv(
    OUT_COMPLETENESS_MATRIX,
    sep="\t",
    index=False,
    float_format="%.6f"
)


cluster_count_matrix = (
    full_grid
    .pivot(
        index="genome",
        columns="module",
        values="n_clusters_present"
    )
    .reset_index()
)


cluster_count_matrix.to_csv(
    OUT_CLUSTER_COUNT_MATRIX,
    sep="\t",
    index=False
)


## ================================================================== ##
## 10. Module-level completeness summary
##
## Statistics among PRESENT genomes are kept separate from statistics
## across all 631 genomes.
## ================================================================== ##

present_summary = (
    genome_module
    .groupby(
        [
            "module",
            "module_size_clusters",
        ],
        as_index=False
    )
    .agg(
        n_genomes_present=(
            "genome",
            "nunique"
        ),

        n_genomes_full_module=(
            "is_full_module",
            "sum"
        ),

        mean_completeness_present=(
            "module_completeness",
            "mean"
        ),

        median_completeness_present=(
            "module_completeness",
            "median"
        ),

        min_completeness_present=(
            "module_completeness",
            "min"
        ),

        max_completeness_present=(
            "module_completeness",
            "max"
        ),

        mean_clusters_present=(
            "n_clusters_present",
            "mean"
        ),

        median_clusters_present=(
            "n_clusters_present",
            "median"
        ),

        max_clusters_present=(
            "n_clusters_present",
            "max"
        ),

        total_module_proteins=(
            "n_module_proteins",
            "sum"
        ),

        n_genomes_with_module_paralogues=(
            "has_module_paralogues",
            "sum"
        ),
    )
)


present_summary[
    "pct_all_genomes_present"
] = (
    100
    * present_summary[
        "n_genomes_present"
    ]
    / EXPECTED_GENOMES
)


present_summary[
    "pct_present_genomes_full_module"
] = (
    100
    * present_summary[
        "n_genomes_full_module"
    ]
    /
    present_summary[
        "n_genomes_present"
    ]
)


present_summary[
    "pct_present_genomes_with_paralogues"
] = (
    100
    * present_summary[
        "n_genomes_with_module_paralogues"
    ]
    /
    present_summary[
        "n_genomes_present"
    ]
)


all_summary = (
    full_grid
    .groupby(
        "module",
        as_index=False
    )
    .agg(
        mean_completeness_all_631=(
            "module_completeness",
            "mean"
        ),
        median_completeness_all_631=(
            "module_completeness",
            "median"
        ),
    )
)


module_summary = (
    present_summary
    .merge(
        all_summary,
        on="module",
        how="left",
        validate="one_to_one"
    )
)


## ================================================================== ##
## 11. Cross-check against Stage 07 module_summary.tsv
##
## This validates:
##
##     module size
##     number of genomes containing >=1 member cluster
##
## against the independently generated Stage 07 summary.
## ================================================================== ##

require_columns(
    old_module_summary,
    [
        "module",
        "module_size_clusters",
        "n_genomes_with_module",
    ],
    MODULE_SUMMARY_07.name
)


old_check = (
    old_module_summary[
        [
            "module",
            "module_size_clusters",
            "n_genomes_with_module",
        ]
    ]
    .copy()
)


old_check[
    "module_size_clusters"
] = pd.to_numeric(
    old_check[
        "module_size_clusters"
    ],
    errors="raise"
).astype(int)


old_check[
    "n_genomes_with_module"
] = pd.to_numeric(
    old_check[
        "n_genomes_with_module"
    ],
    errors="raise"
).astype(int)


check = (
    module_summary[
        [
            "module",
            "module_size_clusters",
            "n_genomes_present",
        ]
    ]
    .merge(
        old_check,
        on="module",
        how="outer",
        suffixes=(
            "_new",
            "_old"
        ),
        indicator=True
    )
)


bad_module_summary = check[
    (check["_merge"] != "both")
    |
    (
        check[
            "module_size_clusters_new"
        ]
        !=
        check[
            "module_size_clusters_old"
        ]
    )
    |
    (
        check[
            "n_genomes_present"
        ]
        !=
        check[
            "n_genomes_with_module"
        ]
    )
]


if len(bad_module_summary) > 0:
    fail(
        "Stage 09 module counts disagree with "
        "07_cooccurrence_network/module_summary.tsv:\n"
        + bad_module_summary.to_string(
            index=False
        )
    )


print()
print(
    "Stage 07 module-summary cross-check: PASS"
)


## ================================================================== ##
## 12. Sort and write module completeness summary
## ================================================================== ##

module_summary = (
    module_summary
    .sort_values(
        "module"
    )
    .reset_index(drop=True)
)


module_summary.to_csv(
    OUT_MODULE_SUMMARY,
    sep="\t",
    index=False,
    float_format="%.6f"
)


## ================================================================== ##
## 13. Additional strict QC
## ================================================================== ##

## Completeness must always be between 0 and 1. ##

if (
    (
        full_grid[
            "module_completeness"
        ] < 0
    ).any()
    or
    (
        full_grid[
            "module_completeness"
        ] > 1
    ).any()
):
    fail(
        "Module completeness outside the interval [0,1]."
    )


## Observed module rows may never have zero clusters. ##

if (
    genome_module[
        "n_clusters_present"
    ] < 1
).any():
    fail(
        "Observed genome-module table contains "
        "a zero-presence row."
    )


## n_module_proteins must be >= n distinct clusters. ##

if (
    genome_module[
        "n_module_proteins"
    ]
    <
    genome_module[
        "n_clusters_present"
    ]
).any():
    fail(
        "Protein count is smaller than distinct "
        "cluster count in a genome-module combination."
    )


## ================================================================== ##
## 14. QC file
## ================================================================== ##

n_observed_genome_module_pairs = len(
    genome_module
)

n_full_genome_module_pairs = int(
    genome_module[
        "is_full_module"
    ].sum()
)

n_genome_module_pairs_with_paralogues = int(
    genome_module[
        "has_module_paralogues"
    ].sum()
)

n_cluster_occurrences_with_paralogues = int(
    cluster_presence[
        "has_cluster_paralogue"
    ].sum()
)


qc = pd.DataFrame(
    [
        [
            "input_clustered_proteins",
            n_clustered_proteins,
        ],
        [
            "input_mmseqs_clusters",
            n_mmseqs_clusters,
        ],
        [
            "input_genomes",
            n_all_genomes,
        ],
        [
            "module_assigned_clusters",
            n_network_clusters,
        ],
        [
            "modules",
            n_modules,
        ],
        [
            "proteins_in_modules",
            n_module_proteins,
        ],
        [
            "genomes_with_module_proteins",
            n_module_genomes,
        ],
        [
            "observed_genome_module_pairs",
            n_observed_genome_module_pairs,
        ],
        [
            "full_631x35_grid_rows",
            len(full_grid),
        ],
        [
            "genome_module_pairs_with_all_member_clusters",
            n_full_genome_module_pairs,
        ],
        [
            "genome_module_pairs_with_paralogues",
            n_genome_module_pairs_with_paralogues,
        ],
        [
            "genome_module_cluster_occurrences_with_paralogues",
            n_cluster_occurrences_with_paralogues,
        ],
    ],
    columns=[
        "metric",
        "value",
    ]
)


qc.to_csv(
    OUT_QC,
    sep="\t",
    index=False
)


## ================================================================== ##
## 15. Terminal report
## ================================================================== ##

print()
print("Module completeness summary")
print(
    module_summary[
        [
            "module",
            "module_size_clusters",
            "n_genomes_present",
            "pct_all_genomes_present",
            "mean_completeness_present",
            "median_completeness_present",
            "n_genomes_full_module",
            "pct_present_genomes_full_module",
            "n_genomes_with_module_paralogues",
        ]
    ]
    .to_string(
        index=False
    )
)


print()
print("Overall occurrence QC")
print(
    f"  Proteins assigned to modules:        "
    f"{n_module_proteins:,}"
)
print(
    f"  Observed genome-module pairs:        "
    f"{n_observed_genome_module_pairs:,}"
)
print(
    f"  Genome-module pairs fully present:   "
    f"{n_full_genome_module_pairs:,}"
)
print(
    f"  Genome-module pairs with paralogues: "
    f"{n_genome_module_pairs_with_paralogues:,}"
)
print(
    f"  Cluster occurrences with paralogues: "
    f"{n_cluster_occurrences_with_paralogues:,}"
)


print()
print("=" * 80)
print("SUCCESS")
print("=" * 80)

print(
    f"Protein module membership : "
    f"{OUT_PROTEIN_MEMBERSHIP}"
)
print(
    f"Genome-module-cluster table: "
    f"{OUT_CLUSTER_PRESENCE}"
)
print(
    f"Observed genome-module table: "
    f"{OUT_GENOME_MODULE}"
)
print(
    f"Complete genome-module grid: "
    f"{OUT_GENOME_MODULE_ALL}"
)
print(
    f"Completeness matrix: "
    f"{OUT_COMPLETENESS_MATRIX}"
)
print(
    f"Cluster-count matrix: "
    f"{OUT_CLUSTER_COUNT_MATRIX}"
)
print(
    f"Module completeness summary: "
    f"{OUT_MODULE_SUMMARY}"
)
print(
    f"QC summary: "
    f"{OUT_QC}"
)
