#!/usr/bin/env python3

from pathlib import Path
from itertools import combinations
import sys

import numpy as np
import pandas as pd


## ================================================================== ##
## STAGE 13A
##
## Summarize physical colocalization of protein-family clusters
## belonging to the same MCL-defined co-occurrence module.
##
## Primary analytical unit:
##
##   one genome
##   x one module
##   x one unordered pair of distinct member clusters
##
## Each such observation means that BOTH cluster families occur
## somewhere in that genome.
##
## Questions:
##
##   - Do they occur on the same assembled contig?
##   - If yes, are any copies within 5 / 10 / 20 kb?
##   - Are any copies adjacent?
##   - What is the closest observed separation?
##   - What is the orientation of the closest same-contig pair?
##
## IMPORTANT:
##
## Different-contig observations are NOT interpreted as physically
## distant. Their true genomic separation is unknown because the
## assembly is fragmented.
##
## Distances at 5/10/20 kb refer to INTERGENIC GAP between CDS
## boundaries, as defined in Stage 12A.
## ================================================================== ##


## ================================================================== ##
## Paths
## ================================================================== ##

HERE = Path(__file__).resolve().parent
WORKFLOW = HERE.parent.parent


PAIR_OBSERVATIONS = (
    WORKFLOW
    / "12_neighborhood_extraction"
    / "module_member_pair_observations.tsv"
)


PROTEIN_MODULE_MEMBERSHIP = (
    WORKFLOW
    / "09_module_occurrence"
    / "protein_module_membership.tsv"
)


MODULE_MEMBERSHIP = (
    WORKFLOW
    / "07_cooccurrence_network"
    / "module_membership.tsv"
)


CLUSTER_SUMMARY = (
    WORKFLOW
    / "06_clustering"
    / "cluster_summary.tsv"
)


## ================================================================== ##
## Outputs
## ================================================================== ##

OUT_PAIR_SUMMARY = (
    HERE
    / "module_cluster_pair_colocalization.tsv"
)


OUT_MODULE_SUMMARY = (
    HERE
    / "module_colocalization_summary.tsv"
)


OUT_QC = (
    HERE
    / "module_colocalization_qc.tsv"
)


## ================================================================== ##
## Expected upstream values
## ================================================================== ##

EXPECTED_MODULE_PROTEINS = 10_537
EXPECTED_MODULE_CLUSTERS = 156
EXPECTED_MODULES = 35

EXPECTED_PAIR_OBSERVATIONS = 34_037

## From the authoritative 35-module membership:
##
## sum choose(n_clusters_in_module, 2)
##
## = 560 possible unordered within-module cluster pairs.
##
EXPECTED_POSSIBLE_CLUSTER_PAIRS = 560


## ================================================================== ##
## Utilities
## ================================================================== ##

def fail(message):

    print(
        f"\nERROR: {message}",
        file=sys.stderr,
    )

    sys.exit(1)


def require_columns(
    df,
    columns,
    source,
):

    missing = (
        set(columns)
        -
        set(df.columns)
    )

    if missing:

        fail(
            f"{source} is missing required "
            f"column(s):\n"
            +
            ", ".join(
                sorted(missing)
            )
        )


def pct(
    numerator,
    denominator,
):

    if denominator == 0:
        return np.nan

    return (
        100.0
        *
        numerator
        /
        denominator
    )


def numeric_median(
    values,
):

    values = pd.to_numeric(
        values,
        errors="coerce",
    ).dropna()

    if len(values) == 0:
        return np.nan

    return float(
        values.median()
    )


def numeric_quantile(
    values,
    q,
):

    values = pd.to_numeric(
        values,
        errors="coerce",
    ).dropna()

    if len(values) == 0:
        return np.nan

    return float(
        values.quantile(q)
    )


## ================================================================== ##
## Start
## ================================================================== ##

print("=" * 80)
print("STAGE 13A - MODULE COLOCALIZATION")
print("=" * 80)


## ================================================================== ##
## 1. Read authoritative MCL module membership
## ================================================================== ##

print()
print(
    "Reading MCL module membership..."
)


modules = pd.read_csv(
    MODULE_MEMBERSHIP,
    sep="\t",
    dtype=str,
    keep_default_na=False,
)


require_columns(
    modules,
    [
        "cluster",
        "module",
    ],
    MODULE_MEMBERSHIP.name,
)


modules = modules[
    [
        "cluster",
        "module",
    ]
].copy()


if (
    modules[
        "cluster"
    ]
    .duplicated()
    .any()
):

    fail(
        "A cluster is assigned to more than "
        "one MCL module."
    )


if (
    modules[
        "cluster"
    ]
    .nunique()
    !=
    EXPECTED_MODULE_CLUSTERS
):

    fail(
        f"Expected "
        f"{EXPECTED_MODULE_CLUSTERS} "
        f"module clusters but found "
        f"{modules['cluster'].nunique()}."
    )


if (
    modules[
        "module"
    ]
    .nunique()
    !=
    EXPECTED_MODULES
):

    fail(
        f"Expected "
        f"{EXPECTED_MODULES} modules but found "
        f"{modules['module'].nunique()}."
    )


print(
    f"  Module clusters: "
    f"{modules['cluster'].nunique():,}"
)

print(
    f"  Modules:         "
    f"{modules['module'].nunique():,}"
)


## ================================================================== ##
## 2. Construct ALL possible within-module cluster pairs
##
## This is important because two clusters can belong to the same MCL
## module without ever being observed together in the same genome.
##
## Such a pair should remain in the final table with:
##
##     n_genomes_both = 0
##
## rather than disappearing from the analysis.
## ================================================================== ##

possible_pair_rows = []


for module, group in modules.groupby(
    "module",
    sort=True,
):

    clusters = sorted(
        group[
            "cluster"
        ].unique()
    )


    for (
        cluster_a,
        cluster_b,
    ) in combinations(
        clusters,
        2,
    ):

        possible_pair_rows.append(
            {
                "module":
                    module,

                "cluster_a":
                    cluster_a,

                "cluster_b":
                    cluster_b,
            }
        )


possible_pairs = pd.DataFrame(
    possible_pair_rows
)


if (
    len(
        possible_pairs
    )
    !=
    EXPECTED_POSSIBLE_CLUSTER_PAIRS
):

    fail(
        f"Expected "
        f"{EXPECTED_POSSIBLE_CLUSTER_PAIRS} "
        f"possible within-module cluster pairs "
        f"but found "
        f"{len(possible_pairs)}."
    )


print(
    f"  Possible within-module "
    f"cluster pairs: "
    f"{len(possible_pairs):,}"
)


## ================================================================== ##
## 3. Read protein-level module membership
##
## We independently reconstruct every genome x module x cluster
## presence and therefore every expected genome-level co-present
## cluster pair.
## ================================================================== ##

print()
print(
    "Reading protein-level module membership..."
)


proteins = pd.read_csv(
    PROTEIN_MODULE_MEMBERSHIP,
    sep="\t",
    dtype=str,
    keep_default_na=False,
)


require_columns(
    proteins,
    [
        "genome",
        "protein_id",
        "cluster",
        "module",
    ],
    PROTEIN_MODULE_MEMBERSHIP.name,
)


if len(proteins) != EXPECTED_MODULE_PROTEINS:

    fail(
        f"Expected "
        f"{EXPECTED_MODULE_PROTEINS:,} "
        f"module proteins but found "
        f"{len(proteins):,}."
    )


if (
    proteins[
        [
            "genome",
            "protein_id",
        ]
    ]
    .duplicated()
    .any()
):

    fail(
        "Duplicate genome + protein_id "
        "keys in protein module membership."
    )


if (
    proteins[
        "cluster"
    ]
    .nunique()
    !=
    EXPECTED_MODULE_CLUSTERS
):

    fail(
        "Protein membership does not recover "
        "all 156 module clusters."
    )


if (
    proteins[
        "module"
    ]
    .nunique()
    !=
    EXPECTED_MODULES
):

    fail(
        "Protein membership does not recover "
        "all 35 modules."
    )


## ================================================================== ##
## Verify protein-level module labels against authoritative
## cluster -> module membership.
## ================================================================== ##

protein_module_check = proteins[
    [
        "cluster",
        "module",
    ]
].drop_duplicates()


protein_module_check = (
    protein_module_check
    .merge(
        modules.rename(
            columns={
                "module":
                    "expected_module",
            }
        ),
        on="cluster",
        how="left",
        validate="one_to_one",
    )
)


bad_module_labels = protein_module_check[
    protein_module_check[
        "module"
    ]
    !=
    protein_module_check[
        "expected_module"
    ]
]


if len(bad_module_labels) > 0:

    fail(
        "Protein-level module labels disagree "
        "with authoritative module membership."
    )


## ================================================================== ##
## Genome x module x cluster presence
## ================================================================== ##

presence = (
    proteins[
        [
            "genome",
            "module",
            "cluster",
        ]
    ]
    .drop_duplicates()
)


## Protein multiplicity per genome/module/cluster is retained for QC. ##

protein_counts = (
    proteins
    .groupby(
        [
            "genome",
            "module",
            "cluster",
        ],
        as_index=False,
    )
    .size()
    .rename(
        columns={
            "size":
                "expected_n_proteins",
        }
    )
)


## ================================================================== ##
## Independently reconstruct expected genome-level pair observations
## ================================================================== ##

expected_observation_rows = []


for (
    genome,
    module,
), group in presence.groupby(
    [
        "genome",
        "module",
    ],
    sort=False,
):

    clusters = sorted(
        group[
            "cluster"
        ].unique()
    )


    for (
        cluster_a,
        cluster_b,
    ) in combinations(
        clusters,
        2,
    ):

        expected_observation_rows.append(
            {
                "genome":
                    genome,

                "module":
                    module,

                "cluster_a":
                    cluster_a,

                "cluster_b":
                    cluster_b,
            }
        )


expected_observations = pd.DataFrame(
    expected_observation_rows
)


if (
    len(
        expected_observations
    )
    !=
    EXPECTED_PAIR_OBSERVATIONS
):

    fail(
        f"Independent reconstruction expected "
        f"{EXPECTED_PAIR_OBSERVATIONS:,} "
        f"genome/module/cluster-pair observations "
        f"but produced "
        f"{len(expected_observations):,}."
    )


print(
    f"  Independently reconstructed "
    f"pair observations: "
    f"{len(expected_observations):,}"
)


## ================================================================== ##
## 4. Read Stage-12A physical pair observations
## ================================================================== ##

print()
print(
    "Reading Stage-12A physical pair observations..."
)


obs = pd.read_csv(
    PAIR_OBSERVATIONS,
    sep="\t",
    dtype=str,
    keep_default_na=False,
)


required_pair_columns = [
    "genome",
    "module",
    "cluster_a",
    "cluster_b",

    "n_proteins_cluster_a",
    "n_proteins_cluster_b",

    "n_total_protein_pairs",
    "n_same_contig_protein_pairs",

    "any_same_contig",

    "n_adjacent_protein_pairs",
    "any_adjacent",

    "n_within_5kb_protein_pairs",
    "any_within_5kb",

    "n_within_10kb_protein_pairs",
    "any_within_10kb",

    "n_within_20kb_protein_pairs",
    "any_within_20kb",

    "closest_gene_rank_distance",
    "closest_genes_between",

    "closest_intergenic_gap_bp",
    "closest_cds_overlap_nt",
    "closest_midpoint_distance_bp",

    "closest_orientation_class",
]


require_columns(
    obs,
    required_pair_columns,
    PAIR_OBSERVATIONS.name,
)


if len(obs) != EXPECTED_PAIR_OBSERVATIONS:

    fail(
        f"Expected "
        f"{EXPECTED_PAIR_OBSERVATIONS:,} "
        f"Stage-12A pair observations but found "
        f"{len(obs):,}."
    )


key_columns = [
    "genome",
    "module",
    "cluster_a",
    "cluster_b",
]


if (
    obs[
        key_columns
    ]
    .duplicated()
    .any()
):

    fail(
        "Duplicate genome/module/cluster-pair "
        "observations detected."
    )


## Cluster pairs should already be stored deterministically. ##

if not (
    obs[
        "cluster_a"
    ]
    <
    obs[
        "cluster_b"
    ]
).all():

    fail(
        "One or more Stage-12A cluster pairs "
        "are not stored as cluster_a < cluster_b."
    )


## ================================================================== ##
## 5. Strict key-set comparison:
##
## Stage 09 presence reconstruction
##
##       versus
##
## Stage 12A physical pair observations
## ================================================================== ##

expected_index = pd.MultiIndex.from_frame(
    expected_observations[
        key_columns
    ]
)


observed_index = pd.MultiIndex.from_frame(
    obs[
        key_columns
    ]
)


missing_from_stage12 = (
    expected_index
    .difference(
        observed_index
    )
)


unexpected_in_stage12 = (
    observed_index
    .difference(
        expected_index
    )
)


if len(missing_from_stage12) > 0:

    fail(
        f"{len(missing_from_stage12):,} "
        f"expected pair observations are missing "
        f"from Stage 12A."
    )


if len(unexpected_in_stage12) > 0:

    fail(
        f"{len(unexpected_in_stage12):,} "
        f"unexpected pair observations occur "
        f"in Stage 12A."
    )


## ================================================================== ##
## 6. Convert pair-count / flag columns to numeric
## ================================================================== ##

integer_columns = [
    "n_proteins_cluster_a",
    "n_proteins_cluster_b",

    "n_total_protein_pairs",
    "n_same_contig_protein_pairs",

    "any_same_contig",

    "n_adjacent_protein_pairs",
    "any_adjacent",

    "n_within_5kb_protein_pairs",
    "any_within_5kb",

    "n_within_10kb_protein_pairs",
    "any_within_10kb",

    "n_within_20kb_protein_pairs",
    "any_within_20kb",
]


for column in integer_columns:

    obs[
        column
    ] = pd.to_numeric(
        obs[
            column
        ],
        errors="raise",
    ).astype(int)


for column in [
    "closest_gene_rank_distance",
    "closest_genes_between",
    "closest_intergenic_gap_bp",
    "closest_cds_overlap_nt",
    "closest_midpoint_distance_bp",
]:

    obs[
        column
    ] = pd.to_numeric(
        obs[
            column
        ],
        errors="coerce",
    )


## ================================================================== ##
## 7. Protein multiplicity cross-check
## ================================================================== ##

counts_a = protein_counts.rename(
    columns={
        "cluster":
            "cluster_a",

        "expected_n_proteins":
            "expected_n_proteins_cluster_a",
    }
)


counts_b = protein_counts.rename(
    columns={
        "cluster":
            "cluster_b",

        "expected_n_proteins":
            "expected_n_proteins_cluster_b",
    }
)


obs_check = (
    obs
    .merge(
        counts_a,
        on=[
            "genome",
            "module",
            "cluster_a",
        ],
        how="left",
        validate="many_to_one",
    )
    .merge(
        counts_b,
        on=[
            "genome",
            "module",
            "cluster_b",
        ],
        how="left",
        validate="many_to_one",
    )
)


if (
    obs_check[
        "expected_n_proteins_cluster_a"
    ]
    .isna()
    .any()
    or
    obs_check[
        "expected_n_proteins_cluster_b"
    ]
    .isna()
    .any()
):

    fail(
        "Could not recover protein multiplicities "
        "for one or more pair observations."
    )


if not (
    obs_check[
        "n_proteins_cluster_a"
    ]
    ==
    obs_check[
        "expected_n_proteins_cluster_a"
    ]
).all():

    fail(
        "Cluster-A protein multiplicity disagrees "
        "with Stage 09."
    )


if not (
    obs_check[
        "n_proteins_cluster_b"
    ]
    ==
    obs_check[
        "expected_n_proteins_cluster_b"
    ]
).all():

    fail(
        "Cluster-B protein multiplicity disagrees "
        "with Stage 09."
    )


expected_total_pairs = (
    obs_check[
        "n_proteins_cluster_a"
    ]
    *
    obs_check[
        "n_proteins_cluster_b"
    ]
)


if not (
    obs_check[
        "n_total_protein_pairs"
    ]
    ==
    expected_total_pairs
).all():

    fail(
        "n_total_protein_pairs is inconsistent "
        "with protein multiplicities."
    )


## ================================================================== ##
## 8. Physical-threshold consistency QC
## ================================================================== ##

if not (
    (
        obs[
            "n_adjacent_protein_pairs"
        ]
        <=
        obs[
            "n_same_contig_protein_pairs"
        ]
    )
    &
    (
        obs[
            "n_within_5kb_protein_pairs"
        ]
        <=
        obs[
            "n_within_10kb_protein_pairs"
        ]
    )
    &
    (
        obs[
            "n_within_10kb_protein_pairs"
        ]
        <=
        obs[
            "n_within_20kb_protein_pairs"
        ]
    )
    &
    (
        obs[
            "n_within_20kb_protein_pairs"
        ]
        <=
        obs[
            "n_same_contig_protein_pairs"
        ]
    )
).all():

    fail(
        "Physical pair-count hierarchy failed."
    )


flag_pairs = [
    (
        "n_same_contig_protein_pairs",
        "any_same_contig",
    ),

    (
        "n_adjacent_protein_pairs",
        "any_adjacent",
    ),

    (
        "n_within_5kb_protein_pairs",
        "any_within_5kb",
    ),

    (
        "n_within_10kb_protein_pairs",
        "any_within_10kb",
    ),

    (
        "n_within_20kb_protein_pairs",
        "any_within_20kb",
    ),
]


for count_column, flag_column in flag_pairs:

    expected_flag = (
        obs[
            count_column
        ]
        >
        0
    ).astype(int)


    if not (
        obs[
            flag_column
        ]
        ==
        expected_flag
    ).all():

        fail(
            f"{flag_column} is inconsistent "
            f"with {count_column}."
        )


## ================================================================== ##
## 9. Summarize each possible cluster pair
##
## IMPORTANT:
##
## Percentages are primarily genome-level.
##
## A genome with three paralogues does NOT count three times simply
## because it provides multiple protein-pair combinations.
## ================================================================== ##

print()
print(
    "Summarizing physical association "
    "for all possible module cluster pairs..."
)


obs_groups = {
    key:
        group.copy()

    for key, group
    in obs.groupby(
        [
            "module",
            "cluster_a",
            "cluster_b",
        ],
        sort=False,
    )
}


pair_summary_rows = []


for row in possible_pairs.itertuples(
    index=False
):

    key = (
        row.module,
        row.cluster_a,
        row.cluster_b,
    )


    group = obs_groups.get(
        key
    )


    ## -------------------------------------------------------------- ##
    ## Pair never co-present in any genome
    ## -------------------------------------------------------------- ##

    if group is None:

        pair_summary_rows.append(
            {
                "module":
                    row.module,

                "cluster_a":
                    row.cluster_a,

                "cluster_b":
                    row.cluster_b,

                "n_genomes_both":
                    0,

                "n_genomes_same_contig":
                    0,

                "n_genomes_no_same_contig_evidence":
                    0,

                "n_genomes_within_5kb":
                    0,

                "n_genomes_within_10kb":
                    0,

                "n_genomes_within_20kb":
                    0,

                "n_genomes_adjacent":
                    0,

                "pct_same_contig_of_copresent_genomes":
                    np.nan,

                "pct_no_same_contig_evidence_of_copresent_genomes":
                    np.nan,

                "pct_within_5kb_of_copresent_genomes":
                    np.nan,

                "pct_within_10kb_of_copresent_genomes":
                    np.nan,

                "pct_within_20kb_of_copresent_genomes":
                    np.nan,

                "pct_adjacent_of_copresent_genomes":
                    np.nan,

                "pct_within_5kb_given_same_contig":
                    np.nan,

                "pct_within_10kb_given_same_contig":
                    np.nan,

                "pct_within_20kb_given_same_contig":
                    np.nan,

                "pct_adjacent_given_same_contig":
                    np.nan,

                "sum_protein_pair_combinations":
                    0,

                "sum_same_contig_protein_pair_combinations":
                    0,

                "sum_within_5kb_protein_pair_combinations":
                    0,

                "sum_within_10kb_protein_pair_combinations":
                    0,

                "sum_within_20kb_protein_pair_combinations":
                    0,

                "sum_adjacent_protein_pair_combinations":
                    0,

                "median_closest_intergenic_gap_bp":
                    np.nan,

                "q25_closest_intergenic_gap_bp":
                    np.nan,

                "q75_closest_intergenic_gap_bp":
                    np.nan,

                "median_closest_midpoint_distance_bp":
                    np.nan,

                "median_closest_genes_between":
                    np.nan,

                "n_closest_orientation_codirectional":
                    0,

                "n_closest_orientation_codirectional_right":
                    0,

                "n_closest_orientation_codirectional_left":
                    0,

                "n_closest_orientation_convergent":
                    0,

                "n_closest_orientation_divergent":
                    0,

                "n_closest_orientation_unknown":
                    0,

                "pct_closest_orientation_codirectional":
                    np.nan,

                "pct_closest_orientation_convergent":
                    np.nan,

                "pct_closest_orientation_divergent":
                    np.nan,
            }
        )

        continue


    ## -------------------------------------------------------------- ##
    ## Genome-level evidence
    ## -------------------------------------------------------------- ##

    n_both = len(
        group
    )


    n_same = int(
        group[
            "any_same_contig"
        ]
        .sum()
    )


    n_no_same = (
        n_both
        -
        n_same
    )


    n_5kb = int(
        group[
            "any_within_5kb"
        ]
        .sum()
    )


    n_10kb = int(
        group[
            "any_within_10kb"
        ]
        .sum()
    )


    n_20kb = int(
        group[
            "any_within_20kb"
        ]
        .sum()
    )


    n_adjacent = int(
        group[
            "any_adjacent"
        ]
        .sum()
    )


    ## -------------------------------------------------------------- ##
    ## Closest same-contig occurrence in each genome
    ## -------------------------------------------------------------- ##

    same_contig = group[
        group[
            "any_same_contig"
        ]
        ==
        1
    ].copy()


    orientation_counts = (
        same_contig[
            "closest_orientation_class"
        ]
        .replace(
            "",
            "unknown"
        )
        .value_counts()
        .to_dict()
    )


    n_right = int(
        orientation_counts.get(
            "codirectional_right",
            0,
        )
    )


    n_left = int(
        orientation_counts.get(
            "codirectional_left",
            0,
        )
    )


    n_codirectional = (
        n_right
        +
        n_left
    )


    n_convergent = int(
        orientation_counts.get(
            "convergent",
            0,
        )
    )


    n_divergent = int(
        orientation_counts.get(
            "divergent",
            0,
        )
    )


    n_known_orientation = (
        n_codirectional
        +
        n_convergent
        +
        n_divergent
    )


    n_unknown_orientation = (
        n_same
        -
        n_known_orientation
    )


    pair_summary_rows.append(
        {
            "module":
                row.module,

            "cluster_a":
                row.cluster_a,

            "cluster_b":
                row.cluster_b,

            ## Genome-level denominators ##

            "n_genomes_both":
                n_both,

            "n_genomes_same_contig":
                n_same,

            "n_genomes_no_same_contig_evidence":
                n_no_same,

            "n_genomes_within_5kb":
                n_5kb,

            "n_genomes_within_10kb":
                n_10kb,

            "n_genomes_within_20kb":
                n_20kb,

            "n_genomes_adjacent":
                n_adjacent,


            ## Percent of ALL genomes where both families occur ##

            "pct_same_contig_of_copresent_genomes":
                pct(
                    n_same,
                    n_both,
                ),

            "pct_no_same_contig_evidence_of_copresent_genomes":
                pct(
                    n_no_same,
                    n_both,
                ),

            "pct_within_5kb_of_copresent_genomes":
                pct(
                    n_5kb,
                    n_both,
                ),

            "pct_within_10kb_of_copresent_genomes":
                pct(
                    n_10kb,
                    n_both,
                ),

            "pct_within_20kb_of_copresent_genomes":
                pct(
                    n_20kb,
                    n_both,
                ),

            "pct_adjacent_of_copresent_genomes":
                pct(
                    n_adjacent,
                    n_both,
                ),


            ## Conditional on at least one same-contig occurrence ##

            "pct_within_5kb_given_same_contig":
                pct(
                    n_5kb,
                    n_same,
                ),

            "pct_within_10kb_given_same_contig":
                pct(
                    n_10kb,
                    n_same,
                ),

            "pct_within_20kb_given_same_contig":
                pct(
                    n_20kb,
                    n_same,
                ),

            "pct_adjacent_given_same_contig":
                pct(
                    n_adjacent,
                    n_same,
                ),


            ## Paralog-sensitive protein-pair totals retained
            ## separately from the primary genome-level statistics. ##

            "sum_protein_pair_combinations":
                int(
                    group[
                        "n_total_protein_pairs"
                    ]
                    .sum()
                ),

            "sum_same_contig_protein_pair_combinations":
                int(
                    group[
                        "n_same_contig_protein_pairs"
                    ]
                    .sum()
                ),

            "sum_within_5kb_protein_pair_combinations":
                int(
                    group[
                        "n_within_5kb_protein_pairs"
                    ]
                    .sum()
                ),

            "sum_within_10kb_protein_pair_combinations":
                int(
                    group[
                        "n_within_10kb_protein_pairs"
                    ]
                    .sum()
                ),

            "sum_within_20kb_protein_pair_combinations":
                int(
                    group[
                        "n_within_20kb_protein_pairs"
                    ]
                    .sum()
                ),

            "sum_adjacent_protein_pair_combinations":
                int(
                    group[
                        "n_adjacent_protein_pairs"
                    ]
                    .sum()
                ),


            ## Closest same-contig pair distribution ##

            "median_closest_intergenic_gap_bp":
                numeric_median(
                    same_contig[
                        "closest_intergenic_gap_bp"
                    ]
                ),

            "q25_closest_intergenic_gap_bp":
                numeric_quantile(
                    same_contig[
                        "closest_intergenic_gap_bp"
                    ],
                    0.25,
                ),

            "q75_closest_intergenic_gap_bp":
                numeric_quantile(
                    same_contig[
                        "closest_intergenic_gap_bp"
                    ],
                    0.75,
                ),

            "median_closest_midpoint_distance_bp":
                numeric_median(
                    same_contig[
                        "closest_midpoint_distance_bp"
                    ]
                ),

            "median_closest_genes_between":
                numeric_median(
                    same_contig[
                        "closest_genes_between"
                    ]
                ),


            ## Orientation of closest same-contig pair ##

            "n_closest_orientation_codirectional":
                n_codirectional,

            "n_closest_orientation_codirectional_right":
                n_right,

            "n_closest_orientation_codirectional_left":
                n_left,

            "n_closest_orientation_convergent":
                n_convergent,

            "n_closest_orientation_divergent":
                n_divergent,

            "n_closest_orientation_unknown":
                n_unknown_orientation,

            "pct_closest_orientation_codirectional":
                pct(
                    n_codirectional,
                    n_same,
                ),

            "pct_closest_orientation_convergent":
                pct(
                    n_convergent,
                    n_same,
                ),

            "pct_closest_orientation_divergent":
                pct(
                    n_divergent,
                    n_same,
                ),
        }
    )


pair_summary = pd.DataFrame(
    pair_summary_rows
)


if len(pair_summary) != EXPECTED_POSSIBLE_CLUSTER_PAIRS:

    fail(
        "Final cluster-pair summary does not "
        "contain all possible module pairs."
    )


if (
    pair_summary[
        [
            "module",
            "cluster_a",
            "cluster_b",
        ]
    ]
    .duplicated()
    .any()
):

    fail(
        "Duplicate cluster-pair summary rows."
    )


## ================================================================== ##
## 10. Add cluster-level biological labels
##
## These are descriptive annotations only.
## They do NOT influence the physical colocalization calculations.
## ================================================================== ##

cluster_summary = pd.read_csv(
    CLUSTER_SUMMARY,
    sep="\t",
    dtype=str,
    keep_default_na=False,
)


require_columns(
    cluster_summary,
    [
        "cluster",
    ],
    CLUSTER_SUMMARY.name,
)


if (
    cluster_summary[
        "cluster"
    ]
    .duplicated()
    .any()
):

    fail(
        "Duplicate clusters in cluster_summary.tsv."
    )


useful_cluster_columns = [
    "cluster",

    "representative_protein",

    "n_proteins",
    "n_genomes",

    "mean_heme_count",

    "pct_fegenie_negative",
    "dominant_fegenie_HMM",

    "dominant_signalp_prediction",
    "dominant_deeptmhmm_class",
    "dominant_localization_class",

    "pct_strong_unknown_soluble_periplasmic_like",
]


useful_cluster_columns = [
    column
    for column
    in useful_cluster_columns
    if column
    in cluster_summary.columns
]


cluster_info = cluster_summary[
    useful_cluster_columns
].copy()


cluster_a_info = cluster_info.rename(
    columns={
        column:
            (
                "cluster_a"
                if column == "cluster"
                else
                "cluster_a_" + column
            )

        for column
        in cluster_info.columns
    }
)


cluster_b_info = cluster_info.rename(
    columns={
        column:
            (
                "cluster_b"
                if column == "cluster"
                else
                "cluster_b_" + column
            )

        for column
        in cluster_info.columns
    }
)


pair_summary = (
    pair_summary
    .merge(
        cluster_a_info,
        on="cluster_a",
        how="left",
        validate="many_to_one",
    )
    .merge(
        cluster_b_info,
        on="cluster_b",
        how="left",
        validate="many_to_one",
    )
)


## ================================================================== ##
## 11. Deterministic output ordering
## ================================================================== ##

pair_summary = (
    pair_summary
    .sort_values(
        [
            "module",
            "cluster_a",
            "cluster_b",
        ],
        kind="stable",
    )
    .reset_index(
        drop=True
    )
)


pair_summary.to_csv(
    OUT_PAIR_SUMMARY,
    sep="\t",
    index=False,
    na_rep="",
    float_format="%.6f",
)


## ================================================================== ##
## 12. Module-level summary
##
## Module-level totals use:
##
##   genome x cluster-pair observations
##
## as the denominator.
##
## These are NOT counts of independent genomes because one genome can
## contribute several distinct cluster pairs within a large module.
## ================================================================== ##

print()
print(
    "Building module-level colocalization summary..."
)


module_genome_counts = (
    presence
    .groupby(
        "module"
    )[
        "genome"
    ]
    .nunique()
    .to_dict()
)


module_cluster_counts = (
    modules
    .groupby(
        "module"
    )[
        "cluster"
    ]
    .nunique()
    .to_dict()
)


module_summary_rows = []


for module, group in pair_summary.groupby(
    "module",
    sort=True,
):

    n_clusters = int(
        module_cluster_counts[
            module
        ]
    )


    n_possible_pairs = len(
        group
    )


    n_pairs_copresent = int(
        (
            group[
                "n_genomes_both"
            ]
            >
            0
        )
        .sum()
    )


    n_pairs_never_copresent = (
        n_possible_pairs
        -
        n_pairs_copresent
    )


    n_pairs_same_contig = int(
        (
            group[
                "n_genomes_same_contig"
            ]
            >
            0
        )
        .sum()
    )


    n_pairs_5kb = int(
        (
            group[
                "n_genomes_within_5kb"
            ]
            >
            0
        )
        .sum()
    )


    n_pairs_10kb = int(
        (
            group[
                "n_genomes_within_10kb"
            ]
            >
            0
        )
        .sum()
    )


    n_pairs_20kb = int(
        (
            group[
                "n_genomes_within_20kb"
            ]
            >
            0
        )
        .sum()
    )


    n_pairs_adjacent = int(
        (
            group[
                "n_genomes_adjacent"
            ]
            >
            0
        )
        .sum()
    )


    total_copresent = int(
        group[
            "n_genomes_both"
        ]
        .sum()
    )


    total_same_contig = int(
        group[
            "n_genomes_same_contig"
        ]
        .sum()
    )


    total_5kb = int(
        group[
            "n_genomes_within_5kb"
        ]
        .sum()
    )


    total_10kb = int(
        group[
            "n_genomes_within_10kb"
        ]
        .sum()
    )


    total_20kb = int(
        group[
            "n_genomes_within_20kb"
        ]
        .sum()
    )


    total_adjacent = int(
        group[
            "n_genomes_adjacent"
        ]
        .sum()
    )


    total_codirectional = int(
        group[
            "n_closest_orientation_codirectional"
        ]
        .sum()
    )


    total_convergent = int(
        group[
            "n_closest_orientation_convergent"
        ]
        .sum()
    )


    total_divergent = int(
        group[
            "n_closest_orientation_divergent"
        ]
        .sum()
    )


    module_summary_rows.append(
        {
            "module":
                module,

            "n_clusters":
                n_clusters,

            "n_genomes_with_at_least_one_module_cluster":
                int(
                    module_genome_counts[
                        module
                    ]
                ),

            "n_possible_cluster_pairs":
                n_possible_pairs,

            "n_cluster_pairs_copresent_at_least_once":
                n_pairs_copresent,

            "n_cluster_pairs_never_copresent":
                n_pairs_never_copresent,

            "n_cluster_pairs_with_same_contig_evidence":
                n_pairs_same_contig,

            "n_cluster_pairs_with_5kb_evidence":
                n_pairs_5kb,

            "n_cluster_pairs_with_10kb_evidence":
                n_pairs_10kb,

            "n_cluster_pairs_with_20kb_evidence":
                n_pairs_20kb,

            "n_cluster_pairs_with_adjacency_evidence":
                n_pairs_adjacent,


            ## Genome x cluster-pair observations ##

            "n_genome_cluster_pair_copresence_observations":
                total_copresent,

            "n_genome_cluster_pair_same_contig_observations":
                total_same_contig,

            "n_genome_cluster_pair_within_5kb_observations":
                total_5kb,

            "n_genome_cluster_pair_within_10kb_observations":
                total_10kb,

            "n_genome_cluster_pair_within_20kb_observations":
                total_20kb,

            "n_genome_cluster_pair_adjacent_observations":
                total_adjacent,


            ## Weighted across genome x cluster-pair observations ##

            "weighted_pct_same_contig_of_copresent":
                pct(
                    total_same_contig,
                    total_copresent,
                ),

            "weighted_pct_within_5kb_of_copresent":
                pct(
                    total_5kb,
                    total_copresent,
                ),

            "weighted_pct_within_10kb_of_copresent":
                pct(
                    total_10kb,
                    total_copresent,
                ),

            "weighted_pct_within_20kb_of_copresent":
                pct(
                    total_20kb,
                    total_copresent,
                ),

            "weighted_pct_adjacent_of_copresent":
                pct(
                    total_adjacent,
                    total_copresent,
                ),


            ## Conditional on same-contig evidence ##

            "weighted_pct_within_5kb_given_same_contig":
                pct(
                    total_5kb,
                    total_same_contig,
                ),

            "weighted_pct_within_10kb_given_same_contig":
                pct(
                    total_10kb,
                    total_same_contig,
                ),

            "weighted_pct_within_20kb_given_same_contig":
                pct(
                    total_20kb,
                    total_same_contig,
                ),

            "weighted_pct_adjacent_given_same_contig":
                pct(
                    total_adjacent,
                    total_same_contig,
                ),


            ## Unweighted pair-level summary ##

            "median_pair_pct_same_contig":
                numeric_median(
                    group.loc[
                        group[
                            "n_genomes_both"
                        ]
                        >
                        0,
                        "pct_same_contig_of_copresent_genomes",
                    ]
                ),

            "median_pair_pct_within_10kb":
                numeric_median(
                    group.loc[
                        group[
                            "n_genomes_both"
                        ]
                        >
                        0,
                        "pct_within_10kb_of_copresent_genomes",
                    ]
                ),


            ## Orientation of closest same-contig pair ##

            "n_closest_orientation_codirectional":
                total_codirectional,

            "n_closest_orientation_convergent":
                total_convergent,

            "n_closest_orientation_divergent":
                total_divergent,

            "pct_closest_orientation_codirectional":
                pct(
                    total_codirectional,
                    total_same_contig,
                ),

            "pct_closest_orientation_convergent":
                pct(
                    total_convergent,
                    total_same_contig,
                ),

            "pct_closest_orientation_divergent":
                pct(
                    total_divergent,
                    total_same_contig,
                ),
        }
    )


module_summary = pd.DataFrame(
    module_summary_rows
)


if len(module_summary) != EXPECTED_MODULES:

    fail(
        f"Expected "
        f"{EXPECTED_MODULES} module summary rows "
        f"but created "
        f"{len(module_summary)}."
    )


module_summary.to_csv(
    OUT_MODULE_SUMMARY,
    sep="\t",
    index=False,
    na_rep="",
    float_format="%.6f",
)


## ================================================================== ##
## 13. Global cross-checks
## ================================================================== ##

if (
    int(
        pair_summary[
            "n_genomes_both"
        ]
        .sum()
    )
    !=
    EXPECTED_PAIR_OBSERVATIONS
):

    fail(
        "Pair-summary co-presence total does not "
        "recover all Stage-12A observations."
    )


stage12_same_contig = int(
    obs[
        "any_same_contig"
    ]
    .sum()
)


stage12_5kb = int(
    obs[
        "any_within_5kb"
    ]
    .sum()
)


stage12_10kb = int(
    obs[
        "any_within_10kb"
    ]
    .sum()
)


stage12_20kb = int(
    obs[
        "any_within_20kb"
    ]
    .sum()
)


stage12_adjacent = int(
    obs[
        "any_adjacent"
    ]
    .sum()
)


if (
    int(
        pair_summary[
            "n_genomes_same_contig"
        ]
        .sum()
    )
    !=
    stage12_same_contig
):

    fail(
        "Same-contig summary total disagrees "
        "with Stage 12A."
    )


if (
    int(
        pair_summary[
            "n_genomes_within_5kb"
        ]
        .sum()
    )
    !=
    stage12_5kb
):

    fail(
        "5-kb summary total disagrees "
        "with Stage 12A."
    )


if (
    int(
        pair_summary[
            "n_genomes_within_10kb"
        ]
        .sum()
    )
    !=
    stage12_10kb
):

    fail(
        "10-kb summary total disagrees "
        "with Stage 12A."
    )


if (
    int(
        pair_summary[
            "n_genomes_within_20kb"
        ]
        .sum()
    )
    !=
    stage12_20kb
):

    fail(
        "20-kb summary total disagrees "
        "with Stage 12A."
    )


if (
    int(
        pair_summary[
            "n_genomes_adjacent"
        ]
        .sum()
    )
    !=
    stage12_adjacent
):

    fail(
        "Adjacency summary total disagrees "
        "with Stage 12A."
    )


## ================================================================== ##
## 14. QC table
## ================================================================== ##

n_pairs_copresent = int(
    (
        pair_summary[
            "n_genomes_both"
        ]
        >
        0
    )
    .sum()
)


n_pairs_same_contig = int(
    (
        pair_summary[
            "n_genomes_same_contig"
        ]
        >
        0
    )
    .sum()
)


n_pairs_10kb = int(
    (
        pair_summary[
            "n_genomes_within_10kb"
        ]
        >
        0
    )
    .sum()
)


n_pairs_adjacent = int(
    (
        pair_summary[
            "n_genomes_adjacent"
        ]
        >
        0
    )
    .sum()
)


qc = pd.DataFrame(
    [
        [
            "module_proteins",
            len(
                proteins
            ),
        ],

        [
            "module_clusters",
            modules[
                "cluster"
            ]
            .nunique(),
        ],

        [
            "modules",
            modules[
                "module"
            ]
            .nunique(),
        ],

        [
            "possible_within_module_cluster_pairs",
            len(
                pair_summary
            ),
        ],

        [
            "cluster_pairs_copresent_at_least_once",
            n_pairs_copresent,
        ],

        [
            "cluster_pairs_with_same_contig_evidence",
            n_pairs_same_contig,
        ],

        [
            "cluster_pairs_with_10kb_evidence",
            n_pairs_10kb,
        ],

        [
            "cluster_pairs_with_adjacency_evidence",
            n_pairs_adjacent,
        ],

        [
            "expected_genome_cluster_pair_observations_from_stage09",
            len(
                expected_observations
            ),
        ],

        [
            "observed_genome_cluster_pair_observations_stage12",
            len(
                obs
            ),
        ],

        [
            "pair_observation_keyset_matches_stage09",
            1,
        ],

        [
            "protein_multiplicity_crosscheck_pass",
            1,
        ],

        [
            "physical_threshold_hierarchy_pass",
            1,
        ],

        [
            "genome_cluster_pair_same_contig_observations",
            stage12_same_contig,
        ],

        [
            "genome_cluster_pair_within_5kb_observations",
            stage12_5kb,
        ],

        [
            "genome_cluster_pair_within_10kb_observations",
            stage12_10kb,
        ],

        [
            "genome_cluster_pair_within_20kb_observations",
            stage12_20kb,
        ],

        [
            "genome_cluster_pair_adjacent_observations",
            stage12_adjacent,
        ],
    ],
    columns=[
        "metric",
        "value",
    ],
)


qc.to_csv(
    OUT_QC,
    sep="\t",
    index=False,
)


## ================================================================== ##
## 15. Terminal summary
## ================================================================== ##

print()
print(
    "Global module-colocalization summary"
)


print(
    f"  Possible within-module pairs:    "
    f"{len(pair_summary):,}"
)

print(
    f"  Co-present in >=1 genome:        "
    f"{n_pairs_copresent:,}"
)

print(
    f"  Same-contig evidence:            "
    f"{n_pairs_same_contig:,}"
)

print(
    f"  Within-10-kb evidence:           "
    f"{n_pairs_10kb:,}"
)

print(
    f"  Adjacency evidence:              "
    f"{n_pairs_adjacent:,}"
)


print()
print(
    "Genome x cluster-pair observations"
)

print(
    f"  Co-present:                      "
    f"{len(obs):,}"
)

print(
    f"  Same contig:                     "
    f"{stage12_same_contig:,} "
    f"("
    f"{pct(stage12_same_contig, len(obs)):.2f}%"
    f")"
)

print(
    f"  Within 5 kb:                     "
    f"{stage12_5kb:,} "
    f"("
    f"{pct(stage12_5kb, len(obs)):.2f}%"
    f")"
)

print(
    f"  Within 10 kb:                    "
    f"{stage12_10kb:,} "
    f"("
    f"{pct(stage12_10kb, len(obs)):.2f}%"
    f")"
)

print(
    f"  Within 20 kb:                    "
    f"{stage12_20kb:,} "
    f"("
    f"{pct(stage12_20kb, len(obs)):.2f}%"
    f")"
)

print(
    f"  Adjacent:                        "
    f"{stage12_adjacent:,} "
    f"("
    f"{pct(stage12_adjacent, len(obs)):.2f}%"
    f")"
)


if stage12_same_contig > 0:

    print()
    print(
        "Conditional on same-contig evidence"
    )

    print(
        f"  Within 5 kb:                     "
        f"{pct(stage12_5kb, stage12_same_contig):.2f}%"
    )

    print(
        f"  Within 10 kb:                    "
        f"{pct(stage12_10kb, stage12_same_contig):.2f}%"
    )

    print(
        f"  Within 20 kb:                    "
        f"{pct(stage12_20kb, stage12_same_contig):.2f}%"
    )

    print(
        f"  Adjacent:                        "
        f"{pct(stage12_adjacent, stage12_same_contig):.2f}%"
    )


## ================================================================== ##
## 16. Spotlight biologically anchored two-cluster modules
##
## This is terminal reporting only. It does not change any data.
## ================================================================== ##

for spotlight_module in [
    "Module_20",
    "Module_35",
]:

    spotlight = pair_summary[
        pair_summary[
            "module"
        ]
        ==
        spotlight_module
    ]


    if len(spotlight) != 1:
        continue


    r = spotlight.iloc[0]


    print()
    print(
        f"{spotlight_module} spotlight"
    )

    print(
        f"  Pair:                            "
        f"{r['cluster_a']} <-> "
        f"{r['cluster_b']}"
    )

    print(
        f"  Genomes with both:               "
        f"{int(r['n_genomes_both']):,}"
    )

    print(
        f"  Same contig:                     "
        f"{int(r['n_genomes_same_contig']):,} "
        f"("
        f"{r['pct_same_contig_of_copresent_genomes']:.2f}%"
        f")"
    )

    print(
        f"  Within 5 kb:                     "
        f"{int(r['n_genomes_within_5kb']):,}"
    )

    print(
        f"  Within 10 kb:                    "
        f"{int(r['n_genomes_within_10kb']):,}"
    )

    print(
        f"  Within 20 kb:                    "
        f"{int(r['n_genomes_within_20kb']):,}"
    )

    print(
        f"  Adjacent:                        "
        f"{int(r['n_genomes_adjacent']):,}"
    )

    if pd.notna(
        r[
            "median_closest_intergenic_gap_bp"
        ]
    ):

        print(
            f"  Median closest intergenic gap:   "
            f"{r['median_closest_intergenic_gap_bp']:.1f} bp"
        )

        print(
            f"  Median genes between:            "
            f"{r['median_closest_genes_between']:.1f}"
        )


print()
print("=" * 80)
print("SUCCESS")
print("=" * 80)


print(
    f"Cluster-pair summary: "
    f"{OUT_PAIR_SUMMARY}"
)

print(
    f"Module summary:       "
    f"{OUT_MODULE_SUMMARY}"
)

print(
    f"QC:                   "
    f"{OUT_QC}"
)
