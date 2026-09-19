#!/usr/bin/env python3

from pathlib import Path
import sys

import numpy as np
import pandas as pd


## ================================================================== ##
## STAGE 13B - CUMULATIVE MMSEQS2 GENE-WINDOW CONSERVATION
##
## Primary biological unit:
##
##     MMseqs2 protein-family cluster
##
## This complements exact gene-offset conservation.
##
## Exact offsets ask:
##
##     "Is Cluster_B exactly +2 genes from Cluster_A?"
##
## Cumulative windows ask:
##
##     "Is Cluster_B anywhere within the first 2 downstream genes
##      of Cluster_A?"
##
## Windows:
##
##     1, 2, 3, 5, 10, 20 genes
##
## Negative = upstream
## Positive = downstream
##
## after orientation to the focal gene's transcriptional direction.
##
##
## CENSORING LOGIC
## ----------------
##
## Positive:
##     Neighbor is actually observed inside the requested window.
##     This remains informative even if the contig ends before the
##     complete requested window.
##
## Informative negative:
##     Neighbor is not observed AND the complete requested gene
##     window is available.
##
## Censored:
##     Neighbor is not observed AND the complete requested window
##     cannot be observed.
##
##
## IMPORTANT REPORTING PRINCIPLE
## -----------------------------
##
## Percentages are always accompanied by:
##
##     numerator
##     denominator
##     support string (e.g. 7/8)
##
## Therefore:
##
##     1/1 = 100%
##
## remains visibly based on only one informative observation.
##
## No "conserved" threshold is imposed in this stage.
## ================================================================== ##


## ================================================================== ##
## Paths
## ================================================================== ##

HERE = Path(__file__).resolve().parent
WORKFLOW = HERE.parent.parent


NEIGHBORHOODS = (
    WORKFLOW
    / "12_neighborhood_extraction"
    / "observed_neighborhood_genes.tsv"
)


OBSERVABILITY = (
    WORKFLOW
    / "12_neighborhood_extraction"
    / "focal_context_observability.tsv"
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

OUT_WINDOWS = (
    HERE
    / "cluster_neighbor_gene_window_conservation.tsv"
)


OUT_QC = (
    HERE
    / "mmseqs_gene_window_conservation_qc.tsv"
)


## ================================================================== ##
## Analysis parameters
## ================================================================== ##

GENE_RADII = [
    1,
    2,
    3,
    5,
    10,
    20,
]


## ================================================================== ##
## Expected upstream values
## ================================================================== ##

EXPECTED_FOCAL_PROTEINS = 10_537
EXPECTED_FOCAL_CLUSTERS = 156
EXPECTED_MODULES = 35
EXPECTED_MMSEQS_CLUSTERS = 1_390


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
            f"{source} is missing required column(s):\n"
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


## ================================================================== ##
## Start
## ================================================================== ##

print("=" * 80)
print("STAGE 13B - CUMULATIVE MMSEQS2 GENE-WINDOW CONSERVATION")
print("=" * 80)


## ================================================================== ##
## 1. Read focal occurrences and gene-space observability
## ================================================================== ##

print()
print("Reading focal observability...")


focals = pd.read_csv(
    OBSERVABILITY,
    sep="\t",
    dtype=str,
    keep_default_na=False,
)


require_columns(
    focals,
    [
        "genome",
        "focal_protein_id",
        "focal_cluster",
        "focal_module",

        "n_oriented_upstream_genes_available",
        "n_oriented_downstream_genes_available",
    ],
    OBSERVABILITY.name,
)


if len(focals) != EXPECTED_FOCAL_PROTEINS:

    fail(
        f"Expected "
        f"{EXPECTED_FOCAL_PROTEINS:,} "
        f"focal proteins but found "
        f"{len(focals):,}."
    )


if (
    focals[
        [
            "genome",
            "focal_protein_id",
        ]
    ]
    .duplicated()
    .any()
):

    fail(
        "Duplicate focal protein keys."
    )


if (
    focals[
        "focal_cluster"
    ]
    .nunique()
    !=
    EXPECTED_FOCAL_CLUSTERS
):

    fail(
        f"Expected "
        f"{EXPECTED_FOCAL_CLUSTERS} "
        f"focal MMseqs2 clusters."
    )


if (
    focals[
        "focal_module"
    ]
    .nunique()
    !=
    EXPECTED_MODULES
):

    fail(
        f"Expected "
        f"{EXPECTED_MODULES} "
        f"MCL modules."
    )


for column in [
    "n_oriented_upstream_genes_available",
    "n_oriented_downstream_genes_available",
]:

    focals[
        column
    ] = pd.to_numeric(
        focals[
            column
        ],
        errors="raise",
    ).astype(int)


focals[
    "focal_key"
] = (
    focals[
        "genome"
    ]
    +
    "\t"
    +
    focals[
        "focal_protein_id"
    ]
)


print(
    f"  Focal proteins: "
    f"{len(focals):,}"
)

print(
    f"  Focal clusters: "
    f"{focals['focal_cluster'].nunique():,}"
)


## ================================================================== ##
## 2. Read observed neighborhoods
## ================================================================== ##

print()
print("Reading observed neighborhoods...")


neigh = pd.read_csv(
    NEIGHBORHOODS,
    sep="\t",
    dtype=str,
    keep_default_na=False,
)


require_columns(
    neigh,
    [
        "genome",

        "focal_protein_id",
        "focal_cluster",
        "focal_module",

        "neighbor_protein_id",
        "neighbor_cluster",

        "oriented_gene_offset",

        "is_focal",
    ],
    NEIGHBORHOODS.name,
)


for column in [
    "oriented_gene_offset",
    "is_focal",
]:

    neigh[
        column
    ] = pd.to_numeric(
        neigh[
            column
        ],
        errors="raise",
    ).astype(int)


neigh[
    "focal_key"
] = (
    neigh[
        "genome"
    ]
    +
    "\t"
    +
    neigh[
        "focal_protein_id"
    ]
)


## ================================================================== ##
## 3. Keep MMseqs2-defined neighboring proteins within +/-20 genes
##
## The focal self-row is removed.
##
## Nearby paralogues belonging to the SAME cluster as the focal are
## intentionally retained.
## ================================================================== ##

mmseqs_neighbors = neigh[
    (
        neigh[
            "is_focal"
        ]
        ==
        0
    )
    &
    (
        neigh[
            "neighbor_cluster"
        ]
        .str.strip()
        !=
        ""
    )
    &
    (
        neigh[
            "oriented_gene_offset"
        ]
        !=
        0
    )
    &
    (
        neigh[
            "oriented_gene_offset"
        ]
        .abs()
        <=
        max(
            GENE_RADII
        )
    )
].copy()


mmseqs_neighbors[
    "direction"
] = np.where(
    mmseqs_neighbors[
        "oriented_gene_offset"
    ]
    <
    0,
    "upstream",
    "downstream",
)


mmseqs_neighbors[
    "absolute_gene_offset"
] = (
    mmseqs_neighbors[
        "oriented_gene_offset"
    ]
    .abs()
)


print(
    f"  MMseqs2 neighbor observations "
    f"within +/-20 genes: "
    f"{len(mmseqs_neighbors):,}"
)

print(
    f"  Distinct neighboring clusters:  "
    f"{mmseqs_neighbors['neighbor_cluster'].nunique():,}"
)


## ================================================================== ##
## 4. Module lookup
## ================================================================== ##

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


if (
    modules[
        "cluster"
    ]
    .duplicated()
    .any()
):

    fail(
        "A cluster occurs more than once in "
        "module_membership.tsv."
    )


module_lookup = dict(
    zip(
        modules[
            "cluster"
        ],
        modules[
            "module"
        ],
    )
)


## ================================================================== ##
## 5. MMseqs2 cluster metadata
##
## Cluster ID remains primary identity.
## Metadata are interpretation/display layers only.
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
    .nunique()
    !=
    EXPECTED_MMSEQS_CLUSTERS
):

    fail(
        f"Expected "
        f"{EXPECTED_MMSEQS_CLUSTERS:,} "
        f"MMseqs2 clusters but found "
        f"{cluster_summary['cluster'].nunique():,}."
    )


metadata_columns = [
    "cluster",

    "representative_protein",

    "n_proteins",
    "n_genomes",

    "mean_heme_count",

    "dominant_fegenie_HMM",

    "pct_fegenie_negative",

    "dominant_signalp_prediction",
    "dominant_deeptmhmm_class",
    "dominant_localization_class",

    "pct_strong_unknown_soluble_periplasmic_like",
]


metadata_columns = [
    column

    for column
    in metadata_columns

    if column
    in cluster_summary.columns
]


cluster_info = cluster_summary[
    metadata_columns
].copy()


## ================================================================== ##
## 6. Candidate focal-family / neighbor-family associations
##
## At least one actual observation within +/-20 genes is required
## for an association to enter this table.
##
## Once an association exists, BOTH upstream and downstream windows
## are evaluated, including zero-positive directions.
## ================================================================== ##

candidate_pairs = (
    mmseqs_neighbors[
        [
            "focal_cluster",
            "focal_module",
            "neighbor_cluster",
        ]
    ]
    .drop_duplicates()
    .sort_values(
        [
            "focal_module",
            "focal_cluster",
            "neighbor_cluster",
        ],
        kind="stable",
    )
    .reset_index(
        drop=True
    )
)


print()
print(
    f"Observed focal/neighbor family associations: "
    f"{len(candidate_pairs):,}"
)


## ================================================================== ##
## Pre-index data for efficient repeated use
## ================================================================== ##

focal_groups = {
    cluster:
        group.copy()

    for cluster, group
    in focals.groupby(
        "focal_cluster",
        sort=False,
    )
}


hit_groups = {
    key:
        group.copy()

    for key, group
    in mmseqs_neighbors.groupby(
        [
            "focal_cluster",
            "neighbor_cluster",
        ],
        sort=False,
    )
}


## ================================================================== ##
## 7. Cumulative gene-window conservation
## ================================================================== ##

print()
print(
    "Calculating cumulative directional gene-window conservation..."
)


rows = []


for candidate in candidate_pairs.itertuples(
    index=False
):

    focal_cluster = (
        candidate.focal_cluster
    )

    focal_module = (
        candidate.focal_module
    )

    neighbor_cluster = (
        candidate.neighbor_cluster
    )


    neighbor_module = module_lookup.get(
        neighbor_cluster,
        "",
    )


    focal_base = focal_groups[
        focal_cluster
    ].copy()


    pair_hits = hit_groups[
        (
            focal_cluster,
            neighbor_cluster,
        )
    ].copy()


    for direction in [
        "upstream",
        "downstream",
    ]:

        direction_hits = pair_hits[
            pair_hits[
                "direction"
            ]
            ==
            direction
        ].copy()


        if direction == "upstream":

            availability_column = (
                "n_oriented_upstream_genes_available"
            )

        else:

            availability_column = (
                "n_oriented_downstream_genes_available"
            )


        for radius in GENE_RADII:

            focal_subset = (
                focal_base.copy()
            )


            ## ------------------------------------------------------ ##
            ## Actual neighbor-protein observations in this window
            ## ------------------------------------------------------ ##

            hits_in_window = direction_hits[
                direction_hits[
                    "absolute_gene_offset"
                ]
                <=
                radius
            ].copy()


            n_neighbor_protein_hits = len(
                hits_in_window
            )


            positive_keys = set(
                hits_in_window[
                    "focal_key"
                ]
            )


            focal_subset[
                "positive"
            ] = (
                focal_subset[
                    "focal_key"
                ]
                .isin(
                    positive_keys
                )
                .astype(int)
            )


            ## ------------------------------------------------------ ##
            ## Complete gene-window observability
            ## ------------------------------------------------------ ##

            focal_subset[
                "full_window_observable"
            ] = (
                focal_subset[
                    availability_column
                ]
                >=
                radius
            ).astype(int)


            ## ------------------------------------------------------ ##
            ## Positive observations remain informative even if the
            ## full requested radius is not available.
            ## ------------------------------------------------------ ##

            focal_subset[
                "informative"
            ] = (
                (
                    focal_subset[
                        "positive"
                    ]
                    ==
                    1
                )
                |
                (
                    focal_subset[
                        "full_window_observable"
                    ]
                    ==
                    1
                )
            ).astype(int)


            focal_subset[
                "partial_window_positive"
            ] = (
                (
                    focal_subset[
                        "positive"
                    ]
                    ==
                    1
                )
                &
                (
                    focal_subset[
                        "full_window_observable"
                    ]
                    ==
                    0
                )
            ).astype(int)


            ## ------------------------------------------------------ ##
            ## Occurrence-level statistics
            ## ------------------------------------------------------ ##

            n_occ_total = len(
                focal_subset
            )


            n_occ_full = int(
                focal_subset[
                    "full_window_observable"
                ]
                .sum()
            )


            n_occ_positive = int(
                focal_subset[
                    "positive"
                ]
                .sum()
            )


            n_occ_partial_positive = int(
                focal_subset[
                    "partial_window_positive"
                ]
                .sum()
            )


            n_occ_informative = int(
                focal_subset[
                    "informative"
                ]
                .sum()
            )


            n_occ_negative = (
                n_occ_informative
                -
                n_occ_positive
            )


            n_occ_censored = (
                n_occ_total
                -
                n_occ_informative
            )


            ## ------------------------------------------------------ ##
            ## Genome-level statistics
            ##
            ## If a genome has multiple focal paralogues:
            ##
            ## positive genome:
            ##     >=1 focal occurrence is positive
            ##
            ## informative negative genome:
            ##     no positives AND every focal copy has a complete
            ##     requested window
            ##
            ## otherwise:
            ##     censored
            ## ------------------------------------------------------ ##

            genome_status = (
                focal_subset
                .groupby(
                    "genome",
                    as_index=False,
                )
                .agg(
                    n_focal_occurrences=(
                        "focal_key",
                        "size"
                    ),

                    n_full_window_occurrences=(
                        "full_window_observable",
                        "sum"
                    ),

                    n_positive_occurrences=(
                        "positive",
                        "sum"
                    ),

                    n_partial_positive_occurrences=(
                        "partial_window_positive",
                        "sum"
                    ),
                )
            )


            genome_status[
                "positive_genome"
            ] = (
                genome_status[
                    "n_positive_occurrences"
                ]
                >
                0
            )


            genome_status[
                "fully_observable_genome"
            ] = (
                genome_status[
                    "n_full_window_occurrences"
                ]
                ==
                genome_status[
                    "n_focal_occurrences"
                ]
            )


            genome_status[
                "informative_genome"
            ] = (
                genome_status[
                    "positive_genome"
                ]
                |
                genome_status[
                    "fully_observable_genome"
                ]
            )


            n_genomes_total = len(
                genome_status
            )


            n_genomes_fully_observable = int(
                genome_status[
                    "fully_observable_genome"
                ]
                .sum()
            )


            n_genomes_positive = int(
                genome_status[
                    "positive_genome"
                ]
                .sum()
            )


            n_genomes_informative = int(
                genome_status[
                    "informative_genome"
                ]
                .sum()
            )


            n_genomes_negative = (
                n_genomes_informative
                -
                n_genomes_positive
            )


            n_genomes_censored = (
                n_genomes_total
                -
                n_genomes_informative
            )


            n_genomes_partial_positive = int(
                (
                    (
                        genome_status[
                            "positive_genome"
                        ]
                    )
                    &
                    (
                        genome_status[
                            "n_partial_positive_occurrences"
                        ]
                        >
                        0
                    )
                )
                .sum()
            )


            ## ------------------------------------------------------ ##
            ## Store counts immediately beside their percentages.
            ## ------------------------------------------------------ ##

            rows.append(
                {
                    "focal_cluster":
                        focal_cluster,

                    "focal_module":
                        focal_module,

                    "neighbor_cluster":
                        neighbor_cluster,

                    "neighbor_module":
                        neighbor_module,

                    "neighbor_same_mcl_module_as_focal":
                        int(
                            neighbor_module != ""
                            and
                            neighbor_module
                            ==
                            focal_module
                        ),

                    "neighbor_same_mmseqs_cluster_as_focal":
                        int(
                            neighbor_cluster
                            ==
                            focal_cluster
                        ),

                    "direction":
                        direction,

                    "gene_radius":
                        radius,


                    ## ---------------------------------------------- ##
                    ## Raw protein hits
                    ## ---------------------------------------------- ##

                    "n_neighbor_protein_hits_in_window":
                        n_neighbor_protein_hits,


                    ## ---------------------------------------------- ##
                    ## Occurrence level
                    ## ---------------------------------------------- ##

                    "n_focal_occurrences_total":
                        n_occ_total,

                    "n_focal_occurrences_full_window_observable":
                        n_occ_full,

                    "n_informative_focal_occurrences":
                        n_occ_informative,

                    "n_focal_occurrences_with_neighbor":
                        n_occ_positive,

                    "occurrence_support":
                        (
                            f"{n_occ_positive}/"
                            f"{n_occ_informative}"
                        ),

                    "pct_informative_occurrences_with_neighbor":
                        pct(
                            n_occ_positive,
                            n_occ_informative,
                        ),

                    "n_informative_negative_occurrences":
                        n_occ_negative,

                    "n_censored_focal_occurrences":
                        n_occ_censored,

                    "n_positive_occurrences_from_partial_windows":
                        n_occ_partial_positive,


                    ## ---------------------------------------------- ##
                    ## Genome level
                    ## ---------------------------------------------- ##

                    "n_focal_genomes_total":
                        n_genomes_total,

                    "n_genomes_all_focal_copies_full_window":
                        n_genomes_fully_observable,

                    "n_informative_genomes":
                        n_genomes_informative,

                    "n_genomes_with_neighbor":
                        n_genomes_positive,

                    "genome_support":
                        (
                            f"{n_genomes_positive}/"
                            f"{n_genomes_informative}"
                        ),

                    "pct_informative_genomes_with_neighbor":
                        pct(
                            n_genomes_positive,
                            n_genomes_informative,
                        ),

                    "n_informative_negative_genomes":
                        n_genomes_negative,

                    "n_censored_genomes":
                        n_genomes_censored,

                    "n_positive_genomes_from_partial_windows":
                        n_genomes_partial_positive,
                }
            )


windows = pd.DataFrame(
    rows
)


## ================================================================== ##
## 8. Strict structural QC
## ================================================================== ##

expected_rows = (
    len(
        candidate_pairs
    )
    *
    2
    *
    len(
        GENE_RADII
    )
)


if len(windows) != expected_rows:

    fail(
        f"Expected "
        f"{expected_rows:,} "
        f"gene-window rows but produced "
        f"{len(windows):,}."
    )


if (
    windows[
        [
            "focal_cluster",
            "neighbor_cluster",
            "direction",
            "gene_radius",
        ]
    ]
    .duplicated()
    .any()
):

    fail(
        "Duplicate focal/neighbor/direction/"
        "radius rows detected."
    )


## Positive counts must never exceed informative denominators. ##

if (
    windows[
        "n_focal_occurrences_with_neighbor"
    ]
    >
    windows[
        "n_informative_focal_occurrences"
    ]
).any():

    fail(
        "Occurrence positives exceed "
        "informative denominator."
    )


if (
    windows[
        "n_genomes_with_neighbor"
    ]
    >
    windows[
        "n_informative_genomes"
    ]
).any():

    fail(
        "Genome positives exceed "
        "informative denominator."
    )


## ================================================================== ##
## 9. Cumulative-window monotonicity QC
##
## As gene radius increases, the number of positive focal occurrences
## and positive genomes may stay equal or increase, but may never
## decrease.
## ================================================================== ##

print()
print(
    "Checking cumulative-window monotonicity..."
)


for (
    focal_cluster,
    neighbor_cluster,
    direction,
), group in windows.groupby(
    [
        "focal_cluster",
        "neighbor_cluster",
        "direction",
    ],
    sort=False,
):

    group = group.sort_values(
        "gene_radius"
    )


    occurrence_positive = (
        group[
            "n_focal_occurrences_with_neighbor"
        ]
        .to_numpy()
    )


    genome_positive = (
        group[
            "n_genomes_with_neighbor"
        ]
        .to_numpy()
    )


    protein_hits = (
        group[
            "n_neighbor_protein_hits_in_window"
        ]
        .to_numpy()
    )


    if (
        np.diff(
            occurrence_positive
        )
        <
        0
    ).any():

        fail(
            "Occurrence-positive counts decreased "
            f"with increasing radius for "
            f"{focal_cluster} -> "
            f"{neighbor_cluster} "
            f"({direction})."
        )


    if (
        np.diff(
            genome_positive
        )
        <
        0
    ).any():

        fail(
            "Genome-positive counts decreased "
            f"with increasing radius for "
            f"{focal_cluster} -> "
            f"{neighbor_cluster} "
            f"({direction})."
        )


    if (
        np.diff(
            protein_hits
        )
        <
        0
    ).any():

        fail(
            "Raw neighbor-protein hit counts "
            "decreased with increasing radius."
        )


## ================================================================== ##
## 10. Add cluster metadata
## ================================================================== ##

def add_metadata(
    table,
):

    focal_info = cluster_info.rename(
        columns={
            column:
                (
                    "focal_cluster"
                    if column == "cluster"
                    else
                    "focal_" + column
                )

            for column
            in cluster_info.columns
        }
    )


    neighbor_info = cluster_info.rename(
        columns={
            column:
                (
                    "neighbor_cluster"
                    if column == "cluster"
                    else
                    "neighbor_" + column
                )

            for column
            in cluster_info.columns
        }
    )


    return (
        table
        .merge(
            focal_info,
            on="focal_cluster",
            how="left",
            validate="many_to_one",
        )
        .merge(
            neighbor_info,
            on="neighbor_cluster",
            how="left",
            validate="many_to_one",
        )
    )


windows = add_metadata(
    windows
)


## ================================================================== ##
## 11. Deterministic ordering
## ================================================================== ##

windows = (
    windows
    .sort_values(
        [
            "focal_module",
            "focal_cluster",
            "neighbor_cluster",
            "direction",
            "gene_radius",
        ],
        kind="stable",
    )
    .reset_index(
        drop=True
    )
)


## ================================================================== ##
## 12. Write output
## ================================================================== ##

windows.to_csv(
    OUT_WINDOWS,
    sep="\t",
    index=False,
    na_rep="",
    float_format="%.6f",
)


## ================================================================== ##
## 13. QC summary
## ================================================================== ##

qc = pd.DataFrame(
    [
        [
            "focal_module_proteins",
            len(
                focals
            ),
        ],

        [
            "focal_mmseqs_clusters",
            focals[
                "focal_cluster"
            ]
            .nunique(),
        ],

        [
            "mcl_modules",
            focals[
                "focal_module"
            ]
            .nunique(),
        ],

        [
            "mmseqs_neighbor_observations_within_20_genes",
            len(
                mmseqs_neighbors
            ),
        ],

        [
            "distinct_neighbor_mmseqs_clusters",
            mmseqs_neighbors[
                "neighbor_cluster"
            ]
            .nunique(),
        ],

        [
            "distinct_focal_neighbor_associations",
            len(
                candidate_pairs
            ),
        ],

        [
            "gene_radii_tested",
            len(
                GENE_RADII
            ),
        ],

        [
            "gene_window_conservation_rows",
            len(
                windows
            ),
        ],

        [
            "positive_occurrence_denominator_qc_pass",
            1,
        ],

        [
            "positive_genome_denominator_qc_pass",
            1,
        ],

        [
            "cumulative_window_monotonicity_qc_pass",
            1,
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
## 14. Terminal summary
## ================================================================== ##

print()
print(
    "Cumulative MMseqs2 gene-window summary"
)


print(
    f"  Focal proteins:                 "
    f"{len(focals):,}"
)

print(
    f"  Focal clusters:                 "
    f"{focals['focal_cluster'].nunique():,}"
)

print(
    f"  MMseqs2 neighbor observations:  "
    f"{len(mmseqs_neighbors):,}"
)

print(
    f"  Neighbor MMseqs2 clusters:      "
    f"{mmseqs_neighbors['neighbor_cluster'].nunique():,}"
)

print(
    f"  Focal-neighbor associations:    "
    f"{len(candidate_pairs):,}"
)

print(
    f"  Gene-window rows:               "
    f"{len(windows):,}"
)


## ================================================================== ##
## 15. Spotlight relationships
## ================================================================== ##

spotlights = [
    (
        "Cluster_00035",
        "Cluster_00048",
    ),

    (
        "Cluster_00206",
        "Cluster_00228",
    ),
]


for (
    focal_cluster,
    neighbor_cluster,
) in spotlights:

    spotlight = windows[
        (
            windows[
                "focal_cluster"
            ]
            ==
            focal_cluster
        )
        &
        (
            windows[
                "neighbor_cluster"
            ]
            ==
            neighbor_cluster
        )
    ].copy()


    if len(
        spotlight
    ) == 0:

        continue


    print()
    print(
        f"Spotlight: "
        f"{focal_cluster} -> "
        f"{neighbor_cluster}"
    )


    display = spotlight[
        [
            "direction",
            "gene_radius",

            "n_focal_occurrences_with_neighbor",
            "n_informative_focal_occurrences",
            "occurrence_support",
            "pct_informative_occurrences_with_neighbor",

            "n_genomes_with_neighbor",
            "n_informative_genomes",
            "genome_support",
            "pct_informative_genomes_with_neighbor",

            "n_censored_focal_occurrences",
            "n_censored_genomes",
        ]
    ].copy()


    print(
        display.to_string(
            index=False
        )
    )


print()
print("=" * 80)
print("SUCCESS")
print("=" * 80)


print(
    f"Gene-window conservation: "
    f"{OUT_WINDOWS}"
)

print(
    f"QC: "
    f"{OUT_QC}"
)
