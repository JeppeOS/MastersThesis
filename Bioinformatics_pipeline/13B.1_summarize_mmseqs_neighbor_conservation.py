#!/usr/bin/env python3

from pathlib import Path
import sys

import numpy as np
import pandas as pd


## ================================================================== ##
## STAGE 13B
##
## MMseqs2-family neighborhood conservation around every focal
## MCL-module protein family.
##
## PRIMARY UNIT:
##
##     MMseqs2 protein-family cluster
##
## GlobDB COG/product annotations are deliberately NOT used to define
## conserved neighborhood units in this stage.
##
## Two complementary analyses are produced:
##
##   1. exact oriented gene-order conservation
##   2. directional physical-window conservation
##
## Focal transcription is normalized:
##
##     negative offsets = upstream
##     positive offsets = downstream
##
## Physical-window denominators are censoring-aware:
##
##     positive hit
##         -> informative even if the contig is truncated farther out
##
##     no hit + complete requested window
##         -> informative negative
##
##     no hit + incomplete requested window
##         -> censored / unknown
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

OUT_GENE_ORDER = (
    HERE
    / "cluster_neighbor_gene_order_conservation.tsv"
)


OUT_BP_WINDOWS = (
    HERE
    / "cluster_neighbor_bp_window_conservation.tsv"
)


OUT_ASSOCIATION_SUMMARY = (
    HERE
    / "cluster_neighbor_association_summary.tsv"
)


OUT_QC = (
    HERE
    / "mmseqs_neighbor_conservation_qc.tsv"
)


## ================================================================== ##
## Locked analysis parameters
## ================================================================== ##

MAX_GENE_OFFSET = 20

BP_RADII = [
    5_000,
    10_000,
    20_000,
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


def radius_label(radius):

    return (
        f"{radius // 1000}kb"
    )


## ================================================================== ##
## Start
## ================================================================== ##

print("=" * 80)
print("STAGE 13B - MMSEQS2 NEIGHBOR CONSERVATION")
print("=" * 80)


## ================================================================== ##
## 1. Read focal occurrence / observability information
## ================================================================== ##

print()
print("Reading focal observability...")


focals = pd.read_csv(
    OBSERVABILITY,
    sep="\t",
    dtype=str,
    keep_default_na=False,
)


required_focal_columns = [
    "genome",
    "focal_protein_id",
    "focal_cluster",
    "focal_module",

    "focal_strand",

    "n_oriented_upstream_genes_available",
    "n_oriented_downstream_genes_available",

    "upstream_5kb_observable",
    "downstream_5kb_observable",

    "upstream_10kb_observable",
    "downstream_10kb_observable",

    "upstream_20kb_observable",
    "downstream_20kb_observable",
]


require_columns(
    focals,
    required_focal_columns,
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
        "Duplicate focal occurrence keys."
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
        f"focal clusters but found "
        f"{focals['focal_cluster'].nunique()}."
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
        f"modules but found "
        f"{focals['focal_module'].nunique()}."
    )


for column in [
    "n_oriented_upstream_genes_available",
    "n_oriented_downstream_genes_available",

    "upstream_5kb_observable",
    "downstream_5kb_observable",

    "upstream_10kb_observable",
    "downstream_10kb_observable",

    "upstream_20kb_observable",
    "downstream_20kb_observable",
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
print("Reading observed neighborhood table...")


neigh = pd.read_csv(
    NEIGHBORHOODS,
    sep="\t",
    dtype=str,
    keep_default_na=False,
)


required_neighbor_columns = [
    "genome",

    "focal_protein_id",
    "focal_cluster",
    "focal_module",

    "neighbor_protein_id",
    "neighbor_cluster",
    "neighbor_module",

    "oriented_gene_offset",
    "oriented_midpoint_offset_bp",

    "intergenic_gap_bp",

    "is_focal",
]


require_columns(
    neigh,
    required_neighbor_columns,
    NEIGHBORHOODS.name,
)


for column in [
    "oriented_gene_offset",
    "intergenic_gap_bp",
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
    "oriented_midpoint_offset_bp"
] = pd.to_numeric(
    neigh[
        "oriented_midpoint_offset_bp"
    ],
    errors="raise",
)


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


## Every neighborhood focal must exist in the focal table. ##

unknown_focal_keys = (
    set(
        neigh[
            "focal_key"
        ]
    )
    -
    set(
        focals[
            "focal_key"
        ]
    )
)


if unknown_focal_keys:

    fail(
        f"{len(unknown_focal_keys):,} "
        f"neighborhood focal keys are absent "
        f"from the observability table."
    )


## ================================================================== ##
## 3. Retain MMseqs2-defined neighboring protein families
##
## Exclude only the focal-self row.
##
## A nearby paralogue belonging to the SAME MMseqs cluster as the
## focal is retained if it occurs at a non-zero gene offset.
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
].copy()


print(
    f"  MMseqs2 neighbor observations: "
    f"{len(mmseqs_neighbors):,}"
)


print(
    f"  Distinct neighbor clusters:    "
    f"{mmseqs_neighbors['neighbor_cluster'].nunique():,}"
)


## ================================================================== ##
## 4. Read cluster/module information
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
        f"MMseqs clusters but found "
        f"{cluster_summary['cluster'].nunique():,}."
    )


## ================================================================== ##
## Cluster metadata propagated as interpretation layers.
##
## Sequence-family identity remains cluster ID.
## ================================================================== ##

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
## 5. Exact gene-order conservation
##
## One output row represents:
##
##     focal_cluster
##     x neighbor_cluster
##     x exact oriented gene offset
##
## Example:
##
##     Cluster_00035
##     Cluster_00048
##     offset = -1
##
## means Cluster_00048 is immediately upstream of the focal,
## after normalizing focal transcription direction.
## ================================================================== ##

print()
print("Calculating exact gene-order conservation...")


gene_order_hits = mmseqs_neighbors[
    (
        mmseqs_neighbors[
            "oriented_gene_offset"
        ]
        !=
        0
    )
    &
    (
        mmseqs_neighbors[
            "oriented_gene_offset"
        ]
        .abs()
        <=
        MAX_GENE_OFFSET
    )
].copy()


## A focal occurrence can have at most one gene at one exact offset.
##
## Therefore focal + offset + neighbor_cluster should not duplicate. ##

if (
    gene_order_hits[
        [
            "focal_key",
            "oriented_gene_offset",
            "neighbor_cluster",
        ]
    ]
    .duplicated()
    .any()
):

    fail(
        "Duplicate focal/offset/neighbor-cluster "
        "gene-order observations detected."
    )


gene_order_rows = []


gene_order_groups = gene_order_hits.groupby(
    [
        "focal_cluster",
        "neighbor_cluster",
        "oriented_gene_offset",
    ],
    sort=True,
)


for (
    focal_cluster,
    neighbor_cluster,
    offset,
), hits in gene_order_groups:

    offset = int(
        offset
    )


    focal_subset = focals[
        focals[
            "focal_cluster"
        ]
        ==
        focal_cluster
    ].copy()


    if offset < 0:

        observable = (
            focal_subset[
                "n_oriented_upstream_genes_available"
            ]
            >=
            abs(
                offset
            )
        )

        direction = (
            "upstream"
        )

    else:

        observable = (
            focal_subset[
                "n_oriented_downstream_genes_available"
            ]
            >=
            abs(
                offset
            )
        )

        direction = (
            "downstream"
        )


    focal_subset[
        "offset_observable"
    ] = observable.astype(int)


    positive_keys = set(
        hits[
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


    ## A positive at the exact offset must necessarily be observable. ##

    impossible_positive = focal_subset[
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
                "offset_observable"
            ]
            ==
            0
        )
    ]


    if len(impossible_positive) > 0:

        fail(
            f"Gene-order positive occurred at "
            f"unobservable offset for "
            f"{focal_cluster} / "
            f"{neighbor_cluster} / "
            f"{offset}."
        )


    n_occ_total = len(
        focal_subset
    )


    n_occ_observable = int(
        focal_subset[
            "offset_observable"
        ]
        .sum()
    )


    n_occ_positive = int(
        focal_subset[
            "positive"
        ]
        .sum()
    )


    ## -------------------------------------------------------------- ##
    ## Genome-level censoring logic
    ##
    ## A genome is positive if ANY focal copy has the relationship.
    ##
    ## A genome is an informative negative only if ALL focal copies
    ## are observable at this offset and NONE is positive.
    ##
    ## Otherwise the genome is censored.
    ## -------------------------------------------------------------- ##

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

            n_observable_occurrences=(
                "offset_observable",
                "sum"
            ),

            n_positive_occurrences=(
                "positive",
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
            "n_observable_occurrences"
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


    n_genomes_censored = (
        n_genomes_total
        -
        n_genomes_informative
    )


    focal_module = (
        focal_subset[
            "focal_module"
        ]
        .iloc[0]
    )


    neighbor_module = module_lookup.get(
        neighbor_cluster,
        "",
    )


    gene_order_rows.append(
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

            "oriented_gene_offset":
                offset,

            "direction":
                direction,

            "absolute_gene_offset":
                abs(
                    offset
                ),

            ## Occurrence level ##

            "n_focal_occurrences_total":
                n_occ_total,

            "n_focal_occurrences_offset_observable":
                n_occ_observable,

            "n_focal_occurrences_with_neighbor_at_offset":
                n_occ_positive,

            "pct_observable_occurrences_with_neighbor_at_offset":
                pct(
                    n_occ_positive,
                    n_occ_observable,
                ),

            ## Genome level ##

            "n_focal_genomes_total":
                n_genomes_total,

            "n_informative_genomes":
                n_genomes_informative,

            "n_censored_genomes":
                n_genomes_censored,

            "n_genomes_with_neighbor_at_offset":
                n_genomes_positive,

            "pct_informative_genomes_with_neighbor_at_offset":
                pct(
                    n_genomes_positive,
                    n_genomes_informative,
                ),
        }
    )


gene_order = pd.DataFrame(
    gene_order_rows
)


## ================================================================== ##
## 6. Physical-window conservation
##
## We evaluate directional windows separately:
##
##     upstream 5 kb
##     downstream 5 kb
##     upstream 10 kb
##     ...
##
## Distance is the intergenic gap between CDS boundaries.
##
## Positive observations remain valid on partially assembled contigs.
## Only negative calls require complete window observability.
## ================================================================== ##

print()
print("Calculating physical-window conservation...")


## Direction based on normalized GENE ORDER rather than midpoint,
## which remains robust even for overlapping CDSs. ##

physical_hits = mmseqs_neighbors[
    mmseqs_neighbors[
        "oriented_gene_offset"
    ]
    !=
    0
].copy()


physical_hits[
    "direction"
] = np.where(
    physical_hits[
        "oriented_gene_offset"
    ]
    <
    0,
    "upstream",
    "downstream",
)


bp_rows = []


## Only analyze focal/neighbor combinations that are actually observed
## somewhere within our Stage-12 20-kb physical envelope. ##

candidate_pairs = (
    physical_hits[
        [
            "focal_cluster",
            "neighbor_cluster",
        ]
    ]
    .drop_duplicates()
    .sort_values(
        [
            "focal_cluster",
            "neighbor_cluster",
        ]
    )
)


for candidate in candidate_pairs.itertuples(
    index=False
):

    focal_cluster = (
        candidate.focal_cluster
    )

    neighbor_cluster = (
        candidate.neighbor_cluster
    )


    focal_subset_base = focals[
        focals[
            "focal_cluster"
        ]
        ==
        focal_cluster
    ].copy()


    focal_module = (
        focal_subset_base[
            "focal_module"
        ]
        .iloc[0]
    )


    neighbor_module = module_lookup.get(
        neighbor_cluster,
        "",
    )


    pair_hits = physical_hits[
        (
            physical_hits[
                "focal_cluster"
            ]
            ==
            focal_cluster
        )
        &
        (
            physical_hits[
                "neighbor_cluster"
            ]
            ==
            neighbor_cluster
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


        for radius in BP_RADII:

            label = radius_label(
                radius
            )


            complete_column = (
                f"{direction}_{label}_observable"
            )


            if complete_column not in focal_subset_base.columns:

                fail(
                    f"Missing observability column "
                    f"{complete_column}."
                )


            focal_subset = (
                focal_subset_base
                .copy()
            )


            ## ------------------------------------------------------ ##
            ## Positive = at least one neighbor-family member exists
            ## within the requested directional intergenic gap.
            ## ------------------------------------------------------ ##

            hits_in_window = direction_hits[
                direction_hits[
                    "intergenic_gap_bp"
                ]
                <=
                radius
            ]


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


            focal_subset[
                "full_window_observable"
            ] = (
                focal_subset[
                    complete_column
                ]
                .astype(int)
            )


            ## ------------------------------------------------------ ##
            ## Informative occurrence:
            ##
            ##     positive
            ##
            ## OR
            ##
            ##     no hit but complete window observed
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


            n_occ_censored = (
                n_occ_total
                -
                n_occ_informative
            )


            ## ------------------------------------------------------ ##
            ## Genome-level censoring logic:
            ##
            ## positive genome:
            ##     any focal copy is positive
            ##
            ## informative negative genome:
            ##     no positive copies
            ##     AND every focal copy has a complete window
            ##
            ## otherwise:
            ##     censored genome
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

                    n_partial_window_positive_occurrences=(
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


            n_genomes_censored = (
                n_genomes_total
                -
                n_genomes_informative
            )


            n_genomes_positive_from_partial = int(
                (
                    (
                        genome_status[
                            "positive_genome"
                        ]
                    )
                    &
                    (
                        genome_status[
                            "n_partial_window_positive_occurrences"
                        ]
                        >
                        0
                    )
                )
                .sum()
            )


            bp_rows.append(
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

                    "radius_bp":
                        radius,

                    "radius_label":
                        label,

                    ## Occurrence level ##

                    "n_focal_occurrences_total":
                        n_occ_total,

                    "n_focal_occurrences_full_window_observable":
                        n_occ_full,

                    "n_informative_focal_occurrences":
                        n_occ_informative,

                    "n_censored_focal_occurrences":
                        n_occ_censored,

                    "n_focal_occurrences_with_neighbor":
                        n_occ_positive,

                    "n_positive_occurrences_from_partial_windows":
                        n_occ_partial_positive,

                    "pct_informative_occurrences_with_neighbor":
                        pct(
                            n_occ_positive,
                            n_occ_informative,
                        ),

                    ## Genome level ##

                    "n_focal_genomes_total":
                        n_genomes_total,

                    "n_informative_genomes":
                        n_genomes_informative,

                    "n_censored_genomes":
                        n_genomes_censored,

                    "n_genomes_with_neighbor":
                        n_genomes_positive,

                    "n_positive_genomes_with_partial_window_hit":
                        n_genomes_positive_from_partial,

                    "pct_informative_genomes_with_neighbor":
                        pct(
                            n_genomes_positive,
                            n_genomes_informative,
                        ),
                }
            )


bp_windows = pd.DataFrame(
    bp_rows
)


## ================================================================== ##
## 7. Add MMseqs2 cluster metadata
##
## Cluster IDs remain the primary identity.
## These annotations are display/interpretation fields only.
## ================================================================== ##

def add_cluster_metadata(
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


    table = (
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


    return table


gene_order = add_cluster_metadata(
    gene_order
)


bp_windows = add_cluster_metadata(
    bp_windows
)


## ================================================================== ##
## 8. Deterministic ordering
## ================================================================== ##

gene_order = (
    gene_order
    .sort_values(
        [
            "focal_module",
            "focal_cluster",
            "oriented_gene_offset",
            "neighbor_cluster",
        ],
        kind="stable",
    )
    .reset_index(
        drop=True
    )
)


bp_windows = (
    bp_windows
    .sort_values(
        [
            "focal_module",
            "focal_cluster",
            "neighbor_cluster",
            "direction",
            "radius_bp",
        ],
        kind="stable",
    )
    .reset_index(
        drop=True
    )
)


## ================================================================== ##
## 9. Write primary tables
## ================================================================== ##

gene_order.to_csv(
    OUT_GENE_ORDER,
    sep="\t",
    index=False,
    na_rep="",
    float_format="%.6f",
)


bp_windows.to_csv(
    OUT_BP_WINDOWS,
    sep="\t",
    index=False,
    na_rep="",
    float_format="%.6f",
)


## ================================================================== ##
## 10. Compact association summary
##
## One row per observed focal-cluster / neighbor-cluster combination.
##
## This provides convenient ranking without imposing any
## "conserved" threshold.
## ================================================================== ##

print()
print("Building compact MMseqs association summary...")


association_keys = (
    pd.concat(
        [
            gene_order[
                [
                    "focal_cluster",
                    "focal_module",
                    "neighbor_cluster",
                    "neighbor_module",
                    "neighbor_same_mcl_module_as_focal",
                    "neighbor_same_mmseqs_cluster_as_focal",
                ]
            ],

            bp_windows[
                [
                    "focal_cluster",
                    "focal_module",
                    "neighbor_cluster",
                    "neighbor_module",
                    "neighbor_same_mcl_module_as_focal",
                    "neighbor_same_mmseqs_cluster_as_focal",
                ]
            ],
        ],
        ignore_index=True,
    )
    .drop_duplicates()
)


summary_rows = []


for row in association_keys.itertuples(
    index=False
):

    focal_cluster = row.focal_cluster
    neighbor_cluster = row.neighbor_cluster


    go = gene_order[
        (
            gene_order[
                "focal_cluster"
            ]
            ==
            focal_cluster
        )
        &
        (
            gene_order[
                "neighbor_cluster"
            ]
            ==
            neighbor_cluster
        )
    ].copy()


    bp = bp_windows[
        (
            bp_windows[
                "focal_cluster"
            ]
            ==
            focal_cluster
        )
        &
        (
            bp_windows[
                "neighbor_cluster"
            ]
            ==
            neighbor_cluster
        )
    ].copy()


    ## Best exact gene-order position by occurrence-level frequency. ##

    if len(go) > 0:

        go_best = (
            go
            .sort_values(
                [
                    "pct_observable_occurrences_with_neighbor_at_offset",
                    "n_focal_occurrences_with_neighbor_at_offset",
                    "absolute_gene_offset",
                ],
                ascending=[
                    False,
                    False,
                    True,
                ],
                kind="stable",
            )
            .iloc[0]
        )


        best_offset = int(
            go_best[
                "oriented_gene_offset"
            ]
        )


        best_offset_pct = float(
            go_best[
                "pct_observable_occurrences_with_neighbor_at_offset"
            ]
        )


        best_offset_n = int(
            go_best[
                "n_focal_occurrences_with_neighbor_at_offset"
            ]
        )

    else:

        best_offset = np.nan
        best_offset_pct = np.nan
        best_offset_n = 0


    def bp_value(
        direction,
        radius,
        column,
    ):

        subset = bp[
            (
                bp[
                    "direction"
                ]
                ==
                direction
            )
            &
            (
                bp[
                    "radius_bp"
                ]
                ==
                radius
            )
        ]


        if len(subset) == 0:
            return np.nan


        return subset.iloc[0][
            column
        ]


    summary_rows.append(
        {
            "focal_cluster":
                focal_cluster,

            "focal_module":
                row.focal_module,

            "neighbor_cluster":
                neighbor_cluster,

            "neighbor_module":
                row.neighbor_module,

            "neighbor_same_mcl_module_as_focal":
                row.neighbor_same_mcl_module_as_focal,

            "neighbor_same_mmseqs_cluster_as_focal":
                row.neighbor_same_mmseqs_cluster_as_focal,

            "best_oriented_gene_offset":
                best_offset,

            "best_offset_positive_occurrences":
                best_offset_n,

            "best_offset_pct_observable_occurrences":
                best_offset_pct,

            "pct_occurrences_upstream_5kb":
                bp_value(
                    "upstream",
                    5_000,
                    "pct_informative_occurrences_with_neighbor",
                ),

            "pct_occurrences_upstream_10kb":
                bp_value(
                    "upstream",
                    10_000,
                    "pct_informative_occurrences_with_neighbor",
                ),

            "pct_occurrences_upstream_20kb":
                bp_value(
                    "upstream",
                    20_000,
                    "pct_informative_occurrences_with_neighbor",
                ),

            "pct_occurrences_downstream_5kb":
                bp_value(
                    "downstream",
                    5_000,
                    "pct_informative_occurrences_with_neighbor",
                ),

            "pct_occurrences_downstream_10kb":
                bp_value(
                    "downstream",
                    10_000,
                    "pct_informative_occurrences_with_neighbor",
                ),

            "pct_occurrences_downstream_20kb":
                bp_value(
                    "downstream",
                    20_000,
                    "pct_informative_occurrences_with_neighbor",
                ),

            "pct_genomes_upstream_10kb":
                bp_value(
                    "upstream",
                    10_000,
                    "pct_informative_genomes_with_neighbor",
                ),

            "pct_genomes_downstream_10kb":
                bp_value(
                    "downstream",
                    10_000,
                    "pct_informative_genomes_with_neighbor",
                ),
        }
    )


association_summary = pd.DataFrame(
    summary_rows
)


association_summary = add_cluster_metadata(
    association_summary
)


association_summary = (
    association_summary
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


association_summary.to_csv(
    OUT_ASSOCIATION_SUMMARY,
    sep="\t",
    index=False,
    na_rep="",
    float_format="%.6f",
)


## ================================================================== ##
## 11. QC
## ================================================================== ##

n_focal_clusters_gene_order = (
    gene_order[
        "focal_cluster"
    ]
    .nunique()
)


n_focal_clusters_bp = (
    bp_windows[
        "focal_cluster"
    ]
    .nunique()
)


## Every positive physical observation must be represented by an
## informative denominator. ##

bad_bp = bp_windows[
    (
        bp_windows[
            "n_focal_occurrences_with_neighbor"
        ]
        >
        bp_windows[
            "n_informative_focal_occurrences"
        ]
    )
]


if len(bad_bp) > 0:

    fail(
        "Positive physical-window counts exceed "
        "informative denominators."
    )


bad_genome_bp = bp_windows[
    (
        bp_windows[
            "n_genomes_with_neighbor"
        ]
        >
        bp_windows[
            "n_informative_genomes"
        ]
    )
]


if len(bad_genome_bp) > 0:

    fail(
        "Positive genome counts exceed "
        "informative genome denominators."
    )


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
            "mmseqs_neighbor_observations",
            len(
                mmseqs_neighbors
            ),
        ],

        [
            "distinct_mmseqs_neighbor_clusters",
            mmseqs_neighbors[
                "neighbor_cluster"
            ]
            .nunique(),
        ],

        [
            "gene_order_conservation_rows",
            len(
                gene_order
            ),
        ],

        [
            "gene_order_focal_clusters_represented",
            n_focal_clusters_gene_order,
        ],

        [
            "bp_window_conservation_rows",
            len(
                bp_windows
            ),
        ],

        [
            "bp_window_focal_clusters_represented",
            n_focal_clusters_bp,
        ],

        [
            "compact_association_rows",
            len(
                association_summary
            ),
        ],

        [
            "physical_window_positive_denominator_qc_pass",
            1,
        ],

        [
            "physical_window_genome_denominator_qc_pass",
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
## 12. Terminal summary
## ================================================================== ##

print()
print("MMseqs2 neighborhood-conservation summary")


print(
    f"  Focal proteins:                 "
    f"{len(focals):,}"
)

print(
    f"  Focal clusters:                 "
    f"{focals['focal_cluster'].nunique():,}"
)

print(
    f"  MMseqs neighbor observations:   "
    f"{len(mmseqs_neighbors):,}"
)

print(
    f"  Neighbor MMseqs clusters:       "
    f"{mmseqs_neighbors['neighbor_cluster'].nunique():,}"
)

print(
    f"  Gene-order summary rows:        "
    f"{len(gene_order):,}"
)

print(
    f"  Physical-window summary rows:   "
    f"{len(bp_windows):,}"
)

print(
    f"  Compact associations:           "
    f"{len(association_summary):,}"
)


## ================================================================== ##
## 13. Spotlight known biologically anchored relationships
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

    print()
    print(
        f"Spotlight: "
        f"{focal_cluster} -> "
        f"{neighbor_cluster}"
    )


    go = gene_order[
        (
            gene_order[
                "focal_cluster"
            ]
            ==
            focal_cluster
        )
        &
        (
            gene_order[
                "neighbor_cluster"
            ]
            ==
            neighbor_cluster
        )
    ]


    if len(go) > 0:

        display_go = (
            go[
                [
                    "oriented_gene_offset",
                    "direction",
                    "n_focal_occurrences_offset_observable",
                    "n_focal_occurrences_with_neighbor_at_offset",
                    "pct_observable_occurrences_with_neighbor_at_offset",
                    "n_informative_genomes",
                    "n_genomes_with_neighbor_at_offset",
                    "pct_informative_genomes_with_neighbor_at_offset",
                ]
            ]
            .sort_values(
                "oriented_gene_offset"
            )
        )


        print()
        print(
            "  Exact oriented gene offsets:"
        )

        print(
            display_go.to_string(
                index=False
            )
        )


    bp = bp_windows[
        (
            bp_windows[
                "focal_cluster"
            ]
            ==
            focal_cluster
        )
        &
        (
            bp_windows[
                "neighbor_cluster"
            ]
            ==
            neighbor_cluster
        )
    ]


    if len(bp) > 0:

        display_bp = (
            bp[
                [
                    "direction",
                    "radius_bp",
                    "radius_label",
                    "n_informative_focal_occurrences",
                    "n_focal_occurrences_with_neighbor",
                    "pct_informative_occurrences_with_neighbor",
                    "n_informative_genomes",
                    "n_genomes_with_neighbor",
                    "pct_informative_genomes_with_neighbor",
                ]
            ]
            .sort_values(
                [
                    "direction",
                    "radius_bp",
                ],
                kind="stable",
            )
        )

        print()
        print(
            "  Physical windows:"
        )

        print(
            display_bp
            .drop(
                columns="radius_bp"
            )
            .to_string(
                index=False
            )
        )


print()
print("=" * 80)
print("SUCCESS")
print("=" * 80)


print(
    f"Gene-order conservation: "
    f"{OUT_GENE_ORDER}"
)

print(
    f"Physical-window conservation: "
    f"{OUT_BP_WINDOWS}"
)

print(
    f"Compact association summary: "
    f"{OUT_ASSOCIATION_SUMMARY}"
)

print(
    f"QC: "
    f"{OUT_QC}"
)