#!/usr/bin/env python3

from pathlib import Path
import sys

import numpy as np
import pandas as pd


## ================================================================== ##
## STAGE 13D - EVIDENCE-PRESERVING CONSENSUS NEIGHBORHOODS
##
## PRIMARY BIOLOGICAL FRAMEWORK:
##
##     focal MCL-associated MMseqs2 protein families
##
## SECONDARY INTERPRETATION:
##
##     GlobDB COG / product
##
##
## IMPORTANT:
##
## This script DOES NOT impose a conserved/not-conserved threshold.
##
## Instead it preserves:
##
##     positive numerator
##     informative denominator
##     percentage
##     censoring
##     exact position
##     cumulative gene-window position
##     physical bp-window position
##
##
## Three outputs:
##
## 1. focal_consensus_position_summary.tsv
##
##    One row for every:
##
##        focal MMseqs cluster x oriented offset -20 ... +20
##
##    Includes:
##
##        dominant MMseqs neighbor family
##        dominant COG
##        product label
##        annotation resolution
##        same-strand support
##
##
## 2. focal_mmseqs_consensus_associations.tsv
##
##    One row for every:
##
##        focal cluster
##        x neighbor MMseqs cluster
##        x direction
##
##    Includes:
##
##        exact best offset
##        smallest cumulative gene radius at which the
##        maximum observed support is reached
##        5/10/20-kb physical support
##
##
## 3. focal_cog_consensus_associations.tsv
##
##    Same idea for COG context.
##
##
## "Saturation radius" means:
##
##     the smallest tested radius at which the final observed
##     numerator is reached.
##
## Example:
##
##     +1 = 2
##     +2 = 6
##     +3 = 7
##     +5 = 7
##
## gives:
##
##     occurrence_saturation_gene_radius = 3
##
## This captures flexible local architecture without hiding the
## original exact-offset observations.
## ================================================================== ##


## ================================================================== ##
## Paths
## ================================================================== ##

HERE = Path(__file__).resolve().parent
WORKFLOW = HERE.parent.parent


OBSERVABILITY = (
    WORKFLOW
    / "12_neighborhood_extraction"
    / "focal_context_observability.tsv"
)


NEIGHBORHOODS = (
    WORKFLOW
    / "12_neighborhood_extraction"
    / "observed_neighborhood_genes.tsv"
)


MMSEQ_EXACT = (
    WORKFLOW
    / "13_neighborhood_conservation"
    / "13B_mmseqs_neighbor_conservation"
    / "cluster_neighbor_gene_order_conservation.tsv"
)


MMSEQ_GENE_WINDOWS = (
    WORKFLOW
    / "13_neighborhood_conservation"
    / "13B_mmseqs_neighbor_conservation"
    / "cluster_neighbor_gene_window_conservation.tsv"
)


MMSEQ_BP_WINDOWS = (
    WORKFLOW
    / "13_neighborhood_conservation"
    / "13B_mmseqs_neighbor_conservation"
    / "cluster_neighbor_bp_window_conservation.tsv"
)


COG_EXACT = (
    WORKFLOW
    / "13_neighborhood_conservation"
    / "13C_annotation_context"
    / "focal_cog_gene_order_context.tsv"
)


COG_WINDOWS = (
    WORKFLOW
    / "13_neighborhood_conservation"
    / "13C_annotation_context"
    / "focal_cog_gene_window_context.tsv"
)


## ================================================================== ##
## Outputs
## ================================================================== ##

OUT_POSITIONS = (
    HERE
    / "focal_consensus_position_summary.tsv"
)


OUT_MMSEQ = (
    HERE
    / "focal_mmseqs_consensus_associations.tsv"
)


OUT_COG = (
    HERE
    / "focal_cog_consensus_associations.tsv"
)


OUT_QC = (
    HERE
    / "consensus_neighborhood_qc.tsv"
)


## ================================================================== ##
## Locked dimensions
## ================================================================== ##

MAX_GENE_OFFSET = 20


GENE_RADII = [
    1,
    2,
    3,
    5,
    10,
    20,
]


BP_RADII = [
    5_000,
    10_000,
    20_000,
]


EXPECTED_FOCAL_PROTEINS = 10_537
EXPECTED_FOCAL_CLUSTERS = 156
EXPECTED_MODULES = 35
EXPECTED_POSITION_ROWS = 156 * 41


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
            f"{source} missing required column(s): "
            +
            ", ".join(
                sorted(missing)
            )
        )


def pct(
    numerator,
    denominator,
):

    numerator = pd.to_numeric(
        numerator,
        errors="coerce",
    )

    denominator = pd.to_numeric(
        denominator,
        errors="coerce",
    )


    return np.where(
        denominator > 0,
        100.0
        *
        numerator
        /
        denominator,
        np.nan,
    )


def support_string(
    numerator,
    denominator,
):

    numerator = (
        pd.Series(
            numerator
        )
        .fillna(
            0
        )
        .astype(int)
        .astype(str)
    )


    denominator = (
        pd.Series(
            denominator
        )
        .fillna(
            0
        )
        .astype(int)
        .astype(str)
    )


    return (
        numerator
        +
        "/"
        +
        denominator
    )


def top_label_summary(
    df,
    group_columns,
    label_column,
    output_prefix,
):

    ## -------------------------------------------------------------- ##
    ## Return:
    ##
    ##     dominant_<prefix>
    ##     n_positions_with_dominant_<prefix>
    ##     n_distinct_<prefix>s
    ##
    ## Ties are resolved lexicographically for deterministic output.
    ## -------------------------------------------------------------- ##

    valid = df[
        df[
            label_column
        ]
        .astype(str)
        .str.strip()
        !=
        ""
    ].copy()


    if len(
        valid
    ) == 0:

        return pd.DataFrame(
            columns=(
                group_columns
                +
                [
                    f"dominant_{output_prefix}",
                    f"n_positions_with_dominant_{output_prefix}",
                    f"n_distinct_{output_prefix}s",
                ]
            )
        )


    counts = (
        valid
        .groupby(
            group_columns
            +
            [
                label_column,
            ],
            as_index=False,
        )
        .size()
        .rename(
            columns={
                "size":
                    "label_count",
            }
        )
    )


    distinct = (
        counts
        .groupby(
            group_columns,
            as_index=False,
        )
        .size()
        .rename(
            columns={
                "size":
                    f"n_distinct_{output_prefix}s",
            }
        )
    )


    counts = counts.sort_values(
        group_columns
        +
        [
            "label_count",
            label_column,
        ],
        ascending=(
            [True]
            *
            len(
                group_columns
            )
            +
            [
                False,
                True,
            ]
        ),
        kind="stable",
    )


    top = (
        counts
        .drop_duplicates(
            group_columns
        )
        [
            group_columns
            +
            [
                label_column,
                "label_count",
            ]
        ]
        .rename(
            columns={
                label_column:
                    f"dominant_{output_prefix}",

                "label_count":
                    f"n_positions_with_dominant_{output_prefix}",
            }
        )
    )


    return top.merge(
        distinct,
        on=group_columns,
        how="left",
        validate="one_to_one",
    )


def first_radius_at_max(
    df,
    group_columns,
    positive_column,
    radius_column,
):

    ## -------------------------------------------------------------- ##
    ## Smallest radius at which the maximum observed numerator is
    ## reached.
    ##
    ## Returns the entire selected row.
    ## -------------------------------------------------------------- ##

    work = df.copy()


    work[
        "_maximum_positive"
    ] = (
        work
        .groupby(
            group_columns
        )[
            positive_column
        ]
        .transform(
            "max"
        )
    )


    candidate = work[
        work[
            positive_column
        ]
        ==
        work[
            "_maximum_positive"
        ]
    ].copy()


    candidate = (
        candidate
        .sort_values(
            group_columns
            +
            [
                radius_column,
            ],
            kind="stable",
        )
        .drop_duplicates(
            group_columns
        )
        .reset_index(
            drop=True
        )
    )


    return candidate


def best_exact_offset(
    df,
    group_columns,
    positive_column,
):

    ## -------------------------------------------------------------- ##
    ## Choose exact offset with:
    ##
    ##   1. largest numerator
    ##   2. smallest absolute offset
    ##   3. signed offset for deterministic final tie break
    ## -------------------------------------------------------------- ##

    work = df.copy()


    return (
        work
        .sort_values(
            group_columns
            +
            [
                positive_column,
                "absolute_gene_offset",
                "oriented_gene_offset",
            ],
            ascending=(
                [True]
                *
                len(
                    group_columns
                )
                +
                [
                    False,
                    True,
                    True,
                ]
            ),
            kind="stable",
        )
        .drop_duplicates(
            group_columns
        )
        .reset_index(
            drop=True
        )
    )


## ================================================================== ##
## Start
## ================================================================== ##

print("=" * 80)
print("STAGE 13D - EVIDENCE-PRESERVING CONSENSUS NEIGHBORHOODS")
print("=" * 80)


## ================================================================== ##
## 1. Read focal occurrence information
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


if len(
    focals
) != EXPECTED_FOCAL_PROTEINS:

    fail(
        f"Expected "
        f"{EXPECTED_FOCAL_PROTEINS:,} "
        f"focal proteins but found "
        f"{len(focals):,}."
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
        "Unexpected focal-cluster count."
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
        "Unexpected module count."
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


if focals[
    "focal_key"
].duplicated().any():

    fail(
        "Duplicate focal keys."
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
## 2. Read complete observed neighborhood rows
## ================================================================== ##

print()
print("Reading observed genomic neighborhoods...")


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

        "oriented_gene_offset",

        "neighbor_protein_id",
        "neighbor_cluster",

        "neighbor_annotation_accepted",
        "neighbor_globdb_cog",
        "neighbor_globdb_product",

        "neighbor_same_strand_as_focal",
        "is_focal",
    ],
    NEIGHBORHOODS.name,
)


for column in [
    "oriented_gene_offset",
    "neighbor_annotation_accepted",
    "neighbor_same_strand_as_focal",
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


neigh[
    "neighbor_cluster"
] = (
    neigh[
        "neighbor_cluster"
    ]
    .astype(str)
    .str.strip()
)


neigh[
    "neighbor_cog"
] = (
    neigh[
        "neighbor_globdb_cog"
    ]
    .astype(str)
    .str.strip()
)


neigh[
    "neighbor_product"
] = (
    neigh[
        "neighbor_globdb_product"
    ]
    .astype(str)
    .str.strip()
)


## Annotation evidence should only exist on accepted GlobDB mappings. ##

bad_annotation = neigh[
    (
        (
            neigh[
                "neighbor_cog"
            ]
            !=
            ""
        )
        |
        (
            neigh[
                "neighbor_product"
            ]
            !=
            ""
        )
    )
    &
    (
        neigh[
            "neighbor_annotation_accepted"
        ]
        !=
        1
    )
]


if len(
    bad_annotation
) > 0:

    fail(
        f"{len(bad_annotation):,} "
        f"GlobDB annotations occur on non-accepted mappings."
    )


## Keep exact gene-order space only. ##

position_rows = neigh[
    neigh[
        "oriented_gene_offset"
    ]
    .abs()
    <=
    MAX_GENE_OFFSET
].copy()


position_rows[
    "has_mmseqs_cluster"
] = (
    position_rows[
        "neighbor_cluster"
    ]
    !=
    ""
).astype(int)


position_rows[
    "has_cog"
] = (
    position_rows[
        "neighbor_cog"
    ]
    !=
    ""
).astype(int)


position_rows[
    "has_product"
] = (
    position_rows[
        "neighbor_product"
    ]
    !=
    ""
).astype(int)


position_rows[
    "without_mmseqs_or_cog"
] = (
    (
        position_rows[
            "neighbor_cluster"
        ]
        ==
        ""
    )
    &
    (
        position_rows[
            "neighbor_cog"
        ]
        ==
        ""
    )
).astype(int)


print(
    f"  Exact-position neighborhood rows: "
    f"{len(position_rows):,}"
)


## ================================================================== ##
## 3. Build complete focal x offset observability grid
##
## This guarantees that positions with zero observable sequence still
## appear explicitly in the consensus output.
## ================================================================== ##

print()
print("Building exact-position observability grid...")


offsets = pd.DataFrame(
    {
        "oriented_gene_offset":
            list(
                range(
                    -MAX_GENE_OFFSET,
                    MAX_GENE_OFFSET + 1,
                )
            )
    }
)


status = focals[
    [
        "focal_key",
        "genome",
        "focal_cluster",
        "focal_module",
        "n_oriented_upstream_genes_available",
        "n_oriented_downstream_genes_available",
    ]
].merge(
    offsets,
    how="cross",
)


status[
    "offset_observable"
] = np.select(
    [
        status[
            "oriented_gene_offset"
        ]
        ==
        0,

        status[
            "oriented_gene_offset"
        ]
        <
        0,

        status[
            "oriented_gene_offset"
        ]
        >
        0,
    ],
    [
        1,

        (
            status[
                "n_oriented_upstream_genes_available"
            ]
            >=
            status[
                "oriented_gene_offset"
            ]
            .abs()
        ).astype(int),

        (
            status[
                "n_oriented_downstream_genes_available"
            ]
            >=
            status[
                "oriented_gene_offset"
            ]
            .abs()
        ).astype(int),
    ],
    default=0,
).astype(int)


## ================================================================== ##
## 4. Validate that every observable exact gene position has exactly
##    one observed neighborhood row.
## ================================================================== ##

if (
    position_rows[
        [
            "focal_key",
            "oriented_gene_offset",
        ]
    ]
    .duplicated()
    .any()
):

    fail(
        "Duplicate focal occurrence / exact-offset "
        "neighborhood rows detected."
    )


expected_keys = status.loc[
    status[
        "offset_observable"
    ]
    ==
    1,
    [
        "focal_key",
        "oriented_gene_offset",
    ],
].copy()


actual_keys = position_rows[
    [
        "focal_key",
        "oriented_gene_offset",
    ]
].copy()


key_check = expected_keys.merge(
    actual_keys,
    on=[
        "focal_key",
        "oriented_gene_offset",
    ],
    how="outer",
    indicator=True,
)


if (
    key_check[
        "_merge"
    ]
    !=
    "both"
).any():

    mismatch = (
        key_check[
            "_merge"
        ]
        !=
        "both"
    ).sum()

    fail(
        f"{mismatch:,} exact-position observability "
        f"keys disagree between Stage 12A and Stage 12B."
    )


## ================================================================== ##
## 5. Base exact-position denominators
## ================================================================== ##

POSITION_GROUP = [
    "focal_cluster",
    "focal_module",
    "oriented_gene_offset",
]


position_base = (
    status
    .groupby(
        POSITION_GROUP,
        as_index=False,
    )
    .agg(
        n_focal_occurrences_total=(
            "focal_key",
            "size",
        ),

        n_observable_occurrences=(
            "offset_observable",
            "sum",
        ),

        n_focal_genomes_total=(
            "genome",
            "nunique",
        ),
    )
)


genome_position_status = (
    status
    .groupby(
        POSITION_GROUP
        +
        [
            "genome",
        ],
        as_index=False,
    )
    .agg(
        n_focal_occurrences=(
            "focal_key",
            "size",
        ),

        n_observable_occurrences=(
            "offset_observable",
            "sum",
        ),
    )
)


genome_position_status[
    "all_focal_copies_observable"
] = (
    genome_position_status[
        "n_focal_occurrences"
    ]
    ==
    genome_position_status[
        "n_observable_occurrences"
    ]
).astype(int)


genome_position_summary = (
    genome_position_status
    .groupby(
        POSITION_GROUP,
        as_index=False,
    )
    .agg(
        n_genomes_all_focal_copies_observable=(
            "all_focal_copies_observable",
            "sum",
        )
    )
)


position_base = position_base.merge(
    genome_position_summary,
    on=POSITION_GROUP,
    how="left",
    validate="one_to_one",
)


## ================================================================== ##
## 6. Position-level annotation coverage
## ================================================================== ##

position_coverage = (
    position_rows
    .groupby(
        POSITION_GROUP,
        as_index=False,
    )
    .agg(
        n_observed_positions=(
            "focal_key",
            "size",
        ),

        n_positions_with_mmseqs_cluster=(
            "has_mmseqs_cluster",
            "sum",
        ),

        n_positions_with_cog=(
            "has_cog",
            "sum",
        ),

        n_positions_with_product=(
            "has_product",
            "sum",
        ),

        n_positions_without_mmseqs_or_cog=(
            "without_mmseqs_or_cog",
            "sum",
        ),

        n_positions_same_strand_as_focal=(
            "neighbor_same_strand_as_focal",
            "sum",
        ),
    )
)


positions = position_base.merge(
    position_coverage,
    on=POSITION_GROUP,
    how="left",
    validate="one_to_one",
)


coverage_count_columns = [
    "n_observed_positions",
    "n_positions_with_mmseqs_cluster",
    "n_positions_with_cog",
    "n_positions_with_product",
    "n_positions_without_mmseqs_or_cog",
    "n_positions_same_strand_as_focal",
]


for column in coverage_count_columns:

    positions[
        column
    ] = (
        positions[
            column
        ]
        .fillna(
            0
        )
        .astype(int)
    )


if (
    positions[
        "n_observed_positions"
    ]
    !=
    positions[
        "n_observable_occurrences"
    ]
).any():

    fail(
        "Observed exact-position row counts disagree "
        "with observability denominators."
    )


positions[
    "pct_observable_positions_with_mmseqs_cluster"
] = pct(
    positions[
        "n_positions_with_mmseqs_cluster"
    ],
    positions[
        "n_observable_occurrences"
    ],
)


positions[
    "pct_observable_positions_with_cog"
] = pct(
    positions[
        "n_positions_with_cog"
    ],
    positions[
        "n_observable_occurrences"
    ],
)


positions[
    "pct_observable_positions_with_product"
] = pct(
    positions[
        "n_positions_with_product"
    ],
    positions[
        "n_observable_occurrences"
    ],
)


positions[
    "pct_observable_positions_without_mmseqs_or_cog"
] = pct(
    positions[
        "n_positions_without_mmseqs_or_cog"
    ],
    positions[
        "n_observable_occurrences"
    ],
)


positions[
    "pct_observable_positions_same_strand_as_focal"
] = pct(
    positions[
        "n_positions_same_strand_as_focal"
    ],
    positions[
        "n_observable_occurrences"
    ],
)


## ================================================================== ##
## 7. Dominant MMseqs2 family at every exact position
## ================================================================== ##

print()
print("Summarizing dominant MMseqs2 families by exact position...")


top_mmseq = top_label_summary(
    position_rows,
    POSITION_GROUP,
    "neighbor_cluster",
    "mmseqs_cluster",
)


positions = positions.merge(
    top_mmseq,
    on=POSITION_GROUP,
    how="left",
    validate="one_to_one",
)


positions[
    "n_positions_with_dominant_mmseqs_cluster"
] = (
    positions[
        "n_positions_with_dominant_mmseqs_cluster"
    ]
    .fillna(
        0
    )
    .astype(int)
)


positions[
    "n_distinct_mmseqs_clusters"
] = (
    positions[
        "n_distinct_mmseqs_clusters"
    ]
    .fillna(
        0
    )
    .astype(int)
)


positions[
    "dominant_mmseqs_support_all_observable"
] = support_string(
    positions[
        "n_positions_with_dominant_mmseqs_cluster"
    ],
    positions[
        "n_observable_occurrences"
    ],
)


positions[
    "pct_all_observable_with_dominant_mmseqs_cluster"
] = pct(
    positions[
        "n_positions_with_dominant_mmseqs_cluster"
    ],
    positions[
        "n_observable_occurrences"
    ],
)


positions[
    "dominant_mmseqs_support_mmseqs_resolved"
] = support_string(
    positions[
        "n_positions_with_dominant_mmseqs_cluster"
    ],
    positions[
        "n_positions_with_mmseqs_cluster"
    ],
)


positions[
    "pct_mmseqs_resolved_with_dominant_mmseqs_cluster"
] = pct(
    positions[
        "n_positions_with_dominant_mmseqs_cluster"
    ],
    positions[
        "n_positions_with_mmseqs_cluster"
    ],
)


## ================================================================== ##
## 8. Dominant COG at every exact position
## ================================================================== ##

print("Summarizing dominant COGs by exact position...")


top_cog = top_label_summary(
    position_rows,
    POSITION_GROUP,
    "neighbor_cog",
    "cog",
)


positions = positions.merge(
    top_cog,
    on=POSITION_GROUP,
    how="left",
    validate="one_to_one",
)


positions[
    "n_positions_with_dominant_cog"
] = (
    positions[
        "n_positions_with_dominant_cog"
    ]
    .fillna(
        0
    )
    .astype(int)
)


positions[
    "n_distinct_cogs"
] = (
    positions[
        "n_distinct_cogs"
    ]
    .fillna(
        0
    )
    .astype(int)
)


positions[
    "dominant_cog_support_all_observable"
] = support_string(
    positions[
        "n_positions_with_dominant_cog"
    ],
    positions[
        "n_observable_occurrences"
    ],
)


positions[
    "pct_all_observable_with_dominant_cog"
] = pct(
    positions[
        "n_positions_with_dominant_cog"
    ],
    positions[
        "n_observable_occurrences"
    ],
)


positions[
    "dominant_cog_support_cog_resolved"
] = support_string(
    positions[
        "n_positions_with_dominant_cog"
    ],
    positions[
        "n_positions_with_cog"
    ],
)


positions[
    "pct_cog_resolved_with_dominant_cog"
] = pct(
    positions[
        "n_positions_with_dominant_cog"
    ],
    positions[
        "n_positions_with_cog"
    ],
)


## ================================================================== ##
## 9. Product label associated specifically with the dominant COG
## ================================================================== ##

print("Assigning product labels to dominant COGs...")


product_counts = (
    position_rows[
        (
            position_rows[
                "neighbor_cog"
            ]
            !=
            ""
        )
        &
        (
            position_rows[
                "neighbor_product"
            ]
            !=
            ""
        )
    ]
    .groupby(
        POSITION_GROUP
        +
        [
            "neighbor_cog",
            "neighbor_product",
        ],
        as_index=False,
    )
    .size()
    .rename(
        columns={
            "size":
                "product_count",
        }
    )
)


dominant_cog_keys = top_cog[
    POSITION_GROUP
    +
    [
        "dominant_cog",
    ]
].rename(
    columns={
        "dominant_cog":
            "neighbor_cog",
    }
)


product_counts = product_counts.merge(
    dominant_cog_keys,
    on=(
        POSITION_GROUP
        +
        [
            "neighbor_cog",
        ]
    ),
    how="inner",
    validate="many_to_one",
)


product_counts = product_counts.sort_values(
    POSITION_GROUP
    +
    [
        "product_count",
        "neighbor_product",
    ],
    ascending=(
        [True]
        *
        len(
            POSITION_GROUP
        )
        +
        [
            False,
            True,
        ]
    ),
    kind="stable",
)


dominant_products = (
    product_counts
    .drop_duplicates(
        POSITION_GROUP
    )
    [
        POSITION_GROUP
        +
        [
            "neighbor_product",
            "product_count",
        ]
    ]
    .rename(
        columns={
            "neighbor_product":
                "dominant_product_for_dominant_cog",

            "product_count":
                "n_positions_with_dominant_product_for_dominant_cog",
        }
    )
)


positions = positions.merge(
    dominant_products,
    on=POSITION_GROUP,
    how="left",
    validate="one_to_one",
)


## ================================================================== ##
## 10. Join genome-level support for the dominant exact MMseqs family
## ================================================================== ##

print()
print("Joining exact MMseqs genome-level support...")


mm_exact = pd.read_csv(
    MMSEQ_EXACT,
    sep="\t",
    dtype=str,
    keep_default_na=False,
)


require_columns(
    mm_exact,
    [
        "focal_cluster",
        "neighbor_cluster",
        "oriented_gene_offset",

        "n_focal_occurrences_with_neighbor_at_offset",
        "n_focal_occurrences_offset_observable",

        "n_informative_genomes",
        "n_genomes_with_neighbor_at_offset",
        "pct_informative_genomes_with_neighbor_at_offset",
    ],
    MMSEQ_EXACT.name,
)


numeric_mm_exact = [
    "oriented_gene_offset",
    "n_focal_occurrences_with_neighbor_at_offset",
    "n_focal_occurrences_offset_observable",
    "n_informative_genomes",
    "n_genomes_with_neighbor_at_offset",
    "pct_informative_genomes_with_neighbor_at_offset",
]


for column in numeric_mm_exact:

    mm_exact[
        column
    ] = pd.to_numeric(
        mm_exact[
            column
        ],
        errors="raise",
    )


mm_support = mm_exact[
    [
        "focal_cluster",
        "oriented_gene_offset",
        "neighbor_cluster",

        "n_informative_genomes",
        "n_genomes_with_neighbor_at_offset",
        "pct_informative_genomes_with_neighbor_at_offset",
    ]
].rename(
    columns={
        "neighbor_cluster":
            "dominant_mmseqs_cluster",

        "n_informative_genomes":
            "dominant_mmseqs_n_informative_genomes",

        "n_genomes_with_neighbor_at_offset":
            "dominant_mmseqs_n_genomes_with_feature",

        "pct_informative_genomes_with_neighbor_at_offset":
            "dominant_mmseqs_pct_informative_genomes_with_feature",
    }
)


positions = positions.merge(
    mm_support,
    on=[
        "focal_cluster",
        "oriented_gene_offset",
        "dominant_mmseqs_cluster",
    ],
    how="left",
    validate="one_to_one",
)


## Focal position 0 is the focal family by definition. ##

is_focal_position = (
    positions[
        "oriented_gene_offset"
    ]
    ==
    0
)


positions.loc[
    is_focal_position,
    "dominant_mmseqs_n_informative_genomes",
] = positions.loc[
    is_focal_position,
    "n_focal_genomes_total",
]


positions.loc[
    is_focal_position,
    "dominant_mmseqs_n_genomes_with_feature",
] = positions.loc[
    is_focal_position,
    "n_focal_genomes_total",
]


positions.loc[
    is_focal_position,
    "dominant_mmseqs_pct_informative_genomes_with_feature",
] = 100.0


## ================================================================== ##
## 11. Join genome-level support for dominant COG
## ================================================================== ##

print("Joining exact COG genome-level support...")


cog_exact = pd.read_csv(
    COG_EXACT,
    sep="\t",
    dtype=str,
    keep_default_na=False,
)


require_columns(
    cog_exact,
    [
        "focal_cluster",
        "neighbor_cog",
        "oriented_gene_offset",

        "n_informative_genomes",
        "n_genomes_with_this_cog",
        "pct_informative_genomes_with_this_cog",
    ],
    COG_EXACT.name,
)


for column in [
    "oriented_gene_offset",
    "n_informative_genomes",
    "n_genomes_with_this_cog",
    "pct_informative_genomes_with_this_cog",
]:

    cog_exact[
        column
    ] = pd.to_numeric(
        cog_exact[
            column
        ],
        errors="raise",
    )


cog_support = cog_exact[
    [
        "focal_cluster",
        "oriented_gene_offset",
        "neighbor_cog",

        "n_informative_genomes",
        "n_genomes_with_this_cog",
        "pct_informative_genomes_with_this_cog",
    ]
].rename(
    columns={
        "neighbor_cog":
            "dominant_cog",

        "n_informative_genomes":
            "dominant_cog_n_informative_genomes",

        "n_genomes_with_this_cog":
            "dominant_cog_n_genomes_with_feature",

        "pct_informative_genomes_with_this_cog":
            "dominant_cog_pct_informative_genomes_with_feature",
    }
)


positions = positions.merge(
    cog_support,
    on=[
        "focal_cluster",
        "oriented_gene_offset",
        "dominant_cog",
    ],
    how="left",
    validate="one_to_one",
)


## ================================================================== ##
## 12. Position direction / display fields
## ================================================================== ##

positions[
    "position_type"
] = np.where(
    positions[
        "oriented_gene_offset"
    ]
    ==
    0,
    "focal",
    "neighbor",
)


positions[
    "direction"
] = np.select(
    [
        positions[
            "oriented_gene_offset"
        ]
        <
        0,

        positions[
            "oriented_gene_offset"
        ]
        >
        0,
    ],
    [
        "upstream",
        "downstream",
    ],
    default="focal",
)


positions[
    "absolute_gene_offset"
] = (
    positions[
        "oriented_gene_offset"
    ]
    .abs()
)


## ================================================================== ##
## 13. MMseq2 directional consensus associations
## ================================================================== ##

print()
print("Building MMseqs2 directional consensus associations...")


mm_gene = pd.read_csv(
    MMSEQ_GENE_WINDOWS,
    sep="\t",
    dtype=str,
    keep_default_na=False,
)


require_columns(
    mm_gene,
    [
        "focal_cluster",
        "focal_module",
        "neighbor_cluster",
        "neighbor_module",

        "neighbor_same_mcl_module_as_focal",
        "neighbor_same_mmseqs_cluster_as_focal",

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
    ],
    MMSEQ_GENE_WINDOWS.name,
)


for column in [
    "neighbor_same_mcl_module_as_focal",
    "neighbor_same_mmseqs_cluster_as_focal",
    "gene_radius",

    "n_focal_occurrences_with_neighbor",
    "n_informative_focal_occurrences",
    "pct_informative_occurrences_with_neighbor",

    "n_genomes_with_neighbor",
    "n_informative_genomes",
    "pct_informative_genomes_with_neighbor",

    "n_censored_focal_occurrences",
    "n_censored_genomes",
]:

    mm_gene[
        column
    ] = pd.to_numeric(
        mm_gene[
            column
        ],
        errors="raise",
    )


MM_GROUP = [
    "focal_cluster",
    "focal_module",
    "neighbor_cluster",
    "neighbor_module",
    "neighbor_same_mcl_module_as_focal",
    "neighbor_same_mmseqs_cluster_as_focal",
    "direction",
]


## Occurrence-level saturation. ##

mm_occ_sat = first_radius_at_max(
    mm_gene,
    MM_GROUP,
    "n_focal_occurrences_with_neighbor",
    "gene_radius",
)


mm_occ_sat = mm_occ_sat[
    MM_GROUP
    +
    [
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
].rename(
    columns={
        "gene_radius":
            "occurrence_saturation_gene_radius",

        "n_focal_occurrences_with_neighbor":
            "occurrence_saturation_n_positive",

        "n_informative_focal_occurrences":
            "occurrence_saturation_n_informative",

        "occurrence_support":
            "occurrence_saturation_support",

        "pct_informative_occurrences_with_neighbor":
            "occurrence_saturation_pct",

        "n_genomes_with_neighbor":
            "n_positive_genomes_at_occurrence_saturation_radius",

        "n_informative_genomes":
            "n_informative_genomes_at_occurrence_saturation_radius",

        "genome_support":
            "genome_support_at_occurrence_saturation_radius",

        "pct_informative_genomes_with_neighbor":
            "genome_pct_at_occurrence_saturation_radius",

        "n_censored_focal_occurrences":
            "n_censored_occurrences_at_saturation",

        "n_censored_genomes":
            "n_censored_genomes_at_saturation",
    }
)


mm_assoc = mm_occ_sat.copy()


mm_assoc[
    "direction_has_gene_window_hits"
] = (
    mm_assoc[
        "occurrence_saturation_n_positive"
    ]
    >
    0
).astype(int)


## Zero-hit direction has no biologically meaningful saturation radius. ##

zero_mm = (
    mm_assoc[
        "direction_has_gene_window_hits"
    ]
    ==
    0
)


mm_assoc.loc[
    zero_mm,
    "occurrence_saturation_gene_radius",
] = np.nan


## Genome-level saturation radius separately. ##

mm_genome_sat = first_radius_at_max(
    mm_gene,
    MM_GROUP,
    "n_genomes_with_neighbor",
    "gene_radius",
)


mm_genome_sat = mm_genome_sat[
    MM_GROUP
    +
    [
        "gene_radius",
        "n_genomes_with_neighbor",
        "n_informative_genomes",
        "genome_support",
        "pct_informative_genomes_with_neighbor",
    ]
].rename(
    columns={
        "gene_radius":
            "genome_saturation_gene_radius",

        "n_genomes_with_neighbor":
            "genome_saturation_n_positive",

        "n_informative_genomes":
            "genome_saturation_n_informative",

        "genome_support":
            "genome_saturation_support",

        "pct_informative_genomes_with_neighbor":
            "genome_saturation_pct",
    }
)


mm_assoc = mm_assoc.merge(
    mm_genome_sat,
    on=MM_GROUP,
    how="left",
    validate="one_to_one",
)


mm_assoc.loc[
    mm_assoc[
        "genome_saturation_n_positive"
    ]
    ==
    0,
    "genome_saturation_gene_radius",
] = np.nan


## ================================================================== ##
## 14. Best exact MMseqs offset
## ================================================================== ##

mm_exact_best = best_exact_offset(
    mm_exact,
    [
        "focal_cluster",
        "neighbor_cluster",
        "direction",
    ],
    "n_focal_occurrences_with_neighbor_at_offset",
)


mm_exact_best = mm_exact_best[
    [
        "focal_cluster",
        "neighbor_cluster",
        "direction",

        "oriented_gene_offset",

        "n_focal_occurrences_with_neighbor_at_offset",
        "n_focal_occurrences_offset_observable",
        "pct_observable_occurrences_with_neighbor_at_offset",

        "n_genomes_with_neighbor_at_offset",
        "n_informative_genomes",
        "pct_informative_genomes_with_neighbor_at_offset",
    ]
].rename(
    columns={
        "oriented_gene_offset":
            "best_exact_gene_offset",

        "n_focal_occurrences_with_neighbor_at_offset":
            "best_exact_n_positive_occurrences",

        "n_focal_occurrences_offset_observable":
            "best_exact_n_observable_occurrences",

        "pct_observable_occurrences_with_neighbor_at_offset":
            "best_exact_pct_occurrences",

        "n_genomes_with_neighbor_at_offset":
            "best_exact_n_positive_genomes",

        "n_informative_genomes":
            "best_exact_n_informative_genomes",

        "pct_informative_genomes_with_neighbor_at_offset":
            "best_exact_pct_genomes",
    }
)


mm_exact_best[
    "best_exact_occurrence_support"
] = support_string(
    mm_exact_best[
        "best_exact_n_positive_occurrences"
    ],
    mm_exact_best[
        "best_exact_n_observable_occurrences"
    ],
)


mm_exact_best[
    "best_exact_genome_support"
] = support_string(
    mm_exact_best[
        "best_exact_n_positive_genomes"
    ],
    mm_exact_best[
        "best_exact_n_informative_genomes"
    ],
)


mm_assoc = mm_assoc.merge(
    mm_exact_best,
    on=[
        "focal_cluster",
        "neighbor_cluster",
        "direction",
    ],
    how="left",
    validate="one_to_one",
)


## ================================================================== ##
## 15. Physical 5/10/20-kb support
## ================================================================== ##

print("Joining physical MMseqs2-window evidence...")


mm_bp = pd.read_csv(
    MMSEQ_BP_WINDOWS,
    sep="\t",
    dtype=str,
    keep_default_na=False,
)


require_columns(
    mm_bp,
    [
        "focal_cluster",
        "neighbor_cluster",
        "direction",
        "radius_bp",

        "n_focal_occurrences_with_neighbor",
        "n_informative_focal_occurrences",
        "pct_informative_occurrences_with_neighbor",

        "n_genomes_with_neighbor",
        "n_informative_genomes",
        "pct_informative_genomes_with_neighbor",
    ],
    MMSEQ_BP_WINDOWS.name,
)


for column in [
    "radius_bp",

    "n_focal_occurrences_with_neighbor",
    "n_informative_focal_occurrences",
    "pct_informative_occurrences_with_neighbor",

    "n_genomes_with_neighbor",
    "n_informative_genomes",
    "pct_informative_genomes_with_neighbor",
]:

    mm_bp[
        column
    ] = pd.to_numeric(
        mm_bp[
            column
        ],
        errors="raise",
    )


for radius in BP_RADII:

    label = (
        f"{radius // 1000}kb"
    )


    tmp = mm_bp[
        mm_bp[
            "radius_bp"
        ]
        ==
        radius
    ][
        [
            "focal_cluster",
            "neighbor_cluster",
            "direction",

            "n_focal_occurrences_with_neighbor",
            "n_informative_focal_occurrences",
            "pct_informative_occurrences_with_neighbor",

            "n_genomes_with_neighbor",
            "n_informative_genomes",
            "pct_informative_genomes_with_neighbor",
        ]
    ].copy()


    tmp[
        f"physical_{label}_occurrence_support"
    ] = support_string(
        tmp[
            "n_focal_occurrences_with_neighbor"
        ],
        tmp[
            "n_informative_focal_occurrences"
        ],
    )


    tmp[
        f"physical_{label}_genome_support"
    ] = support_string(
        tmp[
            "n_genomes_with_neighbor"
        ],
        tmp[
            "n_informative_genomes"
        ],
    )


    tmp = tmp.rename(
        columns={
            "n_focal_occurrences_with_neighbor":
                f"physical_{label}_n_positive_occurrences",

            "n_informative_focal_occurrences":
                f"physical_{label}_n_informative_occurrences",

            "pct_informative_occurrences_with_neighbor":
                f"physical_{label}_pct_occurrences",

            "n_genomes_with_neighbor":
                f"physical_{label}_n_positive_genomes",

            "n_informative_genomes":
                f"physical_{label}_n_informative_genomes",

            "pct_informative_genomes_with_neighbor":
                f"physical_{label}_pct_genomes",
        }
    )


    mm_assoc = mm_assoc.merge(
        tmp,
        on=[
            "focal_cluster",
            "neighbor_cluster",
            "direction",
        ],
        how="left",
        validate="one_to_one",
    )


## ================================================================== ##
## 16. COG directional consensus associations
## ================================================================== ##

print()
print("Building COG directional consensus associations...")


cog_windows = pd.read_csv(
    COG_WINDOWS,
    sep="\t",
    dtype=str,
    keep_default_na=False,
)


require_columns(
    cog_windows,
    [
        "focal_cluster",
        "focal_module",
        "neighbor_cog",

        "direction",
        "gene_radius",

        "n_focal_occurrences_with_cog",
        "n_informative_focal_occurrences",
        "occurrence_support",
        "pct_informative_occurrences_with_cog",

        "n_genomes_with_cog",
        "n_informative_genomes",
        "genome_support",
        "pct_informative_genomes_with_cog",

        "n_censored_focal_occurrences",
        "n_censored_genomes",

        "pct_gene_annotation_coverage_in_fully_observable_windows",

        "dominant_product",
    ],
    COG_WINDOWS.name,
)


for column in [
    "gene_radius",

    "n_focal_occurrences_with_cog",
    "n_informative_focal_occurrences",
    "pct_informative_occurrences_with_cog",

    "n_genomes_with_cog",
    "n_informative_genomes",
    "pct_informative_genomes_with_cog",

    "n_censored_focal_occurrences",
    "n_censored_genomes",

    "pct_gene_annotation_coverage_in_fully_observable_windows",
]:

    cog_windows[
        column
    ] = pd.to_numeric(
        cog_windows[
            column
        ],
        errors="coerce",
    )


COG_GROUP = [
    "focal_cluster",
    "focal_module",
    "neighbor_cog",
    "direction",
]


cog_occ_sat = first_radius_at_max(
    cog_windows,
    COG_GROUP,
    "n_focal_occurrences_with_cog",
    "gene_radius",
)


cog_assoc = cog_occ_sat[
    COG_GROUP
    +
    [
        "gene_radius",

        "n_focal_occurrences_with_cog",
        "n_informative_focal_occurrences",
        "occurrence_support",
        "pct_informative_occurrences_with_cog",

        "n_genomes_with_cog",
        "n_informative_genomes",
        "genome_support",
        "pct_informative_genomes_with_cog",

        "n_censored_focal_occurrences",
        "n_censored_genomes",

        "pct_gene_annotation_coverage_in_fully_observable_windows",

        "dominant_product",
    ]
].rename(
    columns={
        "gene_radius":
            "occurrence_saturation_gene_radius",

        "n_focal_occurrences_with_cog":
            "occurrence_saturation_n_positive",

        "n_informative_focal_occurrences":
            "occurrence_saturation_n_informative",

        "occurrence_support":
            "occurrence_saturation_support",

        "pct_informative_occurrences_with_cog":
            "occurrence_saturation_pct",

        "n_genomes_with_cog":
            "n_positive_genomes_at_occurrence_saturation_radius",

        "n_informative_genomes":
            "n_informative_genomes_at_occurrence_saturation_radius",

        "genome_support":
            "genome_support_at_occurrence_saturation_radius",

        "pct_informative_genomes_with_cog":
            "genome_pct_at_occurrence_saturation_radius",

        "n_censored_focal_occurrences":
            "n_censored_occurrences_at_saturation",

        "n_censored_genomes":
            "n_censored_genomes_at_saturation",

        "pct_gene_annotation_coverage_in_fully_observable_windows":
            "annotation_coverage_pct_at_occurrence_saturation_radius",

        "dominant_product":
            "dominant_product_at_occurrence_saturation_radius",
    }
)


cog_assoc[
    "direction_has_gene_window_hits"
] = (
    cog_assoc[
        "occurrence_saturation_n_positive"
    ]
    >
    0
).astype(int)


cog_assoc.loc[
    cog_assoc[
        "direction_has_gene_window_hits"
    ]
    ==
    0,
    "occurrence_saturation_gene_radius",
] = np.nan


## Genome saturation separately. ##

cog_genome_sat = first_radius_at_max(
    cog_windows,
    COG_GROUP,
    "n_genomes_with_cog",
    "gene_radius",
)


cog_genome_sat = cog_genome_sat[
    COG_GROUP
    +
    [
        "gene_radius",
        "n_genomes_with_cog",
        "n_informative_genomes",
        "genome_support",
        "pct_informative_genomes_with_cog",
    ]
].rename(
    columns={
        "gene_radius":
            "genome_saturation_gene_radius",

        "n_genomes_with_cog":
            "genome_saturation_n_positive",

        "n_informative_genomes":
            "genome_saturation_n_informative",

        "genome_support":
            "genome_saturation_support",

        "pct_informative_genomes_with_cog":
            "genome_saturation_pct",
    }
)


cog_assoc = cog_assoc.merge(
    cog_genome_sat,
    on=COG_GROUP,
    how="left",
    validate="one_to_one",
)


cog_assoc.loc[
    cog_assoc[
        "genome_saturation_n_positive"
    ]
    ==
    0,
    "genome_saturation_gene_radius",
] = np.nan


## ================================================================== ##
## 17. Best exact COG offset
## ================================================================== ##

cog_exact_best = best_exact_offset(
    cog_exact,
    [
        "focal_cluster",
        "neighbor_cog",
        "direction",
    ],
    "n_focal_occurrences_with_this_cog",
)


cog_exact_best = cog_exact_best[
    [
        "focal_cluster",
        "neighbor_cog",
        "direction",

        "oriented_gene_offset",

        "n_focal_occurrences_with_this_cog",
        "n_focal_occurrences_offset_observable",
        "pct_all_observable_with_this_cog",

        "n_genomes_with_this_cog",
        "n_informative_genomes",
        "pct_informative_genomes_with_this_cog",

        "dominant_product",
    ]
].rename(
    columns={
        "oriented_gene_offset":
            "best_exact_gene_offset",

        "n_focal_occurrences_with_this_cog":
            "best_exact_n_positive_occurrences",

        "n_focal_occurrences_offset_observable":
            "best_exact_n_observable_occurrences",

        "pct_all_observable_with_this_cog":
            "best_exact_pct_occurrences",

        "n_genomes_with_this_cog":
            "best_exact_n_positive_genomes",

        "n_informative_genomes":
            "best_exact_n_informative_genomes",

        "pct_informative_genomes_with_this_cog":
            "best_exact_pct_genomes",

        "dominant_product":
            "best_exact_dominant_product",
    }
)


cog_exact_best[
    "best_exact_occurrence_support"
] = support_string(
    cog_exact_best[
        "best_exact_n_positive_occurrences"
    ],
    cog_exact_best[
        "best_exact_n_observable_occurrences"
    ],
)


cog_exact_best[
    "best_exact_genome_support"
] = support_string(
    cog_exact_best[
        "best_exact_n_positive_genomes"
    ],
    cog_exact_best[
        "best_exact_n_informative_genomes"
    ],
)


cog_assoc = cog_assoc.merge(
    cog_exact_best,
    on=[
        "focal_cluster",
        "neighbor_cog",
        "direction",
    ],
    how="left",
    validate="one_to_one",
)


## ================================================================== ##
## 18. QC
## ================================================================== ##

print()
print("Running consensus QC...")


if len(
    positions
) != EXPECTED_POSITION_ROWS:

    fail(
        f"Expected "
        f"{EXPECTED_POSITION_ROWS:,} "
        f"position-summary rows but found "
        f"{len(positions):,}."
    )


position_counts = (
    positions
    .groupby(
        "focal_cluster"
    )
    .size()
)


if (
    position_counts
    !=
    41
).any():

    fail(
        "A focal cluster does not have exactly "
        "41 consensus-position rows."
    )


if (
    mm_assoc[
        "occurrence_saturation_n_positive"
    ]
    >
    mm_assoc[
        "occurrence_saturation_n_informative"
    ]
).any():

    fail(
        "MMseq occurrence numerator exceeds denominator."
    )


if (
    mm_assoc[
        "genome_saturation_n_positive"
    ]
    >
    mm_assoc[
        "genome_saturation_n_informative"
    ]
).any():

    fail(
        "MMseq genome numerator exceeds denominator."
    )


if (
    cog_assoc[
        "occurrence_saturation_n_positive"
    ]
    >
    cog_assoc[
        "occurrence_saturation_n_informative"
    ]
).any():

    fail(
        "COG occurrence numerator exceeds denominator."
    )


if (
    cog_assoc[
        "genome_saturation_n_positive"
    ]
    >
    cog_assoc[
        "genome_saturation_n_informative"
    ]
).any():

    fail(
        "COG genome numerator exceeds denominator."
    )


## ================================================================== ##
## 19. Deterministic sorting
## ================================================================== ##

positions = (
    positions
    .sort_values(
        [
            "focal_module",
            "focal_cluster",
            "oriented_gene_offset",
        ],
        kind="stable",
    )
    .reset_index(
        drop=True
    )
)


mm_assoc = (
    mm_assoc
    .sort_values(
        [
            "focal_module",
            "focal_cluster",
            "neighbor_cluster",
            "direction",
        ],
        kind="stable",
    )
    .reset_index(
        drop=True
    )
)


cog_assoc = (
    cog_assoc
    .sort_values(
        [
            "focal_module",
            "focal_cluster",
            "neighbor_cog",
            "direction",
        ],
        kind="stable",
    )
    .reset_index(
        drop=True
    )
)


## ================================================================== ##
## 20. Write outputs
## ================================================================== ##

positions.to_csv(
    OUT_POSITIONS,
    sep="\t",
    index=False,
    na_rep="",
    float_format="%.6f",
)


mm_assoc.to_csv(
    OUT_MMSEQ,
    sep="\t",
    index=False,
    na_rep="",
    float_format="%.6f",
)


cog_assoc.to_csv(
    OUT_COG,
    sep="\t",
    index=False,
    na_rep="",
    float_format="%.6f",
)


## ================================================================== ##
## 21. QC table
## ================================================================== ##

qc = pd.DataFrame(
    [
        [
            "focal_proteins",
            len(
                focals
            ),
        ],

        [
            "focal_clusters",
            focals[
                "focal_cluster"
            ]
            .nunique(),
        ],

        [
            "modules",
            focals[
                "focal_module"
            ]
            .nunique(),
        ],

        [
            "consensus_position_rows",
            len(
                positions
            ),
        ],

        [
            "expected_consensus_position_rows",
            EXPECTED_POSITION_ROWS,
        ],

        [
            "mmseq_directional_consensus_rows",
            len(
                mm_assoc
            ),
        ],

        [
            "mmseq_directional_rows_with_hits",
            int(
                mm_assoc[
                    "direction_has_gene_window_hits"
                ]
                .sum()
            ),
        ],

        [
            "cog_directional_consensus_rows",
            len(
                cog_assoc
            ),
        ],

        [
            "cog_directional_rows_with_hits",
            int(
                cog_assoc[
                    "direction_has_gene_window_hits"
                ]
                .sum()
            ),
        ],

        [
            "position_observability_qc_pass",
            1,
        ],

        [
            "position_row_count_qc_pass",
            1,
        ],

        [
            "mmseq_denominator_qc_pass",
            1,
        ],

        [
            "cog_denominator_qc_pass",
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
## 22. Terminal summary
## ================================================================== ##

print()
print("Consensus-neighborhood summary")


print(
    f"  Focal proteins:                    "
    f"{len(focals):,}"
)

print(
    f"  Focal clusters:                    "
    f"{focals['focal_cluster'].nunique():,}"
)

print(
    f"  Exact consensus-position rows:     "
    f"{len(positions):,}"
)

print(
    f"  MMseq directional associations:    "
    f"{len(mm_assoc):,}"
)

print(
    f"  MMseq directions with hits:        "
    f"{mm_assoc['direction_has_gene_window_hits'].sum():,}"
)

print(
    f"  COG directional associations:      "
    f"{len(cog_assoc):,}"
)

print(
    f"  COG directions with hits:          "
    f"{cog_assoc['direction_has_gene_window_hits'].sum():,}"
)


## ================================================================== ##
## 23. Spotlight known relationships
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


for focal_cluster, neighbor_cluster in spotlights:

    spot = mm_assoc[
        (
            mm_assoc[
                "focal_cluster"
            ]
            ==
            focal_cluster
        )
        &
        (
            mm_assoc[
                "neighbor_cluster"
            ]
            ==
            neighbor_cluster
        )
    ].copy()


    print()
    print(
        f"MMseq spotlight: "
        f"{focal_cluster} -> {neighbor_cluster}"
    )


    if len(
        spot
    ) == 0:

        print(
            "  Relationship not found."
        )

        continue


    display_columns = [
        "direction",
        "direction_has_gene_window_hits",

        "best_exact_gene_offset",
        "best_exact_occurrence_support",
        "best_exact_pct_occurrences",

        "occurrence_saturation_gene_radius",
        "occurrence_saturation_support",
        "occurrence_saturation_pct",

        "genome_saturation_gene_radius",
        "genome_saturation_support",
        "genome_saturation_pct",

        "physical_5kb_occurrence_support",
        "physical_5kb_pct_occurrences",

        "physical_10kb_occurrence_support",
        "physical_10kb_pct_occurrences",

        "physical_20kb_occurrence_support",
        "physical_20kb_pct_occurrences",
    ]


    print(
        spot[
            display_columns
        ]
        .to_string(
            index=False
        )
    )


print()
print("=" * 80)
print("SUCCESS")
print("=" * 80)


print(
    f"Exact-position consensus: "
    f"{OUT_POSITIONS}"
)

print(
    f"MMseq consensus:          "
    f"{OUT_MMSEQ}"
)

print(
    f"COG consensus:            "
    f"{OUT_COG}"
)

print(
    f"QC:                       "
    f"{OUT_QC}"
)
