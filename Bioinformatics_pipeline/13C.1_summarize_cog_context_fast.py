#!/usr/bin/env python3

from pathlib import Path
import sys

import numpy as np
import pandas as pd


## ================================================================== ##
## STAGE 13C1 - FAST COG ANNOTATION CONTEXT
##
## PRIMARY FRAMEWORK:
##     focal MMseqs2 protein-family cluster
##
## SECONDARY FUNCTIONAL LAYER:
##     GlobDB COG
##
## GlobDB product is retained as a human-readable interpretation field.
##
## This script produces the SAME conceptual analyses as the original
## Stage-13C implementation, but avoids repeated dataframe scans.
##
## Analyses:
##
##   1. exact oriented COG positions
##   2. cumulative directional COG gene windows
##
## Percentages always retain numerator + denominator.
##
## Missing COG annotation means "COG not observed", not proof of
## biological absence.
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

OUT_EXACT = (
    HERE
    / "focal_cog_gene_order_context.tsv"
)

OUT_WINDOWS = (
    HERE
    / "focal_cog_gene_window_context.tsv"
)

OUT_QC = (
    HERE
    / "cog_context_fast_qc.tsv"
)


## ================================================================== ##
## Parameters
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

EXPECTED_FOCAL_PROTEINS = 10_537
EXPECTED_FOCAL_CLUSTERS = 156
EXPECTED_MODULES = 35


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
            f"{source} missing column(s): "
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
        numerator
    )

    denominator = pd.to_numeric(
        denominator
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


def dominant_product_table(
    df,
    group_columns,
):

    nonempty = df[
        df[
            "neighbor_product"
        ]
        .astype(str)
        .str.strip()
        !=
        ""
    ].copy()


    output_columns = (
        group_columns
        +
        [
            "dominant_product",
            "n_hits_with_dominant_product",
            "n_distinct_nonempty_product_labels",
        ]
    )


    if len(
        nonempty
    ) == 0:

        return pd.DataFrame(
            columns=output_columns
        )


    counts = (
        nonempty
        .groupby(
            group_columns
            +
            [
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


    counts = counts.sort_values(
        group_columns
        +
        [
            "product_count",
            "neighbor_product",
        ],
        ascending=(
            [True] * len(group_columns)
            +
            [
                False,
                True,
            ]
        ),
        kind="stable",
    )


    dominant = (
        counts
        .drop_duplicates(
            group_columns
        )
        [
            group_columns
            +
            [
                "neighbor_product",
                "product_count",
            ]
        ]
        .rename(
            columns={
                "neighbor_product":
                    "dominant_product",

                "product_count":
                    "n_hits_with_dominant_product",
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
                    "n_distinct_nonempty_product_labels",
            }
        )
    )


    return dominant.merge(
        distinct,
        on=group_columns,
        how="left",
        validate="one_to_one",
    )


## ================================================================== ##
## Start
## ================================================================== ##

print("=" * 80)
print("STAGE 13C1 - FAST COG ANNOTATION CONTEXT")
print("=" * 80)


## ================================================================== ##
## 1. Focal occurrences
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
        f"focal proteins; found "
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
        "Unexpected MCL-module count."
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
## 2. Neighborhood context
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

        "oriented_gene_offset",
        "is_focal",

        "neighbor_annotation_accepted",
        "neighbor_globdb_cog",
        "neighbor_globdb_product",

        "neighbor_cluster",
    ],
    NEIGHBORHOODS.name,
)


for column in [
    "oriented_gene_offset",
    "is_focal",
    "neighbor_annotation_accepted",
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


## GlobDB evidence must come only from accepted annotation matches. ##

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
        f"GlobDB annotation rows occur on "
        f"non-accepted mappings."
    )


context = neigh[
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
        MAX_GENE_OFFSET
    )
].copy()


context[
    "direction"
] = np.where(
    context[
        "oriented_gene_offset"
    ]
    <
    0,
    "upstream",
    "downstream",
)


context[
    "absolute_gene_offset"
] = (
    context[
        "oriented_gene_offset"
    ]
    .abs()
)


context[
    "has_cog"
] = (
    context[
        "neighbor_cog"
    ]
    !=
    ""
).astype(int)


context[
    "has_product"
] = (
    context[
        "neighbor_product"
    ]
    !=
    ""
).astype(int)


context[
    "has_mmseqs_cluster"
] = (
    context[
        "neighbor_cluster"
    ]
    .astype(str)
    .str.strip()
    !=
    ""
).astype(int)


cog_context = context[
    context[
        "has_cog"
    ]
    ==
    1
].copy()


print(
    f"  Ordinary rows +/-20 genes: "
    f"{len(context):,}"
)

print(
    f"  Rows with COG:              "
    f"{context['has_cog'].sum():,}"
)

print(
    f"  Rows with product:          "
    f"{context['has_product'].sum():,}"
)

print(
    f"  Distinct COGs:              "
    f"{cog_context['neighbor_cog'].nunique():,}"
)


## ================================================================== ##
## 3. Precompute exact-offset observability ONCE
## ================================================================== ##

print()
print(
    "Precomputing exact-offset denominators..."
)


exact_status_parts = []


for direction in [
    "upstream",
    "downstream",
]:

    if direction == "upstream":

        availability_column = (
            "n_oriented_upstream_genes_available"
        )

        sign = -1

    else:

        availability_column = (
            "n_oriented_downstream_genes_available"
        )

        sign = 1


    for distance in range(
        1,
        MAX_GENE_OFFSET + 1,
    ):

        tmp = focals[
            [
                "focal_key",
                "genome",
                "focal_cluster",
                "focal_module",
                availability_column,
            ]
        ].copy()


        tmp[
            "oriented_gene_offset"
        ] = (
            sign
            *
            distance
        )


        tmp[
            "offset_observable"
        ] = (
            tmp[
                availability_column
            ]
            >=
            distance
        ).astype(int)


        exact_status_parts.append(
            tmp[
                [
                    "focal_key",
                    "genome",
                    "focal_cluster",
                    "focal_module",
                    "oriented_gene_offset",
                    "offset_observable",
                ]
            ]
        )


exact_status = pd.concat(
    exact_status_parts,
    ignore_index=True,
)


EXACT_BASE = [
    "focal_cluster",
    "focal_module",
    "oriented_gene_offset",
]


exact_base_stats = (
    exact_status
    .groupby(
        EXACT_BASE,
        as_index=False,
    )
    .agg(
        n_focal_occurrences_total=(
            "focal_key",
            "size",
        ),

        n_focal_occurrences_offset_observable=(
            "offset_observable",
            "sum",
        ),

        n_focal_genomes_total=(
            "genome",
            "nunique",
        ),
    )
)


exact_genome_status = (
    exact_status
    .groupby(
        EXACT_BASE
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


exact_genome_status[
    "fully_observable_genome"
] = (
    exact_genome_status[
        "n_focal_occurrences"
    ]
    ==
    exact_genome_status[
        "n_observable_occurrences"
    ]
).astype(int)


exact_genome_full_counts = (
    exact_genome_status
    .groupby(
        EXACT_BASE,
        as_index=False,
    )
    .agg(
        n_genomes_all_focal_copies_offset_observable=(
            "fully_observable_genome",
            "sum",
        )
    )
)


## ================================================================== ##
## 4. Exact COG-position context - vectorized
## ================================================================== ##

print(
    "Calculating exact COG-position context..."
)


EXACT_TARGET = (
    EXACT_BASE
    +
    [
        "neighbor_cog",
    ]
)


exact_positive = (
    cog_context
    .groupby(
        EXACT_TARGET,
        as_index=False,
    )
    .agg(
        n_focal_occurrences_with_this_cog=(
            "focal_key",
            "nunique",
        )
    )
)


resolved_positions = (
    cog_context
    .groupby(
        [
            "focal_cluster",
            "focal_module",
            "oriented_gene_offset",
        ],
        as_index=False,
    )
    .agg(
        n_positions_with_any_cog_annotation=(
            "focal_key",
            "nunique",
        )
    )
)


exact_positive_genomes = (
    cog_context[
        EXACT_TARGET
        +
        [
            "genome",
        ]
    ]
    .drop_duplicates()
    .merge(
        exact_genome_status[
            EXACT_BASE
            +
            [
                "genome",
                "fully_observable_genome",
            ]
        ],
        on=(
            EXACT_BASE
            +
            [
                "genome",
            ]
        ),
        how="left",
        validate="many_to_one",
    )
)


exact_positive_genomes[
    "positive_genome_not_fully_observable"
] = (
    1
    -
    exact_positive_genomes[
        "fully_observable_genome"
    ]
)


exact_positive_genome_stats = (
    exact_positive_genomes
    .groupby(
        EXACT_TARGET,
        as_index=False,
    )
    .agg(
        n_genomes_with_this_cog=(
            "genome",
            "nunique",
        ),

        n_positive_genomes_not_fully_observable=(
            "positive_genome_not_fully_observable",
            "sum",
        ),
    )
)


exact_products = dominant_product_table(
    cog_context,
    EXACT_TARGET,
)


exact = (
    exact_positive
    .merge(
        exact_base_stats,
        on=EXACT_BASE,
        how="left",
        validate="many_to_one",
    )
    .merge(
        exact_genome_full_counts,
        on=EXACT_BASE,
        how="left",
        validate="many_to_one",
    )
    .merge(
        resolved_positions,
        on=EXACT_BASE,
        how="left",
        validate="many_to_one",
    )
    .merge(
        exact_positive_genome_stats,
        on=EXACT_TARGET,
        how="left",
        validate="one_to_one",
    )
    .merge(
        exact_products,
        on=EXACT_TARGET,
        how="left",
        validate="one_to_one",
    )
)


exact[
    "direction"
] = np.where(
    exact[
        "oriented_gene_offset"
    ]
    <
    0,
    "upstream",
    "downstream",
)


exact[
    "absolute_gene_offset"
] = (
    exact[
        "oriented_gene_offset"
    ]
    .abs()
)


exact[
    "n_informative_genomes"
] = (
    exact[
        "n_genomes_all_focal_copies_offset_observable"
    ]
    +
    exact[
        "n_positive_genomes_not_fully_observable"
    ]
)


exact[
    "support_all_observable"
] = (
    exact[
        "n_focal_occurrences_with_this_cog"
    ]
    .astype(str)
    +
    "/"
    +
    exact[
        "n_focal_occurrences_offset_observable"
    ]
    .astype(str)
)


exact[
    "pct_all_observable_with_this_cog"
] = pct(
    exact[
        "n_focal_occurrences_with_this_cog"
    ],
    exact[
        "n_focal_occurrences_offset_observable"
    ],
)


exact[
    "support_cog_resolved"
] = (
    exact[
        "n_focal_occurrences_with_this_cog"
    ]
    .astype(str)
    +
    "/"
    +
    exact[
        "n_positions_with_any_cog_annotation"
    ]
    .astype(str)
)


exact[
    "pct_cog_resolved_with_this_cog"
] = pct(
    exact[
        "n_focal_occurrences_with_this_cog"
    ],
    exact[
        "n_positions_with_any_cog_annotation"
    ],
)


exact[
    "pct_observable_positions_cog_resolved"
] = pct(
    exact[
        "n_positions_with_any_cog_annotation"
    ],
    exact[
        "n_focal_occurrences_offset_observable"
    ],
)


exact[
    "genome_support"
] = (
    exact[
        "n_genomes_with_this_cog"
    ]
    .astype(str)
    +
    "/"
    +
    exact[
        "n_informative_genomes"
    ]
    .astype(str)
)


exact[
    "pct_informative_genomes_with_this_cog"
] = pct(
    exact[
        "n_genomes_with_this_cog"
    ],
    exact[
        "n_informative_genomes"
    ],
)


if (
    exact[
        "n_focal_occurrences_with_this_cog"
    ]
    >
    exact[
        "n_focal_occurrences_offset_observable"
    ]
).any():

    fail(
        "Exact COG positives exceed observable denominator."
    )


## ================================================================== ##
## 5. Precompute cumulative-window observability ONCE
## ================================================================== ##

print()
print(
    "Precomputing cumulative-window denominators..."
)


window_status_parts = []


for direction in [
    "upstream",
    "downstream",
]:

    if direction == "upstream":

        availability_column = (
            "n_oriented_upstream_genes_available"
        )

    else:

        availability_column = (
            "n_oriented_downstream_genes_available"
        )


    for radius in GENE_RADII:

        tmp = focals[
            [
                "focal_key",
                "genome",
                "focal_cluster",
                "focal_module",
                availability_column,
            ]
        ].copy()


        tmp[
            "direction"
        ] = direction


        tmp[
            "gene_radius"
        ] = radius


        tmp[
            "full_window_observable"
        ] = (
            tmp[
                availability_column
            ]
            >=
            radius
        ).astype(int)


        window_status_parts.append(
            tmp[
                [
                    "focal_key",
                    "genome",
                    "focal_cluster",
                    "focal_module",
                    "direction",
                    "gene_radius",
                    "full_window_observable",
                ]
            ]
        )


window_status = pd.concat(
    window_status_parts,
    ignore_index=True,
)


WINDOW_BASE = [
    "focal_cluster",
    "focal_module",
    "direction",
    "gene_radius",
]


window_base_stats = (
    window_status
    .groupby(
        WINDOW_BASE,
        as_index=False,
    )
    .agg(
        n_focal_occurrences_total=(
            "focal_key",
            "size",
        ),

        n_focal_occurrences_full_window_observable=(
            "full_window_observable",
            "sum",
        ),

        n_focal_genomes_total=(
            "genome",
            "nunique",
        ),
    )
)


window_genome_status = (
    window_status
    .groupby(
        WINDOW_BASE
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

        n_full_window_occurrences=(
            "full_window_observable",
            "sum",
        ),
    )
)


window_genome_status[
    "fully_observable_genome"
] = (
    window_genome_status[
        "n_focal_occurrences"
    ]
    ==
    window_genome_status[
        "n_full_window_occurrences"
    ]
).astype(int)


window_genome_full_counts = (
    window_genome_status
    .groupby(
        WINDOW_BASE,
        as_index=False,
    )
    .agg(
        n_genomes_all_focal_copies_full_window=(
            "fully_observable_genome",
            "sum",
        )
    )
)


## ================================================================== ##
## 6. Expand actual COG hits across cumulative radii
##
## Each observed hit at offset 2 contributes to radii:
##
##     2, 3, 5, 10, 20
##
## but not radius 1.
##
## This replaces thousands of repeated dataframe filters.
## ================================================================== ##

print(
    "Expanding COG hits across cumulative windows..."
)


expanded_cog_parts = []


for radius in GENE_RADII:

    tmp = cog_context[
        cog_context[
            "absolute_gene_offset"
        ]
        <=
        radius
    ][
        [
            "focal_key",
            "genome",
            "focal_cluster",
            "focal_module",
            "neighbor_cog",
            "neighbor_product",
            "direction",
        ]
    ].copy()


    tmp[
        "gene_radius"
    ] = radius


    expanded_cog_parts.append(
        tmp
    )


expanded_cog = pd.concat(
    expanded_cog_parts,
    ignore_index=True,
)


WINDOW_TARGET = (
    WINDOW_BASE
    +
    [
        "neighbor_cog",
    ]
)


## ================================================================== ##
## 7. Positive focal occurrences
## ================================================================== ##

positive_occurrences = (
    expanded_cog[
        WINDOW_TARGET
        +
        [
            "focal_key",
            "genome",
        ]
    ]
    .drop_duplicates()
    .merge(
        window_status[
            [
                "focal_key",
                "focal_cluster",
                "focal_module",
                "direction",
                "gene_radius",
                "full_window_observable",
            ]
        ],
        on=[
            "focal_key",
            "focal_cluster",
            "focal_module",
            "direction",
            "gene_radius",
        ],
        how="left",
        validate="many_to_one",
    )
)


positive_occurrences[
    "positive_from_partial_window"
] = (
    1
    -
    positive_occurrences[
        "full_window_observable"
    ]
)


positive_occurrence_stats = (
    positive_occurrences
    .groupby(
        WINDOW_TARGET,
        as_index=False,
    )
    .agg(
        n_focal_occurrences_with_cog=(
            "focal_key",
            "nunique",
        ),

        n_positive_occurrences_from_partial_windows=(
            "positive_from_partial_window",
            "sum",
        ),
    )
)


## ================================================================== ##
## 8. Positive genomes
## ================================================================== ##

positive_genomes = (
    positive_occurrences[
        WINDOW_TARGET
        +
        [
            "genome",
        ]
    ]
    .drop_duplicates()
    .merge(
        window_genome_status[
            WINDOW_BASE
            +
            [
                "genome",
                "fully_observable_genome",
            ]
        ],
        on=(
            WINDOW_BASE
            +
            [
                "genome",
            ]
        ),
        how="left",
        validate="many_to_one",
    )
)


positive_genomes[
    "positive_genome_not_fully_observable"
] = (
    1
    -
    positive_genomes[
        "fully_observable_genome"
    ]
)


positive_genome_stats = (
    positive_genomes
    .groupby(
        WINDOW_TARGET,
        as_index=False,
    )
    .agg(
        n_genomes_with_cog=(
            "genome",
            "nunique",
        ),

        n_positive_genomes_not_fully_observable=(
            "positive_genome_not_fully_observable",
            "sum",
        ),
    )
)


## ================================================================== ##
## 9. Annotation coverage in fully observable windows
## ================================================================== ##

print(
    "Calculating annotation coverage..."
)


expanded_context_parts = []


for radius in GENE_RADII:

    tmp = context[
        context[
            "absolute_gene_offset"
        ]
        <=
        radius
    ][
        [
            "focal_key",
            "focal_cluster",
            "focal_module",
            "direction",
            "has_cog",
        ]
    ].copy()


    tmp[
        "gene_radius"
    ] = radius


    expanded_context_parts.append(
        tmp
    )


expanded_context = pd.concat(
    expanded_context_parts,
    ignore_index=True,
)


expanded_context = expanded_context.merge(
    window_status[
        [
            "focal_key",
            "focal_cluster",
            "focal_module",
            "direction",
            "gene_radius",
            "full_window_observable",
        ]
    ],
    on=[
        "focal_key",
        "focal_cluster",
        "focal_module",
        "direction",
        "gene_radius",
    ],
    how="left",
    validate="many_to_one",
)


full_context = expanded_context[
    expanded_context[
        "full_window_observable"
    ]
    ==
    1
]


annotation_coverage = (
    full_context
    .groupby(
        WINDOW_BASE,
        as_index=False,
    )
    .agg(
        n_genes_in_fully_observable_windows=(
            "has_cog",
            "size",
        ),

        n_genes_with_cog_annotation_in_fully_observable_windows=(
            "has_cog",
            "sum",
        ),
    )
)


annotation_coverage[
    "pct_gene_annotation_coverage_in_fully_observable_windows"
] = pct(
    annotation_coverage[
        "n_genes_with_cog_annotation_in_fully_observable_windows"
    ],
    annotation_coverage[
        "n_genes_in_fully_observable_windows"
    ],
)


## ================================================================== ##
## 10. Build complete candidate grid
##
## A COG observed anywhere around a focal family is evaluated in:
##
##     both directions
##     x all six gene radii
##
## so zero-positive directions remain explicit.
## ================================================================== ##

candidate_cogs = (
    cog_context[
        [
            "focal_cluster",
            "focal_module",
            "neighbor_cog",
        ]
    ]
    .drop_duplicates()
)


direction_radius_grid = pd.DataFrame(
    [
        {
            "direction":
                direction,

            "gene_radius":
                radius,
        }

        for direction
        in [
            "upstream",
            "downstream",
        ]

        for radius
        in GENE_RADII
    ]
)


window_template = candidate_cogs.merge(
    direction_radius_grid,
    how="cross",
)


window_products = dominant_product_table(
    expanded_cog,
    WINDOW_TARGET,
)


windows = (
    window_template
    .merge(
        window_base_stats,
        on=WINDOW_BASE,
        how="left",
        validate="many_to_one",
    )
    .merge(
        window_genome_full_counts,
        on=WINDOW_BASE,
        how="left",
        validate="many_to_one",
    )
    .merge(
        positive_occurrence_stats,
        on=WINDOW_TARGET,
        how="left",
        validate="one_to_one",
    )
    .merge(
        positive_genome_stats,
        on=WINDOW_TARGET,
        how="left",
        validate="one_to_one",
    )
    .merge(
        annotation_coverage,
        on=WINDOW_BASE,
        how="left",
        validate="many_to_one",
    )
    .merge(
        window_products,
        on=WINDOW_TARGET,
        how="left",
        validate="one_to_one",
    )
)


zero_columns = [
    "n_focal_occurrences_with_cog",
    "n_positive_occurrences_from_partial_windows",
    "n_genomes_with_cog",
    "n_positive_genomes_not_fully_observable",
]


for column in zero_columns:

    windows[
        column
    ] = (
        windows[
            column
        ]
        .fillna(
            0
        )
        .astype(int)
    )


windows[
    "n_informative_focal_occurrences"
] = (
    windows[
        "n_focal_occurrences_full_window_observable"
    ]
    +
    windows[
        "n_positive_occurrences_from_partial_windows"
    ]
)


windows[
    "n_censored_focal_occurrences"
] = (
    windows[
        "n_focal_occurrences_total"
    ]
    -
    windows[
        "n_informative_focal_occurrences"
    ]
)


windows[
    "n_informative_genomes"
] = (
    windows[
        "n_genomes_all_focal_copies_full_window"
    ]
    +
    windows[
        "n_positive_genomes_not_fully_observable"
    ]
)


windows[
    "n_censored_genomes"
] = (
    windows[
        "n_focal_genomes_total"
    ]
    -
    windows[
        "n_informative_genomes"
    ]
)


windows[
    "occurrence_support"
] = (
    windows[
        "n_focal_occurrences_with_cog"
    ]
    .astype(str)
    +
    "/"
    +
    windows[
        "n_informative_focal_occurrences"
    ]
    .astype(str)
)


windows[
    "pct_informative_occurrences_with_cog"
] = pct(
    windows[
        "n_focal_occurrences_with_cog"
    ],
    windows[
        "n_informative_focal_occurrences"
    ],
)


windows[
    "genome_support"
] = (
    windows[
        "n_genomes_with_cog"
    ]
    .astype(str)
    +
    "/"
    +
    windows[
        "n_informative_genomes"
    ]
    .astype(str)
)


windows[
    "pct_informative_genomes_with_cog"
] = pct(
    windows[
        "n_genomes_with_cog"
    ],
    windows[
        "n_informative_genomes"
    ],
)


## ================================================================== ##
## 11. QC
## ================================================================== ##

if (
    windows[
        "n_focal_occurrences_with_cog"
    ]
    >
    windows[
        "n_informative_focal_occurrences"
    ]
).any():

    fail(
        "Window occurrence positives exceed denominator."
    )


if (
    windows[
        "n_genomes_with_cog"
    ]
    >
    windows[
        "n_informative_genomes"
    ]
).any():

    fail(
        "Window genome positives exceed denominator."
    )


print(
    "Checking cumulative monotonicity..."
)


for _, group in windows.groupby(
    [
        "focal_cluster",
        "neighbor_cog",
        "direction",
    ],
    sort=False,
):

    group = group.sort_values(
        "gene_radius"
    )


    if (
        np.diff(
            group[
                "n_focal_occurrences_with_cog"
            ]
            .to_numpy()
        )
        <
        0
    ).any():

        fail(
            "Cumulative COG occurrence count decreased "
            "with increasing radius."
        )


    if (
        np.diff(
            group[
                "n_genomes_with_cog"
            ]
            .to_numpy()
        )
        <
        0
    ).any():

        fail(
            "Cumulative COG genome count decreased "
            "with increasing radius."
        )


## ================================================================== ##
## 12. Sort outputs
## ================================================================== ##

exact = (
    exact
    .sort_values(
        [
            "focal_module",
            "focal_cluster",
            "oriented_gene_offset",
            "neighbor_cog",
        ],
        kind="stable",
    )
    .reset_index(
        drop=True
    )
)


windows = (
    windows
    .sort_values(
        [
            "focal_module",
            "focal_cluster",
            "neighbor_cog",
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
## 13. Write
## ================================================================== ##

exact.to_csv(
    OUT_EXACT,
    sep="\t",
    index=False,
    na_rep="",
    float_format="%.6f",
)


windows.to_csv(
    OUT_WINDOWS,
    sep="\t",
    index=False,
    na_rep="",
    float_format="%.6f",
)


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
            "ordinary_context_rows",
            len(
                context
            ),
        ],

        [
            "rows_with_cog",
            int(
                context[
                    "has_cog"
                ]
                .sum()
            ),
        ],

        [
            "rows_with_product",
            int(
                context[
                    "has_product"
                ]
                .sum()
            ),
        ],

        [
            "distinct_cogs",
            cog_context[
                "neighbor_cog"
            ]
            .nunique(),
        ],

        [
            "exact_cog_rows",
            len(
                exact
            ),
        ],

        [
            "cumulative_cog_window_rows",
            len(
                windows
            ),
        ],

        [
            "exact_denominator_qc_pass",
            1,
        ],

        [
            "window_denominator_qc_pass",
            1,
        ],

        [
            "cumulative_monotonicity_qc_pass",
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
    "Fast COG-context summary"
)

print(
    f"  Focal proteins:             "
    f"{len(focals):,}"
)

print(
    f"  Context rows:               "
    f"{len(context):,}"
)

print(
    f"  COG-annotated rows:         "
    f"{context['has_cog'].sum():,}"
)

print(
    f"  Product-annotated rows:     "
    f"{context['has_product'].sum():,}"
)

print(
    f"  Distinct neighboring COGs:  "
    f"{cog_context['neighbor_cog'].nunique():,}"
)

print(
    f"  Exact COG-position rows:    "
    f"{len(exact):,}"
)

print(
    f"  Cumulative COG-window rows: "
    f"{len(windows):,}"
)


print()
print("=" * 80)
print("SUCCESS")
print("=" * 80)

print(
    f"Exact COG context: "
    f"{OUT_EXACT}"
)

print(
    f"COG-window context: "
    f"{OUT_WINDOWS}"
)

print(
    f"QC: "
    f"{OUT_QC}"
)
