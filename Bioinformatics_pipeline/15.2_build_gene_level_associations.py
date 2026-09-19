#!/usr/bin/env python3

from pathlib import Path
from collections import Counter
import sys

import numpy as np
import pandas as pd


## ================================================================== ##
## STAGE 15B - GENE-LEVEL ASSOCIATIONS
##
## PURPOSE
## -------
##
## Analyse the fresh Stage-15A neighborhoods at ACTUAL-GENE resolution.
##
## PRIMARY UNIT:
##
##     actual focal gene -> actual neighboring gene
##
## SECONDARY UNIT:
##
##     MMseqs2 family -> MMseqs2 family
##
## used only to ask whether gene-level relationships recur across
## different genomes.
##
##
## IMPORTANT
## ---------
##
## This stage DOES NOT collapse individual genes back into clusters.
##
## Every actual gene observation is retained.
##
## Family-level summaries are reversible projections over the
## underlying gene-level observations.
##
##
## BIOLOGICALLY INTERESTING NEIGHBORS
## ----------------------------------
##
## A neighboring gene enters the association analysis if it is:
##
##   * FindMeHemes-positive
##   OR
##   * FeGenie-positive
##   OR
##   * assigned to an MMseqs2 family
##   OR
##   * assigned to an MCL module
##
## Ordinary genes remain preserved in the Stage-15A neighborhood
## dataset and region tables, but do not need to enter every
## association-summary calculation.
##
##
## PARALOGUES
## ----------
##
## Paralogue copies remain separate actual genes.
##
## Family recurrence uses:
##
##   * focal-gene support
##   * genome support
##
## so multiple copies do not silently inflate genome-level support.
##
##
## CENSORING
## ---------
##
## Gene-window recurrence uses focal-specific observed gene extent.
##
## A focal is informative for a directional window if:
##
##   1. that full gene radius is observable
##
## OR
##
##   2. the target neighbor family is already observed within the
##      partial window.
##
## This mirrors the Stage-13 logic:
##
## a positive observation remains informative even if the farther edge
## of the requested window is censored.
##
##
## No new clustering is performed here.
## ================================================================== ##


## ================================================================== ##
## Paths
## ================================================================== ##

HERE = Path(__file__).resolve().parent
WORKFLOW = HERE.parent.parent


STAGE15A = (
    WORKFLOW
    / "15_gene_level_analysis"
    / "15A_resolution_recovery"
)


FOCALS = (
    STAGE15A
    / "focal_gene_catalog.tsv"
)


NEIGHBORHOODS = (
    STAGE15A
    / "focal_gene_neighborhoods.tsv"
)


REGION_GENES = (
    STAGE15A
    / "module_region_genes.tsv"
)


REGIONS = (
    STAGE15A
    / "fixed_module_eligibility_regions.tsv"
)


## ================================================================== ##
## Outputs
## ================================================================== ##

OUT_ACTUAL = (
    HERE
    / "actual_gene_neighbor_associations.tsv"
)


OUT_FOCAL_SUMMARY = (
    HERE
    / "focal_gene_local_context_summary.tsv"
)


OUT_REGION_ARCH = (
    HERE
    / "genome_module_interesting_gene_architectures.tsv"
)


OUT_EXACT = (
    HERE
    / "family_exact_gene_offset_associations.tsv"
)


OUT_WINDOWS = (
    HERE
    / "family_gene_window_associations.tsv"
)


OUT_FAMILY_SUMMARY = (
    HERE
    / "family_gene_association_summary.tsv"
)


OUT_UNCLUSTERED = (
    HERE
    / "unclustered_neighbor_annotation_summary.tsv"
)


OUT_FEGENIE = (
    HERE
    / "fegenie_neighbor_association_summary.tsv"
)


OUT_QC = (
    HERE
    / "gene_level_association_qc.tsv"
)


## ================================================================== ##
## Locked Stage-15A expectations
## ================================================================== ##

EXPECTED_FOCAL_ROWS = 12_653
EXPECTED_UNIQUE_ACTUAL_FOCALS = 11_124

EXPECTED_FRESH_NEIGHBORHOOD_ROWS = 317_000

EXPECTED_REGION_GENE_ROWS = 205_204
EXPECTED_FIXED_REGIONS = 8_568

EXPECTED_MODULES = 35
EXPECTED_MODULE_MEMBER_FOCALS = 10_537


## ================================================================== ##
## Gene-window radii
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
## Utilities
## ================================================================== ##

def fail(message):

    print(
        f"\nERROR: {message}",
        file=sys.stderr,
    )

    sys.exit(1)


def require_file(path):

    if not path.exists():

        fail(
            f"Required file does not exist:\n{path}"
        )


def read_tsv(path):

    require_file(
        path
    )

    return pd.read_csv(
        path,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )


def write_tsv(
    df,
    path,
):

    df.to_csv(
        path,
        sep="\t",
        index=False,
        na_rep="",
        float_format="%.6f",
    )


def clean(value):

    if pd.isna(
        value
    ):

        return ""

    return str(
        value
    ).strip()


def numeric(series):

    return pd.to_numeric(
        series.replace(
            "",
            np.nan,
        ),
        errors="coerce",
    )


def flag(series):

    return (
        pd.to_numeric(
            series.replace(
                "",
                np.nan,
            ),
            errors="coerce",
        )
        ==
        1
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


def support(
    numerator,
    denominator,
):

    return (
        f"{int(numerator)}"
        f"/"
        f"{int(denominator)}"
    )


def distribution_string(values):

    values = [
        clean(
            value
        )

        for value
        in values

        if clean(
            value
        )
        !=
        ""
    ]


    if len(
        values
    ) == 0:

        return ""


    counts = Counter(
        values
    )


    return "; ".join(
        f"{value}:{count}"

        for value, count
        in sorted(
            counts.items(),
            key=lambda item: (
                -item[1],
                item[0],
            ),
        )
    )


## ================================================================== ##
## Read Stage 15A
## ================================================================== ##

print("=" * 80)
print("STAGE 15B - GENE-LEVEL ASSOCIATIONS")
print("=" * 80)


print()
print("Reading Stage-15A resolution-recovery data...")


focals = read_tsv(
    FOCALS
)


neigh = read_tsv(
    NEIGHBORHOODS
)


region_genes = read_tsv(
    REGION_GENES
)


regions = read_tsv(
    REGIONS
)


## ================================================================== ##
## Basic QC
## ================================================================== ##

if len(
    focals
) != EXPECTED_FOCAL_ROWS:

    fail(
        f"Expected {EXPECTED_FOCAL_ROWS:,} focal rows; "
        f"found {len(focals):,}."
    )


if (
    focals[
        [
            "genome",
            "protein_id",
        ]
    ]
    .drop_duplicates()
    .shape[0]
    !=
    EXPECTED_UNIQUE_ACTUAL_FOCALS
):

    fail(
        "Unexpected unique actual focal-gene count."
    )


if len(
    neigh
) != EXPECTED_FRESH_NEIGHBORHOOD_ROWS:

    fail(
        f"Expected {EXPECTED_FRESH_NEIGHBORHOOD_ROWS:,} "
        f"fresh neighborhood rows; found {len(neigh):,}."
    )


if len(
    region_genes
) != EXPECTED_REGION_GENE_ROWS:

    fail(
        "Unexpected module-region gene-row count."
    )


if len(
    regions
) != EXPECTED_FIXED_REGIONS:

    fail(
        "Unexpected fixed-region count."
    )


if focals[
    "module"
].nunique() != EXPECTED_MODULES:

    fail(
        "Unexpected module count."
    )


if int(
    pd.to_numeric(
        focals[
            "focal_module_member"
        ],
        errors="coerce",
    )
    .fillna(
        0
    )
    .sum()
) != EXPECTED_MODULE_MEMBER_FOCALS:

    fail(
        "Unexpected number of original module-member focal genes."
    )


print(
    f"  Focal rows:                "
    f"{len(focals):,}"
)

print(
    f"  Unique actual focals:      "
    f"{focals[['genome','protein_id']].drop_duplicates().shape[0]:,}"
)

print(
    f"  Fresh neighborhood rows:   "
    f"{len(neigh):,}"
)

print(
    f"  Fixed-region gene rows:    "
    f"{len(region_genes):,}"
)

print(
    f"  Fixed regions:             "
    f"{len(regions):,}"
)


## ================================================================== ##
## Convert required numeric fields
## ================================================================== ##

numeric_columns = [
    "focal_start",
    "focal_end",
    "focal_gene_rank",

    "neighbor_start",
    "neighbor_end",
    "neighbor_gene_rank",

    "genomic_gene_offset",
    "oriented_gene_offset",

    "genomic_midpoint_offset_bp",
    "oriented_midpoint_offset_bp",
]


for column in numeric_columns:

    neigh[
        column
    ] = pd.to_numeric(
        neigh[
            column
        ],
        errors="raise",
    )


## ================================================================== ##
## Focal observability in gene units
##
## Maximum number of genes actually observed on each transcription-
## oriented side of each actual focal gene.
## ================================================================== ##

print()
print("Calculating focal-specific directional gene observability...")


extent = neigh[
    [
        "focal_id",
        "oriented_gene_offset",
    ]
].copy()


extent[
    "upstream_extent"
] = np.where(
    extent[
        "oriented_gene_offset"
    ]
    <
    0,
    -
    extent[
        "oriented_gene_offset"
    ],
    0,
)


extent[
    "downstream_extent"
] = np.where(
    extent[
        "oriented_gene_offset"
    ]
    >
    0,
    extent[
        "oriented_gene_offset"
    ],
    0,
)


extent = (
    extent
    .groupby(
        "focal_id",
        as_index=False,
    )
    .agg(
        observable_upstream_genes=(
            "upstream_extent",
            "max",
        ),

        observable_downstream_genes=(
            "downstream_extent",
            "max",
        ),
    )
)


focal_meta = focals.merge(
    extent,
    on="focal_id",
    how="left",
    validate="one_to_one",
)


if focal_meta[
    [
        "observable_upstream_genes",
        "observable_downstream_genes",
    ]
].isna().any().any():

    fail(
        "At least one focal gene lacks gene-observability information."
    )


## ================================================================== ##
## Identify biologically interesting actual neighboring genes
## ================================================================== ##

print()
print("Selecting biologically interesting actual neighboring genes...")


obs = neigh[
    pd.to_numeric(
        neigh[
            "is_focal"
        ],
        errors="coerce",
    )
    !=
    1
].copy()


obs[
    "_neighbor_fmh"
] = flag(
    obs[
        "neighbor_findmehemes_positive"
    ]
).astype(int)


obs[
    "_neighbor_fe"
] = flag(
    obs[
        "neighbor_fegenie_positive"
    ]
).astype(int)


obs[
    "_neighbor_clustered"
] = (
    obs[
        "neighbor_mmseq_cluster"
    ]
    !=
    ""
).astype(int)


obs[
    "_neighbor_mcl"
] = (
    obs[
        "neighbor_mcl_module"
    ]
    !=
    ""
).astype(int)


obs[
    "neighbor_is_biologically_interesting"
] = (
    obs[
        [
            "_neighbor_fmh",
            "_neighbor_fe",
            "_neighbor_clustered",
            "_neighbor_mcl",
        ]
    ]
    .max(
        axis=1
    )
    ==
    1
).astype(int)


interesting = obs[
    obs[
        "neighbor_is_biologically_interesting"
    ]
    ==
    1
].copy()


print(
    f"  Non-self neighborhood rows: "
    f"{len(obs):,}"
)

print(
    f"  Interesting neighbor rows:   "
    f"{len(interesting):,}"
)

print(
    f"  Unique neighboring genes:    "
    f"{interesting[['genome','neighbor_protein_id']].drop_duplicates().shape[0]:,}"
)


## ================================================================== ##
## Derive actual gene-to-gene geometry
## ================================================================== ##

print()
print("Deriving actual gene-pair geometry...")


interesting[
    "direction_from_focal"
] = np.where(
    interesting[
        "oriented_gene_offset"
    ]
    <
    0,
    "upstream",
    "downstream",
)


interesting[
    "genes_between"
] = (
    interesting[
        "genomic_gene_offset"
    ]
    .abs()
    -
    1
).clip(
    lower=0
).astype(int)


## Intergenic distance in physical genomic coordinates. ##

interesting[
    "intergenic_gap_bp"
] = np.select(
    [
        interesting[
            "neighbor_end"
        ]
        <
        interesting[
            "focal_start"
        ],

        interesting[
            "neighbor_start"
        ]
        >
        interesting[
            "focal_end"
        ],
    ],
    [
        (
            interesting[
                "focal_start"
            ]
            -
            interesting[
                "neighbor_end"
            ]
            -
            1
        ),

        (
            interesting[
                "neighbor_start"
            ]
            -
            interesting[
                "focal_end"
            ]
            -
            1
        ),
    ],
    default=0,
)


interesting[
    "intergenic_gap_bp"
] = interesting[
    "intergenic_gap_bp"
].clip(
    lower=0
).astype(int)


interesting[
    "is_adjacent"
] = (
    interesting[
        "genes_between"
    ]
    ==
    0
).astype(int)


interesting[
    "within_5kb_gap"
] = (
    interesting[
        "intergenic_gap_bp"
    ]
    <=
    5_000
).astype(int)


interesting[
    "within_10kb_gap"
] = (
    interesting[
        "intergenic_gap_bp"
    ]
    <=
    10_000
).astype(int)


interesting[
    "within_20kb_gap"
] = (
    interesting[
        "intergenic_gap_bp"
    ]
    <=
    20_000
).astype(int)


## ------------------------------------------------------------------ ##
## Genomic left/right identity
## ------------------------------------------------------------------ ##

focal_is_left = (
    interesting[
        "focal_gene_rank"
    ]
    <
    interesting[
        "neighbor_gene_rank"
    ]
)


interesting[
    "left_protein_id"
] = np.where(
    focal_is_left,
    interesting[
        "focal_protein_id"
    ],
    interesting[
        "neighbor_protein_id"
    ],
)


interesting[
    "right_protein_id"
] = np.where(
    focal_is_left,
    interesting[
        "neighbor_protein_id"
    ],
    interesting[
        "focal_protein_id"
    ],
)


interesting[
    "left_strand"
] = np.where(
    focal_is_left,
    interesting[
        "focal_strand"
    ],
    interesting[
        "neighbor_strand"
    ],
)


interesting[
    "right_strand"
] = np.where(
    focal_is_left,
    interesting[
        "neighbor_strand"
    ],
    interesting[
        "focal_strand"
    ],
)


interesting[
    "gene_pair_orientation"
] = np.select(
    [
        (
            interesting[
                "left_strand"
            ]
            ==
            "+"
        )
        &
        (
            interesting[
                "right_strand"
            ]
            ==
            "+"
        ),

        (
            interesting[
                "left_strand"
            ]
            ==
            "-"
        )
        &
        (
            interesting[
                "right_strand"
            ]
            ==
            "-"
        ),

        (
            interesting[
                "left_strand"
            ]
            ==
            "+"
        )
        &
        (
            interesting[
                "right_strand"
            ]
            ==
            "-"
        ),

        (
            interesting[
                "left_strand"
            ]
            ==
            "-"
        )
        &
        (
            interesting[
                "right_strand"
            ]
            ==
            "+"
        ),
    ],
    [
        "codirectional_right",
        "codirectional_left",
        "convergent",
        "divergent",
    ],
    default="",
)


interesting[
    "broad_orientation"
] = np.where(
    interesting[
        "gene_pair_orientation"
    ].isin(
        [
            "codirectional_right",
            "codirectional_left",
        ]
    ),
    "codirectional",
    interesting[
        "gene_pair_orientation"
    ],
)


## ------------------------------------------------------------------ ##
## Transcription-normalized ACTUAL gene order
## ------------------------------------------------------------------ ##

interesting[
    "transcriptional_gene_order"
] = np.select(
    [
        interesting[
            "gene_pair_orientation"
        ]
        ==
        "codirectional_right",

        interesting[
            "gene_pair_orientation"
        ]
        ==
        "codirectional_left",
    ],
    [
        (
            interesting[
                "left_protein_id"
            ]
            +
            "->"
            +
            interesting[
                "right_protein_id"
            ]
        ),

        (
            interesting[
                "right_protein_id"
            ]
            +
            "->"
            +
            interesting[
                "left_protein_id"
            ]
        ),
    ],
    default="",
)


## ------------------------------------------------------------------ ##
## Transcription-normalized MMseq family order, when both exist
## ------------------------------------------------------------------ ##

focal_cluster = interesting[
    "focal_mmseq_cluster"
]


neighbor_cluster = interesting[
    "neighbor_mmseq_cluster"
]


left_cluster = np.where(
    focal_is_left,
    focal_cluster,
    neighbor_cluster,
)


right_cluster = np.where(
    focal_is_left,
    neighbor_cluster,
    focal_cluster,
)


both_clustered = (
    focal_cluster
    !=
    ""
) & (
    neighbor_cluster
    !=
    ""
)


interesting[
    "transcriptional_family_order"
] = np.where(
    both_clustered
    &
    (
        interesting[
            "gene_pair_orientation"
        ]
        ==
        "codirectional_right"
    ),
    left_cluster
    +
    "->"
    +
    right_cluster,
    np.where(
        both_clustered
        &
        (
            interesting[
                "gene_pair_orientation"
            ]
            ==
            "codirectional_left"
        ),
        right_cluster
        +
        "->"
        +
        left_cluster,
        "",
    ),
)


## ================================================================== ##
## Relationship to focal context module
## ================================================================== ##

interesting[
    "neighbor_module_relationship"
] = np.select(
    [
        interesting[
            "neighbor_mcl_module"
        ]
        ==
        interesting[
            "focal_context_module"
        ],

        (
            interesting[
                "neighbor_mcl_module"
            ]
            !=
            ""
        )
        &
        (
            interesting[
                "neighbor_mcl_module"
            ]
            !=
            interesting[
                "focal_context_module"
            ]
        ),

        (
            interesting[
                "neighbor_mmseq_cluster"
            ]
            !=
            ""
        )
        &
        (
            interesting[
                "neighbor_mcl_module"
            ]
            ==
            ""
        ),
    ],
    [
        "same_context_module",
        "different_mcl_module",
        "clustered_no_mcl_module",
    ],
    default="unclustered_interesting_gene",
)


## ================================================================== ##
## Reciprocal-selected-focal flag
##
## Means that the neighboring actual gene itself qualified as a
## selected focal gene for this same module context.
## ================================================================== ##

interesting[
    "neighbor_is_selected_focal"
] = pd.to_numeric(
    interesting[
        "neighbor_is_selected_focal_for_this_module"
    ],
    errors="coerce",
).fillna(
    0
).astype(int)


## ================================================================== ##
## Actual-gene association table
## ================================================================== ##

actual_columns = [
    "focal_id",
    "focal_region_id",

    "genome",
    "focal_context_module",
    "focal_contig",

    "focal_protein_id",
    "focal_gene_rank",
    "focal_strand",

    "focal_selection_reasons",

    "focal_mmseq_cluster",
    "focal_mcl_module",

    "focal_fegenie_HMMs",
    "focal_number_of_hemes",

    "focal_signalp_prediction",
    "focal_deeptmhmm_class",
    "focal_report_topology",

    "focal_globdb_cog",
    "focal_globdb_product",

    "neighbor_protein_id",
    "neighbor_gene_rank",
    "neighbor_strand",

    "neighbor_mmseq_cluster",
    "neighbor_mcl_module",

    "neighbor_fegenie_positive",
    "neighbor_fegenie_HMMs",
    "neighbor_fegenie_categories",

    "neighbor_findmehemes_positive",
    "neighbor_number_of_hemes",

    "neighbor_signalp_prediction",
    "neighbor_deeptmhmm_class",
    "neighbor_report_topology",

    "neighbor_globdb_cog",
    "neighbor_globdb_product",

    "neighbor_inside_original_fixed_region",
    "neighbor_is_selected_focal",

    "neighbor_module_relationship",

    "genomic_gene_offset",
    "oriented_gene_offset",
    "direction_from_focal",

    "genes_between",
    "intergenic_gap_bp",

    "is_adjacent",
    "within_5kb_gap",
    "within_10kb_gap",
    "within_20kb_gap",

    "gene_pair_orientation",
    "broad_orientation",

    "transcriptional_gene_order",
    "transcriptional_family_order",
]


actual_columns = [
    column

    for column
    in actual_columns

    if column
    in interesting.columns
]


actual = interesting[
    actual_columns
].copy()


actual = actual.sort_values(
    [
        "focal_context_module",
        "genome",
        "focal_protein_id",
        "oriented_gene_offset",
        "neighbor_protein_id",
    ],
    kind="stable",
).reset_index(
    drop=True
)


write_tsv(
    actual,
    OUT_ACTUAL,
)


print()
print("Actual-gene association table")


print(
    f"  Directed actual-gene observations: "
    f"{len(actual):,}"
)

print(
    f"  Unique focal genes represented:    "
    f"{actual['focal_id'].nunique():,}"
)

print(
    f"  Neighbor genes from same module:   "
    f"{int((actual['neighbor_module_relationship'] == 'same_context_module').sum()):,}"
)

print(
    f"  Neighbor genes from other modules: "
    f"{int((actual['neighbor_module_relationship'] == 'different_mcl_module').sum()):,}"
)


## ================================================================== ##
## Focal-gene local context summary
##
## One row per actual focal gene × context module.
## ================================================================== ##

print()
print("Summarizing local context around each actual focal gene...")


summary_rows = []


actual_groups = {
    focal_id:
        group

    for focal_id, group
    in actual.groupby(
        "focal_id",
        sort=False,
    )
}


for focal in focal_meta.itertuples(
    index=False
):

    group = actual_groups.get(
        focal.focal_id
    )


    if group is None:

        group = actual.iloc[
            0:0
        ]


    same_module = (
        group[
            "neighbor_module_relationship"
        ]
        ==
        "same_context_module"
    )


    other_module = (
        group[
            "neighbor_module_relationship"
        ]
        ==
        "different_mcl_module"
    )


    clustered_no_module = (
        group[
            "neighbor_module_relationship"
        ]
        ==
        "clustered_no_mcl_module"
    )


    unclustered = (
        group[
            "neighbor_module_relationship"
        ]
        ==
        "unclustered_interesting_gene"
    )


    neighbor_fe = flag(
        group[
            "neighbor_fegenie_positive"
        ]
    ) if len(group) else pd.Series(
        dtype=bool
    )


    neighbor_fmh = flag(
        group[
            "neighbor_findmehemes_positive"
        ]
    ) if len(group) else pd.Series(
        dtype=bool
    )


    if len(
        group
    ) > 0:

        nearest_index = (
            group[
                "oriented_gene_offset"
            ]
            .abs()
            .idxmin()
        )


        nearest = group.loc[
            nearest_index
        ]


        nearest_id = nearest[
            "neighbor_protein_id"
        ]


        nearest_cluster = nearest[
            "neighbor_mmseq_cluster"
        ]


        nearest_offset = nearest[
            "oriented_gene_offset"
        ]


        nearest_gap = nearest[
            "intergenic_gap_bp"
        ]


    else:

        nearest_id = ""
        nearest_cluster = ""
        nearest_offset = np.nan
        nearest_gap = np.nan


    summary_rows.append(
        {
            "focal_id":
                focal.focal_id,

            "genome":
                focal.genome,

            "module":
                focal.module,

            "protein_id":
                focal.protein_id,

            "focal_reasons":
                focal.focal_reasons,

            "mmseq_cluster":
                focal.mmseq_cluster,

            "mcl_module":
                focal.mcl_module,

            "fegenie_HMMs":
                focal.fegenie_HMMs,

            "number_of_hemes":
                focal.number_of_hemes,

            "report_topology":
                focal.report_topology,

            "observable_upstream_genes":
                focal.observable_upstream_genes,

            "observable_downstream_genes":
                focal.observable_downstream_genes,

            "n_interesting_neighbor_rows":
                len(
                    group
                ),

            "n_unique_interesting_neighbor_genes":
                group[
                    "neighbor_protein_id"
                ].nunique(),

            "n_neighbor_fegenie_positive":
                int(
                    neighbor_fe.sum()
                )
                if len(
                    group
                )
                else
                0,

            "n_neighbor_findmehemes_positive":
                int(
                    neighbor_fmh.sum()
                )
                if len(
                    group
                )
                else
                0,

            "n_same_context_module_neighbors":
                int(
                    same_module.sum()
                ),

            "n_different_mcl_module_neighbors":
                int(
                    other_module.sum()
                ),

            "n_clustered_no_mcl_module_neighbors":
                int(
                    clustered_no_module.sum()
                ),

            "n_unclustered_interesting_neighbors":
                int(
                    unclustered.sum()
                ),

            "n_distinct_neighbor_mmseq_families":
                group.loc[
                    group[
                        "neighbor_mmseq_cluster"
                    ]
                    !=
                    "",
                    "neighbor_mmseq_cluster",
                ]
                .nunique(),

            "nearest_interesting_neighbor":
                nearest_id,

            "nearest_neighbor_mmseq_cluster":
                nearest_cluster,

            "nearest_oriented_gene_offset":
                nearest_offset,

            "nearest_intergenic_gap_bp":
                nearest_gap,
        }
    )


focal_summary = pd.DataFrame(
    summary_rows
)


write_tsv(
    focal_summary,
    OUT_FOCAL_SUMMARY,
)


## ================================================================== ##
## Genome × module × fixed-region interesting-gene architecture
##
## This is an ordered GENE-LEVEL representation, useful for Stage 16.
##
## It does not replace the full region-gene table.
## ================================================================== ##

print()
print("Building fixed-region gene-level architecture summaries...")


rg = region_genes.copy()


rg[
    "_interesting"
] = (
    (
        pd.to_numeric(
            rg[
                "is_automatic_focal"
            ],
            errors="coerce",
        )
        ==
        1
    )
    |
    (
        rg[
            "mmseq_cluster"
        ]
        !=
        ""
    )
    |
    (
        rg[
            "mcl_module"
        ]
        !=
        ""
    )
).astype(int)


rg[
    "gene_rank"
] = pd.to_numeric(
    rg[
        "gene_rank"
    ],
    errors="raise",
)


def gene_label(row):

    parts = []


    cluster = clean(
        row[
            "mmseq_cluster"
        ]
    )


    hmm = clean(
        row[
            "fegenie_HMMs"
        ]
    )


    hemes = clean(
        row[
            "number_of_hemes"
        ]
    )


    cog = clean(
        row[
            "globdb_cog"
        ]
    )


    if cluster != "":

        parts.append(
            cluster
        )


    if hmm != "":

        parts.append(
            f"FeGenie={hmm}"
        )


    if hemes != "":

        parts.append(
            f"hemes={hemes}"
        )


    if (
        cog != ""
        and
        len(
            parts
        )
        <
        3
    ):

        parts.append(
            cog
        )


    if len(
        parts
    ) == 0:

        parts.append(
            "unannotated"
        )


    return (
        f"{row['protein_id']}"
        f"({row['strand']})"
        f"["
        +
        "|".join(
            parts
        )
        +
        "]"
    )


architecture_rows = []


for (
    region_id,
    genome,
    module,
    contig,
), group in rg.groupby(
    [
        "region_id",
        "genome",
        "module",
        "contig",
    ],
    sort=True,
):

    group = group.sort_values(
        "gene_rank"
    )


    interesting_group = group[
        group[
            "_interesting"
        ]
        ==
        1
    ]


    architecture_rows.append(
        {
            "region_id":
                region_id,

            "genome":
                genome,

            "module":
                module,

            "contig":
                contig,

            "region_gene_rank_start":
                int(
                    group[
                        "gene_rank"
                    ].min()
                ),

            "region_gene_rank_end":
                int(
                    group[
                        "gene_rank"
                    ].max()
                ),

            "n_all_genes":
                len(
                    group
                ),

            "n_interesting_genes":
                len(
                    interesting_group
                ),

            "n_module_genes":
                int(
                    pd.to_numeric(
                        group[
                            "focal_module_member"
                        ],
                        errors="coerce",
                    )
                    .fillna(
                        0
                    )
                    .sum()
                ),

            "n_fegenie_positive_genes":
                int(
                    flag(
                        group[
                            "fegenie_positive"
                        ]
                    ).sum()
                ),

            "n_findmehemes_positive_genes":
                int(
                    flag(
                        group[
                            "findmehemes_positive"
                        ]
                    ).sum()
                ),

            "interesting_protein_order":
                " -> ".join(
                    interesting_group[
                        "protein_id"
                    ]
                ),

            "interesting_gene_architecture":
                " -> ".join(
                    gene_label(
                        row
                    )

                    for _, row
                    in interesting_group.iterrows()
                ),
        }
    )


region_arch = pd.DataFrame(
    architecture_rows
)


write_tsv(
    region_arch,
    OUT_REGION_ARCH,
)


## ================================================================== ##
## MMseq family-level recurrence
##
## Secondary projection only.
##
## We require BOTH actual genes to have an MMseq family assignment.
## ================================================================== ##

print()
print("Projecting actual-gene observations onto MMseq family pairs...")


family_obs = interesting[
    (
        interesting[
            "focal_mmseq_cluster"
        ]
        !=
        ""
    )
    &
    (
        interesting[
            "neighbor_mmseq_cluster"
        ]
        !=
        ""
    )
].copy()


print(
    f"  Family-resolvable observations: "
    f"{len(family_obs):,}"
)

print(
    f"  Distinct focal families:        "
    f"{family_obs['focal_mmseq_cluster'].nunique():,}"
)

print(
    f"  Distinct neighbor families:     "
    f"{family_obs['neighbor_mmseq_cluster'].nunique():,}"
)


## ================================================================== ##
## Focal-family denominator cache
## ================================================================== ##

family_focal_meta = focal_meta[
    focal_meta[
        "mmseq_cluster"
    ]
    !=
    ""
].copy()


focal_group_cache = {
    key:
        group.copy()

    for key, group
    in family_focal_meta.groupby(
        [
            "module",
            "mmseq_cluster",
        ],
        sort=False,
    )
}


## ================================================================== ##
## Support helper
##
## Genome-level denominator is deliberately conservative:
##
## Positive genome:
##     >=1 focal copy has the target family.
##
## Negative-but-informative genome:
##     ALL focal copies of that family in that genome have the full
##     requested directional window observable.
##
## Thus one observable paralogue does not falsely make a genome
## informative if another paralogue is censored.
## ================================================================== ##

def support_counts(
    candidate_focals,
    positive_focal_ids,
    direction,
    radius,
):

    if direction == "upstream":

        extent_column = (
            "observable_upstream_genes"
        )


    else:

        extent_column = (
            "observable_downstream_genes"
        )


    positive_focal_ids = set(
        positive_focal_ids
    )


    observable_mask = (
        candidate_focals[
            extent_column
        ]
        >=
        radius
    )


    observable_ids = set(
        candidate_focals.loc[
            observable_mask,
            "focal_id",
        ]
    )


    informative_ids = (
        observable_ids
        |
        positive_focal_ids
    )


    n_positive_focals = len(
        positive_focal_ids
    )


    n_informative_focals = len(
        informative_ids
    )


    ## -------------------------------------------------------------- ##
    ## Conservative genome-level support
    ## -------------------------------------------------------------- ##

    n_positive_genomes = 0
    n_informative_genomes = 0


    for genome, group in candidate_focals.groupby(
        "genome",
        sort=False,
    ):

        ids = set(
            group[
                "focal_id"
            ]
        )


        positive = (
            len(
                ids
                &
                positive_focal_ids
            )
            >
            0
        )


        all_observable = (
            group[
                extent_column
            ]
            >=
            radius
        ).all()


        informative = (
            positive
            or
            all_observable
        )


        if positive:

            n_positive_genomes += 1


        if informative:

            n_informative_genomes += 1


    return {
        "n_positive_focals":
            n_positive_focals,

        "n_informative_focals":
            n_informative_focals,

        "focal_occurrence_support":
            support(
                n_positive_focals,
                n_informative_focals,
            )
            if n_informative_focals
            else
            "",

        "focal_occurrence_pct":
            pct(
                n_positive_focals,
                n_informative_focals,
            ),

        "n_positive_genomes":
            n_positive_genomes,

        "n_informative_genomes":
            n_informative_genomes,

        "genome_support":
            support(
                n_positive_genomes,
                n_informative_genomes,
            )
            if n_informative_genomes
            else
            "",

        "genome_pct":
            pct(
                n_positive_genomes,
                n_informative_genomes,
            ),
    }


## ================================================================== ##
## Exact gene-offset associations
## ================================================================== ##

print()
print("Calculating exact gene-offset recurrence...")


exact_rows = []


exact_group_columns = [
    "focal_context_module",
    "focal_mmseq_cluster",
    "neighbor_mmseq_cluster",
    "direction_from_focal",
    "oriented_gene_offset",
]


for key, group in family_obs.groupby(
    exact_group_columns,
    sort=True,
):

    (
        context_module,
        focal_family,
        neighbor_family,
        direction,
        offset,
    ) = key


    radius = abs(
        int(
            offset
        )
    )


    candidate = focal_group_cache.get(
        (
            context_module,
            focal_family,
        )
    )


    if candidate is None:

        fail(
            f"Missing focal denominator for "
            f"{context_module} / {focal_family}"
        )


    positive_focal_ids = set(
        group[
            "focal_id"
        ]
    )


    counts = support_counts(
        candidate,
        positive_focal_ids,
        direction,
        radius,
    )


    pair_counts = (
        group
        .groupby(
            "focal_id"
        )[
            "neighbor_protein_id"
        ]
        .nunique()
    )


    exact_rows.append(
        {
            "context_module":
                context_module,

            "focal_cluster":
                focal_family,

            "neighbor_cluster":
                neighbor_family,

            "direction":
                direction,

            "oriented_gene_offset":
                int(
                    offset
                ),

            "absolute_gene_offset":
                radius,

            "n_actual_gene_pair_rows":
                len(
                    group
                ),

            "n_distinct_neighbor_genes":
                group[
                    [
                        "genome",
                        "neighbor_protein_id",
                    ]
                ]
                .drop_duplicates()
                .shape[0],

            "n_focals_with_multiple_neighbor_copies":
                int(
                    (
                        pair_counts
                        >
                        1
                    )
                    .sum()
                ),

            **counts,
        }
    )


exact = pd.DataFrame(
    exact_rows
)


write_tsv(
    exact,
    OUT_EXACT,
)


## ================================================================== ##
## Cumulative directional gene-window recurrence
## ================================================================== ##

print("Calculating cumulative directional gene-window recurrence...")


window_rows = []


association_columns = [
    "focal_context_module",
    "focal_mmseq_cluster",
    "neighbor_mmseq_cluster",
    "direction_from_focal",
]


for key, group in family_obs.groupby(
    association_columns,
    sort=True,
):

    (
        context_module,
        focal_family,
        neighbor_family,
        direction,
    ) = key


    candidate = focal_group_cache[
        (
            context_module,
            focal_family,
        )
    ]


    for radius in GENE_RADII:

        within = group[
            group[
                "oriented_gene_offset"
            ].abs()
            <=
            radius
        ]


        positive_focal_ids = set(
            within[
                "focal_id"
            ]
        )


        counts = support_counts(
            candidate,
            positive_focal_ids,
            direction,
            radius,
        )


        copy_counts = (
            within
            .groupby(
                "focal_id"
            )[
                "neighbor_protein_id"
            ]
            .nunique()
        )


        window_rows.append(
            {
                "context_module":
                    context_module,

                "focal_cluster":
                    focal_family,

                "neighbor_cluster":
                    neighbor_family,

                "direction":
                    direction,

                "gene_radius":
                    radius,

                "n_actual_gene_pair_rows":
                    len(
                        within
                    ),

                "n_distinct_neighbor_genes":
                    within[
                        [
                            "genome",
                            "neighbor_protein_id",
                        ]
                    ]
                    .drop_duplicates()
                    .shape[0],

                "n_focals_with_multiple_neighbor_copies":
                    int(
                        (
                            copy_counts
                            >
                            1
                        )
                        .sum()
                    ),

                **counts,
            }
        )


windows = pd.DataFrame(
    window_rows
)


write_tsv(
    windows,
    OUT_WINDOWS,
)


## ================================================================== ##
## Compact family-association summary
## ================================================================== ##

print("Building compact family-association summary...")


family_summary_rows = []


for key, group in family_obs.groupby(
    association_columns,
    sort=True,
):

    (
        context_module,
        focal_family,
        neighbor_family,
        direction,
    ) = key


    exact_group = exact[
        (
            exact[
                "context_module"
            ]
            ==
            context_module
        )
        &
        (
            exact[
                "focal_cluster"
            ]
            ==
            focal_family
        )
        &
        (
            exact[
                "neighbor_cluster"
            ]
            ==
            neighbor_family
        )
        &
        (
            exact[
                "direction"
            ]
            ==
            direction
        )
    ].copy()


    window_group = windows[
        (
            windows[
                "context_module"
            ]
            ==
            context_module
        )
        &
        (
            windows[
                "focal_cluster"
            ]
            ==
            focal_family
        )
        &
        (
            windows[
                "neighbor_cluster"
            ]
            ==
            neighbor_family
        )
        &
        (
            windows[
                "direction"
            ]
            ==
            direction
        )
    ].copy()


    ## -------------------------------------------------------------- ##
    ## Best exact offset:
    ##
    ## priority:
    ##   1. largest focal-gene numerator
    ##   2. largest genome numerator
    ##   3. highest occurrence %
    ##   4. closest absolute offset
    ## -------------------------------------------------------------- ##

    exact_group[
        "_occ_pct"
    ] = pd.to_numeric(
        exact_group[
            "focal_occurrence_pct"
        ],
        errors="coerce",
    ).fillna(
        -1
    )


    best_exact = (
        exact_group
        .sort_values(
            [
                "n_positive_focals",
                "n_positive_genomes",
                "_occ_pct",
                "absolute_gene_offset",
                "oriented_gene_offset",
            ],
            ascending=[
                False,
                False,
                False,
                True,
                True,
            ],
        )
        .iloc[0]
    )


    ## -------------------------------------------------------------- ##
    ## Saturation radius:
    ##
    ## smallest tested radius where the maximum numerator is reached.
    ## -------------------------------------------------------------- ##

    max_focal_n = window_group[
        "n_positive_focals"
    ].max()


    occurrence_saturation = (
        window_group[
            window_group[
                "n_positive_focals"
            ]
            ==
            max_focal_n
        ]
        .sort_values(
            "gene_radius"
        )
        .iloc[0]
    )


    max_genome_n = window_group[
        "n_positive_genomes"
    ].max()


    genome_saturation = (
        window_group[
            window_group[
                "n_positive_genomes"
            ]
            ==
            max_genome_n
        ]
        .sort_values(
            "gene_radius"
        )
        .iloc[0]
    )


    candidate = focal_group_cache[
        (
            context_module,
            focal_family,
        )
    ]


    all_copy_counts = (
        group
        .groupby(
            "focal_id"
        )[
            "neighbor_protein_id"
        ]
        .nunique()
    )


    ## Raw physical proximity observations.
    ##
    ## These are NOT censor-adjusted physical-window denominators.
    ## They simply describe observed gene-level pairs.
    ## -------------------------------------------------------------- ##

    within5 = group[
        group[
            "intergenic_gap_bp"
        ]
        <=
        5_000
    ]


    within10 = group[
        group[
            "intergenic_gap_bp"
        ]
        <=
        10_000
    ]


    within20 = group[
        group[
            "intergenic_gap_bp"
        ]
        <=
        20_000
    ]


    family_summary_rows.append(
        {
            "context_module":
                context_module,

            "focal_cluster":
                focal_family,

            "neighbor_cluster":
                neighbor_family,

            "direction":
                direction,

            "same_mmseq_family":
                int(
                    focal_family
                    ==
                    neighbor_family
                ),

            "n_focal_occurrences_total":
                len(
                    candidate
                ),

            "n_focal_genomes_total":
                candidate[
                    "genome"
                ].nunique(),

            "n_actual_gene_pair_rows":
                len(
                    group
                ),

            "n_focal_genes_with_neighbor":
                group[
                    "focal_id"
                ].nunique(),

            "n_genomes_with_neighbor":
                group[
                    "genome"
                ].nunique(),

            "n_distinct_actual_neighbor_genes":
                group[
                    [
                        "genome",
                        "neighbor_protein_id",
                    ]
                ]
                .drop_duplicates()
                .shape[0],

            "n_focals_with_multiple_neighbor_copies":
                int(
                    (
                        all_copy_counts
                        >
                        1
                    )
                    .sum()
                ),

            "maximum_neighbor_copies_per_focal":
                int(
                    all_copy_counts.max()
                ),

            "best_exact_gene_offset":
                int(
                    best_exact[
                        "oriented_gene_offset"
                    ]
                ),

            "best_exact_occurrence_support":
                best_exact[
                    "focal_occurrence_support"
                ],

            "best_exact_genome_support":
                best_exact[
                    "genome_support"
                ],

            "occurrence_saturation_gene_radius":
                int(
                    occurrence_saturation[
                        "gene_radius"
                    ]
                ),

            "occurrence_saturation_support":
                occurrence_saturation[
                    "focal_occurrence_support"
                ],

            "genome_saturation_gene_radius":
                int(
                    genome_saturation[
                        "gene_radius"
                    ]
                ),

            "genome_saturation_support":
                genome_saturation[
                    "genome_support"
                ],

            "n_focal_genes_with_neighbor_within_5kb_observed":
                within5[
                    "focal_id"
                ].nunique(),

            "n_focal_genes_with_neighbor_within_10kb_observed":
                within10[
                    "focal_id"
                ].nunique(),

            "n_focal_genes_with_neighbor_within_20kb_observed":
                within20[
                    "focal_id"
                ].nunique(),

            "n_genomes_with_neighbor_within_5kb_observed":
                within5[
                    "genome"
                ].nunique(),

            "n_genomes_with_neighbor_within_10kb_observed":
                within10[
                    "genome"
                ].nunique(),

            "n_genomes_with_neighbor_within_20kb_observed":
                within20[
                    "genome"
                ].nunique(),

            "median_absolute_gene_offset":
                group[
                    "oriented_gene_offset"
                ]
                .abs()
                .median(),

            "median_intergenic_gap_bp":
                group[
                    "intergenic_gap_bp"
                ].median(),

            "dominant_pair_orientation":
                (
                    group[
                        "broad_orientation"
                    ]
                    .value_counts()
                    .index[0]
                    if len(
                        group
                    )
                    else
                    ""
                ),

            "orientation_distribution":
                distribution_string(
                    group[
                        "broad_orientation"
                    ]
                ),

            "transcriptional_family_order_distribution":
                distribution_string(
                    group[
                        "transcriptional_family_order"
                    ]
                ),
        }
    )


family_summary = pd.DataFrame(
    family_summary_rows
)


write_tsv(
    family_summary,
    OUT_FAMILY_SUMMARY,
)


## ================================================================== ##
## FeGenie association summary
##
## This explicitly surfaces the FeGenie evidence reintroduced in 15A.
## ================================================================== ##

print()
print("Summarizing FeGenie-positive neighboring genes...")


fe_obs = interesting[
    flag(
        interesting[
            "neighbor_fegenie_positive"
        ]
    )
].copy()


fe_rows = []


for key, group in fe_obs.groupby(
    [
        "focal_context_module",
        "focal_mmseq_cluster",
        "neighbor_fegenie_HMMs",
        "direction_from_focal",
    ],
    sort=True,
):

    (
        context_module,
        focal_family,
        neighbor_hmm,
        direction,
    ) = key


    fe_rows.append(
        {
            "context_module":
                context_module,

            "focal_cluster":
                focal_family,

            "neighbor_fegenie_HMM":
                neighbor_hmm,

            "direction":
                direction,

            "n_actual_gene_pair_rows":
                len(
                    group
                ),

            "n_focal_genes":
                group[
                    "focal_id"
                ].nunique(),

            "n_genomes":
                group[
                    "genome"
                ].nunique(),

            "n_neighbor_genes":
                group[
                    [
                        "genome",
                        "neighbor_protein_id",
                    ]
                ]
                .drop_duplicates()
                .shape[0],

            "neighbor_mmseq_cluster_distribution":
                distribution_string(
                    group[
                        "neighbor_mmseq_cluster"
                    ]
                ),

            "neighbor_mcl_module_distribution":
                distribution_string(
                    group[
                        "neighbor_mcl_module"
                    ]
                ),

            "median_absolute_gene_offset":
                group[
                    "oriented_gene_offset"
                ]
                .abs()
                .median(),

            "median_intergenic_gap_bp":
                group[
                    "intergenic_gap_bp"
                ].median(),
        }
    )


fegenie_summary = pd.DataFrame(
    fe_rows
)


write_tsv(
    fegenie_summary,
    OUT_FEGENIE,
)


## ================================================================== ##
## Unclustered interesting neighbors
##
## IMPORTANT:
##
## This is ANNOTATION-LEVEL recurrence, not a homology-family analysis.
##
## If an interesting neighboring gene has no MMseq family, we are not
## allowed to imply that two such genes are homologues merely because
## they share a COG, FeGenie HMM, heme count or topology.
## ================================================================== ##

print("Summarizing unclustered interesting genes at annotation level...")


unclustered = interesting[
    interesting[
        "neighbor_mmseq_cluster"
    ]
    ==
    ""
].copy()


def annotation_key(row):

    hmm = clean(
        row[
            "neighbor_fegenie_HMMs"
        ]
    )


    cog = clean(
        row[
            "neighbor_globdb_cog"
        ]
    )


    hemes = clean(
        row[
            "neighbor_number_of_hemes"
        ]
    )


    topology = clean(
        row[
            "neighbor_report_topology"
        ]
    )


    if hmm != "":

        return (
            "FeGenie:"
            +
            hmm
        )


    if cog != "":

        return (
            "COG:"
            +
            cog
        )


    return (
        "Unclustered_FMH:"
        f"hemes={hemes if hemes != '' else 'NA'};"
        f"topology={topology if topology != '' else 'NA'}"
    )


if len(
    unclustered
) > 0:

    unclustered[
        "annotation_key"
    ] = unclustered.apply(
        annotation_key,
        axis=1,
    )


    unclustered_summary = (
        unclustered
        .groupby(
            [
                "focal_context_module",
                "focal_mmseq_cluster",
                "direction_from_focal",
                "annotation_key",
            ],
            as_index=False,
        )
        .agg(
            n_actual_gene_rows=(
                "neighbor_protein_id",
                "size",
            ),

            n_actual_neighbor_genes=(
                "neighbor_protein_id",
                "nunique",
            ),

            n_focal_genes=(
                "focal_id",
                "nunique",
            ),

            n_genomes=(
                "genome",
                "nunique",
            ),

            median_absolute_gene_offset=(
                "oriented_gene_offset",
                lambda values:
                    values.abs().median(),
            ),

            median_intergenic_gap_bp=(
                "intergenic_gap_bp",
                "median",
            ),
        )
        .rename(
            columns={
                "focal_context_module":
                    "context_module",

                "focal_mmseq_cluster":
                    "focal_cluster",

                "direction_from_focal":
                    "direction",
            }
        )
    )


else:

    unclustered_summary = pd.DataFrame(
        columns=[
            "context_module",
            "focal_cluster",
            "direction",
            "annotation_key",
            "n_actual_gene_rows",
            "n_actual_neighbor_genes",
            "n_focal_genes",
            "n_genomes",
            "median_absolute_gene_offset",
            "median_intergenic_gap_bp",
        ]
    )


write_tsv(
    unclustered_summary,
    OUT_UNCLUSTERED,
)


## ================================================================== ##
## QC
## ================================================================== ##

print()
print("Running Stage-15B QC...")


## Every actual association must be non-self. ##

if (
    actual[
        "focal_protein_id"
    ]
    ==
    actual[
        "neighbor_protein_id"
    ]
).any():

    fail(
        "Self-pair detected in actual-gene association table."
    )


## Every association must be same contig. ##

if (
    interesting[
        "focal_contig"
    ]
    !=
    interesting[
        "neighbor_contig"
    ]
).any():

    fail(
        "Cross-contig actual-gene association detected."
    )


## Every exact-family row must have a nonzero offset. ##

if (
    exact[
        "oriented_gene_offset"
    ]
    ==
    0
).any():

    fail(
        "Zero-offset non-self family association detected."
    )


## Window numerator must never exceed denominator. ##

if (
    windows[
        "n_positive_focals"
    ]
    >
    windows[
        "n_informative_focals"
    ]
).any():

    fail(
        "Focal occurrence numerator exceeds denominator."
    )


if (
    windows[
        "n_positive_genomes"
    ]
    >
    windows[
        "n_informative_genomes"
    ]
).any():

    fail(
        "Genome numerator exceeds denominator."
    )


## Cumulative numerator must be monotonic with increasing radius. ##

monotonic_failures = 0


for _, group in windows.groupby(
    [
        "context_module",
        "focal_cluster",
        "neighbor_cluster",
        "direction",
    ],
    sort=False,
):

    group = group.sort_values(
        "gene_radius"
    )


    focal_values = group[
        "n_positive_focals"
    ].to_numpy()


    genome_values = group[
        "n_positive_genomes"
    ].to_numpy()


    if np.any(
        np.diff(
            focal_values
        )
        <
        0
    ):

        monotonic_failures += 1


    if np.any(
        np.diff(
            genome_values
        )
        <
        0
    ):

        monotonic_failures += 1


if monotonic_failures != 0:

    fail(
        f"Gene-window monotonicity failed "
        f"{monotonic_failures} times."
    )


qc = pd.DataFrame(
    [
        [
            "stage15A_focal_rows",
            len(
                focals
            ),
            EXPECTED_FOCAL_ROWS,
            int(
                len(
                    focals
                )
                ==
                EXPECTED_FOCAL_ROWS
            ),
        ],

        [
            "stage15A_unique_actual_focals",
            focals[
                [
                    "genome",
                    "protein_id",
                ]
            ]
            .drop_duplicates()
            .shape[0],
            EXPECTED_UNIQUE_ACTUAL_FOCALS,
            int(
                focals[
                    [
                        "genome",
                        "protein_id",
                    ]
                ]
                .drop_duplicates()
                .shape[0]
                ==
                EXPECTED_UNIQUE_ACTUAL_FOCALS
            ),
        ],

        [
            "fresh_neighborhood_rows",
            len(
                neigh
            ),
            EXPECTED_FRESH_NEIGHBORHOOD_ROWS,
            int(
                len(
                    neigh
                )
                ==
                EXPECTED_FRESH_NEIGHBORHOOD_ROWS
            ),
        ],

        [
            "interesting_actual_gene_association_rows",
            len(
                actual
            ),
            "",
            1,
        ],

        [
            "actual_association_focals",
            actual[
                "focal_id"
            ].nunique(),
            "",
            1,
        ],

        [
            "family_resolvable_gene_associations",
            len(
                family_obs
            ),
            "",
            1,
        ],

        [
            "exact_family_association_rows",
            len(
                exact
            ),
            "",
            1,
        ],

        [
            "gene_window_association_rows",
            len(
                windows
            ),
            "",
            1,
        ],

        [
            "compact_family_associations",
            len(
                family_summary
            ),
            "",
            1,
        ],

        [
            "fegenie_neighbor_rows",
            len(
                fe_obs
            ),
            "",
            1,
        ],

        [
            "unclustered_interesting_neighbor_rows",
            len(
                unclustered
            ),
            "",
            1,
        ],

        [
            "fixed_region_architectures",
            len(
                region_arch
            ),
            EXPECTED_FIXED_REGIONS,
            int(
                len(
                    region_arch
                )
                ==
                EXPECTED_FIXED_REGIONS
            ),
        ],

        [
            "gene_window_monotonicity_failures",
            monotonic_failures,
            0,
            int(
                monotonic_failures
                ==
                0
            ),
        ],
    ],
    columns=[
        "metric",
        "observed",
        "expected",
        "pass",
    ],
)


write_tsv(
    qc,
    OUT_QC,
)


## ================================================================== ##
## Terminal summary
## ================================================================== ##

print()
print("Gene-level association summary")


print(
    f"  Actual focal genes:                    "
    f"{len(focals):,}"
)

print(
    f"  Interesting actual gene-pair rows:     "
    f"{len(actual):,}"
)

print(
    f"  Family-resolvable pair rows:           "
    f"{len(family_obs):,}"
)

print(
    f"  Exact family-offset rows:              "
    f"{len(exact):,}"
)

print(
    f"  Cumulative family-window rows:         "
    f"{len(windows):,}"
)

print(
    f"  Compact family associations:           "
    f"{len(family_summary):,}"
)

print(
    f"  FeGenie-positive neighbor rows:        "
    f"{len(fe_obs):,}"
)

print(
    f"  Unclustered interesting neighbor rows: "
    f"{len(unclustered):,}"
)

print(
    f"  Fixed-region architectures:            "
    f"{len(region_arch):,}"
)

print(
    f"  Window monotonicity failures:          "
    f"{monotonic_failures:,}"
)


## ================================================================== ##
## Module 20 / Module 35 spotlights
## ================================================================== ##

for spotlight in [
    "Module_20",
    "Module_35",
]:

    print()
    print("=" * 80)
    print(
        f"{spotlight} GENE-LEVEL ASSOCIATION SPOTLIGHT"
    )
    print("=" * 80)


    spot_actual = actual[
        actual[
            "focal_context_module"
        ]
        ==
        spotlight
    ]


    spot_family = family_summary[
        family_summary[
            "context_module"
        ]
        ==
        spotlight
    ].copy()


    print(
        f"Actual gene-pair rows:       "
        f"{len(spot_actual):,}"
    )

    print(
        f"Actual focal genes:          "
        f"{spot_actual['focal_id'].nunique():,}"
    )

    print(
        f"Distinct neighbor genes:     "
        f"{spot_actual[['genome','neighbor_protein_id']].drop_duplicates().shape[0]:,}"
    )

    print(
        f"Same-module neighbor rows:   "
        f"{int((spot_actual['neighbor_module_relationship'] == 'same_context_module').sum()):,}"
    )

    print(
        f"Other-module neighbor rows:  "
        f"{int((spot_actual['neighbor_module_relationship'] == 'different_mcl_module').sum()):,}"
    )

    print(
        f"Family associations:         "
        f"{len(spot_family):,}"
    )


    if len(
        spot_family
    ) > 0:

        spot_family[
            "_rank"
        ] = (
            spot_family[
                "n_genomes_with_neighbor"
            ]
        )


        show_columns = [
            "focal_cluster",
            "neighbor_cluster",
            "direction",

            "n_focal_genes_with_neighbor",
            "n_genomes_with_neighbor",

            "best_exact_gene_offset",
            "best_exact_occurrence_support",

            "occurrence_saturation_gene_radius",
            "occurrence_saturation_support",

            "genome_saturation_gene_radius",
            "genome_saturation_support",

            "median_intergenic_gap_bp",

            "dominant_pair_orientation",
            "transcriptional_family_order_distribution",
        ]


        print()
        print(
            spot_family
            .sort_values(
                [
                    "_rank",
                    "n_focal_genes_with_neighbor",
                    "focal_cluster",
                    "neighbor_cluster",
                ],
                ascending=[
                    False,
                    False,
                    True,
                    True,
                ],
            )[
                show_columns
            ]
            .head(
                30
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
    f"Actual gene associations:    "
    f"{OUT_ACTUAL}"
)

print(
    f"Focal-gene summaries:        "
    f"{OUT_FOCAL_SUMMARY}"
)

print(
    f"Region architectures:        "
    f"{OUT_REGION_ARCH}"
)

print(
    f"Exact family offsets:        "
    f"{OUT_EXACT}"
)

print(
    f"Family gene windows:         "
    f"{OUT_WINDOWS}"
)

print(
    f"Family association summary:  "
    f"{OUT_FAMILY_SUMMARY}"
)

print(
    f"FeGenie neighbor summary:    "
    f"{OUT_FEGENIE}"
)

print(
    f"Unclustered annotations:     "
    f"{OUT_UNCLUSTERED}"
)

print(
    f"QC:                          "
    f"{OUT_QC}"
)
