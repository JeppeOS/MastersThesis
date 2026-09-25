#!/usr/bin/env python3

from pathlib import Path
import sys

import numpy as np
import pandas as pd


## ================================================================== ##
## STAGE 66
##
## POSITIONAL ARCHITECTURE SUMMARY FOR CLUSTER_00121
##
## For every C121_local family:
##
##   - quantify presence across 19 GlobDB + 4 MAG loci;
##   - select the occurrence nearest the focal gene in each locus;
##   - summarize relative gene position;
##   - summarize relative strand orientation;
##   - identify shared recurrent families;
##   - inspect recurrent unknown families.
##
## No biological family is promoted or removed automatically.
## ================================================================== ##


PROJECT = (
    Path.home()
    / "methanotrophs"
    / "methanotroph_project"
    / "jeppe"
    / "metagenome"
)

WORK = (
    PROJECT
    / "comparative_analysis"
    / "Cluster_00121"
)

LOCAL = (
    WORK
    / "local_clustering"
)


GENES = (
    WORK
    / "Cluster_00121_combined_neighborhood_genes_local_families.tsv"
)

SUMMARY = (
    LOCAL
    / "C121_local_family_summary.tsv"
)


OUT_STATS = (
    LOCAL
    / "C121_family_positional_statistics.tsv"
)

OUT_SHARED = (
    LOCAL
    / "C121_shared_recurrent_architecture_families.tsv"
)

OUT_UNKNOWN = (
    LOCAL
    / "C121_shared_recurrent_unknown_architecture_families.tsv"
)


EXPECTED_REGIONS = 23
EXPECTED_GLOBDB = 19
EXPECTED_MAG = 4


def fail(message):

    print(
        f"ERROR: {message}",
        file=sys.stderr,
    )

    sys.exit(1)


def clean(value):

    if value is None or pd.isna(value):
        return ""

    return str(value).strip()


for path in [GENES, SUMMARY]:

    if not path.is_file():
        fail(f"Missing input: {path}")


## ================================================================== ##
## Read data
## ================================================================== ##

genes = pd.read_csv(
    GENES,
    sep="\t",
    dtype=str,
    keep_default_na=False,
)

summary = pd.read_csv(
    SUMMARY,
    sep="\t",
    dtype=str,
    keep_default_na=False,
)


required = {
    "focal_id",
    "source_dataset",
    "local_family_id",
    "oriented_gene_offset",
    "oriented_midpoint_offset_bp",
    "oriented_neighbor_strand",
    "is_focal",
}

missing = required - set(genes.columns)

if missing:

    fail(
        "Gene table missing columns: "
        + ", ".join(sorted(missing))
    )


if genes["focal_id"].nunique() != EXPECTED_REGIONS:

    fail(
        f"Expected {EXPECTED_REGIONS} regions; "
        f"found {genes['focal_id'].nunique()}."
    )


genes["oriented_gene_offset"] = pd.to_numeric(
    genes["oriented_gene_offset"],
    errors="raise",
)

genes["oriented_midpoint_offset_bp"] = pd.to_numeric(
    genes["oriented_midpoint_offset_bp"],
    errors="raise",
)

genes["is_focal"] = pd.to_numeric(
    genes["is_focal"],
    errors="coerce",
).fillna(0).astype(int)


## ================================================================== ##
## Choose nearest occurrence of each family in each focal region
##
## This prevents a locus with two paralogues from contributing twice
## to positional statistics.
## ================================================================== ##

nearest = (
    genes
    .assign(
        abs_gene_offset=lambda x:
            x["oriented_gene_offset"].abs(),

        abs_bp_offset=lambda x:
            x["oriented_midpoint_offset_bp"].abs(),
    )
    .sort_values(
        [
            "focal_id",
            "local_family_id",
            "abs_gene_offset",
            "abs_bp_offset",
            "neighbor_protein_id",
        ],
        kind="stable",
    )
    .drop_duplicates(
        [
            "focal_id",
            "local_family_id",
        ]
    )
    .reset_index(drop=True)
)


## ================================================================== ##
## Family positional statistics
## ================================================================== ##

rows = []


for family, group in nearest.groupby(
    "local_family_id",
    sort=False,
):

    globdb = group[
        group["source_dataset"] == "GlobDB"
    ]

    mag = group[
        group["source_dataset"] == "Metagenome"
    ]


    offsets = group[
        "oriented_gene_offset"
    ].astype(float)

    bp_offsets = group[
        "oriented_midpoint_offset_bp"
    ].astype(float)


    strands = group[
        "oriented_neighbor_strand"
    ].map(clean)


    strand_counts = (
        strands[
            strands != ""
        ]
        .value_counts()
    )


    if len(strand_counts) > 0:

        dominant_strand = (
            strand_counts.index[0]
        )

        dominant_strand_fraction = (
            strand_counts.iloc[0]
            /
            strand_counts.sum()
        )

    else:

        dominant_strand = ""

        dominant_strand_fraction = np.nan


    exact_offset_counts = (
        offsets
        .astype(int)
        .value_counts()
    )


    dominant_gene_offset = int(
        exact_offset_counts.index[0]
    )

    dominant_gene_offset_fraction = (
        exact_offset_counts.iloc[0]
        /
        len(group)
    )


    rows.append(
        {
            "local_family_id":
                family,

            "n_regions":
                group["focal_id"].nunique(),

            "n_GlobDB_regions":
                globdb["focal_id"].nunique(),

            "n_metagenome_regions":
                mag["focal_id"].nunique(),

            "prevalence_all_23":
                group["focal_id"].nunique()
                /
                EXPECTED_REGIONS,

            "prevalence_GlobDB_19":
                globdb["focal_id"].nunique()
                /
                EXPECTED_GLOBDB,

            "prevalence_metagenome_4":
                mag["focal_id"].nunique()
                /
                EXPECTED_MAG,

            "median_oriented_gene_offset":
                float(
                    offsets.median()
                ),

            "min_oriented_gene_offset":
                int(
                    offsets.min()
                ),

            "max_oriented_gene_offset":
                int(
                    offsets.max()
                ),

            "dominant_exact_gene_offset":
                dominant_gene_offset,

            "dominant_exact_gene_offset_fraction":
                dominant_gene_offset_fraction,

            "median_oriented_midpoint_bp":
                float(
                    bp_offsets.median()
                ),

            "min_oriented_midpoint_bp":
                float(
                    bp_offsets.min()
                ),

            "max_oriented_midpoint_bp":
                float(
                    bp_offsets.max()
                ),

            "dominant_oriented_strand":
                dominant_strand,

            "dominant_oriented_strand_fraction":
                dominant_strand_fraction,

            "n_forward_relative_to_focal":
                int(
                    (
                        strands == "+"
                    ).sum()
                ),

            "n_reverse_relative_to_focal":
                int(
                    (
                        strands == "-"
                    ).sum()
                ),

            "is_focal_family":
                int(
                    group[
                        "is_focal"
                    ].max()
                    ==
                    1
                ),
        }
    )


stats = pd.DataFrame(rows)


## ================================================================== ##
## Add annotation/family summaries
## ================================================================== ##

stats = stats.merge(
    summary,
    on="local_family_id",
    how="left",
    validate="one_to_one",
    suffixes=(
        "_position",
        "_family",
    ),
)


if stats[
    "n_proteins"
].isna().any():

    fail(
        "Could not map all positional families "
        "to Stage-65 family summary."
    )


## ================================================================== ##
## Simple positional-consistency descriptors
##
## These are descriptive labels only.
## ================================================================== ##

stats[
    "position_span_genes"
] = (
    stats[
        "max_oriented_gene_offset"
    ]
    -
    stats[
        "min_oriented_gene_offset"
    ]
)


stats[
    "high_strand_consistency"
] = (
    pd.to_numeric(
        stats[
            "dominant_oriented_strand_fraction"
        ],
        errors="coerce",
    )
    >=
    0.80
).astype(int)


stats[
    "high_exact_position_consistency"
] = (
    pd.to_numeric(
        stats[
            "dominant_exact_gene_offset_fraction"
        ],
        errors="coerce",
    )
    >=
    0.50
).astype(int)


stats[
    "shared_between_datasets"
] = (
    (
        pd.to_numeric(
            stats[
                "n_GlobDB_regions_position"
            ],
            errors="coerce",
        )
        >
        0
    )
    &
    (
        pd.to_numeric(
            stats[
                "n_metagenome_regions_position"
            ],
            errors="coerce",
        )
        >
        0
    )
).astype(int)


## ================================================================== ##
## Sort by biological usefulness for architecture inspection
## ================================================================== ##

stats = (
    stats
    .sort_values(
        [
            "is_focal_family",
            "shared_between_datasets",
            "n_regions",
            "dominant_exact_gene_offset_fraction",
            "dominant_oriented_strand_fraction",
            "local_family_id",
        ],
        ascending=[
            False,
            False,
            False,
            False,
            False,
            True,
        ],
        kind="stable",
    )
    .reset_index(drop=True)
)


stats.to_csv(
    OUT_STATS,
    sep="\t",
    index=False,
    na_rep="",
    float_format="%.6f",
)


## ================================================================== ##
## Shared recurrent families
## ================================================================== ##

shared = stats[
    (
        stats[
            "shared_between_datasets"
        ]
        ==
        1
    )
    &
    (
        pd.to_numeric(
            stats[
                "n_regions"
            ],
            errors="coerce",
        )
        >=
        2
    )
].copy()


shared.to_csv(
    OUT_SHARED,
    sep="\t",
    index=False,
    na_rep="",
    float_format="%.6f",
)


## ================================================================== ##
## Shared recurrent UNKNOWN families
## ================================================================== ##

unknown = shared[
    pd.to_numeric(
        shared[
            "functionally_unannotated"
        ],
        errors="coerce",
    )
    .fillna(0)
    .astype(int)
    ==
    1
].copy()


unknown.to_csv(
    OUT_UNKNOWN,
    sep="\t",
    index=False,
    na_rep="",
    float_format="%.6f",
)


## ================================================================== ##
## Terminal summary
## ================================================================== ##

print("=" * 80)

print(
    "STAGE 66 - CLUSTER_00121 POSITIONAL ARCHITECTURE"
)

print("=" * 80)

print()

print(
    f"Families analysed:                 "
    f"{len(stats)}"
)

print(
    f"Shared recurrent families:         "
    f"{len(shared)}"
)

print(
    f"Shared recurrent unknown families: "
    f"{len(unknown)}"
)


print()
print(
    "Most prevalent shared families"
)

print(
    shared[
        [
            "local_family_id",
            "n_regions",
            "n_GlobDB_regions_position",
            "n_metagenome_regions_position",
            "median_oriented_gene_offset",
            "position_span_genes",
            "dominant_oriented_strand",
            "dominant_oriented_strand_fraction",
            "dominant_cog",
            "dominant_informative_product",
        ]
    ]
    .head(20)
    .to_string(
        index=False
    )
)


print()
print(
    "Shared recurrent UNKNOWN families"
)

if len(
    unknown
) == 0:

    print(
        "  None."
    )

else:

    print(
        unknown[
            [
                "local_family_id",
                "n_regions",
                "n_GlobDB_regions_position",
                "n_metagenome_regions_position",
                "median_oriented_gene_offset",
                "position_span_genes",
                "dominant_exact_gene_offset",
                "dominant_exact_gene_offset_fraction",
                "dominant_oriented_strand",
                "dominant_oriented_strand_fraction",
            ]
        ]
        .to_string(
            index=False
        )
    )


print()
print("Outputs:")

print(
    f"  {OUT_STATS}"
)

print(
    f"  {OUT_SHARED}"
)

print(
    f"  {OUT_UNKNOWN}"
)
