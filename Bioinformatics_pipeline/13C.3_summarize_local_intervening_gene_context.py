#!/usr/bin/env python3

from pathlib import Path
import sys

import numpy as np
import pandas as pd


## ================================================================== ##
## STAGE 13C2 - LOCAL INTERVENING GENE CONTEXT
##
## Authoritative local version of the intervening-gene analysis.
##
## PRIMARY FRAMEWORK:
##
##     MCL-associated MMseqs2 protein-family pairs
##
## LOCALITY REQUIREMENT:
##
##     closest same-contig pair must be within <=20 kb
##     using the Stage-12A intergenic-gap definition.
##
## We additionally retain:
##
##     adjacency
##     <=5 kb
##     <=10 kb
##     <=20 kb
##
## This avoids treating widely separated proteins on the same contig
## as members of a local genomic neighborhood.
##
## MMseqs2 remains the primary sequence-family identity.
## GlobDB COG/product are secondary annotation evidence.
##
## Unannotated Prodigal genes are retained.
##
## The unrestricted Stage-13C2 output is NOT deleted; it remains an
## archival same-contig analysis.
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


GENE_CATALOG = (
    WORKFLOW
    / "11_genome_annotations"
    / "annotated_gene_catalog.tsv"
)


OUT_GENES = (
    HERE
    / "local_module_pair_intervening_genes.tsv"
)


OUT_OBSERVATIONS = (
    HERE
    / "local_module_pair_intervening_observations.tsv"
)


OUT_PATTERNS = (
    HERE
    / "local_module_pair_intervening_pattern_summary.tsv"
)


OUT_PAIR_SUMMARY = (
    HERE
    / "local_module_pair_intervening_pair_summary.tsv"
)


OUT_QC = (
    HERE
    / "local_intervening_gene_context_qc.tsv"
)


## ================================================================== ##
## Expected upstream values
## ================================================================== ##

EXPECTED_PAIR_OBSERVATIONS = 34_037
EXPECTED_WITHIN_20KB = 2_757


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


def detect_column(
    df,
    candidates,
    label,
    required=True,
):

    for column in candidates:

        if column in df.columns:

            return column


    if required:

        fail(
            f"Could not identify {label} column.\n"
            f"Tried: "
            +
            ", ".join(
                candidates
            )
            +
            "\nAvailable columns:\n"
            +
            ", ".join(
                df.columns
            )
        )


    return None


def clean(value):

    if pd.isna(value):

        return ""

    return str(
        value
    ).strip()


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
print("STAGE 13C2 - LOCAL INTERVENING GENE CONTEXT")
print("=" * 80)


## ================================================================== ##
## 1. Read Stage-12A pair observations
## ================================================================== ##

print()
print("Reading module-family pair observations...")


pairs = pd.read_csv(
    PAIR_OBSERVATIONS,
    sep="\t",
    dtype=str,
    keep_default_na=False,
)


require_columns(
    pairs,
    [
        "genome",
        "module",
        "cluster_a",
        "cluster_b",

        "any_same_contig",
        "any_adjacent",
        "any_within_5kb",
        "any_within_10kb",
        "any_within_20kb",

        "closest_contig",
        "closest_genes_between",
        "closest_intergenic_gap_bp",

        "closest_left_protein",
        "closest_left_cluster",
        "closest_left_strand",

        "closest_right_protein",
        "closest_right_cluster",
        "closest_right_strand",

        "closest_orientation_class",
    ],
    PAIR_OBSERVATIONS.name,
)


if len(
    pairs
) != EXPECTED_PAIR_OBSERVATIONS:

    fail(
        f"Expected "
        f"{EXPECTED_PAIR_OBSERVATIONS:,} "
        f"pair observations but found "
        f"{len(pairs):,}."
    )


for column in [
    "any_same_contig",
    "any_adjacent",
    "any_within_5kb",
    "any_within_10kb",
    "any_within_20kb",
]:

    pairs[
        column
    ] = pd.to_numeric(
        pairs[
            column
        ],
        errors="raise",
    ).astype(int)


pairs[
    "closest_genes_between"
] = pd.to_numeric(
    pairs[
        "closest_genes_between"
    ].replace(
        "",
        np.nan,
    ),
    errors="coerce",
)


pairs[
    "closest_intergenic_gap_bp"
] = pd.to_numeric(
    pairs[
        "closest_intergenic_gap_bp"
    ].replace(
        "",
        np.nan,
    ),
    errors="coerce",
)


## ================================================================== ##
## 2. Restrict to authoritative LOCAL pair observations
##
## any_within_20kb already uses the Stage-12A intergenic-gap
## definition, so this reproduces the established physical criterion.
## ================================================================== ##

local = pairs[
    pairs[
        "any_within_20kb"
    ]
    ==
    1
].copy()


if len(
    local
) != EXPECTED_WITHIN_20KB:

    fail(
        f"Expected "
        f"{EXPECTED_WITHIN_20KB:,} "
        f"within-20-kb observations but found "
        f"{len(local):,}."
    )


## Every <=20-kb observation must have same-contig evidence. ##

if (
    local[
        "any_same_contig"
    ]
    !=
    1
).any():

    fail(
        "A <=20-kb pair lacks same-contig evidence."
    )


## Validate the Stage-12 flag against the actual closest gap. ##

if (
    local[
        "closest_intergenic_gap_bp"
    ]
    >
    20_000
).any():

    fail(
        "A pair flagged within 20 kb has "
        "closest_intergenic_gap_bp > 20,000."
    )


## Physical hierarchy must hold. ##

if (
    local[
        "any_within_5kb"
    ]
    >
    local[
        "any_within_10kb"
    ]
).any():

    fail(
        "5-kb flag occurs without 10-kb flag."
    )


if (
    local[
        "any_within_10kb"
    ]
    >
    local[
        "any_within_20kb"
    ]
).any():

    fail(
        "10-kb flag occurs without 20-kb flag."
    )


nonadjacent = local[
    local[
        "closest_genes_between"
    ]
    >
    0
].copy()


adjacent = local[
    local[
        "closest_genes_between"
    ]
    ==
    0
].copy()


print(
    f"  All pair observations:         "
    f"{len(pairs):,}"
)

print(
    f"  Local <=20-kb observations:    "
    f"{len(local):,}"
)

print(
    f"  Local <=10-kb observations:    "
    f"{local['any_within_10kb'].sum():,}"
)

print(
    f"  Local <=5-kb observations:     "
    f"{local['any_within_5kb'].sum():,}"
)

print(
    f"  Local adjacent observations:   "
    f"{len(adjacent):,}"
)

print(
    f"  Local non-adjacent observations:"
    f" {len(nonadjacent):,}"
)


## ================================================================== ##
## 3. Read master Prodigal gene catalogue
## ================================================================== ##

print()
print("Reading annotated Prodigal gene catalogue...")


catalog = pd.read_csv(
    GENE_CATALOG,
    sep="\t",
    dtype=str,
    keep_default_na=False,
)


require_columns(
    catalog,
    [
        "genome",
        "protein_id",
        "contig",
        "contig_gene_rank",
    ],
    GENE_CATALOG.name,
)


COG_COLUMN = detect_column(
    catalog,
    [
        "globdb_cog",
        "cog",
        "COG",
        "cog_id",
        "COG_id",
    ],
    "COG",
)


PRODUCT_COLUMN = detect_column(
    catalog,
    [
        "globdb_product",
        "product",
        "annotation_product",
    ],
    "product",
)


GENE_COLUMN = detect_column(
    catalog,
    [
        "globdb_gene",
        "gene",
        "gene_name",
    ],
    "gene",
    required=False,
)


CLUSTER_COLUMN = detect_column(
    catalog,
    [
        "cluster",
        "mmseqs_cluster",
    ],
    "MMseqs cluster",
    required=False,
)


MODULE_COLUMN = detect_column(
    catalog,
    [
        "module",
        "mcl_module",
    ],
    "module",
    required=False,
)


STRAND_COLUMN = detect_column(
    catalog,
    [
        "strand",
    ],
    "strand",
    required=False,
)


ANNOTATION_ACCEPTED_COLUMN = detect_column(
    catalog,
    [
        "annotation_accepted",
    ],
    "annotation accepted",
    required=False,
)


ANNOTATION_MATCH_COLUMN = detect_column(
    catalog,
    [
        "annotation_match_type",
    ],
    "annotation match type",
    required=False,
)


catalog[
    "contig_gene_rank"
] = pd.to_numeric(
    catalog[
        "contig_gene_rank"
    ],
    errors="raise",
).astype(int)


print(
    f"  Catalogue genes: "
    f"{len(catalog):,}"
)

print(
    f"  COG column:      "
    f"{COG_COLUMN}"
)

print(
    f"  Product column:  "
    f"{PRODUCT_COLUMN}"
)


## ================================================================== ##
## 4. Compact catalogue
## ================================================================== ##

keep_columns = [
    "genome",
    "protein_id",
    "contig",
    "contig_gene_rank",
    COG_COLUMN,
    PRODUCT_COLUMN,
]


for optional in [
    GENE_COLUMN,
    CLUSTER_COLUMN,
    MODULE_COLUMN,
    STRAND_COLUMN,
    ANNOTATION_ACCEPTED_COLUMN,
    ANNOTATION_MATCH_COLUMN,
]:

    if (
        optional is not None
        and
        optional not in keep_columns
    ):

        keep_columns.append(
            optional
        )


catalog_small = catalog[
    keep_columns
].copy()


rename_map = {
    COG_COLUMN:
        "cog",

    PRODUCT_COLUMN:
        "product",
}


optional_rename = [
    (
        GENE_COLUMN,
        "gene",
    ),
    (
        CLUSTER_COLUMN,
        "cluster",
    ),
    (
        MODULE_COLUMN,
        "module",
    ),
    (
        STRAND_COLUMN,
        "strand",
    ),
    (
        ANNOTATION_ACCEPTED_COLUMN,
        "annotation_accepted",
    ),
    (
        ANNOTATION_MATCH_COLUMN,
        "annotation_match_type",
    ),
]


for old, new in optional_rename:

    if old is not None:

        rename_map[
            old
        ] = new


catalog_small = catalog_small.rename(
    columns=rename_map
)


for column in [
    "gene",
    "cluster",
    "module",
    "strand",
    "annotation_accepted",
    "annotation_match_type",
]:

    if column not in catalog_small.columns:

        catalog_small[
            column
        ] = ""


## ================================================================== ##
## 5. Recover left/right ranks
## ================================================================== ##

protein_ranks = catalog_small[
    [
        "genome",
        "protein_id",
        "contig",
        "contig_gene_rank",
    ]
].copy()


left_lookup = protein_ranks.rename(
    columns={
        "protein_id":
            "closest_left_protein",

        "contig":
            "catalog_left_contig",

        "contig_gene_rank":
            "left_rank",
    }
)


right_lookup = protein_ranks.rename(
    columns={
        "protein_id":
            "closest_right_protein",

        "contig":
            "catalog_right_contig",

        "contig_gene_rank":
            "right_rank",
    }
)


nonadjacent = (
    nonadjacent
    .merge(
        left_lookup,
        on=[
            "genome",
            "closest_left_protein",
        ],
        how="left",
        validate="many_to_one",
    )
    .merge(
        right_lookup,
        on=[
            "genome",
            "closest_right_protein",
        ],
        how="left",
        validate="many_to_one",
    )
)


if (
    nonadjacent[
        "left_rank"
    ]
    .isna()
    .any()
    or
    nonadjacent[
        "right_rank"
    ]
    .isna()
    .any()
):

    fail(
        "Could not recover catalogue ranks "
        "for all local non-adjacent pair proteins."
    )


nonadjacent[
    "left_rank"
] = nonadjacent[
    "left_rank"
].astype(int)


nonadjacent[
    "right_rank"
] = nonadjacent[
    "right_rank"
].astype(int)


if (
    nonadjacent[
        "catalog_left_contig"
    ]
    !=
    nonadjacent[
        "closest_contig"
    ]
).any():

    fail(
        "Closest-left contig mismatch."
    )


if (
    nonadjacent[
        "catalog_right_contig"
    ]
    !=
    nonadjacent[
        "closest_contig"
    ]
).any():

    fail(
        "Closest-right contig mismatch."
    )


calculated_between = (
    nonadjacent[
        "right_rank"
    ]
    -
    nonadjacent[
        "left_rank"
    ]
    -
    1
)


if (
    calculated_between
    !=
    nonadjacent[
        "closest_genes_between"
    ]
).any():

    fail(
        "Catalogue gene ranks disagree with "
        "Stage-12 closest_genes_between."
    )


## ================================================================== ##
## 6. Pre-index only contigs actually needed
##
## This avoids unnecessarily building groups for the entire genome
## catalogue when only local pair contigs matter.
## ================================================================== ##

print()
print("Preparing local contigs...")


needed_contigs = (
    nonadjacent[
        [
            "genome",
            "closest_contig",
        ]
    ]
    .drop_duplicates()
    .rename(
        columns={
            "closest_contig":
                "contig",
        }
    )
)


catalog_needed = catalog_small.merge(
    needed_contigs,
    on=[
        "genome",
        "contig",
    ],
    how="inner",
    validate="many_to_one",
)


catalog_groups = {
    key:
        group
        .sort_values(
            "contig_gene_rank"
        )
        .copy()

    for key, group
    in catalog_needed.groupby(
        [
            "genome",
            "contig",
        ],
        sort=False,
    )
}


print(
    f"  Local contigs required: "
    f"{len(needed_contigs):,}"
)


## ================================================================== ##
## 7. Extract intervening genes
## ================================================================== ##

print()
print("Extracting local intervening genes...")


gene_rows = []
observation_rows = []


for pair in nonadjacent.itertuples(
    index=False
):

    key = (
        pair.genome,
        pair.closest_contig,
    )


    if key not in catalog_groups:

        fail(
            f"Missing catalogue contig: "
            f"{pair.genome} / "
            f"{pair.closest_contig}"
        )


    contig = catalog_groups[
        key
    ]


    between = contig[
        (
            contig[
                "contig_gene_rank"
            ]
            >
            pair.left_rank
        )
        &
        (
            contig[
                "contig_gene_rank"
            ]
            <
            pair.right_rank
        )
    ].copy()


    expected = int(
        pair.closest_genes_between
    )


    if len(
        between
    ) != expected:

        fail(
            f"Expected {expected} intervening genes "
            f"but found {len(between)} for "
            f"{pair.genome} / "
            f"{pair.cluster_a} <-> {pair.cluster_b}."
        )


    between = between.sort_values(
        "contig_gene_rank"
    )


    functional_pattern = []
    cog_pattern = []
    mmseqs_pattern = []
    product_pattern = []


    for position, gene in enumerate(
        between.itertuples(
            index=False
        ),
        start=1,
    ):

        g = gene._asdict()


        cog = clean(
            g.get(
                "cog",
                "",
            )
        )


        product = clean(
            g.get(
                "product",
                "",
            )
        )


        mmseqs_cluster = clean(
            g.get(
                "cluster",
                "",
            )
        )


        module = clean(
            g.get(
                "module",
                "",
            )
        )


        gene_name = clean(
            g.get(
                "gene",
                "",
            )
        )


        strand = clean(
            g.get(
                "strand",
                "",
            )
        )


        annotation_accepted = clean(
            g.get(
                "annotation_accepted",
                "",
            )
        )


        annotation_match_type = clean(
            g.get(
                "annotation_match_type",
                "",
            )
        )


        if mmseqs_cluster != "":

            functional_identity = (
                f"MMSEQ:{mmseqs_cluster}"
            )

        elif cog != "":

            functional_identity = (
                f"COG:{cog}"
            )

        else:

            functional_identity = (
                "UNANNOTATED"
            )


        functional_pattern.append(
            functional_identity
        )


        cog_pattern.append(
            cog
            if cog != ""
            else
            "NA"
        )


        mmseqs_pattern.append(
            mmseqs_cluster
            if mmseqs_cluster != ""
            else
            "NA"
        )


        product_pattern.append(
            product
            if product != ""
            else
            "NA"
        )


        gene_rows.append(
            {
                "genome":
                    pair.genome,

                "module":
                    pair.module,

                "cluster_a":
                    pair.cluster_a,

                "cluster_b":
                    pair.cluster_b,

                "any_adjacent":
                    pair.any_adjacent,

                "any_within_5kb":
                    pair.any_within_5kb,

                "any_within_10kb":
                    pair.any_within_10kb,

                "any_within_20kb":
                    pair.any_within_20kb,

                "closest_intergenic_gap_bp":
                    pair.closest_intergenic_gap_bp,

                "closest_left_cluster":
                    pair.closest_left_cluster,

                "closest_left_strand":
                    pair.closest_left_strand,

                "closest_right_cluster":
                    pair.closest_right_cluster,

                "closest_right_strand":
                    pair.closest_right_strand,

                "closest_orientation_class":
                    pair.closest_orientation_class,

                "n_intervening_genes":
                    expected,

                "intervening_position_left_to_right":
                    position,

                "intervening_protein_id":
                    g[
                        "protein_id"
                    ],

                "intervening_gene_rank":
                    g[
                        "contig_gene_rank"
                    ],

                "intervening_strand":
                    strand,

                "intervening_mmseqs_cluster":
                    mmseqs_cluster,

                "intervening_module":
                    module,

                "intervening_cog":
                    cog,

                "intervening_gene":
                    gene_name,

                "intervening_product":
                    product,

                "annotation_accepted":
                    annotation_accepted,

                "annotation_match_type":
                    annotation_match_type,

                "functional_identity":
                    functional_identity,
            }
        )


    observation_rows.append(
        {
            "genome":
                pair.genome,

            "module":
                pair.module,

            "cluster_a":
                pair.cluster_a,

            "cluster_b":
                pair.cluster_b,

            "any_adjacent":
                pair.any_adjacent,

            "any_within_5kb":
                pair.any_within_5kb,

            "any_within_10kb":
                pair.any_within_10kb,

            "any_within_20kb":
                pair.any_within_20kb,

            "closest_intergenic_gap_bp":
                pair.closest_intergenic_gap_bp,

            "closest_left_cluster":
                pair.closest_left_cluster,

            "closest_right_cluster":
                pair.closest_right_cluster,

            "closest_orientation_class":
                pair.closest_orientation_class,

            "n_intervening_genes":
                expected,

            "intervening_functional_pattern":
                " | ".join(
                    functional_pattern
                ),

            "intervening_mmseqs_pattern":
                " | ".join(
                    mmseqs_pattern
                ),

            "intervening_cog_pattern":
                " | ".join(
                    cog_pattern
                ),

            "intervening_product_pattern":
                " | ".join(
                    product_pattern
                ),
        }
    )


genes = pd.DataFrame(
    gene_rows
)


observations = pd.DataFrame(
    observation_rows
)


## ================================================================== ##
## 8. Pair-level local denominator summary
##
## This is important for interpretation:
##
##     1/1 = 100%
##     7/7 = 100%
##     44/44 = 100%
##
## remain visibly different evidence strengths.
## ================================================================== ##

pair_summary = (
    local
    .groupby(
        [
            "module",
            "cluster_a",
            "cluster_b",
        ],
        as_index=False,
    )
    .agg(
        n_local_within_20kb=(
            "any_within_20kb",
            "sum",
        ),

        n_local_within_10kb=(
            "any_within_10kb",
            "sum",
        ),

        n_local_within_5kb=(
            "any_within_5kb",
            "sum",
        ),

        n_local_adjacent=(
            "any_adjacent",
            "sum",
        ),
    )
)


nonadj_counts = (
    nonadjacent
    .groupby(
        [
            "module",
            "cluster_a",
            "cluster_b",
        ],
        as_index=False,
    )
    .size()
    .rename(
        columns={
            "size":
                "n_local_nonadjacent",
        }
    )
)


pair_summary = pair_summary.merge(
    nonadj_counts,
    on=[
        "module",
        "cluster_a",
        "cluster_b",
    ],
    how="left",
    validate="one_to_one",
)


pair_summary[
    "n_local_nonadjacent"
] = (
    pair_summary[
        "n_local_nonadjacent"
    ]
    .fillna(
        0
    )
    .astype(int)
)


pair_summary[
    "adjacency_support"
] = (
    pair_summary[
        "n_local_adjacent"
    ]
    .astype(str)
    +
    "/"
    +
    pair_summary[
        "n_local_within_20kb"
    ]
    .astype(str)
)


pair_summary[
    "pct_local_observations_adjacent"
] = (
    100.0
    *
    pair_summary[
        "n_local_adjacent"
    ]
    /
    pair_summary[
        "n_local_within_20kb"
    ]
)


pair_summary[
    "within_5kb_support"
] = (
    pair_summary[
        "n_local_within_5kb"
    ]
    .astype(str)
    +
    "/"
    +
    pair_summary[
        "n_local_within_20kb"
    ]
    .astype(str)
)


pair_summary[
    "pct_local_observations_within_5kb"
] = (
    100.0
    *
    pair_summary[
        "n_local_within_5kb"
    ]
    /
    pair_summary[
        "n_local_within_20kb"
    ]
)


## ================================================================== ##
## 9. Pattern support
##
## Denominators:
##
##     all local <=20-kb observations
##     local non-adjacent observations
##
## A pattern based on n=1 therefore remains visibly n=1.
## ================================================================== ##

print()
print("Summarizing local intervening patterns...")


if len(
    observations
) > 0:

    patterns = (
        observations
        .groupby(
            [
                "module",
                "cluster_a",
                "cluster_b",
                "n_intervening_genes",
                "intervening_functional_pattern",
                "intervening_cog_pattern",
            ],
            as_index=False,
        )
        .size()
        .rename(
            columns={
                "size":
                    "n_local_observations_with_pattern",
            }
        )
    )


    patterns = patterns.merge(
        pair_summary[
            [
                "module",
                "cluster_a",
                "cluster_b",
                "n_local_within_20kb",
                "n_local_nonadjacent",
            ]
        ],
        on=[
            "module",
            "cluster_a",
            "cluster_b",
        ],
        how="left",
        validate="many_to_one",
    )


    patterns[
        "support_among_local_nonadjacent"
    ] = (
        patterns[
            "n_local_observations_with_pattern"
        ]
        .astype(str)
        +
        "/"
        +
        patterns[
            "n_local_nonadjacent"
        ]
        .astype(str)
    )


    patterns[
        "pct_among_local_nonadjacent"
    ] = (
        100.0
        *
        patterns[
            "n_local_observations_with_pattern"
        ]
        /
        patterns[
            "n_local_nonadjacent"
        ]
    )


    patterns[
        "support_among_all_local"
    ] = (
        patterns[
            "n_local_observations_with_pattern"
        ]
        .astype(str)
        +
        "/"
        +
        patterns[
            "n_local_within_20kb"
        ]
        .astype(str)
    )


    patterns[
        "pct_among_all_local"
    ] = (
        100.0
        *
        patterns[
            "n_local_observations_with_pattern"
        ]
        /
        patterns[
            "n_local_within_20kb"
        ]
    )

else:

    patterns = pd.DataFrame()


## ================================================================== ##
## 10. Sort + write
## ================================================================== ##

if len(
    genes
) > 0:

    genes = (
        genes
        .sort_values(
            [
                "module",
                "cluster_a",
                "cluster_b",
                "genome",
                "intervening_position_left_to_right",
            ],
            kind="stable",
        )
        .reset_index(
            drop=True
        )
    )


if len(
    observations
) > 0:

    observations = (
        observations
        .sort_values(
            [
                "module",
                "cluster_a",
                "cluster_b",
                "genome",
            ],
            kind="stable",
        )
        .reset_index(
            drop=True
        )
    )


if len(
    patterns
) > 0:

    patterns = (
        patterns
        .sort_values(
            [
                "module",
                "cluster_a",
                "cluster_b",
                "n_local_observations_with_pattern",
                "intervening_functional_pattern",
            ],
            ascending=[
                True,
                True,
                True,
                False,
                True,
            ],
            kind="stable",
        )
        .reset_index(
            drop=True
        )
    )


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


genes.to_csv(
    OUT_GENES,
    sep="\t",
    index=False,
    na_rep="",
)


observations.to_csv(
    OUT_OBSERVATIONS,
    sep="\t",
    index=False,
    na_rep="",
)


patterns.to_csv(
    OUT_PATTERNS,
    sep="\t",
    index=False,
    na_rep="",
    float_format="%.6f",
)


pair_summary.to_csv(
    OUT_PAIR_SUMMARY,
    sep="\t",
    index=False,
    na_rep="",
    float_format="%.6f",
)


## ================================================================== ##
## 11. QC
## ================================================================== ##

expected_gene_rows = int(
    nonadjacent[
        "closest_genes_between"
    ]
    .sum()
)


if len(
    genes
) != expected_gene_rows:

    fail(
        f"Expected "
        f"{expected_gene_rows:,} "
        f"intervening genes but extracted "
        f"{len(genes):,}."
    )


qc = pd.DataFrame(
    [
        [
            "all_pair_observations",
            len(
                pairs
            ),
        ],

        [
            "local_within_20kb_pair_observations",
            len(
                local
            ),
        ],

        [
            "local_within_10kb_pair_observations",
            int(
                local[
                    "any_within_10kb"
                ]
                .sum()
            ),
        ],

        [
            "local_within_5kb_pair_observations",
            int(
                local[
                    "any_within_5kb"
                ]
                .sum()
            ),
        ],

        [
            "local_adjacent_pair_observations",
            len(
                adjacent
            ),
        ],

        [
            "local_nonadjacent_pair_observations",
            len(
                nonadjacent
            ),
        ],

        [
            "expected_local_intervening_gene_rows",
            expected_gene_rows,
        ],

        [
            "extracted_local_intervening_gene_rows",
            len(
                genes
            ),
        ],

        [
            "local_intervening_genes_with_mmseqs_cluster",
            (
                int(
                    (
                        genes[
                            "intervening_mmseqs_cluster"
                        ]
                        !=
                        ""
                    )
                    .sum()
                )
                if len(
                    genes
                ) > 0
                else
                0
            ),
        ],

        [
            "local_intervening_genes_with_cog",
            (
                int(
                    (
                        genes[
                            "intervening_cog"
                        ]
                        !=
                        ""
                    )
                    .sum()
                )
                if len(
                    genes
                ) > 0
                else
                0
            ),
        ],

        [
            "local_intervening_genes_with_product",
            (
                int(
                    (
                        genes[
                            "intervening_product"
                        ]
                        !=
                        ""
                    )
                    .sum()
                )
                if len(
                    genes
                ) > 0
                else
                0
            ),
        ],

        [
            "local_intervening_patterns",
            len(
                patterns
            ),
        ],

        [
            "gene_count_qc_pass",
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
print("Local intervening-gene summary")


print(
    f"  <=20-kb pair observations:      "
    f"{len(local):,}"
)

print(
    f"  <=10-kb pair observations:      "
    f"{local['any_within_10kb'].sum():,}"
)

print(
    f"  <=5-kb pair observations:       "
    f"{local['any_within_5kb'].sum():,}"
)

print(
    f"  Adjacent pair observations:     "
    f"{len(adjacent):,}"
)

print(
    f"  Non-adjacent local observations:"
    f" {len(nonadjacent):,}"
)

print(
    f"  Intervening genes extracted:    "
    f"{len(genes):,}"
)

print(
    f"  With MMseqs2 family:            "
    f"{(genes['intervening_mmseqs_cluster'] != '').sum():,}"
)

print(
    f"  With COG:                       "
    f"{(genes['intervening_cog'] != '').sum():,}"
)

print(
    f"  With product:                   "
    f"{(genes['intervening_product'] != '').sum():,}"
)

print(
    f"  Distinct local patterns:        "
    f"{len(patterns):,}"
)


## ================================================================== ##
## 13. Spotlight Module 20 / 35
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


for cluster_x, cluster_y in spotlights:

    spot = genes[
        (
            (
                genes[
                    "cluster_a"
                ]
                ==
                cluster_x
            )
            &
            (
                genes[
                    "cluster_b"
                ]
                ==
                cluster_y
            )
        )
        |
        (
            (
                genes[
                    "cluster_a"
                ]
                ==
                cluster_y
            )
            &
            (
                genes[
                    "cluster_b"
                ]
                ==
                cluster_x
            )
        )
    ].copy()


    print()
    print(
        f"Spotlight: "
        f"{cluster_x} <-> {cluster_y}"
    )


    if len(
        spot
    ) == 0:

        print(
            "  No intervening genes."
        )

        continue


    print(
        spot[
            [
                "genome",
                "closest_left_cluster",
                "closest_right_cluster",
                "closest_intergenic_gap_bp",
                "n_intervening_genes",
                "intervening_position_left_to_right",
                "intervening_protein_id",
                "annotation_accepted",
                "annotation_match_type",
                "intervening_mmseqs_cluster",
                "intervening_cog",
                "intervening_gene",
                "intervening_product",
            ]
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
    f"Local intervening genes:  "
    f"{OUT_GENES}"
)

print(
    f"Local observations:       "
    f"{OUT_OBSERVATIONS}"
)

print(
    f"Local pattern summary:    "
    f"{OUT_PATTERNS}"
)

print(
    f"Local pair summary:       "
    f"{OUT_PAIR_SUMMARY}"
)

print(
    f"QC:                       "
    f"{OUT_QC}"
)
