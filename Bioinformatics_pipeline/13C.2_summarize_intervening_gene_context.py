#!/usr/bin/env python3

from pathlib import Path
from collections import Counter
import sys

import numpy as np
import pandas as pd


## ================================================================== ##
## STAGE 13C2 - INTERVENING GENE CONTEXT
##
## Purpose:
##
## Identify the actual Prodigal genes lying between the closest
## same-contig occurrences of two MMseqs2 families belonging to the
## same MCL module.
##
## Example:
##
##     Cluster_A -> X -> Cluster_B
##
## This stage does NOT redefine family relationships.
##
## Primary family identity remains:
##
##     MMseqs2 cluster
##
## Annotation fields:
##
##     COG      = functional grouping / evidence
##     product  = human-readable interpretation
##
## Unannotated Prodigal genes are retained.
##
##
## Outputs:
##
## 1. one row per actual intervening gene
## 2. one row per non-adjacent pair observation
## 3. recurring intervening-pattern summary
##
## Percentages always retain numerator and denominator.
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
    / "module_pair_intervening_genes.tsv"
)


OUT_OBSERVATIONS = (
    HERE
    / "module_pair_intervening_observations.tsv"
)


OUT_PATTERNS = (
    HERE
    / "module_pair_intervening_pattern_summary.tsv"
)


OUT_QC = (
    HERE
    / "intervening_gene_context_qc.tsv"
)


## ================================================================== ##
## Expected upstream values
## ================================================================== ##

EXPECTED_PAIR_OBSERVATIONS = 34_037
EXPECTED_SAME_CONTIG = 8_080


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
print("STAGE 13C2 - INTERVENING GENE CONTEXT")
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

        "closest_contig",
        "closest_genes_between",

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


pairs[
    "any_same_contig"
] = pd.to_numeric(
    pairs[
        "any_same_contig"
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


same_contig = pairs[
    pairs[
        "any_same_contig"
    ]
    ==
    1
].copy()


if len(
    same_contig
) != EXPECTED_SAME_CONTIG:

    fail(
        f"Expected "
        f"{EXPECTED_SAME_CONTIG:,} "
        f"same-contig pair observations but found "
        f"{len(same_contig):,}."
    )


nonadjacent = same_contig[
    same_contig[
        "closest_genes_between"
    ]
    >
    0
].copy()


print(
    f"  Pair observations:                "
    f"{len(pairs):,}"
)

print(
    f"  Same-contig observations:         "
    f"{len(same_contig):,}"
)

print(
    f"  Same-contig with intervening gene:"
    f" {len(nonadjacent):,}"
)


## ================================================================== ##
## 2. Read master Prodigal gene catalogue
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
        "neighbor_globdb_cog",
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
        "neighbor_globdb_product",
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
## 3. Build compact catalogue
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


if GENE_COLUMN is not None:

    rename_map[
        GENE_COLUMN
    ] = "gene"


if CLUSTER_COLUMN is not None:

    rename_map[
        CLUSTER_COLUMN
    ] = "cluster"


if MODULE_COLUMN is not None:

    rename_map[
        MODULE_COLUMN
    ] = "module"


if STRAND_COLUMN is not None:

    rename_map[
        STRAND_COLUMN
    ] = "strand"


catalog_small = catalog_small.rename(
    columns=rename_map
)


for optional in [
    "gene",
    "cluster",
    "module",
    "strand",
]:

    if optional not in catalog_small.columns:

        catalog_small[
            optional
        ] = ""


## ================================================================== ##
## 4. Recover ranks of the two closest family proteins
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
        "Could not recover catalogue ranks for "
        "all closest pair proteins."
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


## ================================================================== ##
## 5. Strict coordinate/rank QC
## ================================================================== ##

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
        "Closest-left catalogue contig mismatch."
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
        "Closest-right catalogue contig mismatch."
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
        "Stage-12 closest_genes_between disagrees "
        "with catalogue gene ranks."
    )


## ================================================================== ##
## 6. Pre-index catalogue by genome + contig
## ================================================================== ##

print()
print("Extracting intervening genes...")


catalog_groups = {
    key:
        group
        .sort_values(
            "contig_gene_rank"
        )
        .copy()

    for key, group
    in catalog_small.groupby(
        [
            "genome",
            "contig",
        ],
        sort=False,
    )
}


## ================================================================== ##
## 7. Extract each intervening gene
## ================================================================== ##

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
            f"Catalogue contig missing: "
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
            f"Expected {expected} intervening genes but "
            f"found {len(between)} for "
            f"{pair.genome} / "
            f"{pair.module} / "
            f"{pair.cluster_a} <-> {pair.cluster_b}."
        )


    between = between.sort_values(
        "contig_gene_rank"
    )


    cog_pattern = []
    product_pattern = []
    cluster_pattern = []
    functional_pattern = []


    for position, gene in enumerate(
        between.itertuples(
            index=False
        ),
        start=1,
    ):

        gene_dict = gene._asdict()


        cog = clean(
            gene_dict.get(
                "cog",
                "",
            )
        )


        product = clean(
            gene_dict.get(
                "product",
                "",
            )
        )


        cluster = clean(
            gene_dict.get(
                "cluster",
                "",
            )
        )


        module = clean(
            gene_dict.get(
                "module",
                "",
            )
        )


        gene_name = clean(
            gene_dict.get(
                "gene",
                "",
            )
        )


        strand = clean(
            gene_dict.get(
                "strand",
                "",
            )
        )


        ## Functional pattern prioritizes sequence family,
        ## then COG, then explicitly unannotated. ##

        if cluster != "":

            functional_label = (
                f"MMSEQ:{cluster}"
            )

        elif cog != "":

            functional_label = (
                f"COG:{cog}"
            )

        else:

            functional_label = (
                "UNANNOTATED"
            )


        cog_pattern.append(
            cog
            if cog != ""
            else
            "NA"
        )


        product_pattern.append(
            product
            if product != ""
            else
            "NA"
        )


        cluster_pattern.append(
            cluster
            if cluster != ""
            else
            "NA"
        )


        functional_pattern.append(
            functional_label
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

                "closest_left_protein":
                    pair.closest_left_protein,

                "closest_left_cluster":
                    pair.closest_left_cluster,

                "closest_left_strand":
                    pair.closest_left_strand,

                "closest_right_protein":
                    pair.closest_right_protein,

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
                    gene_dict[
                        "protein_id"
                    ],

                "intervening_gene_rank":
                    gene_dict[
                        "contig_gene_rank"
                    ],

                "intervening_strand":
                    strand,

                "intervening_mmseqs_cluster":
                    cluster,

                "intervening_module":
                    module,

                "intervening_cog":
                    cog,

                "intervening_gene":
                    gene_name,

                "intervening_product":
                    product,

                "functional_identity":
                    functional_label,
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

            "closest_left_cluster":
                pair.closest_left_cluster,

            "closest_right_cluster":
                pair.closest_right_cluster,

            "closest_orientation_class":
                pair.closest_orientation_class,

            "n_intervening_genes":
                expected,

            "intervening_mmseqs_pattern":
                " | ".join(
                    cluster_pattern
                ),

            "intervening_cog_pattern":
                " | ".join(
                    cog_pattern
                ),

            "intervening_product_pattern":
                " | ".join(
                    product_pattern
                ),

            "intervening_functional_pattern":
                " | ".join(
                    functional_pattern
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
## 8. Pattern-level support
##
## Important:
##
## Two denominators are retained:
##
## 1. all same-contig observations of this family pair
## 2. non-adjacent same-contig observations of this family pair
##
## Thus a pattern occurring 1/1 remains visibly based on n=1.
## ================================================================== ##

print()
print("Summarizing recurring intervening patterns...")


same_pair_counts = (
    same_contig
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
                "n_same_contig_pair_observations",
        }
    )
)


nonadj_pair_counts = (
    observations
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
                "n_nonadjacent_pair_observations",
        }
    )
)


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
                "n_pair_observations_with_pattern",
        }
    )
)


patterns = (
    patterns
    .merge(
        same_pair_counts,
        on=[
            "module",
            "cluster_a",
            "cluster_b",
        ],
        how="left",
        validate="many_to_one",
    )
    .merge(
        nonadj_pair_counts,
        on=[
            "module",
            "cluster_a",
            "cluster_b",
        ],
        how="left",
        validate="many_to_one",
    )
)


patterns[
    "support_among_nonadjacent"
] = (
    patterns[
        "n_pair_observations_with_pattern"
    ]
    .astype(str)
    +
    "/"
    +
    patterns[
        "n_nonadjacent_pair_observations"
    ]
    .astype(str)
)


patterns[
    "pct_among_nonadjacent"
] = (
    100.0
    *
    patterns[
        "n_pair_observations_with_pattern"
    ]
    /
    patterns[
        "n_nonadjacent_pair_observations"
    ]
)


patterns[
    "support_among_all_same_contig"
] = (
    patterns[
        "n_pair_observations_with_pattern"
    ]
    .astype(str)
    +
    "/"
    +
    patterns[
        "n_same_contig_pair_observations"
    ]
    .astype(str)
)


patterns[
    "pct_among_all_same_contig"
] = (
    100.0
    *
    patterns[
        "n_pair_observations_with_pattern"
    ]
    /
    patterns[
        "n_same_contig_pair_observations"
    ]
)


## ================================================================== ##
## 9. Sort and write
## ================================================================== ##

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


patterns = (
    patterns
    .sort_values(
        [
            "module",
            "cluster_a",
            "cluster_b",
            "n_pair_observations_with_pattern",
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


## ================================================================== ##
## 10. QC
## ================================================================== ##

if len(
    genes
) != int(
    nonadjacent[
        "closest_genes_between"
    ]
    .sum()
):

    fail(
        "Total extracted intervening-gene count "
        "does not equal Stage-12 expected total."
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
            "same_contig_pair_observations",
            len(
                same_contig
            ),
        ],

        [
            "nonadjacent_same_contig_pair_observations",
            len(
                nonadjacent
            ),
        ],

        [
            "expected_intervening_gene_rows",
            int(
                nonadjacent[
                    "closest_genes_between"
                ]
                .sum()
            ),
        ],

        [
            "extracted_intervening_gene_rows",
            len(
                genes
            ),
        ],

        [
            "intervening_genes_with_mmseqs_cluster",
            int(
                (
                    genes[
                        "intervening_mmseqs_cluster"
                    ]
                    !=
                    ""
                )
                .sum()
            ),
        ],

        [
            "intervening_genes_with_cog",
            int(
                (
                    genes[
                        "intervening_cog"
                    ]
                    !=
                    ""
                )
                .sum()
            ),
        ],

        [
            "intervening_genes_with_product",
            int(
                (
                    genes[
                        "intervening_product"
                    ]
                    !=
                    ""
                )
                .sum()
            ),
        ],

        [
            "distinct_intervening_patterns",
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
## 11. Terminal summary
## ================================================================== ##

print()
print("Intervening-gene summary")

print(
    f"  Same-contig pair observations:       "
    f"{len(same_contig):,}"
)

print(
    f"  Non-adjacent pair observations:      "
    f"{len(nonadjacent):,}"
)

print(
    f"  Intervening genes extracted:         "
    f"{len(genes):,}"
)

print(
    f"  Intervening genes with MMseqs family:"
    f" {(genes['intervening_mmseqs_cluster'] != '').sum():,}"
)

print(
    f"  Intervening genes with COG:          "
    f"{(genes['intervening_cog'] != '').sum():,}"
)

print(
    f"  Distinct intervening patterns:       "
    f"{len(patterns):,}"
)


## ================================================================== ##
## 12. Spotlight Module 20 / Module 35
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
    ]


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
                "n_intervening_genes",
                "intervening_position_left_to_right",
                "intervening_protein_id",
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
    f"Intervening genes:        "
    f"{OUT_GENES}"
)

print(
    f"Pair observations:        "
    f"{OUT_OBSERVATIONS}"
)

print(
    f"Pattern summary:          "
    f"{OUT_PATTERNS}"
)

print(
    f"QC:                       "
    f"{OUT_QC}"
)
