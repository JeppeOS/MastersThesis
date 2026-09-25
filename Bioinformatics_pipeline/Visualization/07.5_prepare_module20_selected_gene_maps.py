#!/usr/bin/env python3

## ================================================================== ##
## MODULE 20 - PREPARE FINAL REPRESENTATIVE GENE-MAP DATA
##
## Purpose
## -------
## Create the plotting dataset for the provisional 12 representative
## Module-20 focal regions.
##
## Selection strategy:
##
##   - cover all six P1 Methylobacter architecture classes
##   - include a second example of the dominant P1 architecture
##   - include all three MtrA loci
##   - include two FeGenie-negative contrast architectures
##
## This is a REPRESENTATIVE architecture set, not a frequency-weighted
## sample of the 74 Module-20 regions.
## ================================================================== ##


from pathlib import Path

import pandas as pd


## ================================================================== ##
## Paths
## ================================================================== ##

ROOT = (
    Path.home()
    / "methanotrophs"
    / "methanotroph_project"
    / "jeppe"
    / "genome_analysis_workflow"
)


ARCH_DIR = (
    ROOT
    / "16_visualization"
    / "16A3_module20_representative_selection"
)


GENE_FILE = (
    ARCH_DIR
    / "Module_20_gene_map_data_enriched.tsv"
)


MATRIX_FILE = (
    ARCH_DIR
    / "Module_20_region_architecture_matrix.tsv"
)


OUT_MANIFEST = (
    ARCH_DIR
    / "Module_20_selected_regions.tsv"
)


OUT_GENES = (
    ARCH_DIR
    / "Module_20_selected_gene_map_data.tsv"
)


OUT_QC = (
    ARCH_DIR
    / "Module_20_selected_gene_map_qc.tsv"
)


## ================================================================== ##
## Explicit representative selection
##
## IMPORTANT:
##
## Focal protein ID is part of the key because some genomes contain
## multiple Cluster_00035 loci with different architectures.
## ================================================================== ##

SELECTION = [

    {
        "selection_order": 1,
        "genome": "GCA_003158415",
        "focal_protein_id": "GCA_003158415_000000000196_4",
        "selection_class": "P1",
        "selection_reason": "Dominant Methylobacter_C architecture; Cluster00048 + MtrB local004 + L003",
    },

    {
        "selection_order": 2,
        "genome": "GCA_012960135",
        "focal_protein_id": "GCA_012960135_000000000005_3",
        "selection_class": "architecture_diversity",
        "selection_reason": (
            "Recurrent Methyloprofundus accessory-rich architecture; "
            "MtoA + MtrB local009 + Cluster00166 + Cluster00209 + "
            "four-heme NapC-like family + multicopper oxidase"
        ),
    },

    {
        "selection_order": 3,
        "genome": "MOTU40_052632",
        "focal_protein_id": "MOTU40_052632_000000000219_2",
        "selection_class": "P1",
        "selection_reason": "P1 architecture 02; Cluster00048 + MtrB local004 without L003",
    },

    {
        "selection_order": 4,
        "genome": "GCF_003994235",
        "focal_protein_id": "GCF_003994235_000000000094_32",
        "selection_class": "P1",
        "selection_reason": "Methylobacter_A oryzae; MtrB local010 without Cluster00048",
    },

    {
        "selection_order": 5,
        "genome": "GCF_015476545",
        "focal_protein_id": "GCF_015476545_000000000055_2",
        "selection_class": "P1",
        "selection_reason": "MtrB local010 plus three local015 paralogues",
    },

    {
        "selection_order": 6,
        "genome": "GCA_903873045",
        "focal_protein_id": "GCA_903873045_000000000273_1",
        "selection_class": "P1",
        "selection_reason": "Singleton P1 architecture 05; MtrB local004 + L003 without Cluster00048",
    },

    {
        "selection_order": 7,
        "genome": "GCF_000427625",
        "focal_protein_id": "GCF_000427625_000000000001_412",
        "selection_class": "P1",
        "selection_reason": "Singleton P1 architecture 06; Methylobacter luteus; four local015 paralogues",
    },

    {
        "selection_order": 8,
        "genome": "GCA_006844585",
        "focal_protein_id": "GCA_006844585_000000000276_33",
        "selection_class": "MtrA",
        "selection_reason": "MtrA focal; simple MtrA-MtrB local010 architecture",
    },

    {
        "selection_order": 9,
        "genome": "GCF_030717995",
        "focal_protein_id": "GCF_030717995_000000000001_885",
        "selection_class": "MtrA",
        "selection_reason": "MtrA focal; MtrB local009 + Cluster00166 + four-heme NapC-like family",
    },

    {
        "selection_order": 10,
        "genome": "GCF_000421465",
        "focal_protein_id": "GCF_000421465_000000000002_392",
        "selection_class": "MtrA",
        "selection_reason": "MtrA focal; Methylohalobius crimeensis; unique MtrB local424",
    },

    {
        "selection_order": 11,
        "genome": "GCA_003584895",
        "focal_protein_id": "GCA_003584895_000000000033_2",
        "selection_class": "fegenie_negative_methylobacter_contrast",
        "selection_reason": (
        "FeGenie-negative Cluster_00035 homolog lacking adjacent MtrB "
        "and the major Module-20 partner families; also selected as a "
        "cross-figure genome carrying the canonical Module-10 Cyc2 architecture."
        ),
    },

    {
        "selection_order": 12,
        "genome": "GCA_022600935",
        "focal_protein_id": "GCA_022600935_000000000179_2",
        "selection_class": "FeGenie_negative",
        "selection_reason": "Non-Methylobacter negative contrast with Cluster00166 + four-heme NapC-like context",
    },

]


selection = pd.DataFrame(
    SELECTION
)


## ================================================================== ##
## Read data
## ================================================================== ##

genes = pd.read_csv(
    GENE_FILE,
    sep="\t",
    dtype=str,
    keep_default_na=False
)


matrix = pd.read_csv(
    MATRIX_FILE,
    sep="\t",
    dtype=str,
    keep_default_na=False
)


## ================================================================== ##
## Validate unique selected keys
## ================================================================== ##

if selection.duplicated(
    [
        "genome",
        "focal_protein_id"
    ]
).any():

    raise SystemExit(
        "ERROR: duplicate selected focal-region key."
    )


if len(selection) != 12:

    raise SystemExit(
        f"ERROR: expected 12 selections; found {len(selection)}."
    )


## ================================================================== ##
## Attach architecture metadata
## ================================================================== ##

selected_manifest = selection.merge(
    matrix,
    on=[
        "genome",
        "focal_protein_id"
    ],
    how="left",
    validate="one_to_one"
)


missing = selected_manifest[
    "plot_region_id"
].eq(
    ""
) | selected_manifest[
    "plot_region_id"
].isna()


if missing.any():

    print(
        selected_manifest.loc[
            missing,
            [
                "genome",
                "focal_protein_id"
            ]
        ].to_string(
            index=False
        )
    )

    raise SystemExit(
        "ERROR: one or more selected loci were not found in architecture matrix."
    )


## ================================================================== ##
## Pull selected gene-map regions
## ================================================================== ##

selected_region_ids = set(
    selected_manifest[
        "plot_region_id"
    ]
)


selected_genes = genes[
    genes[
        "plot_region_id"
    ].isin(
        selected_region_ids
    )
].copy()


## ================================================================== ##
## Attach selection metadata to every plotted gene
## ================================================================== ##

selection_meta = selected_manifest[
    [
        "plot_region_id",
        "selection_order",
        "selection_class",
        "selection_reason",
    ]
]


selected_genes = selected_genes.merge(
    selection_meta,
    on="plot_region_id",
    how="left",
    validate="many_to_one"
)


selected_genes[
    "selection_order"
] = pd.to_numeric(
    selected_genes[
        "selection_order"
    ],
    errors="raise"
)


selected_manifest[
    "selection_order"
] = pd.to_numeric(
    selected_manifest[
        "selection_order"
    ],
    errors="raise"
)


## ================================================================== ##
## Plot display label
##
## Keep taxonomy as the visible track label.
## Protein/genome IDs remain in the underlying table.
## ================================================================== ##

selected_genes[
    "track_label"
] = selected_genes[
    "taxonomy_display"
]


## ================================================================== ##
## Initial biologically meaningful classes
##
## This does NOT yet decide every final figure color.
##
## It establishes the features that definitely deserve prominence.
## ================================================================== ##

def gene_class(row):

    if str(
        row[
            "plot_is_focal"
        ]
    ) == "1":

        if row[
            "focal_call"
        ] == "MtrA":

            return "focal_MtrA"

        if row[
            "focal_call"
        ] == "MtoA":

            return "focal_MtoA"

        return "focal_FeGenie_negative"


    if str(
        row.get(
            "is_mtrb",
            "0"
        )
    ) == "1":

        return "MtrB"


    cluster = str(
        row[
            "neighbor_mmseq_cluster"
        ]
    )


    if cluster == "Cluster_00048":

        return "CytC551_552_like"


    if cluster == "Cluster_00166":

        return "Cluster00166"


    if cluster == "Cluster_00209":

        return "CytC553_like"


    local_family = str(
        row[
            "local_family_id"
        ]
    )


    if local_family == "M20_local_007":

        return "four_heme_NapC_like"


    if local_family == "M20_local_014":

        return "multicopper_oxidase"


    if local_family == "M20_local_015":

        return "local015"


    return "background_gene"


selected_genes[
    "figure_class"
] = selected_genes.apply(
    gene_class,
    axis=1
)


## ================================================================== ##
## Initial labels
##
## Every definitely highlighted biological component gets a label.
## Final recurrent unknown-family labels can be added after inspection.
## ================================================================== ##

def gene_label(row):

    gene_class_value = row[
        "figure_class"
    ]


    if gene_class_value == "focal_MtoA":
        return "MtoA"


    if gene_class_value == "focal_MtrA":
        return "MtrA"


    if gene_class_value == "focal_FeGenie_negative":
        return "Cluster00035 homolog"


    if gene_class_value == "MtrB":
        return "MtrB"


    if gene_class_value == "CytC551_552_like":
        return "CytC551/552-like"


    if gene_class_value == "Cluster00166":
        return "2-heme cyt.c"


    if gene_class_value == "CytC553_like":
        return "CytC553-like"


    if gene_class_value == "four_heme_NapC_like":
        return "4-heme NapC-like"


    if gene_class_value == "multicopper_oxidase":
        return "Multicopper oxidase"


    if gene_class_value == "local015":
        return "M20-local015"


    return ""


selected_genes[
    "figure_label"
] = selected_genes.apply(
    gene_label,
    axis=1
)


## ================================================================== ##
## Sort
## ================================================================== ##

selected_manifest = selected_manifest.sort_values(
    "selection_order",
    kind="stable"
)


selected_genes = selected_genes.sort_values(
    [
        "selection_order",
        "plot_start_bp"
    ],
    kind="stable"
)


## ================================================================== ##
## Write
## ================================================================== ##

selected_manifest.to_csv(
    OUT_MANIFEST,
    sep="\t",
    index=False
)


selected_genes.to_csv(
    OUT_GENES,
    sep="\t",
    index=False
)


## ================================================================== ##
## QC
## ================================================================== ##

qc_rows = [

    {
        "metric":
            "selected_regions",
        "value":
            selected_manifest[
                "plot_region_id"
            ].nunique()
    },

    {
        "metric":
            "selected_genomes",
        "value":
            selected_manifest[
                "genome"
            ].nunique()
    },

    {
        "metric":
            "selected_gene_rows",
        "value":
            len(
                selected_genes
            )
    },

    {
        "metric":
            "P1_regions",
        "value":
            int(
                (
                    selected_manifest[
                        "selection_class"
                    ]
                    ==
                    "P1"
                ).sum()
            )
    },

    {
        "metric":
            "MtrA_regions",
        "value":
            int(
                (
                    selected_manifest[
                        "selection_class"
                    ]
                    ==
                    "MtrA"
                ).sum()
            )
    },

    {
        "metric":
            "negative_regions",
        "value":
            int(
                (
                    selected_manifest[
                        "selection_class"
                    ]
                    ==
                    "FeGenie_negative"
                ).sum()
            )
    },

    {
        "metric":
            "selected_regions_with_MtrB",
        "value":
            int(
                pd.to_numeric(
                    selected_manifest[
                        "has_adjacent_mtrb"
                    ],
                    errors="coerce"
                )
                .fillna(
                    0
                )
                .sum()
            )
    },

    {
        "metric":
            "selected_regions_with_Cluster00048",
        "value":
            int(
                pd.to_numeric(
                    selected_manifest[
                        "has_cluster00048"
                    ],
                    errors="coerce"
                )
                .fillna(
                    0
                )
                .sum()
            )
    },

]


pd.DataFrame(
    qc_rows
).to_csv(
    OUT_QC,
    sep="\t",
    index=False
)


## ================================================================== ##
## Console report
## ================================================================== ##

print(
    "=" * 110
)

print(
    "MODULE 20 - PROVISIONAL FINAL REPRESENTATIVES"
)

print(
    "=" * 110
)


show = [
    "selection_order",
    "taxonomy_display",
    "genome",
    "focal_call",
    "mtrb_local_families",
    "has_cluster00048",
    "has_cluster00166",
    "has_cluster00209",
    "selection_reason",
]


print(
    selected_manifest[
        show
    ]
    .to_string(
        index=False
    )
)


print()
print(
    "=" * 110
)

print(
    "SELECTED FIGURE CLASSES"
)

print(
    "=" * 110
)


print(
    selected_genes[
        "figure_class"
    ]
    .value_counts()
    .to_string()
)


print()
print(
    "=" * 110
)

print(
    "OUTPUTS"
)

print(
    "=" * 110
)


print(
    OUT_MANIFEST
)

print(
    OUT_GENES
)

print(
    OUT_QC
)


print()
print(
    "SUCCESS"
)
