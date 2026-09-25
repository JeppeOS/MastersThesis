#!/usr/bin/env python3


## ================================================================== ##
## STAGE 16A3 - MODULE 10 REPRESENTATIVE GENE-MAP SELECTION
##
## Purpose
## -------
##
## Lock the 12 representative Cluster_00050 / Cyc2 focal regions used
## for the final Module-10 gene-map figure.
##
## Selection logic:
##
##   1. Deliberately retain several canonical Methylobacter 4/5 loci
##      to demonstrate repeated conservation of the dominant
##      Methylobacter Cyc2-centered architecture.
##
##   2. Include alternative 4/5 architectures differing in gene order
##      or Module-10 partner-family composition.
##
##   3. Represent both locally resolved Cluster_00222 variants.
##
##   4. Represent distinct 3/5 architectures.
##
##   5. Represent distinct 2/5 architectures, including both
##      Cluster_00064 local-family variants and the unique
##      Cyc2 + Cluster_00069 architecture.
##
## This script performs no biological reclassification.
## It only selects focal regions and attaches plot order / rationale.
## ================================================================== ##


from pathlib import Path

import pandas as pd


## ================================================================== ##
## Paths
## ================================================================== ##

WORKFLOW = Path(
    "~/methanotrophs/methanotroph_project/jeppe/"
    "genome_analysis_workflow"
).expanduser()


INPUT_GENES = (
    WORKFLOW
    / "16_visualization"
    / "16A2_local_architecture_pipeline"
    / "runs"
    / "Module_10"
    / "Module_10_gene_map_data_with_local_families.tsv"
)


INPUT_CANDIDATES = (
    WORKFLOW
    / "16_visualization"
    / "16A3_module10_representative_selection"
    / "Module_10_representative_candidate_table.tsv"
)


OUTPUT_DIR = (
    WORKFLOW
    / "16_visualization"
    / "16A3_module10_representative_selection"
)


OUTPUT_REGIONS = (
    OUTPUT_DIR
    / "Module_10_selected_regions.tsv"
)


OUTPUT_GENES = (
    OUTPUT_DIR
    / "Module_10_selected_gene_map_data.tsv"
)


OUTPUT_QC = (
    OUTPUT_DIR
    / "Module_10_selected_gene_map_qc.tsv"
)


## ================================================================== ##
## Locked representative set
##
## plot_order = top-to-bottom order in the final figure.
## ================================================================== ##

SELECTED = [

    {
        "plot_order": 1,
        "genome": "GCA_023229605",
        "selection_class": "canonical_methylobacter_4of5",
        "selection_reason":
            "Canonical Methylobacter 4/5 architecture; clean representative."
    },

    {
        "plot_order": 2,
        "genome": "GCF_022788635",
        "selection_class": "canonical_methylobacter_4of5",
        "selection_reason":
            "Second canonical Methylobacter 4/5 representative from a distinct Methylobacter lineage."
    },

    {
        "plot_order": 3,
        "genome": "SPIREOTU_01747602",
        "selection_class": "canonical_methylobacter_4of5",
        "selection_reason":
            "Canonical Methylobacter 4/5 core with additional recurrent flanking context."
    },

    {
        "plot_order": 4,
        "genome": "GCA_003584895",
        "selection_class": "canonical_methylobacter_4of5",
        "selection_reason":
            "Canonical Methylobacter 4/5 core with a distinct recurrent upstream context."
    },

    {
        "plot_order": 5,
        "genome": "GCF_039904985",
        "selection_class": "reordered_4of5",
        "selection_reason":
            "Same 4/5 Module-10 family composition as canonical architecture but altered gene order."
    },

    {
        "plot_order": 6,
        "genome": "GCA_039795615",
        "selection_class": "alternative_4of5_C00222_A",
        "selection_reason":
            "Alternative 4/5 architecture containing Cluster_00222 local variant M10_local_031."
    },

    {
        "plot_order": 7,
        "genome": "MOTU40_080091",
        "selection_class": "alternative_4of5_C00222_B",
        "selection_reason":
            "Alternative 4/5 architecture containing distinct Cluster_00222 local variant M10_local_077."
    },

    {
        "plot_order": 8,
        "genome": "GCF_018734325",
        "selection_class": "three_of_five_C00064",
        "selection_reason":
            "Unique 3/5 Cyc2 + Cluster_00091 + Cluster_00064 architecture."
    },

    {
        "plot_order": 9,
        "genome": "GCA_041158225",
        "selection_class": "three_of_five_C00069",
        "selection_reason":
            "3/5 Cyc2 + Cluster_00091 + Cluster_00069 architecture with distinct recurrent context."
    },

    {
        "plot_order": 10,
        "genome": "GCA_003584865",
        "selection_class": "methylobacter_two_of_five",
        "selection_reason":
            "Methylobacter 2/5 contrast with Cluster_00064 local variant M10_local_002."
    },

    {
        "plot_order": 11,
        "genome": "GCF_006175985",
        "selection_class": "two_of_five_C00064_variant_B",
        "selection_reason":
            "Methylotetracoccus 2/5 architecture containing Cluster_00064 local variant M10_local_007."
    },

    {
        "plot_order": 12,
        "genome": "MOTU40_108429",
        "selection_class": "two_of_five_C00069",
        "selection_reason":
            "Unique Cyc2 + Cluster_00069 2/5 architecture."
    },

]


selected = pd.DataFrame(
    SELECTED
)


## ================================================================== ##
## Helpers
## ================================================================== ##

def fail(message):

    raise RuntimeError(
        message
    )


## ================================================================== ##
## Input validation
## ================================================================== ##

for path in [
    INPUT_GENES,
    INPUT_CANDIDATES,
]:

    if not path.exists():

        fail(
            f"Required input does not exist:\n{path}"
        )


if len(
    selected
) != 12:

    fail(
        "Expected exactly 12 selected genomes."
    )


if selected[
    "genome"
].duplicated().any():

    fail(
        "Selected genome list contains duplicates."
    )


## ================================================================== ##
## Read inputs
## ================================================================== ##

genes = pd.read_csv(
    INPUT_GENES,
    sep="\t",
    dtype=str,
    keep_default_na=False
)


candidates = pd.read_csv(
    INPUT_CANDIDATES,
    sep="\t",
    dtype=str,
    keep_default_na=False
)


## ================================================================== ##
## Validate selected genomes against architecture catalogue
## ================================================================== ##

candidate_selected = candidates.loc[
    candidates[
        "genome"
    ].isin(
        selected[
            "genome"
        ]
    )
].copy()


missing = sorted(
    set(
        selected[
            "genome"
        ]
    )
    -
    set(
        candidate_selected[
            "genome"
        ]
    )
)


if missing:

    fail(
        "Selected genomes missing from Module-10 candidate table:\n"
        +
        "\n".join(
            missing
        )
    )


## ------------------------------------------------------------------ ##
## One focal region per selected genome is required here.
## ------------------------------------------------------------------ ##

region_counts = (

    candidate_selected

    .groupby(
        "genome"
    )[
        "plot_region_id"
    ]

    .nunique()

)


bad = region_counts.loc[
    region_counts != 1
]


if not bad.empty:

    print(
        bad
    )

    fail(
        "A selected genome contains more than one Module-10 focal "
        "region. Selection must then be made by plot_region_id."
    )


## ================================================================== ##
## Build selected-region table
## ================================================================== ##

region_columns = [

    "genome",
    "plot_region_id",

    "taxonomy_display",
    "priority_taxon",

    "n_module_clusters_10kb",

    "global_architecture_pattern_id",
    "local_architecture_id",
    "local_architecture_n_regions",

    "module_local_family_order",
    "recurrent_local_family_order",

]


selected_regions = (

    candidate_selected[
        region_columns
    ]

    .drop_duplicates()

    .merge(
        selected,
        on="genome",
        how="inner",
        validate="one_to_one"
    )

    .sort_values(
        "plot_order"
    )

    .reset_index(
        drop=True
    )

)


if len(
    selected_regions
) != 12:

    fail(
        f"Expected 12 selected regions, found {len(selected_regions)}."
    )


## ================================================================== ##
## Plot labels
## ================================================================== ##

selected_regions[
    "track_label"
] = selected_regions[
    "taxonomy_display"
]


## ------------------------------------------------------------------ ##
## Protect against duplicated display labels.
## ------------------------------------------------------------------ ##

duplicate_labels = selected_regions[
    "track_label"
].duplicated(
    keep=False
)


selected_regions.loc[
    duplicate_labels,
    "track_label"
] = (

    selected_regions.loc[
        duplicate_labels,
        "track_label"
    ]

    +
    " ["

    +
    selected_regions.loc[
        duplicate_labels,
        "genome"
    ]

    +
    "]"

)


## ================================================================== ##
## Subset gene-map data
## ================================================================== ##

selected_genes = genes.loc[
    genes[
        "plot_region_id"
    ].isin(
        selected_regions[
            "plot_region_id"
        ]
    )
].copy()


selected_genes = selected_genes.merge(

    selected_regions[
        [
            "plot_region_id",

            "plot_order",
            "track_label",

            "selection_class",
            "selection_reason",

            "global_architecture_pattern_id",
            "local_architecture_id",
            "n_module_clusters_10kb",
        ]
    ],

    on="plot_region_id",

    how="left",

    validate="many_to_one"

)


## ================================================================== ##
## QC
## ================================================================== ##

if selected_genes[
    "plot_region_id"
].nunique() != 12:

    fail(
        "Selected gene table does not contain exactly 12 focal regions."
    )


if selected_genes[
    "genome"
].nunique() != 12:

    fail(
        "Selected gene table does not contain exactly 12 genomes."
    )


if selected_genes[
    [
        "plot_region_id",
        "protein_id"
    ]
].duplicated().any():

    fail(
        "Duplicate plot_region_id + protein_id rows detected."
    )


focal_counts = (

    selected_genes.loc[
        selected_genes[
            "display_class"
        ]
        ==
        "focal_cluster00050"
    ]

    .groupby(
        "plot_region_id"
    )

    .size()

)


if len(
    focal_counts
) != 12:

    fail(
        "Not every selected region contains a focal_cluster00050 row."
    )


if not (
    focal_counts
    ==
    1
).all():

    fail(
        "Expected exactly one focal Cluster_00050 gene per region."
    )


## ================================================================== ##
## Final ordering
## ================================================================== ##

selected_genes = selected_genes.sort_values(

    [
        "plot_order",
        "plot_start_bp",
        "plot_end_bp",
        "protein_id",
    ]

).reset_index(
    drop=True
)


## ================================================================== ##
## Write outputs
## ================================================================== ##

selected_regions.to_csv(
    OUTPUT_REGIONS,
    sep="\t",
    index=False
)


selected_genes.to_csv(
    OUTPUT_GENES,
    sep="\t",
    index=False
)


qc = pd.DataFrame(
    [
        {
            "metric": "n_selected_regions",
            "value": selected_genes[
                "plot_region_id"
            ].nunique()
        },

        {
            "metric": "n_selected_genomes",
            "value": selected_genes[
                "genome"
            ].nunique()
        },

        {
            "metric": "n_selected_gene_rows",
            "value": len(
                selected_genes
            )
        },

        {
            "metric": "n_unique_selected_genes",
            "value": selected_genes[
                [
                    "genome",
                    "protein_id"
                ]
            ].drop_duplicates().shape[0]
        },

        {
            "metric": "n_focal_cyc2_rows",
            "value": int(
                (
                    selected_genes[
                        "display_class"
                    ]
                    ==
                    "focal_cluster00050"
                ).sum()
            )
        },
    ]
)


qc.to_csv(
    OUTPUT_QC,
    sep="\t",
    index=False
)


## ================================================================== ##
## Report
## ================================================================== ##

print(
    "=" * 120
)

print(
    "LOCKED MODULE-10 REPRESENTATIVES"
)

print(
    "=" * 120
)


print(

    selected_regions[
        [
            "plot_order",
            "genome",
            "track_label",

            "n_module_clusters_10kb",

            "global_architecture_pattern_id",
            "local_architecture_id",

            "selection_class",
        ]
    ]

    .to_string(
        index=False
    )

)


print()

print(
    "=" * 120
)

print(
    "QC"
)

print(
    "=" * 120
)


print(
    qc.to_string(
        index=False
    )
)


print()

print(
    f"Selected regions:\n{OUTPUT_REGIONS}"
)

print()

print(
    f"Selected genes:\n{OUTPUT_GENES}"
)

print()

print(
    f"QC:\n{OUTPUT_QC}"
)
