#!/usr/bin/env python3


## ================================================================== ##
## STAGE 16A3 - MODULE 10 LOCAL-ARCHITECTURE INSPECTION
##
## Purpose
## -------
##
## Integrate:
##
##   1. exhaustive Module-10 local-family analysis
##   2. deduplicated global Module-10 cluster architecture
##
## in order to understand:
##
##   - which M10_local_* families correspond to each of the five
##     global Module-10 MMseq families;
##
##   - whether global Module-10 families split into distinct local
##     sequence-family variants;
##
##   - whether regions having the same 4/5, 3/5 or 2/5 global
##     architecture actually contain different local-family
##     arrangements;
##
##   - which recurrent non-module families occur at conserved
##     positions around Cyc2.
##
##
## IMPORTANT
## ---------
##
## This script is DESCRIPTIVE.
##
## It does not automatically choose final representative genomes.
##
## Representative selection follows biological hierarchy:
##
##   global Module-10 completeness
##       >
##   Methylobacter priority
##       >
##   distinct local-family architecture
##       >
##   taxonomic diversity
##       >
##   technical completeness
##
## No opaque composite score is used.
## ================================================================== ##


from pathlib import Path
from collections import Counter

import numpy as np
import pandas as pd


## ================================================================== ##
## Paths
## ================================================================== ##

WORKFLOW = Path(
    "~/methanotrophs/methanotroph_project/jeppe/"
    "genome_analysis_workflow"
).expanduser()


LOCAL_DATA = (
    WORKFLOW
    / "16_visualization"
    / "16A2_local_architecture_pipeline"
    / "runs"
    / "Module_10"
    / "Module_10_gene_map_data_with_local_families.tsv"
)


GLOBAL_ARCHITECTURE = (
    WORKFLOW
    / "16_visualization"
    / "16A_cluster00050_cyc2_gene_map_framework"
    / "Cluster_00050_region_architecture_deduplicated.tsv"
)


OUTPUT_DIR = (
    WORKFLOW
    / "16_visualization"
    / "16A3_module10_representative_selection"
)


OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)


MODULE_CLUSTERS = [

    "Cluster_00050",
    "Cluster_00064",
    "Cluster_00069",
    "Cluster_00091",
    "Cluster_00222",

]


## ================================================================== ##
## Outputs
## ================================================================== ##

OUT_CATALOG = (
    OUTPUT_DIR
    / "Module_10_local_family_catalog.tsv"
)


OUT_MODULE_MAP = (
    OUTPUT_DIR
    / "Module_10_module_cluster_local_family_map.tsv"
)


OUT_REGION_ARCH = (
    OUTPUT_DIR
    / "Module_10_region_local_architectures.tsv"
)


OUT_PATTERN_TABLE = (
    OUTPUT_DIR
    / "Module_10_local_architecture_patterns.tsv"
)


OUT_CANDIDATES = (
    OUTPUT_DIR
    / "Module_10_representative_candidate_table.tsv"
)


## ================================================================== ##
## Helpers
## ================================================================== ##

def fail(message):

    raise RuntimeError(
        message
    )


def clean(series):

    return (
        series
        .fillna("")
        .astype(str)
        .str.strip()
    )


def unique_join(values):

    vals = sorted(
        {
            str(x).strip()
            for x in values
            if pd.notna(x)
            and str(x).strip()
        }
    )

    return ";".join(
        vals
    )


def dominant(values):

    vals = [

        str(x).strip()

        for x in values

        if pd.notna(x)
        and str(x).strip()

    ]


    if not vals:

        return ""


    return Counter(
        vals
    ).most_common(
        1
    )[0][0]


def numeric(series):

    return pd.to_numeric(
        series,
        errors="coerce"
    )


## ================================================================== ##
## Read data
## ================================================================== ##

if not LOCAL_DATA.exists():

    fail(
        f"Missing local-family data:\n{LOCAL_DATA}"
    )


if not GLOBAL_ARCHITECTURE.exists():

    fail(
        f"Missing global architecture table:\n{GLOBAL_ARCHITECTURE}"
    )


genes = pd.read_csv(
    LOCAL_DATA,
    sep="\t",
    dtype=str,
    keep_default_na=False
)


arch = pd.read_csv(
    GLOBAL_ARCHITECTURE,
    sep="\t",
    dtype=str,
    keep_default_na=False
)


## ================================================================== ##
## Required columns
## ================================================================== ##

required_genes = {

    "genome",
    "protein_id",
    "plot_region_id",

    "taxonomy_display",

    "plot_start_bp",
    "plot_end_bp",
    "plot_strand",

    "local_family_id",

    "neighbor_mmseq_cluster",

    "neighbor_globdb_cog",
    "neighbor_globdb_product",

    "neighbor_fegenie_HMMs",
    "neighbor_number_of_hemes",

    "display_class",

}


missing = (
    required_genes
    -
    set(
        genes.columns
    )
)


if missing:

    fail(
        "Missing columns from generalized local-family data:\n"
        +
        "\n".join(
            sorted(
                missing
            )
        )
    )


required_arch = {

    "plot_region_id",
    "genome",

    "n_module_clusters_10kb",
    "clusters_10kb",

    "architecture_pattern_id",

    "priority_taxon",

}


missing = (
    required_arch
    -
    set(
        arch.columns
    )
)


if missing:

    fail(
        "Missing columns from global architecture table:\n"
        +
        "\n".join(
            sorted(
                missing
            )
        )
    )


## ================================================================== ##
## Numeric positions
## ================================================================== ##

genes[
    "plot_start_bp"
] = numeric(
    genes[
        "plot_start_bp"
    ]
)


genes[
    "plot_end_bp"
] = numeric(
    genes[
        "plot_end_bp"
    ]
)


genes[
    "midpoint_bp"
] = (
    genes[
        "plot_start_bp"
    ]
    +
    genes[
        "plot_end_bp"
    ]
) / 2.0


## ================================================================== ##
## Basic QC
## ================================================================== ##

print(
    "=" * 110
)

print(
    "BASIC QC"
)

print(
    "=" * 110
)


print(
    f"Regions:        {genes['plot_region_id'].nunique()}"
)

print(
    f"Genomes:        {genes['genome'].nunique()}"
)

print(
    f"Unique genes:   "
    f"{genes[['genome','protein_id']].drop_duplicates().shape[0]}"
)

print(
    f"Local families: {genes['local_family_id'].nunique()}"
)


if genes[
    "plot_region_id"
].nunique() != 53:

    fail(
        "Expected 53 Module-10 focal regions."
    )


if (
    genes[
        [
            "plot_region_id",
            "protein_id"
        ]
    ].duplicated().any()
):

    fail(
        "Duplicate plot_region_id + protein_id observations found."
    )


## ================================================================== ##
## Local family catalogue
## ================================================================== ##

catalog_rows = []


for family, x in genes.groupby(
    "local_family_id"
):

    strands = [

        s
        for s in x[
            "plot_strand"
        ]
        if s in {
            "+",
            "-"
        }

    ]


    if strands:

        strand_counts = Counter(
            strands
        )

        dominant_strand, dominant_n = strand_counts.most_common(
            1
        )[0]

        strand_fraction = (
            dominant_n
            /
            len(
                strands
            )
        )

    else:

        dominant_strand = ""

        strand_fraction = np.nan


    catalog_rows.append(
        {

            "local_family_id":
                family,

            "n_proteins":
                x[
                    [
                        "genome",
                        "protein_id"
                    ]
                ].drop_duplicates().shape[0],

            "n_regions":
                x[
                    "plot_region_id"
                ].nunique(),

            "n_genomes":
                x[
                    "genome"
                ].nunique(),

            "median_midpoint_bp":
                x[
                    "midpoint_bp"
                ].median(),

            "min_midpoint_bp":
                x[
                    "midpoint_bp"
                ].min(),

            "max_midpoint_bp":
                x[
                    "midpoint_bp"
                ].max(),

            "position_span_bp":
                (
                    x[
                        "midpoint_bp"
                    ].max()
                    -
                    x[
                        "midpoint_bp"
                    ].min()
                ),

            "dominant_strand":
                dominant_strand,

            "strand_conservation_fraction":
                strand_fraction,

            "global_mmseq_clusters":
                unique_join(
                    x[
                        "neighbor_mmseq_cluster"
                    ]
                ),

            "display_classes":
                unique_join(
                    x[
                        "display_class"
                    ]
                ),

            "dominant_cog":
                dominant(
                    x[
                        "neighbor_globdb_cog"
                    ]
                ),

            "dominant_product":
                dominant(
                    x[
                        "neighbor_globdb_product"
                    ]
                ),

            "dominant_fegenie_HMM":
                dominant(
                    x[
                        "neighbor_fegenie_HMMs"
                    ]
                ),

            "heme_counts":
                unique_join(
                    x[
                        "neighbor_number_of_hemes"
                    ]
                ),

        }
    )


catalog = pd.DataFrame(
    catalog_rows
).sort_values(

    [
        "n_regions",
        "position_span_bp",
        "local_family_id",
    ],

    ascending=[
        False,
        True,
        True,
    ]

)


catalog.to_csv(
    OUT_CATALOG,
    sep="\t",
    index=False
)


## ================================================================== ##
## Map global Module-10 clusters -> local families
## ================================================================== ##

module_genes = genes.loc[
    genes[
        "neighbor_mmseq_cluster"
    ].isin(
        MODULE_CLUSTERS
    )
].copy()


map_rows = []


for (
    cluster,
    family
), x in module_genes.groupby(

    [
        "neighbor_mmseq_cluster",
        "local_family_id",
    ]

):


    map_rows.append(
        {

            "global_mmseq_cluster":
                cluster,

            "local_family_id":
                family,

            "n_proteins":
                x[
                    [
                        "genome",
                        "protein_id"
                    ]
                ].drop_duplicates().shape[0],

            "n_regions":
                x[
                    "plot_region_id"
                ].nunique(),

            "n_genomes":
                x[
                    "genome"
                ].nunique(),

            "median_midpoint_bp":
                x[
                    "midpoint_bp"
                ].median(),

            "position_span_bp":
                (
                    x[
                        "midpoint_bp"
                    ].max()
                    -
                    x[
                        "midpoint_bp"
                    ].min()
                ),

            "dominant_strand":
                dominant(
                    x[
                        "plot_strand"
                    ]
                ),

            "dominant_cog":
                dominant(
                    x[
                        "neighbor_globdb_cog"
                    ]
                ),

            "dominant_product":
                dominant(
                    x[
                        "neighbor_globdb_product"
                    ]
                ),

            "dominant_fegenie_HMM":
                dominant(
                    x[
                        "neighbor_fegenie_HMMs"
                    ]
                ),

            "heme_counts":
                unique_join(
                    x[
                        "neighbor_number_of_hemes"
                    ]
                ),

        }
    )


module_map = pd.DataFrame(
    map_rows
).sort_values(

    [
        "global_mmseq_cluster",
        "n_regions",
        "local_family_id",
    ],

    ascending=[
        True,
        False,
        True,
    ]

)


module_map.to_csv(
    OUT_MODULE_MAP,
    sep="\t",
    index=False
)


## ================================================================== ##
## Region-level local architecture
##
## Preserve ALL copies and order them by normalized genomic position.
## ================================================================== ##

region_rows = []


for region_id, x in genes.groupby(
    "plot_region_id"
):


    x = x.sort_values(
        [
            "midpoint_bp",
            "protein_id",
        ]
    )


    meta = arch.loc[
        arch[
            "plot_region_id"
        ]
        ==
        region_id
    ]


    if len(
        meta
    ) != 1:

        fail(
            f"Expected one global architecture row for {region_id}; "
            f"found {len(meta)}."
        )


    meta = meta.iloc[0]


    module_x = x.loc[
        x[
            "neighbor_mmseq_cluster"
        ].isin(
            MODULE_CLUSTERS
        )
    ].copy()


    module_elements = []


    for _, row in module_x.iterrows():

        module_elements.append(

            (
                f"{row['neighbor_mmseq_cluster']}"
                f":{row['local_family_id']}"
                f"@{row['midpoint_bp']:.1f}"
                f"{row['plot_strand']}"
            )

        )


    ## -------------------------------------------------------------- ##
    ## Sequence-family composition only:
    ##
    ## coordinates removed so loci with the same set/order of local
    ## sequence families group together despite small spacing changes.
    ## -------------------------------------------------------------- ##

    module_family_order = [

        (
            f"{row['neighbor_mmseq_cluster']}"
            f":{row['local_family_id']}"
        )

        for _, row in module_x.iterrows()

    ]


    ## -------------------------------------------------------------- ##
    ## Recurrent context:
    ##
    ## descriptive families occurring in >=2 focal regions.
    ##
    ## This is NOT a biological threshold for inclusion in the actual
    ## gene map; it is only a convenient descriptive summary.
    ## -------------------------------------------------------------- ##

    recurrent_families = set(

        catalog.loc[
            catalog[
                "n_regions"
            ]
            >=
            2,
            "local_family_id"
        ]

    )


    recurrent_x = x.loc[
        x[
            "local_family_id"
        ].isin(
            recurrent_families
        )
    ]


    recurrent_order = [

        row[
            "local_family_id"
        ]

        for _, row in recurrent_x.iterrows()

    ]


    region_rows.append(
        {

            "plot_region_id":
                region_id,

            "genome":
                meta[
                    "genome"
                ],

            "taxonomy_display":
                x[
                    "taxonomy_display"
                ].iloc[0],

            "priority_taxon":
                meta[
                    "priority_taxon"
                ],

            "n_module_clusters_10kb":
                int(
                    meta[
                        "n_module_clusters_10kb"
                    ]
                ),

            "global_architecture_pattern_id":
                meta[
                    "architecture_pattern_id"
                ],

            "global_clusters_10kb":
                meta[
                    "clusters_10kb"
                ],

            "n_module_gene_copies":
                len(
                    module_x
                ),

            "module_local_family_order":
                " > ".join(
                    module_family_order
                ),

            "module_local_family_position_signature":
                " > ".join(
                    module_elements
                ),

            "recurrent_local_family_order":
                " > ".join(
                    recurrent_order
                ),

            "n_recurrent_local_families":
                len(
                    recurrent_x
                ),

        }
    )


region_arch = pd.DataFrame(
    region_rows
)


## ================================================================== ##
## Define exact module-local-family patterns
## ================================================================== ##

pattern_lookup = (

    region_arch

    .groupby(
        [
            "global_architecture_pattern_id",
            "module_local_family_order",
        ],
        dropna=False
    )

    .agg(

        n_regions=(
            "plot_region_id",
            "nunique"
        ),

        n_genomes=(
            "genome",
            "nunique"
        ),

        n_methylobacter=(
            "priority_taxon",
            lambda x:
            int(
                (
                    x
                    ==
                    "Methylobacter"
                ).sum()
            )
        ),

        taxa=(
            "taxonomy_display",
            unique_join
        ),

    )

    .reset_index()

    .sort_values(
        [
            "n_regions",
            "global_architecture_pattern_id",
        ],
        ascending=[
            False,
            True,
        ]
    )

    .reset_index(
        drop=True
    )

)


pattern_lookup[
    "local_architecture_id"
] = [

    f"M10_arch_{i:02d}"

    for i in range(
        1,
        len(
            pattern_lookup
        )
        +
        1
    )

]


region_arch = region_arch.merge(

    pattern_lookup[
        [
            "global_architecture_pattern_id",
            "module_local_family_order",
            "local_architecture_id",
            "n_regions",
        ]
    ].rename(
        columns={
            "n_regions":
                "local_architecture_n_regions"
        }
    ),

    on=[
        "global_architecture_pattern_id",
        "module_local_family_order",
    ],

    how="left",

    validate="many_to_one"

)


region_arch.to_csv(
    OUT_REGION_ARCH,
    sep="\t",
    index=False
)


pattern_lookup.to_csv(
    OUT_PATTERN_TABLE,
    sep="\t",
    index=False
)


## ================================================================== ##
## Candidate table
## ================================================================== ##

candidate = region_arch.copy()


candidate[
    "priority_taxon_order"
] = np.where(

    candidate[
        "priority_taxon"
    ]
    ==
    "Methylobacter",

    0,

    1

)


candidate = candidate.sort_values(

    [
        "n_module_clusters_10kb",
        "priority_taxon_order",

        "local_architecture_n_regions",

        "taxonomy_display",
        "genome",
    ],

    ascending=[
        False,
        True,

        True,

        True,
        True,
    ]

).reset_index(
    drop=True
)


candidate[
    "candidate_order"
] = np.arange(
    1,
    len(
        candidate
    )
    +
    1
)


candidate.to_csv(
    OUT_CANDIDATES,
    sep="\t",
    index=False
)


## ================================================================== ##
## Console report
## ================================================================== ##

print()

print(
    "=" * 110
)

print(
    "GLOBAL MODULE CLUSTERS -> LOCAL SEQUENCE FAMILIES"
)

print(
    "=" * 110
)


print(
    module_map.to_string(
        index=False
    )
)


print()

print(
    "=" * 110
)

print(
    "RECURRENT LOCAL FAMILIES - >=2 REGIONS"
)

print(
    "=" * 110
)


print(

    catalog.loc[
        catalog[
            "n_regions"
        ]
        >=
        2,

        [
            "local_family_id",

            "n_regions",
            "n_genomes",

            "median_midpoint_bp",
            "position_span_bp",

            "dominant_strand",
            "strand_conservation_fraction",

            "global_mmseq_clusters",

            "dominant_cog",
            "dominant_product",
            "dominant_fegenie_HMM",
            "heme_counts",
        ]

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
    "LOCAL ARCHITECTURE PATTERNS WITHIN GLOBAL ARCHITECTURES"
)

print(
    "=" * 110
)


print(
    pattern_lookup.to_string(
        index=False
    )
)


print()

print(
    "=" * 110
)

print(
    "4/5 METHYLOBACTER LOCI"
)

print(
    "=" * 110
)


methylobacter_4 = candidate.loc[

    (
        candidate[
            "n_module_clusters_10kb"
        ]
        ==
        4
    )

    &

    (
        candidate[
            "priority_taxon"
        ]
        ==
        "Methylobacter"
    )

]


if methylobacter_4.empty:

    print(
        "[none]"
    )

else:

    print(

        methylobacter_4[
            [
                "genome",
                "taxonomy_display",

                "global_architecture_pattern_id",
                "local_architecture_id",
                "local_architecture_n_regions",

                "module_local_family_order",
                "recurrent_local_family_order",
            ]
        ].to_string(
            index=False
        )

    )


print()

print(
    "=" * 110
)

print(
    "ALL 4/5 LOCI"
)

print(
    "=" * 110
)


print(

    candidate.loc[
        candidate[
            "n_module_clusters_10kb"
        ]
        ==
        4,

        [
            "genome",
            "taxonomy_display",
            "priority_taxon",

            "global_architecture_pattern_id",

            "local_architecture_id",
            "local_architecture_n_regions",

            "module_local_family_order",
            "recurrent_local_family_order",
        ]

    ].to_string(
        index=False
    )

)


print()

print(
    "=" * 110
)

print(
    "3/5 LOCI"
)

print(
    "=" * 110
)


print(

    candidate.loc[
        candidate[
            "n_module_clusters_10kb"
        ]
        ==
        3,

        [
            "genome",
            "taxonomy_display",
            "priority_taxon",

            "global_architecture_pattern_id",

            "local_architecture_id",

            "module_local_family_order",
            "recurrent_local_family_order",
        ]

    ].to_string(
        index=False
    )

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


for path in [

    OUT_CATALOG,
    OUT_MODULE_MAP,
    OUT_REGION_ARCH,
    OUT_PATTERN_TABLE,
    OUT_CANDIDATES,

]:

    print(
        path
    )
