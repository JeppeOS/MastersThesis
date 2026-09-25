#!/usr/bin/env python3

## ================================================================== ##
## STAGE 16A3 - MODULE 20 REGION ARCHITECTURE MATRIX
##
## Purpose
## -------
## Convert the complete Module-20 local-family analysis into one
## biologically interpretable row per Cluster_00035 focal region.
##
## This script integrates:
##
##   - actual focal-region identity
##   - authoritative MtoA / MtrA / FeGenie-negative focal call
##   - MtrB_TIGR03509 recovered from integrated master annotations
##   - MtrB local sequence-family identity
##   - MtrB adjacency / orientation / physical spacing
##   - Cluster_00048
##   - Cluster_00166
##   - Cluster_00209
##   - recurrent local families
##   - taxonomy
##   - architecture signatures
##   - representative-selection priority metadata
##
## IMPORTANT:
##
## This script DOES NOT select the final <=12 regions.
##
## It produces the auditable candidate table from which those regions
## will be chosen.
##
##
## Agreed biological priority:
##
##   1. Methylobacter/Crenothrix + MtoA/MtrA + adjacent MtrB
##
##   2. Other taxa + MtoA/MtrA + adjacent MtrB
##
##   3. Methylobacter/Crenothrix without the canonical positive
##      MtoA/MtrA-MtrB architecture
##
##   plus:
##
##      reserve approximately two FeGenie-negative Cluster_00035
##      regions as explicit biological contrasts.
##
##
## Within priority groups, later representative selection will favour:
##
##   - architecture diversity
##   - MtrA representation
##   - MtrB local-family diversity
##   - taxonomic diversity
##   - non-redundant accessory architectures
##
## ================================================================== ##


from pathlib import Path
import re

import numpy as np
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


PLOT_FILE = (
    ROOT
    / "16_visualization"
    / "16A_module20_gene_map_framework"
    / "Module_20_gene_map_data.tsv"
)


LOCAL_RUN = (
    ROOT
    / "16_visualization"
    / "16A2_local_architecture_pipeline"
    / "runs"
    / "Module_20"
)


MASTER_FILE = (
    ROOT
    / "05_integrated_annotations"
    / "MASTER_PROTEIN_ANNOTATIONS.tsv"
)


OUTPUT_DIR = (
    ROOT
    / "16_visualization"
    / "16A3_module20_representative_selection"
)


OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)


OUT_MATRIX = (
    OUTPUT_DIR
    / "Module_20_region_architecture_matrix.tsv"
)


OUT_PATTERNS = (
    OUTPUT_DIR
    / "Module_20_architecture_patterns.tsv"
)


OUT_COUNTS = (
    OUTPUT_DIR
    / "Module_20_region_local_family_counts.tsv"
)


OUT_FAMILY_CATALOG = (
    OUTPUT_DIR
    / "Module_20_recurrent_local_family_catalog.tsv"
)


OUT_CANDIDATES = (
    OUTPUT_DIR
    / "Module_20_representative_candidate_table.tsv"
)


OUT_ENRICHED_GENES = (
    OUTPUT_DIR
    / "Module_20_gene_map_data_enriched.tsv"
)


OUT_QC = (
    OUTPUT_DIR
    / "Module_20_architecture_matrix_qc.tsv"
)


## ================================================================== ##
## Constants
## ================================================================== ##

EXPECTED_REGIONS = 74

EXPECTED_GENOMES = 71

EXPECTED_GENES = 965


MTRB_HMM = "MtrB_TIGR03509"


TAXA_OF_SPECIAL_INTEREST = {
    "Methylobacter",
    "Crenothrix",
}


## ------------------------------------------------------------------ ##
## Known global families
## ------------------------------------------------------------------ ##

GLOBAL_CLUSTER_FOCAL = "Cluster_00035"

GLOBAL_CLUSTER_PARTNER = "Cluster_00048"

GLOBAL_CLUSTER_166 = "Cluster_00166"

GLOBAL_CLUSTER_209 = "Cluster_00209"


## ------------------------------------------------------------------ ##
## Local families of current biological interest
##
## These are not assumed to be universally functional components.
## They are simply carried into the compact architecture matrix.
## ------------------------------------------------------------------ ##

LOCAL_FAMILY_CORE_PARTNER = "M20_local_002"

LOCAL_FAMILY_UNKNOWN_COMMON = "M20_local_003"

LOCAL_FAMILY_FOUR_HEME_NAPC = "M20_local_007"

LOCAL_FAMILY_CLUSTER166 = "M20_local_008"

LOCAL_FAMILY_CLUSTER209 = "M20_local_013"

LOCAL_FAMILY_MULTICOPPER = "M20_local_014"

LOCAL_FAMILY_PARALOG_BLOCK = "M20_local_015"


## ------------------------------------------------------------------ ##
## Minimum recurrence for family catalog / architecture metadata
##
## This is descriptive only.
##
## No family is removed from the underlying data.
## ------------------------------------------------------------------ ##

RECURRENT_MIN_REGIONS = 3


## ================================================================== ##
## Helpers
## ================================================================== ##

def fail(message):

    raise SystemExit(
        f"ERROR: {message}"
    )


def clean(value):

    if pd.isna(value):
        return ""

    return str(
        value
    ).strip()


def split_annotation(value):

    value = clean(
        value
    )

    if not value:
        return []

    return [
        token.strip()
        for token in re.split(
            r"[;,|]",
            value
        )
        if token.strip()
    ]


def contains_hmm(
    value,
    target
):

    return target in set(
        split_annotation(
            value
        )
    )


def dominant_nonempty(
    series
):

    x = (
        series
        .fillna(
            ""
        )
        .astype(
            str
        )
        .str.strip()
    )


    x = x[
        x != ""
    ]


    if x.empty:
        return ""


    return (
        x
        .value_counts()
        .index[
            0
        ]
    )


def joined_unique(
    series
):

    values = sorted(
        {
            clean(
                x
            )
            for x in series
            if clean(
                x
            )
        }
    )


    return ";".join(
        values
    )


def bool_int(
    value
):

    return int(
        bool(
            value
        )
    )


def classify_interest_taxon(
    value
):

    value = clean(
        value
    )


    if value.startswith(
        "Methylobacter"
    ):

        return "Methylobacter"


    if value == "Crenothrix":

        return "Crenothrix"


    return "Other"


## ================================================================== ##
## Locate local-family membership
## ================================================================== ##

membership_files = sorted(
    LOCAL_RUN.glob(
        "*local_family_membership*.tsv"
    )
)


if len(
    membership_files
) != 1:

    fail(
        "Expected exactly one Module-20 local-family membership table; "
        f"found {len(membership_files)}."
    )


MEMBERSHIP_FILE = membership_files[
    0
]


## ================================================================== ##
## Read input
## ================================================================== ##

print(
    "=" * 110
)

print(
    "STAGE 16A3 - MODULE 20 REGION ARCHITECTURE MATRIX"
)

print(
    "=" * 110
)


print()
print(
    "Reading Module-20 gene-map data..."
)


plot = pd.read_csv(
    PLOT_FILE,
    sep="\t",
    dtype=str,
    keep_default_na=False
)


print(
    f"  rows: "
    f"{len(plot):,}"
)


print()
print(
    "Reading local-family membership..."
)


membership = pd.read_csv(
    MEMBERSHIP_FILE,
    sep="\t",
    dtype=str,
    keep_default_na=False
)


print(
    f"  rows: "
    f"{len(membership):,}"
)


print()
print(
    "Reading integrated annotation master..."
)


master = pd.read_csv(
    MASTER_FILE,
    sep="\t",
    dtype=str,
    keep_default_na=False
)


## ================================================================== ##
## Validate plot structure
## ================================================================== ##

required_plot = {

    "genome",
    "protein_id",

    "plot_region_id",
    "focal_protein_id",
    "focal_call",

    "plot_is_focal",

    "gene_midpoint",
    "oriented_gene_offset",
    "distance_to_focal_cds_bp",
    "plot_strand",

    "neighbor_mmseq_cluster",

    "neighbor_fegenie_HMMs",
    "neighbor_number_of_hemes",

    "neighbor_globdb_cog",
    "neighbor_globdb_product",

    "taxonomy_family",
    "taxonomy_genus",
    "taxonomy_species",
    "taxonomy_display",
    "priority_taxon",

}


missing = (
    required_plot
    -
    set(
        plot.columns
    )
)


if missing:

    fail(
        "Plot table missing columns: "
        +
        ", ".join(
            sorted(
                missing
            )
        )
    )


if len(
    plot
) != EXPECTED_GENES:

    fail(
        f"Expected {EXPECTED_GENES} Module-20 genes; "
        f"found {len(plot)}."
    )


if (
    plot[
        "plot_region_id"
    ].nunique()
    !=
    EXPECTED_REGIONS
):

    fail(
        f"Expected {EXPECTED_REGIONS} focal regions."
    )


if (
    plot[
        "genome"
    ].nunique()
    !=
    EXPECTED_GENOMES
):

    fail(
        f"Expected {EXPECTED_GENOMES} focal genomes."
    )


## ================================================================== ##
## Validate membership
## ================================================================== ##

required_membership = {

    "genome",
    "protein_id",
    "local_family_id",

}


missing = (
    required_membership
    -
    set(
        membership.columns
    )
)


if missing:

    fail(
        "Membership table missing columns: "
        +
        ", ".join(
            sorted(
                missing
            )
        )
    )


membership_map = (
    membership[
        [
            "genome",
            "protein_id",
            "local_family_id"
        ]
    ]
    .drop_duplicates()
)


if membership_map.duplicated(
    [
        "genome",
        "protein_id"
    ]
).any():

    fail(
        "A protein belongs to more than one local family."
    )


## ================================================================== ##
## Attach local-family identity to actual genes
## ================================================================== ##

genes = plot.merge(
    membership_map,
    on=[
        "genome",
        "protein_id"
    ],
    how="left",
    validate="one_to_one"
)


if genes[
    "local_family_id"
].isna().any():

    fail(
        "Some Module-20 genes lack local-family assignments."
    )


## ================================================================== ##
## Supplement missing FeGenie annotations from integrated master
##
## Important:
##
## Many ordinary neighborhood genes do not occur in the integrated
## candidate master. That is expected.
##
## We use the master only as an annotation supplement.
## ================================================================== ##

master_required = {

    "genome",
    "protein_id",
    "fegenie_HMMs",

}


missing = (
    master_required
    -
    set(
        master.columns
    )
)


if missing:

    fail(
        "Integrated master missing columns: "
        +
        ", ".join(
            sorted(
                missing
            )
        )
    )


master_keep_columns = [
    c
    for c in [
        "genome",
        "protein_id",

        "fegenie_positive",
        "fegenie_HMMs",
        "fegenie_categories",

        "findmehemes_positive",
        "number_of_hemes",

        "signalp_prediction",
        "deeptmhmm_class",
        "report_topology",

        "mmseq_cluster",
        "mcl_module",

        "globdb_cog",
        "globdb_gene",
        "globdb_product",
    ]
    if c in master.columns
]


master_keep = master[
    master_keep_columns
].copy()


if master_keep.duplicated(
    [
        "genome",
        "protein_id"
    ]
).any():

    fail(
        "Integrated master contains duplicate genome + protein_id keys."
    )


master_keep[
    "integrated_master_match"
] = 1


master_keep = master_keep.rename(
    columns={
        c:
            f"master_{c}"

        for c in master_keep.columns

        if c not in {
            "genome",
            "protein_id"
        }
    }
)


genes = genes.merge(
    master_keep,
    on=[
        "genome",
        "protein_id"
    ],
    how="left",
    validate="one_to_one"
)


genes[
    "integrated_master_match"
] = (
    genes[
        "master_integrated_master_match"
    ]
    .fillna(
        0
    )
    .astype(
        int
    )
)


## ================================================================== ##
## Effective FeGenie annotation
##
## Prefer integrated master when present.
## Fall back to Stage-15A neighborhood annotation otherwise.
## ================================================================== ##

master_hmm = (
    genes[
        "master_fegenie_HMMs"
    ]
    .fillna(
        ""
    )
)


stage15_hmm = (
    genes[
        "neighbor_fegenie_HMMs"
    ]
    .fillna(
        ""
    )
)


genes[
    "effective_fegenie_HMMs"
] = np.where(
    master_hmm.ne(
        ""
    ),
    master_hmm,
    stage15_hmm
)


genes[
    "is_mtrb"
] = genes[
    "effective_fegenie_HMMs"
].map(
    lambda value:
        int(
            contains_hmm(
                value,
                MTRB_HMM
            )
        )
)


## ================================================================== ##
## Numeric architecture fields
## ================================================================== ##

for column in [
    "gene_midpoint",
    "oriented_gene_offset",
    "distance_to_focal_cds_bp",
]:

    genes[
        column
    ] = pd.to_numeric(
        genes[
            column
        ],
        errors="coerce"
    )


## ================================================================== ##
## Write enriched gene-map table
##
## This will later become the preferred input for the final Module-20
## plot because MtrB annotation has now been restored.
## ================================================================== ##

genes.to_csv(
    OUT_ENRICHED_GENES,
    sep="\t",
    index=False
)


## ================================================================== ##
## Family prevalence
## ================================================================== ##

family_prevalence = (
    genes
    .groupby(
        "local_family_id"
    )
    .agg(

        n_proteins=(
            "protein_id",
            "size"
        ),

        n_regions=(
            "plot_region_id",
            "nunique"
        ),

        n_genomes=(
            "genome",
            "nunique"
        ),

    )
    .reset_index()
)


family_prevalence_map = (
    family_prevalence
    .set_index(
        "local_family_id"
    )[
        "n_regions"
    ]
    .to_dict()
)


recurrent_families = set(
    family_prevalence.loc[
        family_prevalence[
            "n_regions"
        ]
        >=
        RECURRENT_MIN_REGIONS,
        "local_family_id"
    ]
)


## ================================================================== ##
## Recurrent family annotation catalog
## ================================================================== ##

family_catalog_rows = []


for family, group in genes.groupby(
    "local_family_id",
    sort=True
):

    n_regions = group[
        "plot_region_id"
    ].nunique()


    if n_regions < RECURRENT_MIN_REGIONS:
        continue


    mids = group[
        "gene_midpoint"
    ].dropna()


    strands = (
        group[
            "plot_strand"
        ]
        .replace(
            "",
            pd.NA
        )
        .dropna()
    )


    if len(
        strands
    ):

        dominant_strand = (
            strands
            .value_counts()
            .index[
                0
            ]
        )

        strand_conservation = (
            strands
            .eq(
                dominant_strand
            )
            .mean()
        )

    else:

        dominant_strand = ""

        strand_conservation = np.nan


    family_catalog_rows.append(
        {

            "local_family_id":
                family,

            "n_proteins":
                len(
                    group
                ),

            "n_regions":
                n_regions,

            "n_genomes":
                group[
                    "genome"
                ].nunique(),

            "median_midpoint_bp":
                mids.median()
                if len(
                    mids
                )
                else np.nan,

            "position_span_bp":
                (
                    mids.max()
                    -
                    mids.min()
                )
                if len(
                    mids
                )
                else np.nan,

            "dominant_strand":
                dominant_strand,

            "strand_conservation_fraction":
                strand_conservation,

            "n_mtrb":
                int(
                    group[
                        "is_mtrb"
                    ].sum()
                ),

            "dominant_global_cluster":
                dominant_nonempty(
                    group[
                        "neighbor_mmseq_cluster"
                    ]
                ),

            "dominant_fegenie_HMM":
                dominant_nonempty(
                    group[
                        "effective_fegenie_HMMs"
                    ]
                ),

            "dominant_heme_count":
                dominant_nonempty(
                    group[
                        "neighbor_number_of_hemes"
                    ]
                ),

            "dominant_cog":
                dominant_nonempty(
                    group[
                        "neighbor_globdb_cog"
                    ]
                ),

            "dominant_product":
                dominant_nonempty(
                    group[
                        "neighbor_globdb_product"
                    ]
                ),

        }
    )


family_catalog = pd.DataFrame(
    family_catalog_rows
).sort_values(
    [
        "n_regions",
        "local_family_id"
    ],
    ascending=[
        False,
        True
    ],
    kind="stable"
)


family_catalog.to_csv(
    OUT_FAMILY_CATALOG,
    sep="\t",
    index=False
)


## ================================================================== ##
## Full region x local-family count matrix
##
## Counts are retained rather than converting immediately to binary,
## because paralogues can be biologically informative.
## ================================================================== ##

family_counts = (
    genes
    .groupby(
        [
            "plot_region_id",
            "local_family_id"
        ]
    )
    .size()
    .unstack(
        fill_value=0
    )
)


family_counts = family_counts.reindex(
    sorted(
        family_counts.columns
    ),
    axis=1
)


family_counts = (
    family_counts
    .reset_index()
)


family_counts.to_csv(
    OUT_COUNTS,
    sep="\t",
    index=False
)


## ================================================================== ##
## Base one-row-per-region table
## ================================================================== ##

focal_rows = genes[
    genes[
        "plot_is_focal"
    ]
    .astype(
        str
    )
    ==
    "1"
].copy()


if (
    focal_rows[
        "plot_region_id"
    ].duplicated().any()
):

    fail(
        "More than one focal gene is present in at least one region."
    )


if len(
    focal_rows
) != EXPECTED_REGIONS:

    fail(
        f"Expected {EXPECTED_REGIONS} focal rows; "
        f"found {len(focal_rows)}."
    )


region_base = focal_rows[
    [
        "plot_region_id",

        "genome",
        "focal_protein_id",
        "focal_call",

        "taxonomy_family",
        "taxonomy_genus",
        "taxonomy_species",
        "taxonomy_display",
        "priority_taxon",
    ]
].copy()


region_base[
    "interest_taxon"
] = region_base[
    "taxonomy_genus"
].map(
    classify_interest_taxon
)


region_base[
    "is_interest_taxon"
] = region_base[
    "interest_taxon"
].isin(
    TAXA_OF_SPECIAL_INTEREST
).astype(
    int
)


## ================================================================== ##
## Region-level helper functions
## ================================================================== ##

def family_present(
    group,
    family
):

    return int(
        (
            group[
                "local_family_id"
            ]
            ==
            family
        ).any()
    )


def cluster_present(
    group,
    cluster
):

    return int(
        (
            group[
                "neighbor_mmseq_cluster"
            ]
            ==
            cluster
        ).any()
    )


def family_count(
    group,
    family
):

    return int(
        (
            group[
                "local_family_id"
            ]
            ==
            family
        ).sum()
    )


## ================================================================== ##
## Build architecture row for each focal region
## ================================================================== ##

architecture_rows = []


for region_id, group in genes.groupby(
    "plot_region_id",
    sort=True
):

    mtrb = group[
        group[
            "is_mtrb"
        ]
        ==
        1
    ].copy()


    ## -------------------------------------------------------------- ##
    ## MtrB architecture
    ## -------------------------------------------------------------- ##

    if mtrb.empty:

        n_mtrb = 0

        mtrb_families = ""

        mtrb_offsets = ""

        mtrb_median_midpoint = np.nan

        mtrb_min_gap = np.nan

        has_adjacent_mtrb = 0

        all_mtrb_downstream = 0

        all_mtrb_codirectional = 0


    else:

        n_mtrb = len(
            mtrb
        )


        mtrb_families = joined_unique(
            mtrb[
                "local_family_id"
            ]
        )


        offsets = (
            mtrb[
                "oriented_gene_offset"
            ]
            .dropna()
            .astype(
                int
            )
        )


        mtrb_offsets = ";".join(
            str(
                x
            )
            for x in sorted(
                offsets.tolist()
            )
        )


        mtrb_median_midpoint = float(
            mtrb[
                "gene_midpoint"
            ]
            .median()
        )


        mtrb_min_gap = float(
            mtrb[
                "distance_to_focal_cds_bp"
            ]
            .min()
        )


        has_adjacent_mtrb = int(
            offsets.abs().eq(
                1
            ).any()
        )


        all_mtrb_downstream = int(
            (
                mtrb[
                    "gene_midpoint"
                ]
                >
                0
            ).all()
        )


        all_mtrb_codirectional = int(
            (
                mtrb[
                    "plot_strand"
                ]
                ==
                "+"
            ).all()
        )


    ## -------------------------------------------------------------- ##
    ## Recurrent family string
    ##
    ## Retain only families present in >=3 regions for compact
    ## descriptive architecture.
    ## -------------------------------------------------------------- ##

    recurrent_here = sorted(
        set(
            group.loc[
                group[
                    "local_family_id"
                ].isin(
                    recurrent_families
                ),
                "local_family_id"
            ]
        )
    )


    recurrent_family_string = ";".join(
        recurrent_here
    )


    architecture_rows.append(
        {

            "plot_region_id":
                region_id,


            ## ------------------------------------------------------ ##
            ## MtrB
            ## ------------------------------------------------------ ##

            "n_mtrb":
                n_mtrb,

            "has_mtrb":
                int(
                    n_mtrb
                    >
                    0
                ),

            "has_adjacent_mtrb":
                has_adjacent_mtrb,

            "mtrb_local_families":
                mtrb_families,

            "mtrb_oriented_gene_offsets":
                mtrb_offsets,

            "mtrb_median_midpoint_bp":
                mtrb_median_midpoint,

            "mtrb_min_intergenic_gap_bp":
                mtrb_min_gap,

            "mtrb_all_downstream":
                all_mtrb_downstream,

            "mtrb_all_codirectional":
                all_mtrb_codirectional,


            ## ------------------------------------------------------ ##
            ## Global protein families
            ## ------------------------------------------------------ ##

            "has_cluster00048":
                cluster_present(
                    group,
                    GLOBAL_CLUSTER_PARTNER
                ),

            "has_cluster00166":
                cluster_present(
                    group,
                    GLOBAL_CLUSTER_166
                ),

            "has_cluster00209":
                cluster_present(
                    group,
                    GLOBAL_CLUSTER_209
                ),


            ## ------------------------------------------------------ ##
            ## Selected local architecture families
            ## ------------------------------------------------------ ##

            "has_local003":
                family_present(
                    group,
                    LOCAL_FAMILY_UNKNOWN_COMMON
                ),

            "has_local007_four_heme_napc":
                family_present(
                    group,
                    LOCAL_FAMILY_FOUR_HEME_NAPC
                ),

            "has_local008_cluster166_like":
                family_present(
                    group,
                    LOCAL_FAMILY_CLUSTER166
                ),

            "has_local013_cluster209":
                family_present(
                    group,
                    LOCAL_FAMILY_CLUSTER209
                ),

            "has_local014_multicopper":
                family_present(
                    group,
                    LOCAL_FAMILY_MULTICOPPER
                ),

            "n_local015":
                family_count(
                    group,
                    LOCAL_FAMILY_PARALOG_BLOCK
                ),


            ## ------------------------------------------------------ ##
            ## Complete recurrent-family architecture
            ## ------------------------------------------------------ ##

            "recurrent_local_families":
                recurrent_family_string,

            "n_recurrent_local_families":
                len(
                    recurrent_here
                ),

        }
    )


architecture = pd.DataFrame(
    architecture_rows
)


## ================================================================== ##
## Combine with taxonomy/focal metadata
## ================================================================== ##

matrix = region_base.merge(
    architecture,
    on="plot_region_id",
    how="left",
    validate="one_to_one"
)


## ================================================================== ##
## Architecture signatures
##
## CORE:
##
##     intentionally ignores exact MtrB sequence subfamily.
##
## DETAILED:
##
##     retains MtrB subfamily and selected recurrent accessory families.
##
## ================================================================== ##

matrix[
    "core_architecture_signature"
] = (

    "focal="
    +
    matrix[
        "focal_call"
    ]

    +
    "|MtrB="
    +
    matrix[
        "has_adjacent_mtrb"
    ].astype(
        str
    )

    +
    "|C00048="
    +
    matrix[
        "has_cluster00048"
    ].astype(
        str
    )

    +
    "|C00166="
    +
    matrix[
        "has_cluster00166"
    ].astype(
        str
    )

    +
    "|C00209="
    +
    matrix[
        "has_cluster00209"
    ].astype(
        str
    )

)


matrix[
    "detailed_architecture_signature"
] = (

    matrix[
        "core_architecture_signature"
    ]

    +
    "|MtrBfam="
    +
    matrix[
        "mtrb_local_families"
    ].replace(
        "",
        "none"
    )

    +
    "|L003="
    +
    matrix[
        "has_local003"
    ].astype(
        str
    )

    +
    "|L007="
    +
    matrix[
        "has_local007_four_heme_napc"
    ].astype(
        str
    )

    +
    "|L014="
    +
    matrix[
        "has_local014_multicopper"
    ].astype(
        str
    )

    +
    "|L015n="
    +
    matrix[
        "n_local015"
    ].astype(
        str
    )

)


## ================================================================== ##
## Pattern prevalence
## ================================================================== ##

core_counts = (
    matrix[
        "core_architecture_signature"
    ]
    .value_counts()
    .to_dict()
)


detailed_counts = (
    matrix[
        "detailed_architecture_signature"
    ]
    .value_counts()
    .to_dict()
)


matrix[
    "core_architecture_n_regions"
] = matrix[
    "core_architecture_signature"
].map(
    core_counts
)


matrix[
    "detailed_architecture_n_regions"
] = matrix[
    "detailed_architecture_signature"
].map(
    detailed_counts
)


## ================================================================== ##
## Selection-priority metadata
##
## IMPORTANT:
##
## This is a category, NOT a final score.
##
## The final representative set must remain auditable and should not
## be determined by a single opaque numeric ranking.
## ================================================================== ##

matrix[
    "is_positive_focal"
] = matrix[
    "focal_call"
].isin(
    [
        "MtoA",
        "MtrA"
    ]
).astype(
    int
)


matrix[
    "has_canonical_mtrb_architecture"
] = (
    (
        matrix[
            "is_positive_focal"
        ]
        ==
        1
    )
    &
    (
        matrix[
            "has_adjacent_mtrb"
        ]
        ==
        1
    )
).astype(
    int
)


def priority_group(
    row
):

    special_taxon = (
        int(
            row[
                "is_interest_taxon"
            ]
        )
        ==
        1
    )


    canonical = (
        int(
            row[
                "has_canonical_mtrb_architecture"
            ]
        )
        ==
        1
    )


    if (
        special_taxon
        and
        canonical
    ):

        return "P1_taxon_plus_MtoA_MtrA_MtrB"


    if canonical:

        return "P2_MtoA_MtrA_MtrB"


    if special_taxon:

        return "P3_taxon_without_canonical_MtrB"


    return "P4_other"


matrix[
    "selection_priority_group"
] = matrix.apply(
    priority_group,
    axis=1
)


priority_order = {

    "P1_taxon_plus_MtoA_MtrA_MtrB":
        1,

    "P2_MtoA_MtrA_MtrB":
        2,

    "P3_taxon_without_canonical_MtrB":
        3,

    "P4_other":
        4,

}


matrix[
    "selection_priority_rank"
] = matrix[
    "selection_priority_group"
].map(
    priority_order
)


## ================================================================== ##
## Explicit contrast / diversity flags
##
## These flags will guide the final <=12 selection.
## ================================================================== ##

matrix[
    "negative_contrast_candidate"
] = (
    matrix[
        "focal_call"
    ]
    ==
    "FeGenie_negative"
).astype(
    int
)


matrix[
    "mtra_candidate"
] = (
    matrix[
        "focal_call"
    ]
    ==
    "MtrA"
).astype(
    int
)


matrix[
    "rare_detailed_architecture"
] = (
    matrix[
        "detailed_architecture_n_regions"
    ]
    <=
    2
).astype(
    int
)


## ================================================================== ##
## Architecture pattern summary
## ================================================================== ##

patterns = (
    matrix
    .groupby(
        [
            "core_architecture_signature",
            "detailed_architecture_signature"
        ],
        dropna=False
    )
    .agg(

        n_regions=(
            "plot_region_id",
            "size"
        ),

        n_genomes=(
            "genome",
            "nunique"
        ),

        n_Methylobacter=(
            "interest_taxon",
            lambda x:
                int(
                    (
                        x
                        ==
                        "Methylobacter"
                    ).sum()
                )
        ),

        n_Crenothrix=(
            "interest_taxon",
            lambda x:
                int(
                    (
                        x
                        ==
                        "Crenothrix"
                    ).sum()
                )
        ),

        n_MtoA=(
            "focal_call",
            lambda x:
                int(
                    (
                        x
                        ==
                        "MtoA"
                    ).sum()
                )
        ),

        n_MtrA=(
            "focal_call",
            lambda x:
                int(
                    (
                        x
                        ==
                        "MtrA"
                    ).sum()
                )
        ),

        n_FeGenie_negative=(
            "focal_call",
            lambda x:
                int(
                    (
                        x
                        ==
                        "FeGenie_negative"
                    ).sum()
                )
        ),

        genomes=(
            "genome",
            lambda x:
                ";".join(
                    sorted(
                        set(
                            x
                        )
                    )
                )
        ),

    )
    .reset_index()
    .sort_values(
        [
            "n_regions",
            "core_architecture_signature"
        ],
        ascending=[
            False,
            True
        ],
        kind="stable"
    )
)


patterns.to_csv(
    OUT_PATTERNS,
    sep="\t",
    index=False
)


## ================================================================== ##
## Candidate table
##
## Sort only for convenience.
##
## Do NOT interpret this as automatic final selection.
## ================================================================== ##

candidate_table = matrix.sort_values(

    [
        "selection_priority_rank",

        ## Ensure MtrA is visible near top of relevant groups.
        "mtra_candidate",

        ## Prefer architecture diversity within equal priority.
        "detailed_architecture_n_regions",

        "interest_taxon",
        "taxonomy_genus",
        "taxonomy_species",

        "genome",
        "focal_protein_id",
    ],

    ascending=[
        True,
        False,
        True,
        True,
        True,
        True,
        True,
        True,
    ],

    kind="stable"

)


candidate_table.to_csv(
    OUT_CANDIDATES,
    sep="\t",
    index=False
)


matrix.to_csv(
    OUT_MATRIX,
    sep="\t",
    index=False
)


## ================================================================== ##
## QC
## ================================================================== ##

qc = pd.DataFrame(
    [
        {
            "metric":
                "plot_gene_rows",
            "value":
                len(
                    genes
                )
        },

        {
            "metric":
                "focal_regions",
            "value":
                matrix[
                    "plot_region_id"
                ].nunique()
        },

        {
            "metric":
                "focal_genomes",
            "value":
                matrix[
                    "genome"
                ].nunique()
        },

        {
            "metric":
                "integrated_master_matches",
            "value":
                int(
                    genes[
                        "integrated_master_match"
                    ].sum()
                )
        },

        {
            "metric":
                "MtrB_genes",
            "value":
                int(
                    genes[
                        "is_mtrb"
                    ].sum()
                )
        },

        {
            "metric":
                "regions_with_MtrB",
            "value":
                int(
                    matrix[
                        "has_mtrb"
                    ].sum()
                )
        },

        {
            "metric":
                "regions_with_adjacent_MtrB",
            "value":
                int(
                    matrix[
                        "has_adjacent_mtrb"
                    ].sum()
                )
        },

        {
            "metric":
                "regions_with_cluster00048",
            "value":
                int(
                    matrix[
                        "has_cluster00048"
                    ].sum()
                )
        },

        {
            "metric":
                "regions_with_cluster00166",
            "value":
                int(
                    matrix[
                        "has_cluster00166"
                    ].sum()
                )
        },

        {
            "metric":
                "regions_with_cluster00209",
            "value":
                int(
                    matrix[
                        "has_cluster00209"
                    ].sum()
                )
        },

        {
            "metric":
                "core_architecture_patterns",
            "value":
                matrix[
                    "core_architecture_signature"
                ].nunique()
        },

        {
            "metric":
                "detailed_architecture_patterns",
            "value":
                matrix[
                    "detailed_architecture_signature"
                ].nunique()
        },

        {
            "metric":
                "negative_contrast_candidates",
            "value":
                int(
                    matrix[
                        "negative_contrast_candidate"
                    ].sum()
                )
        },

        {
            "metric":
                "MtrA_candidates",
            "value":
                int(
                    matrix[
                        "mtra_candidate"
                    ].sum()
                )
        },

    ]
)


qc.to_csv(
    OUT_QC,
    sep="\t",
    index=False
)


## ================================================================== ##
## Console interpretation summary
## ================================================================== ##

print()
print(
    "=" * 110
)

print(
    "MODULE 20 ARCHITECTURE SUMMARY"
)

print(
    "=" * 110
)


print(
    f"Focal regions:          "
    f"{len(matrix)}"
)


print(
    f"Genomes:                "
    f"{matrix['genome'].nunique()}"
)


print(
    f"Local families:         "
    f"{genes['local_family_id'].nunique()}"
)


print(
    f"Recurrent families >=3: "
    f"{len(recurrent_families)}"
)


print()
print(
    "Focal calls:"
)


print(
    matrix[
        "focal_call"
    ]
    .value_counts()
    .to_string()
)


print()
print(
    "MtrB:"
)


print(
    f"  regions containing MtrB: "
    f"{matrix['has_mtrb'].sum()}/"
    f"{len(matrix)}"
)


print(
    f"  adjacent MtrB:           "
    f"{matrix['has_adjacent_mtrb'].sum()}/"
    f"{len(matrix)}"
)


print()
print(
    "Key family combinations:"
)


print(
    f"  Cluster_00048: "
    f"{matrix['has_cluster00048'].sum()}"
)


print(
    f"  Cluster_00166: "
    f"{matrix['has_cluster00166'].sum()}"
)


print(
    f"  Cluster_00209: "
    f"{matrix['has_cluster00209'].sum()}"
)


triple = (
    (
        matrix[
            "has_cluster00048"
        ]
        ==
        1
    )
    &
    (
        matrix[
            "has_adjacent_mtrb"
        ]
        ==
        1
    )
)


print()
print(
    "Cluster_00048 -> focal -> MtrB candidate architecture:"
)


print(
    f"  {triple.sum()}/{len(matrix)} regions"
)


print()
print(
    "Selection priority groups:"
)


print(
    matrix[
        "selection_priority_group"
    ]
    .value_counts()
    .to_string()
)


print()
print(
    "Taxon x focal/MtrB status:"
)


print(
    pd.crosstab(

        matrix[
            "interest_taxon"
        ],

        [
            matrix[
                "focal_call"
            ],

            matrix[
                "has_adjacent_mtrb"
            ],
        ]

    ).to_string()
)


print()
print(
    "=" * 110
)

print(
    "MtrA REGIONS"
)

print(
    "=" * 110
)


mtra_show = [
    "genome",
    "taxonomy_display",
    "focal_protein_id",

    "mtrb_local_families",

    "has_cluster00048",
    "has_cluster00166",
    "has_cluster00209",

    "core_architecture_signature",
    "detailed_architecture_signature",
]


print(
    matrix[
        matrix[
            "focal_call"
        ]
        ==
        "MtrA"
    ][
        mtra_show
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
    "FEGENIE-NEGATIVE CONTRAST REGIONS"
)

print(
    "=" * 110
)


negative_show = [
    "genome",
    "taxonomy_display",
    "focal_protein_id",

    "interest_taxon",

    "has_adjacent_mtrb",
    "has_cluster00048",
    "has_cluster00166",
    "has_cluster00209",

    "recurrent_local_families",

    "detailed_architecture_n_regions",
]


print(
    matrix[
        matrix[
            "negative_contrast_candidate"
        ]
        ==
        1
    ][
        negative_show
    ]
    .sort_values(
        [
            "interest_taxon",
            "detailed_architecture_n_regions",
            "genome",
        ]
    )
    .to_string(
        index=False
    )
)


print()
print(
    "=" * 110
)

print(
    "MOST COMMON DETAILED ARCHITECTURES"
)

print(
    "=" * 110
)


architecture_counts = (
    matrix[
        "detailed_architecture_signature"
    ]
    .value_counts()
    .rename_axis(
        "architecture"
    )
    .reset_index(
        name="n_regions"
    )
)


print(
    architecture_counts
    .head(
        20
    )
    .to_string(
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
    OUT_MATRIX,
    OUT_PATTERNS,
    OUT_COUNTS,
    OUT_FAMILY_CATALOG,
    OUT_CANDIDATES,
    OUT_ENRICHED_GENES,
    OUT_QC,
]:

    print(
        path
    )


print()
print(
    "SUCCESS"
)
