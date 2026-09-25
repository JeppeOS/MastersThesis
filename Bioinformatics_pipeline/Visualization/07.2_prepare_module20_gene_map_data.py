#!/usr/bin/env python3

## ================================================================== ##
## STAGE 16A - MODULE 20 COMPLETE GENE-MAP FRAMEWORK
##
## Purpose
## -------
## Construct complete, focal-gene-centered neighborhoods for every
## Cluster_00035 protein belonging to Module 20.
##
## Biological observational unit:
##
##     genome + actual Cluster_00035 protein
##
## NOT:
##
##     Stage-15A focal catalogue row
##
## Stage 15A may contain the same actual protein as a focal in several
## module-associated contexts. Those duplicate contexts are collapsed
## here to the actual biological protein/locus.
##
##
## This script:
##
##   1. identifies the 74 unique Cluster_00035 proteins
##
##   2. resolves their canonical Module_20 gene-centered neighborhoods
##
##   3. retains complete Stage-15A context for architecture summaries
##
##   4. identifies MtrB directly using the authoritative FeGenie HMM:
##
##          MtrB_TIGR03509
##
##   5. calculates MtrB proximity at <=5, <=10 and <=20 kb BEFORE
##      cropping the plotting/local-clustering window
##
##   6. records MtrB gene offset, adjacency, order, orientation and
##      physical intergenic gap
##
##   7. normalizes all focal regions so Cluster_00035 points left->right
##
##   8. retains every actual CDS intersecting +/-10 kb for downstream
##      local-family clustering and visualization
##
##   9. attaches authoritative GlobDB r226 taxonomy
##
##  10. preserves whether each Cluster_00035 protein is called:
##
##          MtoA
##          MtrA
##          FeGenie_negative
##
##  11. tracks Cluster_00048 and known recurrent global families
##
##
## No representative genomes are selected here.
##
## Selection of <=12 representative regions occurs only AFTER:
##
##     local-family clustering
##     exhaustive architecture profiling
##     MtrB architecture characterization
##
## ================================================================== ##


from pathlib import Path
import re
import sys

import numpy as np
import pandas as pd


## ================================================================== ##
## Paths
## ================================================================== ##

WORKFLOW = Path(
    "~/methanotrophs/methanotroph_project/jeppe/"
    "genome_analysis_workflow"
).expanduser()


STAGE15A = (
    WORKFLOW
    / "15_gene_level_analysis"
    / "15A_resolution_recovery"
)


FOCAL_FILE = (
    STAGE15A
    / "focal_gene_catalog.tsv"
)


NEIGHBORHOOD_FILE = (
    STAGE15A
    / "focal_gene_neighborhoods.tsv"
)


TAXONOMY_FILE = (
    WORKFLOW
    / "16_visualization"
    / "globdb_r226_taxonomy.tsv"
)


OUTPUT_DIR = (
    WORKFLOW
    / "16_visualization"
    / "16A_module20_gene_map_framework"
)


OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)


OUT_DATA = (
    OUTPUT_DIR
    / "Module_20_gene_map_data.tsv"
)


OUT_REGIONS = (
    OUTPUT_DIR
    / "Module_20_focal_regions.tsv"
)


OUT_QC = (
    OUTPUT_DIR
    / "Module_20_gene_map_qc.tsv"
)


## ================================================================== ##
## Constants
## ================================================================== ##

MODULE = "Module_20"

FOCAL_CLUSTER = "Cluster_00035"


## ------------------------------------------------------------------ ##
## Established Module-20 co-occurring family
## ------------------------------------------------------------------ ##

MODULE_PARTNER = "Cluster_00048"


## ------------------------------------------------------------------ ##
## Previously observed recurrent accessory families
##
## These are retained as metadata only.
##
## The new local architecture pipeline will determine independently
## whether they belong to conserved subarchitectures.
## ------------------------------------------------------------------ ##

KNOWN_ACCESSORY_FAMILIES = {
    "Cluster_00166",
    "Cluster_00209",
}


## ------------------------------------------------------------------ ##
## Authoritative FeGenie MtrB HMM
## ------------------------------------------------------------------ ##

MTRB_FEGENIE_HMM = "MtrB_TIGR03509"


## ------------------------------------------------------------------ ##
## Plot / local-family window
## ------------------------------------------------------------------ ##

PLOT_WINDOW_BP = 10_000


## ------------------------------------------------------------------ ##
## Regression expectations from validated upstream analysis
## ------------------------------------------------------------------ ##

EXPECTED_UNIQUE_FOCALS = 74

EXPECTED_FOCAL_GENOMES = 71


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

    return str(value).strip()


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


def unique_nonempty(series):

    return sorted(
        {
            clean(x)
            for x in series
            if clean(x)
        }
    )


def first_nonempty(series):

    for value in series:

        value = clean(
            value
        )

        if value:
            return value

    return ""


def numeric_flag(value):

    return clean(
        value
    ).lower() in {
        "1",
        "true",
        "t",
        "yes",
    }


def parse_taxonomic_rank(
    taxonomy,
    prefix
):

    for field in clean(
        taxonomy
    ).split(
        ";"
    ):

        field = field.strip()

        if field.startswith(
            prefix
        ):

            return field[
                len(
                    prefix
                ):
            ].strip()

    return ""


def flip_strand(
    strand
):

    strand = clean(
        strand
    )

    if strand == "+":
        return "-"

    if strand == "-":
        return "+"

    return strand


def has_fegenie_hmm(
    value,
    target_hmm
):

    annotations = set(
        split_annotation(
            value
        )
    )

    return int(
        target_hmm
        in
        annotations
    )


def classify_focal_call(
    series
):

    calls = set()

    for value in series:

        calls.update(
            split_annotation(
                value
            )
        )


    has_mtoa = (
        "MtoA"
        in
        calls
    )


    has_mtra = (
        "MtrA"
        in
        calls
    )


    if (
        has_mtoa
        and
        has_mtra
    ):

        fail(
            "A single Cluster_00035 protein has conflicting "
            "MtoA and MtrA FeGenie annotations."
        )


    if has_mtra:
        return "MtrA"


    if has_mtoa:
        return "MtoA"


    return "FeGenie_negative"


def make_taxonomy_display(
    row
):

    genus = clean(
        row[
            "taxonomy_genus"
        ]
    )


    species = clean(
        row[
            "taxonomy_species"
        ]
    )


    genome = clean(
        row[
            "genome"
        ]
    )


    ## -------------------------------------------------------------- ##
    ## GTDB species strings usually already contain genus + species,
    ## e.g.:
    ##
    ##     Methylobacter_C sp003158415
    ##
    ## -------------------------------------------------------------- ##

    if (
        species
        and
        species != genus
    ):

        return species


    if genus:

        return (
            genus
            +
            " "
            +
            genome
        )


    return genome


def classify_priority_taxon(
    genus
):

    genus = clean(
        genus
    )


    ## -------------------------------------------------------------- ##
    ## GTDB-suffixed genera such as Methylobacter_A,
    ## Methylobacter_C etc. are intentionally grouped together.
    ## -------------------------------------------------------------- ##

    if genus.startswith(
        "Methylobacter"
    ):

        return "Methylobacter"


    if genus == "Crenothrix":

        return "Crenothrix"


    if genus == "Methylomonas":

        return "Methylomonas"


    if genus == "Methylococcus":

        return "Methylococcus"


    if genus == "Methylovulum":

        return "Methylovulum"


    return "Other"


def interval_gap(
    neighbor_start,
    neighbor_end,
    focal_start,
    focal_end
):

    ## -------------------------------------------------------------- ##
    ## Physical edge-to-edge distance between CDS intervals.
    ##
    ## Overlap => 0 bp.
    ## -------------------------------------------------------------- ##

    if neighbor_end < focal_start:

        return float(
            focal_start
            -
            neighbor_end
        )


    if neighbor_start > focal_end:

        return float(
            neighbor_start
            -
            focal_end
        )


    return 0.0


## ================================================================== ##
## Validate inputs
## ================================================================== ##

for path in [
    FOCAL_FILE,
    NEIGHBORHOOD_FILE,
    TAXONOMY_FILE,
]:

    if not path.exists():

        fail(
            f"Required input does not exist: {path}"
        )


print(
    "=" * 110
)

print(
    "STAGE 16A - MODULE 20 COMPLETE GENE-MAP FRAMEWORK"
)

print(
    "=" * 110
)


## ================================================================== ##
## 1. Read authoritative Stage-15A focal catalogue
## ================================================================== ##

print()
print(
    "Reading Stage-15A focal catalogue..."
)


focals = pd.read_csv(
    FOCAL_FILE,
    sep="\t",
    dtype=str,
    keep_default_na=False
)


required_focal_columns = {

    "genome",
    "module",
    "protein_id",
    "contig",
    "start",
    "end",
    "strand",

    "mmseq_cluster",
    "mcl_module",

    "fegenie_HMMs",
    "number_of_hemes",

}


missing = (
    required_focal_columns
    -
    set(
        focals.columns
    )
)


if missing:

    fail(
        "Missing focal columns: "
        +
        ", ".join(
            sorted(
                missing
            )
        )
    )


## ================================================================== ##
## 2. Identify actual Module-20 Cluster_00035 proteins
##
## IMPORTANT:
##
## The focal catalogue contains duplicate rows because one actual
## protein can be selected in several module-associated contexts.
##
## Collapse by:
##
##     genome + protein_id
##
## ================================================================== ##

target_focal_rows = focals[
    (
        focals[
            "mmseq_cluster"
        ]
        ==
        FOCAL_CLUSTER
    )
    &
    (
        focals[
            "mcl_module"
        ]
        ==
        MODULE
    )
].copy()


print(
    f"  Raw matching focal rows: "
    f"{len(target_focal_rows):,}"
)


focal_records = []


for (
    genome,
    protein_id
), group in target_focal_rows.groupby(
    [
        "genome",
        "protein_id"
    ],
    sort=True
):


    ## -------------------------------------------------------------- ##
    ## Structural identity must be consistent across duplicate
    ## Stage-15A contexts.
    ## -------------------------------------------------------------- ##

    for column in [
        "contig",
        "start",
        "end",
        "strand",
        "mmseq_cluster",
        "mcl_module",
        "number_of_hemes",
    ]:

        values = unique_nonempty(
            group[
                column
            ]
        )


        if len(
            values
        ) > 1:

            fail(
                f"Conflicting {column} values for "
                f"{genome} / {protein_id}: "
                f"{values}"
            )


    record = {

        "genome":
            genome,

        "focal_protein_id":
            protein_id,

        "focal_contig":
            first_nonempty(
                group[
                    "contig"
                ]
            ),

        "focal_start":
            int(
                float(
                    first_nonempty(
                        group[
                            "start"
                        ]
                    )
                )
            ),

        "focal_end":
            int(
                float(
                    first_nonempty(
                        group[
                            "end"
                        ]
                    )
                )
            ),

        "focal_strand":
            first_nonempty(
                group[
                    "strand"
                ]
            ),

        "focal_mmseq_cluster":
            FOCAL_CLUSTER,

        "focal_mcl_module":
            MODULE,

        "focal_call":
            classify_focal_call(
                group[
                    "fegenie_HMMs"
                ]
            ),

        "focal_number_of_hemes":
            first_nonempty(
                group[
                    "number_of_hemes"
                ]
            ),

        "n_stage15a_focal_context_rows":
            len(
                group
            ),

    }


    record[
        "plot_region_id"
    ] = (
        genome
        +
        "|"
        +
        protein_id
    )


    focal_records.append(
        record
    )


focal_table = pd.DataFrame(
    focal_records
)


n_focals = len(
    focal_table
)


n_focal_genomes = focal_table[
    "genome"
].nunique()


print(
    f"  Unique actual focal proteins: "
    f"{n_focals:,}"
)


print(
    f"  Unique focal genomes:         "
    f"{n_focal_genomes:,}"
)


if (
    n_focals
    !=
    EXPECTED_UNIQUE_FOCALS
):

    fail(
        f"Expected {EXPECTED_UNIQUE_FOCALS} unique "
        f"Cluster_00035 proteins; found {n_focals}."
    )


if (
    n_focal_genomes
    !=
    EXPECTED_FOCAL_GENOMES
):

    fail(
        f"Expected {EXPECTED_FOCAL_GENOMES} genomes; "
        f"found {n_focal_genomes}."
    )


## ================================================================== ##
## 3. Unique focal FeGenie calls
## ================================================================== ##

print()
print(
    "Unique focal FeGenie calls:"
)


print(
    focal_table[
        "focal_call"
    ]
    .value_counts()
    .to_string()
)


## ================================================================== ##
## 4. Read fresh Stage-15A neighborhoods
## ================================================================== ##

print()
print(
    "Reading fresh Stage-15A neighborhoods..."
)


neighbors = pd.read_csv(
    NEIGHBORHOOD_FILE,
    sep="\t",
    dtype=str,
    keep_default_na=False
)


required_neighbor_columns = {

    "genome",

    "neighbor_contig",
    "neighbor_gene_rank",
    "neighbor_protein_id",
    "neighbor_start",
    "neighbor_end",
    "neighbor_strand",

    "neighbor_globdb_cog",
    "neighbor_globdb_gene",
    "neighbor_globdb_product",

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

    "focal_context_module",
    "focal_protein_id",
    "focal_contig",
    "focal_start",
    "focal_end",
    "focal_strand",

    "focal_mmseq_cluster",
    "focal_mcl_module",

    "genomic_gene_offset",
    "oriented_gene_offset",

    "genomic_midpoint_offset_bp",
    "oriented_midpoint_offset_bp",

    "within_20_genes",
    "within_20kb",

    "is_focal",

}


missing = (
    required_neighbor_columns
    -
    set(
        neighbors.columns
    )
)


if missing:

    fail(
        "Missing neighborhood columns: "
        +
        ", ".join(
            sorted(
                missing
            )
        )
    )


print(
    f"  Full neighborhood rows: "
    f"{len(neighbors):,}"
)


## ================================================================== ##
## 5. Select canonical Module-20 focal contexts
##
## Require:
##
##     focal_context_module == Module_20
##     focal_mmseq_cluster == Cluster_00035
##     focal_mcl_module == Module_20
##
## This avoids duplicate copies generated when the same actual protein
## was selected in another module context.
## ================================================================== ##

target_proteins = set(
    focal_table[
        "focal_protein_id"
    ]
)


full_context = neighbors[
    (
        neighbors[
            "focal_context_module"
        ]
        ==
        MODULE
    )
    &
    (
        neighbors[
            "focal_mmseq_cluster"
        ]
        ==
        FOCAL_CLUSTER
    )
    &
    (
        neighbors[
            "focal_mcl_module"
        ]
        ==
        MODULE
    )
    &
    (
        neighbors[
            "focal_protein_id"
        ].isin(
            target_proteins
        )
    )
].copy()


print(
    f"  Canonical Module-20 neighborhood rows: "
    f"{len(full_context):,}"
)


observed_focals = set(
    full_context[
        "focal_protein_id"
    ]
)


missing_focals = sorted(
    target_proteins
    -
    observed_focals
)


if missing_focals:

    fail(
        f"{len(missing_focals)} target proteins lack a canonical "
        "Module_20 neighborhood. First examples: "
        +
        ", ".join(
            missing_focals[:10]
        )
    )


print(
    f"  Focal proteins represented: "
    f"{len(observed_focals):,}/{n_focals:,}"
)


## ================================================================== ##
## 6. One actual gene observation per focal region
## ================================================================== ##

duplicate_key = [
    "genome",
    "focal_protein_id",
    "neighbor_protein_id",
]


duplicates = (
    full_context
    .groupby(
        duplicate_key
    )
    .size()
)


duplicates = duplicates[
    duplicates > 1
]


if len(
    duplicates
) > 0:

    fail(
        f"Found {len(duplicates)} duplicated "
        "focal-region + neighbor-gene combinations."
    )


## ================================================================== ##
## 7. Numeric columns
## ================================================================== ##

for column in [
    "neighbor_start",
    "neighbor_end",
    "focal_start",
    "focal_end",
    "genomic_gene_offset",
    "oriented_gene_offset",
    "genomic_midpoint_offset_bp",
    "oriented_midpoint_offset_bp",
]:

    full_context[
        column
    ] = pd.to_numeric(
        full_context[
            column
        ],
        errors="raise"
    )


## ================================================================== ##
## 8. Normalize orientation around focal protein
##
## Every focal protein is displayed pointing left -> right.
##
## Focal midpoint becomes x = 0.
## ================================================================== ##

full_context[
    "focal_midpoint"
] = (
    full_context[
        "focal_start"
    ]
    +
    full_context[
        "focal_end"
    ]
) / 2


plus_mask = (
    full_context[
        "focal_strand"
    ]
    ==
    "+"
)


minus_mask = (
    full_context[
        "focal_strand"
    ]
    ==
    "-"
)


if not (
    plus_mask
    |
    minus_mask
).all():

    bad = sorted(
        set(
            full_context.loc[
                ~(
                    plus_mask
                    |
                    minus_mask
                ),
                "focal_strand"
            ]
        )
    )

    fail(
        f"Unexpected focal strand values: {bad}"
    )


full_context[
    "plot_start_bp"
] = np.nan


full_context[
    "plot_end_bp"
] = np.nan


full_context[
    "plot_strand"
] = ""


## ------------------------------------------------------------------ ##
## Forward focal
## ------------------------------------------------------------------ ##

full_context.loc[
    plus_mask,
    "plot_start_bp"
] = (
    full_context.loc[
        plus_mask,
        "neighbor_start"
    ]
    -
    full_context.loc[
        plus_mask,
        "focal_midpoint"
    ]
)


full_context.loc[
    plus_mask,
    "plot_end_bp"
] = (
    full_context.loc[
        plus_mask,
        "neighbor_end"
    ]
    -
    full_context.loc[
        plus_mask,
        "focal_midpoint"
    ]
)


full_context.loc[
    plus_mask,
    "plot_strand"
] = full_context.loc[
    plus_mask,
    "neighbor_strand"
]


## ------------------------------------------------------------------ ##
## Reverse focal
## ------------------------------------------------------------------ ##

full_context.loc[
    minus_mask,
    "plot_start_bp"
] = (
    full_context.loc[
        minus_mask,
        "focal_midpoint"
    ]
    -
    full_context.loc[
        minus_mask,
        "neighbor_end"
    ]
)


full_context.loc[
    minus_mask,
    "plot_end_bp"
] = (
    full_context.loc[
        minus_mask,
        "focal_midpoint"
    ]
    -
    full_context.loc[
        minus_mask,
        "neighbor_start"
    ]
)


full_context.loc[
    minus_mask,
    "plot_strand"
] = full_context.loc[
    minus_mask,
    "neighbor_strand"
].map(
    flip_strand
)


full_context[
    "gene_midpoint"
] = (
    full_context[
        "plot_start_bp"
    ]
    +
    full_context[
        "plot_end_bp"
    ]
) / 2


## ================================================================== ##
## 9. Validate our orientation against Stage-15A oriented midpoint
##
## These should agree apart from tiny floating-point differences.
## ================================================================== ##

offset_difference = (
    full_context[
        "gene_midpoint"
    ]
    -
    full_context[
        "oriented_midpoint_offset_bp"
    ]
).abs()


max_offset_difference = float(
    offset_difference.max()
)


print()
print(
    f"Maximum difference from Stage-15A oriented midpoint: "
    f"{max_offset_difference:.6f} bp"
)


if max_offset_difference > 0.001:

    fail(
        "Recomputed normalized coordinates disagree with "
        "Stage-15A oriented midpoint offsets."
    )


## ================================================================== ##
## 10. Region identity
## ================================================================== ##

full_context[
    "plot_region_id"
] = (
    full_context[
        "genome"
    ]
    +
    "|"
    +
    full_context[
        "focal_protein_id"
    ]
)


if (
    full_context[
        "plot_region_id"
    ].nunique()
    !=
    EXPECTED_UNIQUE_FOCALS
):

    fail(
        "Canonical context does not contain exactly 74 regions."
    )


## ================================================================== ##
## 11. Identify actual focal row
## ================================================================== ##

full_context[
    "plot_is_focal"
] = (
    full_context[
        "neighbor_protein_id"
    ]
    ==
    full_context[
        "focal_protein_id"
    ]
).astype(
    int
)


focal_counts = (
    full_context
    .groupby(
        "plot_region_id"
    )[
        "plot_is_focal"
    ]
    .sum()
)


bad_focal_counts = focal_counts[
    focal_counts != 1
]


if len(
    bad_focal_counts
) > 0:

    fail(
        f"{len(bad_focal_counts)} regions do not contain "
        "exactly one focal CDS."
    )


## ================================================================== ##
## 12. Attach unique focal call
## ================================================================== ##

focal_call_map = (
    focal_table
    .set_index(
        [
            "genome",
            "focal_protein_id"
        ]
    )[
        "focal_call"
    ]
    .to_dict()
)


full_context[
    "focal_call"
] = [

    focal_call_map[
        (
            genome,
            focal_protein
        )
    ]

    for genome, focal_protein in zip(

        full_context[
            "genome"
        ],

        full_context[
            "focal_protein_id"
        ],

    )

]


## ================================================================== ##
## 13. Identify MtrB using authoritative FeGenie HMM
##
## Primary MtrB definition:
##
##     neighbor_fegenie_HMMs contains exact token:
##
##         MtrB_TIGR03509
##
## ================================================================== ##

full_context[
    "neighbor_is_mtrb_fegenie"
] = full_context[
    "neighbor_fegenie_HMMs"
].map(
    lambda x:
        has_fegenie_hmm(
            x,
            MTRB_FEGENIE_HMM
        )
)


print()
print(
    "FeGenie MtrB identification:"
)


print(
    f"  MtrB gene observations in complete context: "
    f"{full_context['neighbor_is_mtrb_fegenie'].sum():,}"
)


## ================================================================== ##
## 14. Physical gene-to-focal distance
##
## Distinguish:
##
##     midpoint distance
##     edge-to-edge CDS gap
##
## This is useful for prioritizing genuinely tight MtrA/MtoA-MtrB
## architectures.
## ================================================================== ##

full_context[
    "distance_to_focal_cds_bp"
] = full_context.apply(

    lambda row:
        interval_gap(

            neighbor_start=
                row[
                    "plot_start_bp"
                ],

            neighbor_end=
                row[
                    "plot_end_bp"
                ],

            focal_start=
                -(
                    row[
                        "focal_end"
                    ]
                    -
                    row[
                        "focal_start"
                    ]
                )
                / 2,

            focal_end=
                (
                    row[
                        "focal_end"
                    ]
                    -
                    row[
                        "focal_start"
                    ]
                )
                / 2,

        ),

    axis=1

)


## ================================================================== ##
## 15. Read authoritative r226 taxonomy
##
## IMPORTANT:
##
## globdb_r226_taxonomy.tsv is HEADERLESS.
## ================================================================== ##

print()
print(
    "Reading authoritative r226 taxonomy..."
)


taxonomy = pd.read_csv(

    TAXONOMY_FILE,

    sep="\t",

    header=None,

    names=[
        "genome",
        "taxonomy"
    ],

    usecols=[
        0,
        1
    ],

    dtype=str,

    keep_default_na=False

)


print(
    f"  Taxonomy rows: "
    f"{len(taxonomy):,}"
)


if taxonomy[
    "genome"
].duplicated().any():

    duplicated = taxonomy.loc[

        taxonomy[
            "genome"
        ].duplicated(
            keep=False
        ),

        "genome"

    ].unique()


    fail(
        "Taxonomy contains duplicate genome IDs. "
        "First examples: "
        +
        ", ".join(
            duplicated[:10]
        )
    )


taxonomy[
    "taxonomy_family"
] = taxonomy[
    "taxonomy"
].map(
    lambda x:
        parse_taxonomic_rank(
            x,
            "f__"
        )
)


taxonomy[
    "taxonomy_genus"
] = taxonomy[
    "taxonomy"
].map(
    lambda x:
        parse_taxonomic_rank(
            x,
            "g__"
        )
)


taxonomy[
    "taxonomy_species"
] = taxonomy[
    "taxonomy"
].map(
    lambda x:
        parse_taxonomic_rank(
            x,
            "s__"
        )
)


taxonomy[
    "taxonomy_display"
] = taxonomy.apply(
    make_taxonomy_display,
    axis=1
)


taxonomy_keep = taxonomy[
    [
        "genome",
        "taxonomy",
        "taxonomy_family",
        "taxonomy_genus",
        "taxonomy_species",
        "taxonomy_display",
    ]
]


## ================================================================== ##
## 16. Attach taxonomy
## ================================================================== ##

full_context = full_context.merge(

    taxonomy_keep,

    on="genome",

    how="left",

    validate="many_to_one"

)


missing_taxonomy = full_context[
    "taxonomy"
].isna()


if missing_taxonomy.any():

    missing_genomes = sorted(
        full_context.loc[
            missing_taxonomy,
            "genome"
        ].unique()
    )


    fail(
        f"Missing r226 taxonomy for "
        f"{len(missing_genomes)} genomes: "
        +
        ", ".join(
            missing_genomes[:10]
        )
    )


## ================================================================== ##
## 17. Taxonomic priority metadata
##
## This does NOT select anything.
## ================================================================== ##

full_context[
    "priority_taxon"
] = full_context[
    "taxonomy_genus"
].map(
    classify_priority_taxon
)


taxon_rank = {

    "Methylobacter":
        1,

    "Crenothrix":
        1,

    ## Ranks 2 and 3 are architecture-dependent:
    ##
    ##   2. MtoA/MtrA + nearby MtrB
    ##   3. MtrA
    ##
    ## They are assigned during final representative selection.

    "Methylomonas":
        4,

    "Methylococcus":
        5,

    "Methylovulum":
        6,

    "Other":
        99,

}


full_context[
    "taxon_priority_rank"
] = full_context[
    "priority_taxon"
].map(
    taxon_rank
)


## ================================================================== ##
## 18. Build FULL-context region summaries
##
## IMPORTANT:
##
## Do this BEFORE cropping to +/-10 kb.
##
## Therefore MtrB at 10-20 kb remains visible to representative
## selection even though it will not appear in the final +/-10 kb
## gene-map input.
## ================================================================== ##

region_rows = []


for region_id, group in full_context.groupby(
    "plot_region_id",
    sort=True
):

    focal = group[
        group[
            "plot_is_focal"
        ]
        ==
        1
    ].iloc[
        0
    ]


    mtrb = group[
        group[
            "neighbor_is_mtrb_fegenie"
        ]
        ==
        1
    ].copy()


    partner = group[
        group[
            "neighbor_mmseq_cluster"
        ]
        ==
        MODULE_PARTNER
    ].copy()


    accessory_166 = group[
        group[
            "neighbor_mmseq_cluster"
        ]
        ==
        "Cluster_00166"
    ].copy()


    accessory_209 = group[
        group[
            "neighbor_mmseq_cluster"
        ]
        ==
        "Cluster_00209"
    ].copy()


    ## -------------------------------------------------------------- ##
    ## Closest-row helper
    ## -------------------------------------------------------------- ##

    def closest_row(
        subset
    ):

        if subset.empty:
            return None


        return subset.loc[
            subset[
                "gene_midpoint"
            ]
            .abs()
            .idxmin()
        ]


    closest_mtrb = closest_row(
        mtrb
    )


    closest_partner = closest_row(
        partner
    )


    ## -------------------------------------------------------------- ##
    ## MtrB architecture
    ## -------------------------------------------------------------- ##

    if closest_mtrb is None:

        closest_mtrb_midpoint_bp = np.nan

        closest_mtrb_gap_bp = np.nan

        closest_mtrb_gene_offset = np.nan

        closest_mtrb_order = ""

        closest_mtrb_relative_strand = ""

        closest_mtrb_adjacent = 0


    else:

        closest_mtrb_midpoint_bp = abs(
            float(
                closest_mtrb[
                    "gene_midpoint"
                ]
            )
        )


        closest_mtrb_gap_bp = float(
            closest_mtrb[
                "distance_to_focal_cds_bp"
            ]
        )


        closest_mtrb_gene_offset = int(
            closest_mtrb[
                "oriented_gene_offset"
            ]
        )


        if (
            closest_mtrb[
                "gene_midpoint"
            ]
            <
            0
        ):

            closest_mtrb_order = (
                "MtrB->focal"
            )

        elif (
            closest_mtrb[
                "gene_midpoint"
            ]
            >
            0
        ):

            closest_mtrb_order = (
                "focal->MtrB"
            )

        else:

            closest_mtrb_order = (
                "same_midpoint"
            )


        closest_mtrb_relative_strand = clean(
            closest_mtrb[
                "plot_strand"
            ]
        )


        closest_mtrb_adjacent = int(
            abs(
                closest_mtrb_gene_offset
            )
            ==
            1
        )


    ## -------------------------------------------------------------- ##
    ## Cluster_00048 partner
    ## -------------------------------------------------------------- ##

    if closest_partner is None:

        closest_partner_midpoint_bp = np.nan

    else:

        closest_partner_midpoint_bp = abs(
            float(
                closest_partner[
                    "gene_midpoint"
                ]
            )
        )


    def closest_abs_midpoint(
        subset
    ):

        if subset.empty:
            return np.nan


        return float(
            subset[
                "gene_midpoint"
            ]
            .abs()
            .min()
        )


    region_rows.append(
        {

            "plot_region_id":
                region_id,

            "genome":
                focal[
                    "genome"
                ],

            "focal_protein_id":
                focal[
                    "focal_protein_id"
                ],

            "focal_call":
                focal[
                    "focal_call"
                ],

            "taxonomy_family":
                focal[
                    "taxonomy_family"
                ],

            "taxonomy_genus":
                focal[
                    "taxonomy_genus"
                ],

            "taxonomy_species":
                focal[
                    "taxonomy_species"
                ],

            "taxonomy_display":
                focal[
                    "taxonomy_display"
                ],

            "priority_taxon":
                focal[
                    "priority_taxon"
                ],

            "taxon_priority_rank":
                focal[
                    "taxon_priority_rank"
                ],


            ## ------------------------------------------------------ ##
            ## Full Stage-15A context
            ## ------------------------------------------------------ ##

            "n_genes_full_context":
                len(
                    group
                ),


            ## ------------------------------------------------------ ##
            ## MtrB
            ## ------------------------------------------------------ ##

            "n_mtrb":
                len(
                    mtrb
                ),

            "closest_mtrb_midpoint_bp":
                closest_mtrb_midpoint_bp,

            "closest_mtrb_intergenic_gap_bp":
                closest_mtrb_gap_bp,

            "closest_mtrb_oriented_gene_offset":
                closest_mtrb_gene_offset,

            "closest_mtrb_adjacent":
                closest_mtrb_adjacent,

            "closest_mtrb_order":
                closest_mtrb_order,

            "closest_mtrb_relative_strand":
                closest_mtrb_relative_strand,


            ## ------------------------------------------------------ ##
            ## Known Module-20 families
            ## ------------------------------------------------------ ##

            "n_cluster00048":
                len(
                    partner
                ),

            "closest_cluster00048_midpoint_bp":
                closest_partner_midpoint_bp,

            "n_cluster00166":
                len(
                    accessory_166
                ),

            "closest_cluster00166_midpoint_bp":
                closest_abs_midpoint(
                    accessory_166
                ),

            "n_cluster00209":
                len(
                    accessory_209
                ),

            "closest_cluster00209_midpoint_bp":
                closest_abs_midpoint(
                    accessory_209
                ),

        }
    )


region_summary = pd.DataFrame(
    region_rows
)


## ================================================================== ##
## 19. Explicit physical-distance flags
##
## Keep 5 / 10 / 20 kb separately rather than choosing an arbitrary
## definition of "nearby".
## ================================================================== ##

for threshold in [
    5_000,
    10_000,
    20_000,
]:

    suffix = (
        f"{threshold // 1000}kb"
    )


    region_summary[
        f"mtrb_within_{suffix}"
    ] = (
        region_summary[
            "closest_mtrb_midpoint_bp"
        ]
        <=
        threshold
    ).fillna(
        False
    ).astype(
        int
    )


    region_summary[
        f"cluster00048_within_{suffix}"
    ] = (
        region_summary[
            "closest_cluster00048_midpoint_bp"
        ]
        <=
        threshold
    ).fillna(
        False
    ).astype(
        int
    )


## ================================================================== ##
## 20. Convenience MtrB architecture classification
##
## This is descriptive metadata only.
##
## No biological interpretation is imposed.
## ================================================================== ##

def classify_mtrb_proximity(
    row
):

    if int(
        row[
            "n_mtrb"
        ]
    ) == 0:

        return "no_MtrB"


    if int(
        row[
            "closest_mtrb_adjacent"
        ]
    ) == 1:

        return "adjacent_MtrB"


    distance = row[
        "closest_mtrb_midpoint_bp"
    ]


    if pd.isna(
        distance
    ):

        return "MtrB_distance_unknown"


    if distance <= 5_000:

        return "MtrB_within_5kb"


    if distance <= 10_000:

        return "MtrB_within_10kb"


    if distance <= 20_000:

        return "MtrB_within_20kb"


    return "MtrB_beyond_20kb"


region_summary[
    "mtrb_proximity_class"
] = region_summary.apply(
    classify_mtrb_proximity,
    axis=1
)


region_summary[
    "focal_mtrb_architecture"
] = (
    region_summary[
        "focal_call"
    ]
    +
    "|"
    +
    region_summary[
        "mtrb_proximity_class"
    ]
)


## ================================================================== ##
## 21. Crop to +/-10 kb for local-family analysis and visualization
##
## Retain a CDS if ANY part intersects the plotting interval.
## ================================================================== ##

plot_data = full_context[
    (
        full_context[
            "plot_end_bp"
        ]
        >=
        -PLOT_WINDOW_BP
    )
    &
    (
        full_context[
            "plot_start_bp"
        ]
        <=
        PLOT_WINDOW_BP
    )
].copy()


print()
print(
    f"Genes intersecting +/-{PLOT_WINDOW_BP:,} bp: "
    f"{len(plot_data):,}"
)


if (
    plot_data[
        "plot_region_id"
    ].nunique()
    !=
    EXPECTED_UNIQUE_FOCALS
):

    fail(
        "At least one focal region disappeared after +/-10 kb crop."
    )


## ================================================================== ##
## 22. Confirm focal exactly once after crop
## ================================================================== ##

cropped_focal_counts = (
    plot_data
    .groupby(
        "plot_region_id"
    )[
        "plot_is_focal"
    ]
    .sum()
)


if (
    cropped_focal_counts
    !=
    1
).any():

    fail(
        "A cropped plot region does not contain exactly one focal gene."
    )


## ================================================================== ##
## 23. Visualization-level classes
##
## These are intentionally preliminary.
##
## Final architecture classes/colors will be decided only after the
## exhaustive Module-20 local architecture analysis.
## ================================================================== ##

def display_class(
    row
):

    if int(
        row[
            "plot_is_focal"
        ]
    ) == 1:

        return "focal_cluster00035"


    ## -------------------------------------------------------------- ##
    ## MtrB gets explicit prominence because it is a major biological
    ## criterion for representative-region selection.
    ## -------------------------------------------------------------- ##

    if int(
        row[
            "neighbor_is_mtrb_fegenie"
        ]
    ) == 1:

        return "mtrb"


    if clean(
        row[
            "neighbor_mmseq_cluster"
        ]
    ) == MODULE_PARTNER:

        return "module20_partner"


    if numeric_flag(
        row[
            "neighbor_fegenie_positive"
        ]
    ):

        return "fegenie_other"


    if numeric_flag(
        row[
            "neighbor_findmehemes_positive"
        ]
    ):

        return "findmehemes_other"


    if clean(
        row[
            "neighbor_mmseq_cluster"
        ]
    ):

        return "other_mmseq_family"


    return "background_gene"


plot_data[
    "display_class"
] = plot_data.apply(
    display_class,
    axis=1
)


## ================================================================== ##
## 24. Preliminary labels
##
## Do not over-label here.
##
## Local architecture analysis will determine which recurrent families
## deserve final figure labels.
## ================================================================== ##

def display_label(
    row
):

    if int(
        row[
            "plot_is_focal"
        ]
    ) == 1:

        return clean(
            row[
                "focal_call"
            ]
        )


    if int(
        row[
            "neighbor_is_mtrb_fegenie"
        ]
    ) == 1:

        return "MtrB"


    cluster = clean(
        row[
            "neighbor_mmseq_cluster"
        ]
    )


    if cluster == MODULE_PARTNER:

        return "Cluster_00048"


    hmm = clean(
        row[
            "neighbor_fegenie_HMMs"
        ]
    )


    if hmm:

        return hmm


    return ""


plot_data[
    "display_label"
] = plot_data.apply(
    display_label,
    axis=1
)


## ================================================================== ##
## 25. Generic protein ID required by local architecture pipeline
## ================================================================== ##

plot_data[
    "protein_id"
] = plot_data[
    "neighbor_protein_id"
]


## ================================================================== ##
## 26. Add number of genes actually present in +/-10 kb
## ================================================================== ##

window_counts = (
    plot_data
    .groupby(
        "plot_region_id"
    )
    .size()
    .rename(
        "n_genes_in_10kb_window"
    )
)


region_summary = region_summary.merge(

    window_counts,

    left_on=
        "plot_region_id",

    right_index=True,

    how="left",

    validate="one_to_one"

)


## ================================================================== ##
## 27. Sort
## ================================================================== ##

plot_data = plot_data.sort_values(

    [
        "taxonomy_family",
        "taxonomy_genus",
        "taxonomy_species",

        "genome",
        "focal_protein_id",

        "plot_start_bp",
        "neighbor_gene_rank",
    ],

    kind="stable"

)


region_summary = region_summary.sort_values(

    [
        "taxon_priority_rank",
        "taxonomy_genus",
        "taxonomy_species",

        "genome",
        "focal_protein_id",
    ],

    kind="stable"

)


## ================================================================== ##
## 28. Write outputs
## ================================================================== ##

plot_data.to_csv(
    OUT_DATA,
    sep="\t",
    index=False
)


region_summary.to_csv(
    OUT_REGIONS,
    sep="\t",
    index=False
)


## ================================================================== ##
## 29. QC output
## ================================================================== ##

qc_rows = [

    {
        "metric":
            "raw_cluster00035_focal_rows",
        "value":
            len(
                target_focal_rows
            )
    },

    {
        "metric":
            "unique_cluster00035_focals",
        "value":
            n_focals
    },

    {
        "metric":
            "unique_cluster00035_genomes",
        "value":
            n_focal_genomes
    },

    {
        "metric":
            "canonical_full_context_rows",
        "value":
            len(
                full_context
            )
    },

    {
        "metric":
            "plot_regions",
        "value":
            plot_data[
                "plot_region_id"
            ].nunique()
    },

    {
        "metric":
            "plot_gene_rows",
        "value":
            len(
                plot_data
            )
    },

    {
        "metric":
            "unique_plot_proteins",
        "value":
            plot_data[
                "protein_id"
            ].nunique()
    },

    {
        "metric":
            "mtrb_gene_observations_full_context",
        "value":
            int(
                full_context[
                    "neighbor_is_mtrb_fegenie"
                ].sum()
            )
    },

    {
        "metric":
            "regions_with_mtrb",
        "value":
            int(
                (
                    region_summary[
                        "n_mtrb"
                    ]
                    >
                    0
                ).sum()
            )
    },

    {
        "metric":
            "regions_with_adjacent_mtrb",
        "value":
            int(
                region_summary[
                    "closest_mtrb_adjacent"
                ].sum()
            )
    },

    {
        "metric":
            "regions_with_mtrb_within_5kb",
        "value":
            int(
                region_summary[
                    "mtrb_within_5kb"
                ].sum()
            )
    },

    {
        "metric":
            "regions_with_mtrb_within_10kb",
        "value":
            int(
                region_summary[
                    "mtrb_within_10kb"
                ].sum()
            )
    },

    {
        "metric":
            "regions_with_mtrb_within_20kb",
        "value":
            int(
                region_summary[
                    "mtrb_within_20kb"
                ].sum()
            )
    },

    {
        "metric":
            "regions_with_cluster00048",
        "value":
            int(
                (
                    region_summary[
                        "n_cluster00048"
                    ]
                    >
                    0
                ).sum()
            )
    },

    {
        "metric":
            "max_normalized_offset_validation_difference_bp",
        "value":
            max_offset_difference
    },

]


for focal_call, n in (
    focal_table[
        "focal_call"
    ]
    .value_counts()
    .items()
):

    qc_rows.append(
        {

            "metric":
                f"focal_call_{focal_call}",

            "value":
                int(
                    n
                )

        }
    )


pd.DataFrame(
    qc_rows
).to_csv(
    OUT_QC,
    sep="\t",
    index=False
)


## ================================================================== ##
## 30. Console summary
## ================================================================== ##

print()
print(
    "=" * 110
)

print(
    "MODULE 20 PREPARATION SUMMARY"
)

print(
    "=" * 110
)


print(
    f"Unique Cluster_00035 focals: "
    f"{n_focals}"
)


print(
    f"Focal genomes:                "
    f"{n_focal_genomes}"
)


print(
    f"Plot regions:                 "
    f"{plot_data['plot_region_id'].nunique()}"
)


print(
    f"Full-context gene rows:       "
    f"{len(full_context):,}"
)


print(
    f"Gene rows in +/-10 kb:        "
    f"{len(plot_data):,}"
)


print(
    f"Unique +/-10 kb proteins:     "
    f"{plot_data['protein_id'].nunique():,}"
)


print()
print(
    "Focal calls:"
)


print(
    focal_table[
        "focal_call"
    ]
    .value_counts()
    .to_string()
)


print()
print(
    "Priority taxa among focal regions:"
)


print(
    region_summary[
        "priority_taxon"
    ]
    .value_counts()
    .to_string()
)


print()
print(
    f"FeGenie MtrB ({MTRB_FEGENIE_HMM}):"
)


print(
    f"  MtrB observations:          "
    f"{full_context['neighbor_is_mtrb_fegenie'].sum()}"
)


print(
    f"  regions containing MtrB:    "
    f"{(region_summary['n_mtrb'] > 0).sum()}"
)


print(
    f"  adjacent to focal:          "
    f"{region_summary['closest_mtrb_adjacent'].sum()}"
)


print(
    f"  within 5 kb:                "
    f"{region_summary['mtrb_within_5kb'].sum()}"
)


print(
    f"  within 10 kb:               "
    f"{region_summary['mtrb_within_10kb'].sum()}"
)


print(
    f"  within 20 kb:               "
    f"{region_summary['mtrb_within_20kb'].sum()}"
)


print()
print(
    "Focal-call x MtrB proximity:"
)


print(
    pd.crosstab(

        region_summary[
            "focal_call"
        ],

        region_summary[
            "mtrb_proximity_class"
        ]

    ).to_string()
)


print()
print(
    "Known Cluster_00048 partner:"
)


print(
    f"  regions containing partner: "
    f"{(region_summary['n_cluster00048'] > 0).sum()}"
)


print(
    f"  within 5 kb:                "
    f"{region_summary['cluster00048_within_5kb'].sum()}"
)


print(
    f"  within 10 kb:               "
    f"{region_summary['cluster00048_within_10kb'].sum()}"
)


print(
    f"  within 20 kb:               "
    f"{region_summary['cluster00048_within_20kb'].sum()}"
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
    OUT_DATA
)


print(
    OUT_REGIONS
)


print(
    OUT_QC
)


print()
print(
    "SUCCESS"
)