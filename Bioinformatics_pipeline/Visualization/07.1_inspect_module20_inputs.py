#!/usr/bin/env python3

## ================================================================== ##
## STAGE 16A - MODULE 20 INPUT AUDIT
##
## Purpose
## -------
## Inspect the authoritative Stage-15A tables before constructing
## complete Cluster_00035-centered Module-20 neighborhoods.
##
## This script performs NO filtering for representative genomes and
## NO biological classification beyond identifying the focal
## Cluster_00035 proteins.
##
## It reports:
##
##   - exact schemas of the focal and neighborhood tables
##   - number of Cluster_00035 Module-20 focal proteins/genomes
##   - FeGenie MtoA / MtrA / negative distribution
##   - possible focal-region ID columns
##   - possible coordinate / contig / censoring columns
##   - where MtrB occurs in the neighborhood table
##   - taxonomy matches for the taxa prioritized for visualization
##
## After this audit, the production Module-20 preparation script can
## use explicit authoritative column names rather than guessed ones.
## ================================================================== ##


from pathlib import Path
from collections import Counter
import sys

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


OUT_FILE = (
    WORKFLOW
    / "16_visualization"
    / "16A_module20_gene_map_framework"
    / "Module_20_input_audit.txt"
)


## ================================================================== ##
## Constants
## ================================================================== ##

MODULE = "Module_20"

FOCAL_CLUSTER = "Cluster_00035"


PRIORITY_TAXA = [
    "Methylobacter",
    "Crenothrix",
    "Methylomonas",
    "Methylococcus",
    "Methylovulum",
]


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


def heading(text):

    print()
    print("=" * 110)
    print(text)
    print("=" * 110)


def print_columns(df):

    for index, column in enumerate(
        df.columns,
        start=1
    ):

        print(
            f"{index:3d}. {column}"
        )


def columns_containing(
    columns,
    keywords
):

    result = []

    for column in columns:

        lower = column.lower()

        if any(
            keyword.lower() in lower
            for keyword in keywords
        ):

            result.append(
                column
            )

    return result


## ================================================================== ##
## Validate files
## ================================================================== ##

for path in [
    FOCAL_FILE,
    NEIGHBORHOOD_FILE,
    TAXONOMY_FILE,
]:

    if not path.exists():

        fail(
            f"Required file not found: {path}"
        )


## ================================================================== ##
## Read focal catalogue
## ================================================================== ##

heading(
    "1. FOCAL GENE CATALOGUE"
)


focals = pd.read_csv(
    FOCAL_FILE,
    sep="\t",
    dtype=str,
    keep_default_na=False
)


print(
    f"File: {FOCAL_FILE}"
)


print(
    f"Rows: {len(focals):,}"
)


print(
    f"Columns: {len(focals.columns):,}"
)


print()
print(
    "COLUMN NAMES"
)


print_columns(
    focals
)


## ================================================================== ##
## Required focal fields
## ================================================================== ##

required_focal = [
    "genome",
    "protein_id",
    "mmseq_cluster",
    "mcl_module",
]


missing = [
    column
    for column in required_focal
    if column not in focals.columns
]


if missing:

    fail(
        "Focal catalogue is missing expected columns: "
        +
        ", ".join(
            missing
        )
    )


## ================================================================== ##
## Select authoritative Cluster_00035 focal proteins
## ================================================================== ##

module20_focals = focals[
    (
        focals[
            "mcl_module"
        ]
        ==
        MODULE
    )
    &
    (
        focals[
            "mmseq_cluster"
        ]
        ==
        FOCAL_CLUSTER
    )
].copy()


heading(
    "2. MODULE 20 / CLUSTER_00035 FOCALS"
)


print(
    f"Rows:             {len(module20_focals):,}"
)


print(
    f"Unique proteins:  "
    f"{module20_focals['protein_id'].nunique():,}"
)


print(
    f"Unique genomes:   "
    f"{module20_focals['genome'].nunique():,}"
)


## ------------------------------------------------------------------ ##
## Check whether genome + protein is unique
## ------------------------------------------------------------------ ##

duplicate_focals = (
    module20_focals
    .groupby(
        [
            "genome",
            "protein_id"
        ]
    )
    .size()
)


duplicate_focals = duplicate_focals[
    duplicate_focals > 1
]


print(
    f"Duplicate genome+protein focal rows: "
    f"{len(duplicate_focals):,}"
)


## ================================================================== ##
## FeGenie focal calls
## ================================================================== ##

heading(
    "3. CLUSTER_00035 FEGENIE CALLS"
)


if "fegenie_HMMs" in module20_focals.columns:

    focal_calls = (
        module20_focals[
            "fegenie_HMMs"
        ]
        .map(
            clean
        )
        .replace(
            "",
            "FeGenie_negative"
        )
        .value_counts(
            dropna=False
        )
    )


    print(
        focal_calls.to_string()
    )


else:

    print(
        "No fegenie_HMMs column found."
    )


## ================================================================== ##
## Useful focal fields
## ================================================================== ##

heading(
    "4. FIRST CLUSTER_00035 FOCAL ROWS"
)


show_columns = [

    column

    for column in [
        "genome",
        "protein_id",
        "focal_reasons",
        "mmseq_cluster",
        "mcl_module",
        "fegenie_HMMs",
        "number_of_hemes",
        "signalp_prediction",
        "deeptmhmm_class",
        "report_topology",
        "globdb_cog",
        "globdb_product",
    ]

    if column in module20_focals.columns

]


print(
    module20_focals[
        show_columns
    ]
    .head(
        10
    )
    .to_string(
        index=False
    )
)


## ================================================================== ##
## Read neighborhood schema
## ================================================================== ##

heading(
    "5. FRESH FOCAL-NEIGHBORHOOD TABLE"
)


neighborhoods = pd.read_csv(
    NEIGHBORHOOD_FILE,
    sep="\t",
    dtype=str,
    keep_default_na=False
)


print(
    f"File: {NEIGHBORHOOD_FILE}"
)


print(
    f"Rows: {len(neighborhoods):,}"
)


print(
    f"Columns: {len(neighborhoods.columns):,}"
)


print()
print(
    "COLUMN NAMES"
)


print_columns(
    neighborhoods
)


## ================================================================== ##
## Candidate schema fields
## ================================================================== ##

heading(
    "6. CANDIDATE STRUCTURAL COLUMNS"
)


groups = {

    "FOCAL-related":
        [
            "focal"
        ],

    "PROTEIN / GENE-related":
        [
            "protein",
            "gene"
        ],

    "CONTIG-related":
        [
            "contig",
            "seqid",
            "sequence"
        ],

    "COORDINATE-related":
        [
            "start",
            "end",
            "coord",
            "midpoint",
            "distance",
            "offset"
        ],

    "STRAND-related":
        [
            "strand"
        ],

    "EDGE / CENSOR-related":
        [
            "edge",
            "censor",
            "trunc",
            "boundary"
        ],

    "ANNOTATION-related":
        [
            "fegenie",
            "globdb",
            "cog",
            "product",
            "mmseq"
        ],

}


for label, keywords in groups.items():

    print()
    print(
        f"{label}:"
    )

    matches = columns_containing(
        neighborhoods.columns,
        keywords
    )

    if matches:

        for column in matches:
            print(
                f"  {column}"
            )

    else:

        print(
            "  [none]"
        )


## ================================================================== ##
## Find columns containing the selected focal IDs
## ================================================================== ##

heading(
    "7. WHICH NEIGHBORHOOD COLUMNS CONTAIN FOCAL PROTEIN IDs?"
)


focal_ids = set(
    module20_focals[
        "protein_id"
    ]
)


first_focal_id = (
    module20_focals[
        "protein_id"
    ]
    .iloc[0]
    if len(
        module20_focals
    ) > 0
    else ""
)


columns_with_focal_ids = []


for column in neighborhoods.columns:

    values = neighborhoods[
        column
    ]


    if values.isin(
        focal_ids
    ).any():

        n_matches = int(
            values.isin(
                focal_ids
            ).sum()
        )

        columns_with_focal_ids.append(
            (
                column,
                n_matches
            )
        )


if columns_with_focal_ids:

    for column, n_matches in columns_with_focal_ids:

        print(
            f"{column}: {n_matches:,} rows matching "
            "a Cluster_00035 focal protein ID"
        )


else:

    print(
        "No neighborhood column contained exact Cluster_00035 "
        "protein IDs."
    )


## ================================================================== ##
## Show rows involving the first focal protein
## ================================================================== ##

heading(
    "8. EXAMPLE NEIGHBORHOOD ROWS FOR ONE CLUSTER_00035 FOCAL"
)


print(
    f"Example focal protein: {first_focal_id}"
)


if first_focal_id:

    row_mask = pd.Series(
        False,
        index=neighborhoods.index
    )


    for column in neighborhoods.columns:

        row_mask |= (
            neighborhoods[
                column
            ]
            ==
            first_focal_id
        )


    example_rows = neighborhoods[
        row_mask
    ]


    print(
        f"Rows containing this focal ID somewhere: "
        f"{len(example_rows):,}"
    )


    if len(
        example_rows
    ) > 0:

        print()

        print(
            example_rows
            .head(
                10
            )
            .to_string(
                index=False
            )
        )


## ================================================================== ##
## Search for MtrB annotations
##
## We deliberately search every column here because we do not yet
## assume whether MtrB is represented as FeGenie HMM, product text,
## or another annotation field.
## ================================================================== ##

heading(
    "9. MtrB OCCURRENCES IN FRESH NEIGHBORHOODS"
)


mtrb_columns = []


for column in neighborhoods.columns:

    values = (
        neighborhoods[
            column
        ]
        .astype(
            str
        )
    )


    mask = values.str.contains(
        r"\bMtrB\b",
        case=False,
        regex=True,
        na=False
    )


    n = int(
        mask.sum()
    )


    if n > 0:

        mtrb_columns.append(
            (
                column,
                n
            )
        )


if mtrb_columns:

    for column, n in mtrb_columns:

        print(
            f"{column}: {n:,} MtrB-containing rows"
        )


else:

    print(
        "No literal MtrB annotation found in the neighborhood table."
    )


## ------------------------------------------------------------------ ##
## Show first MtrB-containing rows
## ------------------------------------------------------------------ ##

if mtrb_columns:

    mtrb_mask = pd.Series(
        False,
        index=neighborhoods.index
    )


    for column, _ in mtrb_columns:

        mtrb_mask |= (
            neighborhoods[
                column
            ]
            .astype(
                str
            )
            .str.contains(
                r"\bMtrB\b",
                case=False,
                regex=True,
                na=False
            )
        )


    print()
    print(
        "First MtrB-containing neighborhood rows:"
    )
    print()

    print(
        neighborhoods[
            mtrb_mask
        ]
        .head(
            10
        )
        .to_string(
            index=False
        )
    )


## ================================================================== ##
## Taxonomy
## ================================================================== ##

heading(
    "10. TAXONOMY TABLE"
)


taxonomy = pd.read_csv(
    TAXONOMY_FILE,
    sep="\t",
    dtype=str,
    keep_default_na=False
)


print(
    f"File: {TAXONOMY_FILE}"
)


print(
    f"Rows: {len(taxonomy):,}"
)


print()
print(
    "Columns:"
)


print_columns(
    taxonomy
)


## ================================================================== ##
## Priority-taxon search
##
## Search all fields so we do not assume taxonomy column names.
## ================================================================== ##

heading(
    "11. PRIORITY TAXA PRESENT IN TAXONOMY TABLE"
)


taxonomy_text = (
    taxonomy
    .astype(
        str
    )
    .agg(
        " | ".join,
        axis=1
    )
)


for taxon in PRIORITY_TAXA:

    mask = taxonomy_text.str.contains(
        taxon,
        case=False,
        regex=False
    )


    hits = taxonomy[
        mask
    ]


    print()
    print(
        f"{taxon}: {len(hits):,} taxonomy rows"
    )


    if len(
        hits
    ) > 0:

        print(
            hits
            .head(
                10
            )
            .to_string(
                index=False
            )
        )


## ================================================================== ##
## Priority taxa among Cluster_00035 focal genomes
##
## If a shared genome column can be identified, also count how many
## focal regions belong to each target taxon.
## ================================================================== ##

heading(
    "12. PRIORITY TAXA AMONG CLUSTER_00035 FOCAL GENOMES"
)


taxonomy_genome_candidates = [
    column
    for column in taxonomy.columns
    if column.lower() in {
        "genome",
        "genome_id",
        "accession",
        "assembly",
        "assembly_accession",
    }
]


print(
    "Possible taxonomy genome-key columns: "
    +
    (
        ", ".join(
            taxonomy_genome_candidates
        )
        if taxonomy_genome_candidates
        else
        "[none detected automatically]"
    )
)


for key in taxonomy_genome_candidates:

    matched = taxonomy[
        taxonomy[
            key
        ].isin(
            set(
                module20_focals[
                    "genome"
                ]
            )
        )
    ]


    print()
    print(
        f"Using taxonomy key '{key}': "
        f"{matched[key].nunique():,} focal genomes matched"
    )


    if len(
        matched
    ) == 0:
        continue


    matched_text = (
        matched
        .astype(
            str
        )
        .agg(
            " | ".join,
            axis=1
        )
    )


    for taxon in PRIORITY_TAXA:

        n = int(
            matched_text.str.contains(
                taxon,
                case=False,
                regex=False
            ).sum()
        )


        print(
            f"  {taxon:<15} {n}"
        )


## ================================================================== ##
## Final
## ================================================================== ##

heading(
    "AUDIT COMPLETE"
)


print(
    "The next production script will use these exact schema fields "
    "to build every Cluster_00035-centered Module-20 region before "
    "any representative-genome selection."
)
