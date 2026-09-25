#!/usr/bin/env python3

from pathlib import Path
import re
import sys

import numpy as np
import pandas as pd


## ================================================================== ##
## STAGE 16A - MODULE 35 GENE-MAP DATA PREPARATION
##
## Purpose
## -------
## Prepare a standardized, plot-ready gene table for all genomes
## containing the Module-35 Cyc2 family Cluster_00206.
##
## The resulting coordinates are transcription-normalized so that
## Cluster_00206 points from left to right in every plotted genome.
##
## Every actual Prodigal gene in the requested genomic window is
## retained. Visual importance is encoded separately in display_class.
##
## Plot window:
##
##     focal Cyc2 CDS +/- 10,000 bp
##
## Taxonomy:
##
##     GTDB-style prefixes such as g__, f__, s__ are removed.
##
## The accession/genome ID remains in the data for reproducibility,
## but the intended visible row label is taxonomy_display.
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


OUT_DATA = (
    HERE
    / "Module_35_gene_map_data.tsv"
)


OUT_GENOMES = (
    HERE
    / "Module_35_gene_map_genomes.tsv"
)


OUT_QC = (
    HERE
    / "Module_35_gene_map_qc.tsv"
)


## ================================================================== ##
## Locked biological definition
## ================================================================== ##

TARGET_MODULE = "Module_35"

FOCAL_CLUSTER = "Cluster_00206"

EXPECTED_GENOMES = 8

WINDOW_BP = 10_000


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

    if not path.is_file():

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


def bool01(series):

    return (
        pd.to_numeric(
            series.replace(
                "",
                "0",
            ),
            errors="coerce",
        )
        .fillna(
            0
        )
        .astype(int)
    )


## ================================================================== ##
## Taxonomy parsing
## ================================================================== ##

RANK_PREFIXES = {
    "d": "domain",
    "p": "phylum",
    "c": "class",
    "o": "order",
    "f": "family",
    "g": "genus",
    "s": "species",
}


def strip_rank_prefix(value):

    value = clean(
        value
    )

    return re.sub(
        r"^[A-Za-z]__",
        "",
        value,
    )


def parse_taxonomy_string(taxonomy):

    result = {
        "taxonomy_domain": "",
        "taxonomy_phylum": "",
        "taxonomy_class": "",
        "taxonomy_order": "",
        "taxonomy_family": "",
        "taxonomy_genus": "",
        "taxonomy_species": "",
    }


    for token in clean(
        taxonomy
    ).split(";"):

        token = token.strip()

        match = re.match(
            r"^([dpcofgs])__(.*)$",
            token,
        )


        if match is None:

            continue


        rank_code = match.group(
            1
        )


        value = strip_rank_prefix(
            token
        )


        rank_name = RANK_PREFIXES[
            rank_code
        ]


        result[
            f"taxonomy_{rank_name}"
        ] = value


    return result


def make_taxonomy_display(
    genus,
    species,
):

    genus = clean(
        genus
    )

    species = clean(
        species
    )


    ## -------------------------------------------------------------- ##
    ## Prefer the species-level designation if it already contains
    ## the genus, e.g.
    ##
    ##     Methylobacter tundripaludum
    ##
    ## Otherwise combine genus + species-level value.
    ## -------------------------------------------------------------- ##

    if species != "":

        if (
            genus != ""
            and
            species.lower().startswith(
                genus.lower()
                +
                " "
            )
        ):

            return species


        if genus != "":

            return (
                genus
                +
                " "
                +
                species
            )


        return species


    if genus != "":

        return genus


    return "Unclassified"


## ================================================================== ##
## Read taxonomy
## ================================================================== ##

print("=" * 80)
print("STAGE 16A - MODULE 35 GENE-MAP DATA")
print("=" * 80)


print()
print("Reading taxonomy...")


require_file(
    TAXONOMY_FILE
)


taxonomy_raw = pd.read_csv(
    TAXONOMY_FILE,
    sep="\t",
    header=None,
    names=[
        "genome",
        "taxonomy_string",
    ],
    dtype=str,
    keep_default_na=False,
)


if taxonomy_raw[
    "genome"
].duplicated().any():

    duplicated = (
        taxonomy_raw.loc[
            taxonomy_raw[
                "genome"
            ].duplicated(
                keep=False
            ),
            "genome",
        ]
        .unique()
        .tolist()
    )

    fail(
        "Duplicate genome IDs in taxonomy file. "
        f"Examples: {duplicated[:10]}"
    )


taxonomy_parsed = taxonomy_raw[
    "taxonomy_string"
].apply(
    parse_taxonomy_string
)


taxonomy_parsed = pd.DataFrame(
    taxonomy_parsed.tolist()
)


taxonomy = pd.concat(
    [
        taxonomy_raw.reset_index(
            drop=True
        ),
        taxonomy_parsed.reset_index(
            drop=True
        ),
    ],
    axis=1,
)


taxonomy[
    "taxonomy_display"
] = taxonomy.apply(
    lambda row:
        make_taxonomy_display(
            row[
                "taxonomy_genus"
            ],
            row[
                "taxonomy_species"
            ],
        ),
    axis=1,
)


print(
    f"  Taxonomy entries: "
    f"{len(taxonomy):,}"
)


## ================================================================== ##
## Read Stage 15A
## ================================================================== ##

print()
print("Reading Stage-15A focal catalogue...")


focals = read_tsv(
    FOCAL_FILE
)


print(
    f"  Focal rows: "
    f"{len(focals):,}"
)


## ================================================================== ##
## Identify the eight Module-35 Cyc2 focal genes
## ================================================================== ##

module35_focals = focals[
    (
        focals[
            "module"
        ]
        ==
        TARGET_MODULE
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


if len(
    module35_focals
) != EXPECTED_GENOMES:

    fail(
        f"Expected {EXPECTED_GENOMES} "
        f"{FOCAL_CLUSTER} focal rows in {TARGET_MODULE}; "
        f"found {len(module35_focals)}."
    )


if module35_focals[
    "genome"
].nunique() != EXPECTED_GENOMES:

    fail(
        "The Module-35 Cyc2 focal rows do not represent "
        "eight distinct genomes."
    )


if module35_focals[
    "focal_id"
].duplicated().any():

    fail(
        "Duplicate focal_id among Module-35 Cyc2 focals."
    )


print()
print(
    f"Identified {len(module35_focals)} "
    f"{FOCAL_CLUSTER} / Cyc2 focal genes."
)


## ================================================================== ##
## Read fresh Stage-15A neighborhoods
## ================================================================== ##

print()
print("Reading fresh gene-centered neighborhoods...")


neigh = read_tsv(
    NEIGHBORHOOD_FILE
)


print(
    f"  Fresh neighborhood rows: "
    f"{len(neigh):,}"
)


## ================================================================== ##
## Retain only the eight selected focal neighborhoods
## ================================================================== ##

selected_focal_ids = set(
    module35_focals[
        "focal_id"
    ]
)


plot = neigh[
    neigh[
        "focal_id"
    ].isin(
        selected_focal_ids
    )
].copy()


if plot[
    "focal_id"
].nunique() != EXPECTED_GENOMES:

    fail(
        "Not all eight selected focal neighborhoods were recovered."
    )


## ================================================================== ##
## Numeric coordinates
## ================================================================== ##

coordinate_columns = [
    "focal_start",
    "focal_end",
    "neighbor_start",
    "neighbor_end",
    "neighbor_gene_rank",
]


for column in coordinate_columns:

    if column not in plot.columns:

        fail(
            f"Required neighborhood column missing: {column}"
        )


    plot[
        column
    ] = pd.to_numeric(
        plot[
            column
        ],
        errors="raise",
    )


## ================================================================== ##
## Restrict to focal CDS +/- 10 kb
##
## Genomic interval:
##
##     focal_start - 10 kb
##            through
##     focal_end   + 10 kb
##
## Any CDS intersecting that interval is retained.
## ================================================================== ##

print()
print(
    f"Restricting each neighborhood to focal CDS "
    f"+/- {WINDOW_BP:,} bp..."
)


plot[
    "plot_genomic_left"
] = (
    plot[
        "focal_start"
    ]
    -
    WINDOW_BP
)


plot[
    "plot_genomic_right"
] = (
    plot[
        "focal_end"
    ]
    +
    WINDOW_BP
)


plot = plot[
    (
        plot[
            "neighbor_end"
        ]
        >=
        plot[
            "plot_genomic_left"
        ]
    )
    &
    (
        plot[
            "neighbor_start"
        ]
        <=
        plot[
            "plot_genomic_right"
        ]
    )
].copy()


## ================================================================== ##
## Transcription-normalized plotting coordinates
##
## All focal Cyc2 genes point left -> right.
##
## We align the focal gene midpoint to x = 0.
##
##
## Focal on genomic + strand:
##
##     normalized x = genomic x - focal midpoint
##
##
## Focal on genomic - strand:
##
##     normalized x = focal midpoint - genomic x
##
## and gene direction is flipped.
## ================================================================== ##

print()
print("Normalizing all regions to focal transcription direction...")


plot[
    "focal_midpoint"
] = (
    plot[
        "focal_start"
    ]
    +
    plot[
        "focal_end"
    ]
) / 2.0


focal_plus = (
    plot[
        "focal_strand"
    ]
    ==
    "+"
)


plot[
    "plot_start_bp"
] = np.where(
    focal_plus,

    plot[
        "neighbor_start"
    ]
    -
    plot[
        "focal_midpoint"
    ],

    plot[
        "focal_midpoint"
    ]
    -
    plot[
        "neighbor_end"
    ],
)


plot[
    "plot_end_bp"
] = np.where(
    focal_plus,

    plot[
        "neighbor_end"
    ]
    -
    plot[
        "focal_midpoint"
    ],

    plot[
        "focal_midpoint"
    ]
    -
    plot[
        "neighbor_start"
    ],
)


## -------------------------------------------------------------- ##
## Display strand after normalization.
##
## If focal is +:
##     preserve original gene strand.
##
## If focal is -:
##     flip all strands.
## -------------------------------------------------------------- ##

plot[
    "plot_strand"
] = np.where(
    focal_plus,

    plot[
        "neighbor_strand"
    ],

    np.where(
        plot[
            "neighbor_strand"
        ]
        ==
        "+",
        "-",
        "+",
    ),
)


## ================================================================== ##
## Join taxonomy
## ================================================================== ##

taxonomy_columns = [
    "genome",
    "taxonomy_domain",
    "taxonomy_phylum",
    "taxonomy_class",
    "taxonomy_order",
    "taxonomy_family",
    "taxonomy_genus",
    "taxonomy_species",
    "taxonomy_display",
]


plot = plot.merge(
    taxonomy[
        taxonomy_columns
    ],
    on="genome",
    how="left",
    validate="many_to_one",
)


taxonomy_missing = (
    plot["taxonomy_display"].isna()
    |
    plot["taxonomy_display"].eq("")
)

if taxonomy_missing.any():

    missing_taxonomy = sorted(
        plot.loc[
            taxonomy_missing,
            "genome"
        ].unique()
    )

    fail(
        "Missing taxonomy for plotted genomes: "
        + ", ".join(missing_taxonomy)
    )


## ================================================================== ##
## Clean biologically useful fields
## ================================================================== ##

for column in [
    "neighbor_mmseq_cluster",
    "neighbor_mcl_module",
    "neighbor_fegenie_HMMs",
    "neighbor_globdb_cog",
    "neighbor_globdb_product",
    "neighbor_report_topology",
]:

    if column not in plot.columns:

        plot[
            column
        ] = ""


for column in [
    "neighbor_fegenie_positive",
    "neighbor_findmehemes_positive",
]:

    if column not in plot.columns:

        plot[
            column
        ] = "0"


plot[
    "neighbor_fegenie_positive"
] = bool01(
    plot[
        "neighbor_fegenie_positive"
    ]
)


plot[
    "neighbor_findmehemes_positive"
] = bool01(
    plot[
        "neighbor_findmehemes_positive"
    ]
)


## ================================================================== ##
## Plotting classes
##
## Priority matters.
##
## A Module-35 member should remain a module member even if it is also
## FeGenie-positive.
## ================================================================== ##

is_focal_family = (
    plot[
        "neighbor_mmseq_cluster"
    ]
    ==
    FOCAL_CLUSTER
)


same_module = (
    plot[
        "neighbor_mcl_module"
    ]
    ==
    TARGET_MODULE
)


fegenie_positive = (
    plot[
        "neighbor_fegenie_positive"
    ]
    ==
    1
)


findmehemes_positive = (
    plot[
        "neighbor_findmehemes_positive"
    ]
    ==
    1
)


clustered = (
    plot[
        "neighbor_mmseq_cluster"
    ]
    !=
    ""
)


plot[
    "display_class"
] = np.select(
    [
        is_focal_family,

        same_module,

        fegenie_positive,

        findmehemes_positive,

        clustered,
    ],
    [
        "focal_cyc2",

        "same_module_family",

        "fegenie_other",

        "findmehemes_other",

        "other_mmseq_family",
    ],
    default="background_gene",
)


## ================================================================== ##
## Human-readable labels
##
## Keep plotting labels conservative.
##
## We label:
##
##   Cluster_00206 -> Cyc2
##   Cluster_00228 -> CytC5-like
##
## Other FeGenie hits use FeGenie's HMM label.
##
## Other genes are initially unlabeled unless they are important enough
## to warrant a label in later figure refinement.
## ================================================================== ##

def make_display_label(row):

    cluster = clean(
        row["neighbor_mmseq_cluster"]
    )

    hmm = clean(
        row["neighbor_fegenie_HMMs"]
    )

    product = clean(
        row["neighbor_globdb_product"]
    )


    if cluster == "Cluster_00206":
        return "Cyc2"

    if cluster == "Cluster_00228":
        return "CytC5-like"

    if cluster == "Cluster_00966":
        return "Cyc1"

    if cluster == "Cluster_00176":
        return "CccA-like"

    if cluster == "Cluster_00296":
        return "7-heme cyt. c"

    if cluster == "Cluster_00033":
        return "TrxA-like"

    if hmm != "":
        return hmm

    return ""

plot[
    "display_label"
] = plot.apply(
    make_display_label,
    axis=1,
)


## ================================================================== ##
## Focal flag
## ================================================================== ##

plot[
    "is_plot_focal_gene"
] = (
    plot[
        "neighbor_protein_id"
    ]
    ==
    plot[
        "focal_protein_id"
    ]
).astype(int)


## ================================================================== ##
## Module-35 completeness per genome
## ================================================================== ##

module_clusters_per_genome = (
    plot.loc[
        plot[
            "neighbor_mcl_module"
        ]
        ==
        TARGET_MODULE,
        [
            "genome",
            "neighbor_mmseq_cluster",
        ],
    ]
    .drop_duplicates()
    .groupby(
        "genome"
    )[
        "neighbor_mmseq_cluster"
    ]
    .nunique()
)


plot[
    "module35_clusters_in_plot"
] = plot[
    "genome"
].map(
    module_clusters_per_genome
).fillna(
    0
).astype(int)


plot[
    "module35_status"
] = np.where(
    plot[
        "module35_clusters_in_plot"
    ]
    >=
    2,
    "full",
    "partial",
)


## ================================================================== ##
## Sort rows for plotting
## ================================================================== ##

plot = plot.sort_values(
    [
        "taxonomy_genus",
        "taxonomy_species",
        "genome",
        "plot_start_bp",
        "neighbor_protein_id",
    ],
    kind="stable",
).reset_index(
    drop=True
)


## ================================================================== ##
## Genome-level metadata
## ================================================================== ##

genome_table = (
    plot[
        [
            "genome",
            "taxonomy_family",
            "taxonomy_genus",
            "taxonomy_species",
            "taxonomy_display",
            "module35_status",
            "module35_clusters_in_plot",
        ]
    ]
    .drop_duplicates()
    .sort_values(
        [
            "taxonomy_genus",
            "taxonomy_species",
            "genome",
        ]
    )
    .reset_index(
        drop=True
    )
)


## ================================================================== ##
## QC
## ================================================================== ##

print()
print("Running Module-35 gene-map QC...")


n_genomes = plot[
    "genome"
].nunique()


n_focals = int(
    plot[
        "is_plot_focal_gene"
    ].sum()
)


focals_per_genome = (
    plot[
        plot[
            "is_plot_focal_gene"
        ]
        ==
        1
    ]
    .groupby(
        "genome"
    )
    .size()
)


wrong_focal_direction = plot[
    (
        plot[
            "is_plot_focal_gene"
        ]
        ==
        1
    )
    &
    (
        plot[
            "plot_strand"
        ]
        !=
        "+"
    )
]


duplicated_gene_rows = (
    plot[
        [
            "focal_id",
            "neighbor_protein_id",
        ]
    ]
    .duplicated()
    .sum()
)


qc_rows = [
    {
        "metric":
            "number_of_module35_genomes",

        "observed":
            n_genomes,

        "expected":
            EXPECTED_GENOMES,

        "pass":
            int(
                n_genomes
                ==
                EXPECTED_GENOMES
            ),
    },

    {
        "metric":
            "number_of_plot_focal_genes",

        "observed":
            n_focals,

        "expected":
            EXPECTED_GENOMES,

        "pass":
            int(
                n_focals
                ==
                EXPECTED_GENOMES
            ),
    },

    {
        "metric":
            "genomes_with_exactly_one_focal",

        "observed":
            int(
                (
                    focals_per_genome
                    ==
                    1
                )
                .sum()
            ),

        "expected":
            EXPECTED_GENOMES,

        "pass":
            int(
                (
                    focals_per_genome
                    ==
                    1
                )
                .sum()
                ==
                EXPECTED_GENOMES
            ),
    },

    {
        "metric":
            "focal_genes_pointing_right",

        "observed":
            n_focals
            -
            len(
                wrong_focal_direction
            ),

        "expected":
            EXPECTED_GENOMES,

        "pass":
            int(
                len(
                    wrong_focal_direction
                )
                ==
                0
            ),
    },

    {
        "metric":
            "duplicate_focal_neighbor_rows",

        "observed":
            int(
                duplicated_gene_rows
            ),

        "expected":
            0,

        "pass":
            int(
                duplicated_gene_rows
                ==
                0
            ),
    },

    {
        "metric":
            "taxonomy_complete",

        "observed":
            int(
                plot[
                    "taxonomy_display"
                ]
                .ne(
                    ""
                )
                .all()
            ),

        "expected":
            1,

        "pass":
            int(
                plot[
                    "taxonomy_display"
                ]
                .ne(
                    ""
                )
                .all()
            ),
    },
]


qc = pd.DataFrame(
    qc_rows
)


if (
    qc[
        "pass"
    ]
    !=
    1
).any():

    print()
    print(
        qc.to_string(
            index=False
        )
    )

    fail(
        "At least one Module-35 gene-map QC test failed."
    )


## ================================================================== ##
## Output
## ================================================================== ##

write_tsv(
    plot,
    OUT_DATA,
)


write_tsv(
    genome_table,
    OUT_GENOMES,
)


write_tsv(
    qc,
    OUT_QC,
)


## ================================================================== ##
## Summary
## ================================================================== ##

print()
print("Module 35 gene-map summary")


print(
    f"  Genomes:                 "
    f"{n_genomes}"
)


print(
    f"  Plot genes:              "
    f"{len(plot):,}"
)


print(
    f"  Focal Cyc2 genes:        "
    f"{n_focals}"
)


print(
    f"  Module-35 genes:         "
    f"{int(same_module.loc[plot.index].sum()) if len(plot) else 0}"
)


print()
print("Genome labels")


print(
    genome_table[
        [
            "taxonomy_display",
            "taxonomy_family",
            "module35_status",
        ]
    ]
    .to_string(
        index=False
    )
)


print()
print("Display classes")


print(
    plot[
        "display_class"
    ]
    .value_counts()
    .to_string()
)


print()
print("=" * 80)
print("SUCCESS")
print("=" * 80)


print(
    f"Plot-ready genes: "
    f"{OUT_DATA}"
)


print(
    f"Genome metadata:  "
    f"{OUT_GENOMES}"
)


print(
    f"QC:               "
    f"{OUT_QC}"
)
