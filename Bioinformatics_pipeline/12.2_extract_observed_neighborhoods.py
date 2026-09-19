#!/usr/bin/env python3

from pathlib import Path
from itertools import combinations
import csv
import sys

import numpy as np
import pandas as pd


## ================================================================== ##
## STAGE 12A
##
## Extract complete observed local genomic context around every
## MCL-module protein.
##
## This stage is descriptive only:
##
##   - no conservation thresholds
##   - no consensus neighborhoods
##   - no functional inference
##   - no removal of unannotated genes
##
## A gene is included in a focal neighborhood if it is:
##
##   within +/- 20 genes
##
##        OR
##
##   intersects the focal CDS +/- 20,000 bp
##
## This deliberately preserves a broad raw neighborhood from which
## smaller windows can later be derived without rerunning extraction.
## ================================================================== ##


## ================================================================== ##
## Paths
## ================================================================== ##

HERE = Path(__file__).resolve().parent

WORKFLOW = HERE.parent


GENE_CATALOG = (
    WORKFLOW
    / "11_genome_annotations"
    / "annotated_gene_catalog.tsv"
)


OUT_FOCALS = (
    HERE
    / "focal_module_occurrences.tsv"
)


OUT_NEIGHBORHOODS = (
    HERE
    / "observed_neighborhood_genes.tsv"
)


OUT_MODULE_PAIRS = (
    HERE
    / "module_member_pair_observations.tsv"
)


OUT_QC = (
    HERE
    / "neighborhood_extraction_qc.tsv"
)


## ================================================================== ##
## Locked Stage-12A extraction parameters
## ================================================================== ##

GENE_RADIUS = 20

BP_RADIUS = 20_000

CHUNK_SIZE = 250_000


## ================================================================== ##
## Expected upstream values
## ================================================================== ##

EXPECTED_TOTAL_GENES = 2_002_656

EXPECTED_GENOMES = 631

EXPECTED_MODULE_PROTEINS = 10_537

EXPECTED_MODULES = 35


## ================================================================== ##
## Utilities
## ================================================================== ##

def fail(message):

    print(
        f"\nERROR: {message}",
        file=sys.stderr
    )

    sys.exit(1)


def require_columns(
    available,
    required,
    source,
):

    missing = (
        set(required)
        -
        set(available)
    )

    if missing:

        fail(
            f"{source} is missing required "
            f"column(s):\n"
            +
            ", ".join(
                sorted(missing)
            )
        )


def ordered_unique(values):

    seen = set()

    result = []

    for value in values:

        if value in seen:
            continue

        seen.add(value)

        result.append(value)

    return result


def flip_strand(strand):

    if strand == "+":
        return "-"

    if strand == "-":
        return "+"

    return strand


def cds_relationship(
    a_start,
    a_end,
    b_start,
    b_end,
):

    """
    Compare two inclusive CDS intervals.

    Returns:

        intergenic_gap_bp
            Number of nucleotides lying strictly between the CDSs.
            Zero if they overlap or directly touch.

        overlap_nt
            Number of overlapping nucleotides.
            Zero for non-overlapping genes.
    """

    if a_end < b_start:

        gap = max(
            0,
            b_start
            -
            a_end
            -
            1
        )

        return (
            gap,
            0,
        )


    if b_end < a_start:

        gap = max(
            0,
            a_start
            -
            b_end
            -
            1
        )

        return (
            gap,
            0,
        )


    overlap = (
        min(
            a_end,
            b_end
        )
        -
        max(
            a_start,
            b_start
        )
        +
        1
    )


    return (
        0,
        max(
            0,
            overlap
        ),
    )


def midpoint_offset(
    neighbor_start,
    neighbor_end,
    focal_start,
    focal_end,
):

    """
    Signed midpoint displacement in genomic coordinates.

    Positive = physically to the right on the contig.
    Negative = physically to the left.
    """

    return (
        (
            neighbor_start
            +
            neighbor_end
        )
        -
        (
            focal_start
            +
            focal_end
        )
    ) / 2.0


def orientation_class(
    left_strand,
    right_strand,
):

    pattern = (
        left_strand
        +
        right_strand
    )


    if pattern == "++":

        return (
            pattern,
            "codirectional_right"
        )


    if pattern == "--":

        return (
            pattern,
            "codirectional_left"
        )


    if pattern == "+-":

        return (
            pattern,
            "convergent"
        )


    if pattern == "-+":

        return (
            pattern,
            "divergent"
        )


    return (
        pattern,
        "unknown"
    )


## ================================================================== ##
## Start
## ================================================================== ##

print("=" * 80)

print(
    "STAGE 12A - EXTRACT OBSERVED GENOMIC NEIGHBORHOODS"
)

print("=" * 80)


if not GENE_CATALOG.exists():

    fail(
        f"Cannot find gene catalog:\n"
        f"{GENE_CATALOG}"
    )


## ================================================================== ##
## Inspect available columns
## ================================================================== ##

header = pd.read_csv(
    GENE_CATALOG,
    sep="\t",
    nrows=0
)


available_columns = list(
    header.columns
)


## ================================================================== ##
## Required structural columns
## ================================================================== ##

required_columns = [
    "genome",
    "contig",
    "contig_gene_rank",
    "contig_gene_count",
    "protein_id",
    "start",
    "end",
    "strand",

    "annotation_match_type",
    "annotation_accepted",
    "globdb_cog",
    "globdb_gene",
    "globdb_product",

    "is_exported_heme_cluster_protein",
    "cluster",

    "is_mcl_module_protein",
    "module",
]


require_columns(
    available_columns,
    required_columns,
    GENE_CATALOG.name
)


## ================================================================== ##
## Useful optional evidence columns
##
## These are propagated when they exist in the catalog.
## ================================================================== ##

optional_evidence_columns = [
    "reciprocal_overlap",

    "candidate_source",

    "fegenie_positive",
    "fegenie_HMMs",
    "fegenie_categories",

    "findmehemes_positive",
    "number_of_hemes",

    "signalp_prediction",

    "deeptmhmm_class",
    "deeptmhmm_n_tm_helices",

    "export_evidence",
    "localization_class",
]


optional_evidence_columns = [
    column

    for column
    in optional_evidence_columns

    if column
    in available_columns
]


## ================================================================== ##
## PASS 1
##
## Find all 10,537 module proteins and therefore all contigs that
## contain at least one focal protein.
##
## We do this first so that PASS 2 does not need to keep all
## 2,002,656 genes in memory.
## ================================================================== ##

print()
print(
    "Pass 1: locating module proteins..."
)


pass1_columns = [
    "genome",
    "contig",
    "protein_id",
    "cluster",
    "module",
]


focal_parts = []

total_catalog_rows = 0


for chunk in pd.read_csv(
    GENE_CATALOG,
    sep="\t",
    usecols=pass1_columns,
    dtype=str,
    keep_default_na=False,
    chunksize=CHUNK_SIZE,
):

    total_catalog_rows += len(
        chunk
    )


    focal_chunk = chunk[
        chunk[
            "module"
        ].str.strip()
        != ""
    ].copy()


    if len(
        focal_chunk
    ) > 0:

        focal_parts.append(
            focal_chunk
        )


if (
    total_catalog_rows
    !=
    EXPECTED_TOTAL_GENES
):

    fail(
        f"Expected "
        f"{EXPECTED_TOTAL_GENES:,} "
        f"catalog genes but read "
        f"{total_catalog_rows:,}."
    )


focals_minimal = pd.concat(
    focal_parts,
    ignore_index=True
)


if (
    len(
        focals_minimal
    )
    !=
    EXPECTED_MODULE_PROTEINS
):

    fail(
        f"Expected "
        f"{EXPECTED_MODULE_PROTEINS:,} "
        f"module proteins but found "
        f"{len(focals_minimal):,}."
    )


if (
    focals_minimal[
        [
            "genome",
            "protein_id",
        ]
    ]
    .duplicated()
    .any()
):

    fail(
        "Duplicate module protein keys "
        "detected."
    )


if (
    focals_minimal[
        "module"
    ]
    .nunique()
    !=
    EXPECTED_MODULES
):

    fail(
        f"Expected "
        f"{EXPECTED_MODULES} modules but "
        f"found "
        f"{focals_minimal['module'].nunique()}."
    )


target_contigs = set(
    zip(
        focals_minimal[
            "genome"
        ],

        focals_minimal[
            "contig"
        ],
    )
)


target_contig_keys = {
    genome
    +
    "\t"
    +
    contig

    for (
        genome,
        contig
    )
    in target_contigs
}


print(
    f"  Module proteins: "
    f"{len(focals_minimal):,}"
)

print(
    f"  Modules:         "
    f"{focals_minimal['module'].nunique():,}"
)

print(
    f"  Target contigs:  "
    f"{len(target_contigs):,}"
)


## ================================================================== ##
## PASS 2
##
## Read all genes belonging to those focal-containing contigs.
##
## This preserves ordinary genes, unannotated genes, non-heme genes,
## and proteins outside our MMseqs clustering population.
## ================================================================== ##

print()
print(
    "Pass 2: loading genes from focal-containing contigs..."
)


context_columns = ordered_unique(
    required_columns
    +
    optional_evidence_columns
)


context_parts = []


for chunk in pd.read_csv(
    GENE_CATALOG,
    sep="\t",
    usecols=context_columns,
    dtype=str,
    keep_default_na=False,
    chunksize=CHUNK_SIZE,
):

    key = (
        chunk[
            "genome"
        ]
        +
        "\t"
        +
        chunk[
            "contig"
        ]
    )


    keep = key.isin(
        target_contig_keys
    )


    if keep.any():

        context_parts.append(
            chunk.loc[
                keep
            ].copy()
        )


context = pd.concat(
    context_parts,
    ignore_index=True
)


## ================================================================== ##
## Convert structural coordinates to integers
## ================================================================== ##

for column in [
    "contig_gene_rank",
    "contig_gene_count",
    "start",
    "end",
]:

    context[
        column
    ] = pd.to_numeric(
        context[
            column
        ],
        errors="raise"
    ).astype(
        np.int64
    )


context = (
    context
    .sort_values(
        [
            "genome",
            "contig",
            "contig_gene_rank",
        ],
        kind="stable"
    )
    .reset_index(
        drop=True
    )
)


print(
    f"  Genes retained on target contigs: "
    f"{len(context):,}"
)


## ================================================================== ##
## Strict focal recovery
## ================================================================== ##

context_focals = context[
    context[
        "module"
    ].str.strip()
    != ""
].copy()


if (
    len(
        context_focals
    )
    !=
    EXPECTED_MODULE_PROTEINS
):

    fail(
        "Pass 2 did not recover all module "
        "proteins."
    )


pass1_keys = set(
    zip(
        focals_minimal[
            "genome"
        ],
        focals_minimal[
            "protein_id"
        ],
    )
)


pass2_keys = set(
    zip(
        context_focals[
            "genome"
        ],
        context_focals[
            "protein_id"
        ],
    )
)


if (
    pass1_keys
    !=
    pass2_keys
):

    fail(
        "Module-protein key mismatch between "
        "Pass 1 and Pass 2."
    )


## ================================================================== ##
## Neighborhood output columns
## ================================================================== ##

neighborhood_columns = [
    "genome",

    "focal_protein_id",
    "focal_cluster",
    "focal_module",

    "focal_contig",
    "focal_start",
    "focal_end",
    "focal_strand",
    "focal_gene_rank",

    "neighbor_protein_id",
    "neighbor_start",
    "neighbor_end",
    "neighbor_strand",
    "neighbor_gene_rank",

    "genomic_gene_offset",
    "oriented_gene_offset",

    "genomic_midpoint_offset_bp",
    "oriented_midpoint_offset_bp",

    "intergenic_gap_bp",
    "cds_overlap_nt",

    "within_20_genes",
    "within_20kb",
    "extraction_basis",

    "neighbor_same_strand_as_focal",
    "oriented_neighbor_strand",

    "is_focal",

    "neighbor_same_cluster_as_focal",
    "neighbor_same_module_as_focal",

    "neighbor_annotation_match_type",
    "neighbor_annotation_accepted",

    "neighbor_globdb_cog",
    "neighbor_globdb_gene",
    "neighbor_globdb_product",

    "neighbor_is_exported_heme_cluster_protein",
    "neighbor_cluster",

    "neighbor_is_mcl_module_protein",
    "neighbor_module",
]


for column in optional_evidence_columns:

    neighborhood_columns.append(
        "neighbor_"
        +
        column
    )


## ================================================================== ##
## Focal occurrence output columns
## ================================================================== ##

focal_columns = [
    "genome",

    "focal_protein_id",
    "focal_cluster",
    "focal_module",

    "contig",

    "focal_start",
    "focal_end",
    "focal_strand",

    "contig_gene_rank",
    "contig_gene_count",

    "n_genes_left_available",
    "n_genes_right_available",

    "n_oriented_upstream_genes_available",
    "n_oriented_downstream_genes_available",

    "upstream_20_gene_context_complete",
    "downstream_20_gene_context_complete",

    "n_genes_in_extraction_envelope",
    "n_genes_within_20_gene_window",
    "n_genes_within_20kb_window",

    "n_other_module_proteins_in_envelope",
    "n_same_module_other_proteins_in_envelope",

    "focal_annotation_match_type",
    "focal_annotation_accepted",

    "focal_globdb_cog",
    "focal_globdb_gene",
    "focal_globdb_product",
]


for column in optional_evidence_columns:

    focal_columns.append(
        "focal_"
        +
        column
    )


## ================================================================== ##
## Extract neighborhoods
##
## Output is streamed directly to disk because the final table can
## contain hundreds of thousands of rows.
## ================================================================== ##

print()
print(
    "Extracting observed neighborhoods..."
)


focal_rows = []

neighborhood_row_count = 0

self_row_count = 0

neighborhood_rows_with_product = 0

neighborhood_rows_with_cog = 0


with open(
    OUT_NEIGHBORHOODS,
    "w",
    newline="",
) as handle:

    writer = csv.DictWriter(
        handle,
        fieldnames=neighborhood_columns,
        delimiter="\t",
        lineterminator="\n",
        extrasaction="ignore",
    )


    writer.writeheader()


    grouped = context.groupby(
        [
            "genome",
            "contig",
        ],
        sort=False
    )


    processed_focals = 0


    for (
        genome,
        contig
    ), group in grouped:

        group = (
            group
            .sort_values(
                [
                    "contig_gene_rank",
                    "start",
                    "end",
                    "protein_id",
                ],
                kind="stable"
            )
            .reset_index(
                drop=True
            )
        )


        ranks = (
            group[
                "contig_gene_rank"
            ]
            .to_numpy(
                dtype=np.int64
            )
        )


        starts = (
            group[
                "start"
            ]
            .to_numpy(
                dtype=np.int64
            )
        )


        ends = (
            group[
                "end"
            ]
            .to_numpy(
                dtype=np.int64
            )
        )


        focal_indices = np.flatnonzero(
            group[
                "module"
            ]
            .str.strip()
            .ne("")
            .to_numpy()
        )


        if (
            len(
                focal_indices
            )
            == 0
        ):

            continue


        records = group.to_dict(
            orient="records"
        )


        for focal_index in focal_indices:

            focal = records[
                focal_index
            ]


            focal_rank = int(
                focal[
                    "contig_gene_rank"
                ]
            )

            focal_gene_count = int(
                focal[
                    "contig_gene_count"
                ]
            )

            focal_start = int(
                focal[
                    "start"
                ]
            )

            focal_end = int(
                focal[
                    "end"
                ]
            )

            focal_strand = focal[
                "strand"
            ]


            if focal_strand not in [
                "+",
                "-",
            ]:

                fail(
                    f"Unexpected focal strand "
                    f"{focal_strand!r} for "
                    f"{genome} / "
                    f"{focal['protein_id']}."
                )


            ## ------------------------------------------------------ ##
            ## Extraction masks
            ## ------------------------------------------------------ ##

            genomic_offsets = (
                ranks
                -
                focal_rank
            )


            within_gene_window = (
                np.abs(
                    genomic_offsets
                )
                <=
                GENE_RADIUS
            )


            window_left = (
                focal_start
                -
                BP_RADIUS
            )

            window_right = (
                focal_end
                +
                BP_RADIUS
            )


            within_bp_window = (
                (
                    ends
                    >=
                    window_left
                )
                &
                (
                    starts
                    <=
                    window_right
                )
            )


            keep_mask = (
                within_gene_window
                |
                within_bp_window
            )


            neighbor_indices = np.flatnonzero(
                keep_mask
            )


            ## ------------------------------------------------------ ##
            ## Focal contig-edge awareness in gene space
            ## ------------------------------------------------------ ##

            n_left = (
                focal_rank
                -
                1
            )

            n_right = (
                focal_gene_count
                -
                focal_rank
            )


            if focal_strand == "+":

                n_upstream = (
                    n_left
                )

                n_downstream = (
                    n_right
                )

            else:

                n_upstream = (
                    n_right
                )

                n_downstream = (
                    n_left
                )


            other_module_count = 0

            same_module_other_count = 0


            for neighbor_index in neighbor_indices:

                if (
                    neighbor_index
                    ==
                    focal_index
                ):

                    continue


                neighbor = records[
                    neighbor_index
                ]


                if (
                    neighbor[
                        "module"
                    ].strip()
                    != ""
                ):

                    other_module_count += 1


                    if (
                        neighbor[
                            "module"
                        ]
                        ==
                        focal[
                            "module"
                        ]
                    ):

                        same_module_other_count += 1


            focal_row = {
                "genome":
                    genome,

                "focal_protein_id":
                    focal[
                        "protein_id"
                    ],

                "focal_cluster":
                    focal[
                        "cluster"
                    ],

                "focal_module":
                    focal[
                        "module"
                    ],

                "contig":
                    contig,

                "focal_start":
                    focal_start,

                "focal_end":
                    focal_end,

                "focal_strand":
                    focal_strand,

                "contig_gene_rank":
                    focal_rank,

                "contig_gene_count":
                    focal_gene_count,

                "n_genes_left_available":
                    n_left,

                "n_genes_right_available":
                    n_right,

                "n_oriented_upstream_genes_available":
                    n_upstream,

                "n_oriented_downstream_genes_available":
                    n_downstream,

                "upstream_20_gene_context_complete":
                    int(
                        n_upstream
                        >=
                        GENE_RADIUS
                    ),

                "downstream_20_gene_context_complete":
                    int(
                        n_downstream
                        >=
                        GENE_RADIUS
                    ),

                "n_genes_in_extraction_envelope":
                    int(
                        keep_mask.sum()
                    ),

                "n_genes_within_20_gene_window":
                    int(
                        within_gene_window.sum()
                    ),

                "n_genes_within_20kb_window":
                    int(
                        within_bp_window.sum()
                    ),

                "n_other_module_proteins_in_envelope":
                    other_module_count,

                "n_same_module_other_proteins_in_envelope":
                    same_module_other_count,

                "focal_annotation_match_type":
                    focal[
                        "annotation_match_type"
                    ],

                "focal_annotation_accepted":
                    focal[
                        "annotation_accepted"
                    ],

                "focal_globdb_cog":
                    focal[
                        "globdb_cog"
                    ],

                "focal_globdb_gene":
                    focal[
                        "globdb_gene"
                    ],

                "focal_globdb_product":
                    focal[
                        "globdb_product"
                    ],
            }


            for column in optional_evidence_columns:

                focal_row[
                    "focal_"
                    +
                    column
                ] = focal.get(
                    column,
                    ""
                )


            focal_rows.append(
                focal_row
            )


            ## ------------------------------------------------------ ##
            ## Write every observed neighbor
            ## ------------------------------------------------------ ##

            direction = (
                1
                if
                focal_strand == "+"
                else
                -1
            )


            for neighbor_index in neighbor_indices:

                neighbor = records[
                    neighbor_index
                ]


                neighbor_start = int(
                    neighbor[
                        "start"
                    ]
                )

                neighbor_end = int(
                    neighbor[
                        "end"
                    ]
                )

                neighbor_rank = int(
                    neighbor[
                        "contig_gene_rank"
                    ]
                )


                genomic_offset = (
                    neighbor_rank
                    -
                    focal_rank
                )


                oriented_offset = (
                    genomic_offset
                    *
                    direction
                )


                midpoint_delta = midpoint_offset(
                    neighbor_start,
                    neighbor_end,
                    focal_start,
                    focal_end,
                )


                oriented_midpoint_delta = (
                    midpoint_delta
                    *
                    direction
                )


                (
                    intergenic_gap,
                    overlap_nt,
                ) = cds_relationship(
                    focal_start,
                    focal_end,
                    neighbor_start,
                    neighbor_end,
                )


                in_gene_window = int(
                    abs(
                        genomic_offset
                    )
                    <=
                    GENE_RADIUS
                )


                in_bp_window = int(
                    (
                        neighbor_end
                        >=
                        window_left
                    )
                    and
                    (
                        neighbor_start
                        <=
                        window_right
                    )
                )


                if (
                    in_gene_window
                    and
                    in_bp_window
                ):

                    extraction_basis = (
                        "both"
                    )

                elif in_gene_window:

                    extraction_basis = (
                        "gene_window"
                    )

                elif in_bp_window:

                    extraction_basis = (
                        "bp_window"
                    )

                else:

                    fail(
                        "Internal extraction-mask "
                        "error."
                    )


                neighbor_strand = (
                    neighbor[
                        "strand"
                    ]
                )


                if focal_strand == "+":

                    oriented_neighbor_strand = (
                        neighbor_strand
                    )

                else:

                    oriented_neighbor_strand = (
                        flip_strand(
                            neighbor_strand
                        )
                    )


                is_focal = int(
                    neighbor_index
                    ==
                    focal_index
                )


                if is_focal:

                    self_row_count += 1


                same_cluster = int(
                    focal[
                        "cluster"
                    ].strip()
                    != ""
                    and
                    neighbor[
                        "cluster"
                    ]
                    ==
                    focal[
                        "cluster"
                    ]
                )


                same_module = int(
                    focal[
                        "module"
                    ].strip()
                    != ""
                    and
                    neighbor[
                        "module"
                    ]
                    ==
                    focal[
                        "module"
                    ]
                )


                row = {
                    "genome":
                        genome,

                    "focal_protein_id":
                        focal[
                            "protein_id"
                        ],

                    "focal_cluster":
                        focal[
                            "cluster"
                        ],

                    "focal_module":
                        focal[
                            "module"
                        ],

                    "focal_contig":
                        contig,

                    "focal_start":
                        focal_start,

                    "focal_end":
                        focal_end,

                    "focal_strand":
                        focal_strand,

                    "focal_gene_rank":
                        focal_rank,

                    "neighbor_protein_id":
                        neighbor[
                            "protein_id"
                        ],

                    "neighbor_start":
                        neighbor_start,

                    "neighbor_end":
                        neighbor_end,

                    "neighbor_strand":
                        neighbor_strand,

                    "neighbor_gene_rank":
                        neighbor_rank,

                    "genomic_gene_offset":
                        genomic_offset,

                    "oriented_gene_offset":
                        oriented_offset,

                    "genomic_midpoint_offset_bp":
                        f"{midpoint_delta:.1f}",

                    "oriented_midpoint_offset_bp":
                        f"{oriented_midpoint_delta:.1f}",

                    "intergenic_gap_bp":
                        intergenic_gap,

                    "cds_overlap_nt":
                        overlap_nt,

                    "within_20_genes":
                        in_gene_window,

                    "within_20kb":
                        in_bp_window,

                    "extraction_basis":
                        extraction_basis,

                    "neighbor_same_strand_as_focal":
                        int(
                            neighbor_strand
                            ==
                            focal_strand
                        ),

                    "oriented_neighbor_strand":
                        oriented_neighbor_strand,

                    "is_focal":
                        is_focal,

                    "neighbor_same_cluster_as_focal":
                        same_cluster,

                    "neighbor_same_module_as_focal":
                        same_module,

                    "neighbor_annotation_match_type":
                        neighbor[
                            "annotation_match_type"
                        ],

                    "neighbor_annotation_accepted":
                        neighbor[
                            "annotation_accepted"
                        ],

                    "neighbor_globdb_cog":
                        neighbor[
                            "globdb_cog"
                        ],

                    "neighbor_globdb_gene":
                        neighbor[
                            "globdb_gene"
                        ],

                    "neighbor_globdb_product":
                        neighbor[
                            "globdb_product"
                        ],

                    "neighbor_is_exported_heme_cluster_protein":
                        int(
                            neighbor[
                                "cluster"
                            ].strip()
                            != ""
                        ),

                    "neighbor_cluster":
                        neighbor[
                            "cluster"
                        ],

                    "neighbor_is_mcl_module_protein":
                        int(
                            neighbor[
                                "module"
                            ].strip()
                            != ""
                        ),

                    "neighbor_module":
                        neighbor[
                            "module"
                        ],
                }


                for column in optional_evidence_columns:

                    row[
                        "neighbor_"
                        +
                        column
                    ] = neighbor.get(
                        column,
                        ""
                    )


                writer.writerow(
                    row
                )


                neighborhood_row_count += 1


                if (
                    neighbor[
                        "globdb_product"
                    ].strip()
                    != ""
                ):

                    neighborhood_rows_with_product += 1


                if (
                    neighbor[
                        "globdb_cog"
                    ].strip()
                    != ""
                ):

                    neighborhood_rows_with_cog += 1


            processed_focals += 1


            if (
                processed_focals
                %
                1000
                ==
                0
            ):

                print(
                    f"  Processed "
                    f"{processed_focals:,}/"
                    f"{EXPECTED_MODULE_PROTEINS:,} "
                    f"focal proteins"
                )


## ================================================================== ##
## Focal QC
## ================================================================== ##

if (
    len(
        focal_rows
    )
    !=
    EXPECTED_MODULE_PROTEINS
):

    fail(
        f"Expected "
        f"{EXPECTED_MODULE_PROTEINS:,} "
        f"focal rows but created "
        f"{len(focal_rows):,}."
    )


if (
    self_row_count
    !=
    EXPECTED_MODULE_PROTEINS
):

    fail(
        f"Expected exactly one focal-self "
        f"neighborhood row per focal "
        f"({EXPECTED_MODULE_PROTEINS:,}) "
        f"but found "
        f"{self_row_count:,}."
    )


focal_df = pd.DataFrame(
    focal_rows,
    columns=focal_columns
)


if (
    focal_df[
        [
            "genome",
            "focal_protein_id",
        ]
    ]
    .duplicated()
    .any()
):

    fail(
        "Duplicate focal-protein rows detected."
    )


focal_df = (
    focal_df
    .sort_values(
        [
            "focal_module",
            "focal_cluster",
            "genome",
            "contig",
            "contig_gene_rank",
            "focal_protein_id",
        ],
        kind="stable"
    )
    .reset_index(
        drop=True
    )
)


focal_df.to_csv(
    OUT_FOCALS,
    sep="\t",
    index=False,
    na_rep=""
)


## ================================================================== ##
## MODULE MEMBER PAIR OBSERVATIONS
##
## One row =
##
##   one genome
##   one module
##   one unordered pair of DIFFERENT clusters
##
## where both clusters occur in that genome.
##
## Paralogue occurrences are preserved through counts, while the
## closest same-contig occurrence pair is reported explicitly.
## ================================================================== ##

print()
print(
    "Calculating module-member physical associations..."
)


module_proteins = context[
    context[
        "module"
    ].str.strip()
    != ""
].copy()


pair_rows = []


for (
    genome,
    module
), group in module_proteins.groupby(
    [
        "genome",
        "module",
    ],
    sort=False
):

    clusters_present = sorted(
        group[
            "cluster"
        ]
        .drop_duplicates()
        .tolist()
    )


    if (
        len(
            clusters_present
        )
        <
        2
    ):

        continue


    cluster_records = {
        cluster:
        group[
            group[
                "cluster"
            ]
            ==
            cluster
        ].to_dict(
            orient="records"
        )

        for cluster
        in clusters_present
    }


    for (
        cluster_a,
        cluster_b
    ) in combinations(
        clusters_present,
        2
    ):

        proteins_a = (
            cluster_records[
                cluster_a
            ]
        )

        proteins_b = (
            cluster_records[
                cluster_b
            ]
        )


        n_total_pairs = (
            len(
                proteins_a
            )
            *
            len(
                proteins_b
            )
        )


        same_contig_pairs = []


        for protein_a in proteins_a:

            for protein_b in proteins_b:

                if (
                    protein_a[
                        "contig"
                    ]
                    !=
                    protein_b[
                        "contig"
                    ]
                ):

                    continue


                rank_a = int(
                    protein_a[
                        "contig_gene_rank"
                    ]
                )

                rank_b = int(
                    protein_b[
                        "contig_gene_rank"
                    ]
                )


                rank_distance = abs(
                    rank_a
                    -
                    rank_b
                )


                genes_between = max(
                    0,
                    rank_distance
                    -
                    1
                )


                (
                    gap_bp,
                    overlap_nt,
                ) = cds_relationship(
                    int(
                        protein_a[
                            "start"
                        ]
                    ),
                    int(
                        protein_a[
                            "end"
                        ]
                    ),
                    int(
                        protein_b[
                            "start"
                        ]
                    ),
                    int(
                        protein_b[
                            "end"
                        ]
                    ),
                )


                midpoint_distance = abs(
                    midpoint_offset(
                        int(
                            protein_b[
                                "start"
                            ]
                        ),
                        int(
                            protein_b[
                                "end"
                            ]
                        ),
                        int(
                            protein_a[
                                "start"
                            ]
                        ),
                        int(
                            protein_a[
                                "end"
                            ]
                        ),
                    )
                )


                same_contig_pairs.append(
                    {
                        "protein_a":
                            protein_a,

                        "protein_b":
                            protein_b,

                        "rank_distance":
                            rank_distance,

                        "genes_between":
                            genes_between,

                        "gap_bp":
                            gap_bp,

                        "overlap_nt":
                            overlap_nt,

                        "midpoint_distance":
                            midpoint_distance,
                    }
                )


        n_same_contig_pairs = len(
            same_contig_pairs
        )


        n_adjacent_pairs = sum(
            pair[
                "rank_distance"
            ]
            ==
            1

            for pair
            in same_contig_pairs
        )


        n_within_5kb_pairs = sum(
            pair[
                "gap_bp"
            ]
            <=
            5_000

            for pair
            in same_contig_pairs
        )


        n_within_10kb_pairs = sum(
            pair[
                "gap_bp"
            ]
            <=
            10_000

            for pair
            in same_contig_pairs
        )


        n_within_20kb_pairs = sum(
            pair[
                "gap_bp"
            ]
            <=
            20_000

            for pair
            in same_contig_pairs
        )


        row = {
            "genome":
                genome,

            "module":
                module,

            "cluster_a":
                cluster_a,

            "cluster_b":
                cluster_b,

            "n_proteins_cluster_a":
                len(
                    proteins_a
                ),

            "n_proteins_cluster_b":
                len(
                    proteins_b
                ),

            "n_total_protein_pairs":
                n_total_pairs,

            "n_same_contig_protein_pairs":
                n_same_contig_pairs,

            "same_contig_fraction_of_protein_pairs":
                (
                    n_same_contig_pairs
                    /
                    n_total_pairs
                ),

            "any_same_contig":
                int(
                    n_same_contig_pairs
                    >
                    0
                ),

            "n_adjacent_protein_pairs":
                n_adjacent_pairs,

            "any_adjacent":
                int(
                    n_adjacent_pairs
                    >
                    0
                ),

            "n_within_5kb_protein_pairs":
                n_within_5kb_pairs,

            "any_within_5kb":
                int(
                    n_within_5kb_pairs
                    >
                    0
                ),

            "n_within_10kb_protein_pairs":
                n_within_10kb_pairs,

            "any_within_10kb":
                int(
                    n_within_10kb_pairs
                    >
                    0
                ),

            "n_within_20kb_protein_pairs":
                n_within_20kb_pairs,

            "any_within_20kb":
                int(
                    n_within_20kb_pairs
                    >
                    0
                ),

            "closest_same_contig_protein_a":
                "",

            "closest_same_contig_protein_b":
                "",

            "closest_contig":
                "",

            "closest_gene_rank_distance":
                "",

            "closest_genes_between":
                "",

            "closest_intergenic_gap_bp":
                "",

            "closest_cds_overlap_nt":
                "",

            "closest_midpoint_distance_bp":
                "",

            "closest_left_protein":
                "",

            "closest_left_cluster":
                "",

            "closest_left_strand":
                "",

            "closest_right_protein":
                "",

            "closest_right_cluster":
                "",

            "closest_right_strand":
                "",

            "closest_strand_pattern":
                "",

            "closest_orientation_class":
                "",
        }


        ## ---------------------------------------------------------- ##
        ## Pick closest physical occurrence pair
        ##
        ## Priority:
        ##   1. smallest intergenic gap
        ##   2. smallest midpoint distance
        ##   3. fewest genes between
        ##   4. deterministic protein IDs
        ## ---------------------------------------------------------- ##

        if same_contig_pairs:

            closest = min(
                same_contig_pairs,
                key=lambda x: (
                    x[
                        "gap_bp"
                    ],
                    x[
                        "midpoint_distance"
                    ],
                    x[
                        "genes_between"
                    ],
                    x[
                        "protein_a"
                    ][
                        "protein_id"
                    ],
                    x[
                        "protein_b"
                    ][
                        "protein_id"
                    ],
                )
            )


            protein_a = (
                closest[
                    "protein_a"
                ]
            )

            protein_b = (
                closest[
                    "protein_b"
                ]
            )


            if (
                int(
                    protein_a[
                        "contig_gene_rank"
                    ]
                )
                <
                int(
                    protein_b[
                        "contig_gene_rank"
                    ]
                )
            ):

                left = protein_a
                right = protein_b

            else:

                left = protein_b
                right = protein_a


            (
                strand_pattern,
                relative_orientation,
            ) = orientation_class(
                left[
                    "strand"
                ],
                right[
                    "strand"
                ],
            )


            row.update(
                {
                    "closest_same_contig_protein_a":
                        protein_a[
                            "protein_id"
                        ],

                    "closest_same_contig_protein_b":
                        protein_b[
                            "protein_id"
                        ],

                    "closest_contig":
                        protein_a[
                            "contig"
                        ],

                    "closest_gene_rank_distance":
                        closest[
                            "rank_distance"
                        ],

                    "closest_genes_between":
                        closest[
                            "genes_between"
                        ],

                    "closest_intergenic_gap_bp":
                        closest[
                            "gap_bp"
                        ],

                    "closest_cds_overlap_nt":
                        closest[
                            "overlap_nt"
                        ],

                    "closest_midpoint_distance_bp":
                        f"{closest['midpoint_distance']:.1f}",

                    "closest_left_protein":
                        left[
                            "protein_id"
                        ],

                    "closest_left_cluster":
                        left[
                            "cluster"
                        ],

                    "closest_left_strand":
                        left[
                            "strand"
                        ],

                    "closest_right_protein":
                        right[
                            "protein_id"
                        ],

                    "closest_right_cluster":
                        right[
                            "cluster"
                        ],

                    "closest_right_strand":
                        right[
                            "strand"
                        ],

                    "closest_strand_pattern":
                        strand_pattern,

                    "closest_orientation_class":
                        relative_orientation,
                }
            )


        pair_rows.append(
            row
        )


pair_df = pd.DataFrame(
    pair_rows
)


if len(
    pair_df
) > 0:

    if (
        pair_df[
            [
                "genome",
                "module",
                "cluster_a",
                "cluster_b",
            ]
        ]
        .duplicated()
        .any()
    ):

        fail(
            "Duplicate genome/module/cluster-pair "
            "observations detected."
        )


    pair_df = (
        pair_df
        .sort_values(
            [
                "module",
                "cluster_a",
                "cluster_b",
                "genome",
            ],
            kind="stable"
        )
        .reset_index(
            drop=True
        )
    )


pair_df.to_csv(
    OUT_MODULE_PAIRS,
    sep="\t",
    index=False,
    na_rep="",
    float_format="%.6f"
)


## ================================================================== ##
## Stage 12A QC
## ================================================================== ##

n_focal_genomes = (
    focal_df[
        "genome"
    ]
    .nunique()
)


if (
    n_focal_genomes
    !=
    EXPECTED_GENOMES
):

    fail(
        f"Expected module proteins across "
        f"{EXPECTED_GENOMES} genomes but found "
        f"{n_focal_genomes}."
    )


n_upstream_gene_truncated = int(
    (
        focal_df[
            "upstream_20_gene_context_complete"
        ]
        ==
        0
    )
    .sum()
)


n_downstream_gene_truncated = int(
    (
        focal_df[
            "downstream_20_gene_context_complete"
        ]
        ==
        0
    )
    .sum()
)


n_both_gene_truncated = int(
    (
        (
            focal_df[
                "upstream_20_gene_context_complete"
            ]
            ==
            0
        )
        &
        (
            focal_df[
                "downstream_20_gene_context_complete"
            ]
            ==
            0
        )
    )
    .sum()
)


n_pair_observations = len(
    pair_df
)


if (
    n_pair_observations
    >
    0
):

    n_pair_same_contig = int(
        pair_df[
            "any_same_contig"
        ]
        .sum()
    )


    n_pair_adjacent = int(
        pair_df[
            "any_adjacent"
        ]
        .sum()
    )


    n_pair_10kb = int(
        pair_df[
            "any_within_10kb"
        ]
        .sum()
    )


    n_pair_20kb = int(
        pair_df[
            "any_within_20kb"
        ]
        .sum()
    )

else:

    n_pair_same_contig = 0

    n_pair_adjacent = 0

    n_pair_10kb = 0

    n_pair_20kb = 0


qc = pd.DataFrame(
    [
        [
            "gene_radius",
            GENE_RADIUS,
        ],

        [
            "bp_radius",
            BP_RADIUS,
        ],

        [
            "total_gene_catalog_rows",
            total_catalog_rows,
        ],

        [
            "module_proteins",
            len(
                focal_df
            ),
        ],

        [
            "module_protein_genomes",
            n_focal_genomes,
        ],

        [
            "modules",
            focal_df[
                "focal_module"
            ].nunique(),
        ],

        [
            "focal_contigs",
            len(
                target_contigs
            ),
        ],

        [
            "genes_on_focal_contigs",
            len(
                context
            ),
        ],

        [
            "observed_neighborhood_rows",
            neighborhood_row_count,
        ],

        [
            "focal_self_rows",
            self_row_count,
        ],

        [
            "neighborhood_rows_with_globdb_product",
            neighborhood_rows_with_product,
        ],

        [
            "neighborhood_rows_with_cog",
            neighborhood_rows_with_cog,
        ],

        [
            "focals_without_full_20_gene_upstream_context",
            n_upstream_gene_truncated,
        ],

        [
            "focals_without_full_20_gene_downstream_context",
            n_downstream_gene_truncated,
        ],

        [
            "focals_without_full_20_gene_context_either_side",
            n_both_gene_truncated,
        ],

        [
            "genome_module_cluster_pair_observations",
            n_pair_observations,
        ],

        [
            "cluster_pair_observations_any_same_contig",
            n_pair_same_contig,
        ],

        [
            "cluster_pair_observations_any_adjacent",
            n_pair_adjacent,
        ],

        [
            "cluster_pair_observations_any_within_10kb",
            n_pair_10kb,
        ],

        [
            "cluster_pair_observations_any_within_20kb",
            n_pair_20kb,
        ],
    ],
    columns=[
        "metric",
        "value",
    ]
)


qc.to_csv(
    OUT_QC,
    sep="\t",
    index=False
)


## ================================================================== ##
## Terminal summary
## ================================================================== ##

print()
print(
    "Observed neighborhood summary"
)

print(
    f"  Focal module proteins:       "
    f"{len(focal_df):,}"
)

print(
    f"  Focal-containing contigs:    "
    f"{len(target_contigs):,}"
)

print(
    f"  Genes on focal contigs:      "
    f"{len(context):,}"
)

print(
    f"  Neighborhood rows:           "
    f"{neighborhood_row_count:,}"
)

print(
    f"  Self rows:                   "
    f"{self_row_count:,}"
)


print()
print(
    "Gene-order censoring"
)

print(
    f"  Upstream <20 genes:          "
    f"{n_upstream_gene_truncated:,}"
)

print(
    f"  Downstream <20 genes:        "
    f"{n_downstream_gene_truncated:,}"
)

print(
    f"  Both sides <20 genes:        "
    f"{n_both_gene_truncated:,}"
)


print()
print(
    "Module-member pair observations"
)

print(
    f"  Genome/module/cluster pairs: "
    f"{n_pair_observations:,}"
)

print(
    f"  Any same contig:             "
    f"{n_pair_same_contig:,}"
)

print(
    f"  Any within 20 kb:            "
    f"{n_pair_20kb:,}"
)

print(
    f"  Any within 10 kb:            "
    f"{n_pair_10kb:,}"
)

print(
    f"  Any adjacent:                "
    f"{n_pair_adjacent:,}"
)


print()
print("=" * 80)
print("SUCCESS")
print("=" * 80)


print(
    f"Focal occurrences: "
    f"{OUT_FOCALS}"
)

print(
    f"Observed neighborhoods: "
    f"{OUT_NEIGHBORHOODS}"
)

print(
    f"Module-member pairs: "
    f"{OUT_MODULE_PAIRS}"
)

print(
    f"QC: "
    f"{OUT_QC}"
)
