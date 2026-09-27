#!/usr/bin/env python3

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd


## ================================================================== ##
## STAGE 52
##
## EXTRACT OBSERVED PRIORITY-GENE NEIGHBORHOODS
##
## Purpose
## -------
##
## Extract complete observed genomic context around the 21 priority
## Methylococcales proteins:
##
##   - 6 Cyc2 proteins
##   - 15 proteins with >=5 predicted hemes
##
## A neighboring gene is retained when:
##
##     abs(gene-rank offset) <= 20
##
## OR
##
##     its CDS intersects:
##
##         focal CDS +/- 20,000 bp
##
## Every actual Prodigal gene is retained regardless of annotation.
##
## No conservation threshold is applied.
## No functional interpretation is performed.
## No new focal genes are generated.
##
## Context-observability information from Stage 51 is propagated so
## contig-boundary censoring remains explicit.
## ================================================================== ##


## ================================================================== ##
## Paths
## ================================================================== ##

PROJECT = (
    Path.home()
    / "methanotrophs"
    / "methanotroph_project"
    / "jeppe"
    / "metagenome"
)


WORK_DIR = (
    PROJECT
    / "functional_analysis"
    / "methylococcales_neighborhoods"
)


GENE_CATALOG = (
    WORK_DIR
    / "methylococcales_gene_catalog.tsv"
)


FOCALS = (
    WORK_DIR
    / "priority_focal_genes.tsv"
)


OBSERVABILITY = (
    WORK_DIR
    / "priority_context_observability.tsv"
)


OUT_FOCALS = (
    WORK_DIR
    / "priority_focal_occurrences.tsv"
)


OUT_NEIGHBORHOODS = (
    WORK_DIR
    / "priority_observed_neighborhood_genes.tsv"
)


OUT_UNIQUE_GENES = (
    WORK_DIR
    / "priority_neighborhood_unique_genes.tsv"
)


OUT_QC = (
    WORK_DIR
    / "priority_neighborhood_extraction_qc.tsv"
)


## ================================================================== ##
## Locked extraction parameters
## ================================================================== ##

GENE_RADIUS = 20

BP_RADIUS = 20_000


## ================================================================== ##
## Locked expectations from Stages 50-51
## ================================================================== ##

EXPECTED_GENES = 59_786

EXPECTED_FOCALS = 21

EXPECTED_CYC2 = 6

EXPECTED_GE5_HEME = 15


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
            f"Required file does not exist:\n"
            f"{path}"
        )


def clean(value):

    if pd.isna(value):

        return ""

    return str(
        value
    ).strip()


def flip_strand(strand):

    if strand == "+":

        return "-"

    if strand == "-":

        return "+"

    return strand


def midpoint_offset(
    neighbor_start,
    neighbor_end,
    focal_start,
    focal_end,
):

    ## Signed midpoint displacement in genomic coordinates.
    ##
    ## Positive = physically right of focal.
    ## Negative = physically left of focal.

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


def cds_relationship(
    a_start,
    a_end,
    b_start,
    b_end,
):

    ## Inclusive CDS coordinates.
    ##
    ## Returns:
    ##
    ##     intergenic_gap_bp
    ##     cds_overlap_nt

    if a_end < b_start:

        gap = max(
            0,
            b_start
            -
            a_end
            -
            1,
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
            1,
        )

        return (
            gap,
            0,
        )


    overlap = (
        min(
            a_end,
            b_end,
        )
        -
        max(
            a_start,
            b_start,
        )
        +
        1
    )


    return (
        0,
        max(
            0,
            overlap,
        ),
    )


## ================================================================== ##
## Start
## ================================================================== ##

print("=" * 80)
print("STAGE 52 - EXTRACT PRIORITY OBSERVED NEIGHBORHOODS")
print("=" * 80)


for path in (
    GENE_CATALOG,
    FOCALS,
    OBSERVABILITY,
):

    require_file(
        path
    )


## ================================================================== ##
## 1. Read complete Stage-50 gene catalogue
## ================================================================== ##

print()
print("Reading complete Methylococcales gene catalogue...")


catalog = pd.read_csv(
    GENE_CATALOG,
    sep="\t",
    dtype=str,
    keep_default_na=False,
)


required_catalog = {
    "genome",
    "protein_id",
    "contig",
    "contig_gene_rank",
    "contig_gene_count",
    "start",
    "end",
    "strand",
}


missing = (
    required_catalog
    -
    set(
        catalog.columns
    )
)


if missing:

    fail(
        "Gene catalogue is missing columns:\n"
        +
        "\n".join(
            sorted(
                missing
            )
        )
    )


if len(
    catalog
) != EXPECTED_GENES:

    fail(
        f"Expected {EXPECTED_GENES:,} genes but "
        f"found {len(catalog):,}."
    )


if catalog[
    [
        "genome",
        "protein_id",
    ]
].duplicated().any():

    fail(
        "Duplicate genome + protein_id keys "
        "in complete gene catalogue."
    )


for column in [
    "contig_gene_rank",
    "contig_gene_count",
    "start",
    "end",
]:

    catalog[
        column
    ] = pd.to_numeric(
        catalog[
            column
        ],
        errors="raise",
    ).astype(
        np.int64
    )


catalog = (
    catalog
    .sort_values(
        [
            "genome",
            "contig",
            "contig_gene_rank",
            "start",
            "end",
            "protein_id",
        ],
        kind="stable",
    )
    .reset_index(
        drop=True
    )
)


print(
    f"  Genes: {len(catalog):,}"
)


## ================================================================== ##
## 2. Read priority focal genes
## ================================================================== ##

print()
print("Reading priority focal genes...")


focals = pd.read_csv(
    FOCALS,
    sep="\t",
    dtype=str,
    keep_default_na=False,
)


required_focals = {
    "genome",
    "protein_id",
    "contig",
    "contig_gene_rank",
    "contig_gene_count",
    "start",
    "end",
    "strand",
    "priority_reason",
    "priority_cyc2",
    "priority_ge5_hemes",
    "number_of_hemes",
}


missing = (
    required_focals
    -
    set(
        focals.columns
    )
)


if missing:

    fail(
        "Priority focal table is missing columns:\n"
        +
        "\n".join(
            sorted(
                missing
            )
        )
    )


if len(
    focals
) != EXPECTED_FOCALS:

    fail(
        f"Expected {EXPECTED_FOCALS} focal genes "
        f"but found {len(focals)}."
    )


if focals[
    [
        "genome",
        "protein_id",
    ]
].duplicated().any():

    fail(
        "Duplicate focal genome + protein_id keys."
    )


for column in [
    "contig_gene_rank",
    "contig_gene_count",
    "start",
    "end",
    "number_of_hemes",
    "priority_cyc2",
    "priority_ge5_hemes",
]:

    focals[
        column
    ] = pd.to_numeric(
        focals[
            column
        ],
        errors="raise",
    ).astype(
        np.int64
    )


if int(
    focals[
        "priority_cyc2"
    ].sum()
) != EXPECTED_CYC2:

    fail(
        "Unexpected number of Cyc2 focal genes."
    )


if int(
    focals[
        "priority_ge5_hemes"
    ].sum()
) != EXPECTED_GE5_HEME:

    fail(
        "Unexpected number of >=5-heme focal genes."
    )


## ================================================================== ##
## 3. Read focal observability
## ================================================================== ##

print()
print("Reading Stage-51 observability...")


obs = pd.read_csv(
    OBSERVABILITY,
    sep="\t",
    dtype=str,
    keep_default_na=False,
)


required_obs = {
    "genome",
    "protein_id",
    "contig_length_nt",
    "genomic_left_bp_available",
    "genomic_right_bp_available",
    "oriented_upstream_bp_available",
    "oriented_downstream_bp_available",
    "genomic_left_genes_available",
    "genomic_right_genes_available",
    "oriented_upstream_genes_available",
    "oriented_downstream_genes_available",
    "full_5kb_context_observable",
    "full_10kb_context_observable",
    "full_20kb_context_observable",
    "full_5genes_context_observable",
    "full_10genes_context_observable",
    "full_20genes_context_observable",
    "full_extraction_context_observable",
}


missing = (
    required_obs
    -
    set(
        obs.columns
    )
)


if missing:

    fail(
        "Observability table is missing columns:\n"
        +
        "\n".join(
            sorted(
                missing
            )
        )
    )


if len(
    obs
) != EXPECTED_FOCALS:

    fail(
        f"Expected {EXPECTED_FOCALS} observability rows "
        f"but found {len(obs)}."
    )


if obs[
    [
        "genome",
        "protein_id",
    ]
].duplicated().any():

    fail(
        "Duplicate focal keys in observability table."
    )


## Keep only observability fields not already carried by focal table.

obs_keep = [
    "genome",
    "protein_id",
    "contig_length_nt",

    "genomic_left_bp_available",
    "genomic_right_bp_available",

    "oriented_upstream_bp_available",
    "oriented_downstream_bp_available",

    "genomic_left_genes_available",
    "genomic_right_genes_available",

    "oriented_upstream_genes_available",
    "oriented_downstream_genes_available",

    "upstream_5kb_observable",
    "downstream_5kb_observable",
    "full_5kb_context_observable",

    "upstream_10kb_observable",
    "downstream_10kb_observable",
    "full_10kb_context_observable",

    "upstream_20kb_observable",
    "downstream_20kb_observable",
    "full_20kb_context_observable",

    "upstream_5genes_observable",
    "downstream_5genes_observable",
    "full_5genes_context_observable",

    "upstream_10genes_observable",
    "downstream_10genes_observable",
    "full_10genes_context_observable",

    "upstream_20genes_observable",
    "downstream_20genes_observable",
    "full_20genes_context_observable",

    "full_extraction_context_observable",
]


obs_keep = [
    column

    for column
    in obs_keep

    if column
    in obs.columns
]


focals = focals.merge(
    obs[
        obs_keep
    ],
    on=[
        "genome",
        "protein_id",
    ],
    how="left",
    validate="one_to_one",
)


if focals[
    "contig_length_nt"
].isna().any():

    fail(
        "One or more priority focal genes lack "
        "Stage-51 observability information."
    )


## ================================================================== ##
## 4. Validate focal genes against complete catalogue
## ================================================================== ##

print()
print("Validating focal genes against complete catalogue...")


catalog_keys = set(
    zip(
        catalog[
            "genome"
        ],
        catalog[
            "protein_id"
        ],
    )
)


focal_keys = set(
    zip(
        focals[
            "genome"
        ],
        focals[
            "protein_id"
        ],
    )
)


missing_focal_keys = sorted(
    focal_keys
    -
    catalog_keys
)


if missing_focal_keys:

    fail(
        "Priority focal proteins missing from "
        "complete gene catalogue:\n"
        +
        "\n".join(
            "\t".join(
                key
            )
            for key
            in missing_focal_keys
        )
    )


## ================================================================== ##
## 5. Assign stable focal IDs
## ================================================================== ##

focals = (
    focals
    .sort_values(
        [
            "genome",
            "contig",
            "contig_gene_rank",
            "protein_id",
        ],
        kind="stable",
    )
    .reset_index(
        drop=True
    )
)


focals[
    "focal_id"
] = [
    f"Focal_{index:03d}"

    for index
    in range(
        1,
        len(
            focals
        )
        + 1,
    )
]


focal_id_lookup = {
    (
        row.genome,
        row.protein_id,
    ):
        row.focal_id

    for row
    in focals[
        [
            "genome",
            "protein_id",
            "focal_id",
        ]
    ].itertuples(
        index=False
    )
}


## ================================================================== ##
## 6. Restrict working catalogue to focal-containing contigs
##
## Ordinary genes remain present.
## ================================================================== ##

target_contigs = set(
    zip(
        focals[
            "genome"
        ],
        focals[
            "contig"
        ],
    )
)


context = catalog[
    [
        (
            genome,
            contig
        )
        in target_contigs

        for genome, contig
        in zip(
            catalog[
                "genome"
            ],
            catalog[
                "contig"
            ],
        )
    ]
].copy()


print(
    f"  Focal-containing contigs: "
    f"{len(target_contigs)}"
)

print(
    f"  Genes on focal contigs:   "
    f"{len(context):,}"
)


## ================================================================== ##
## 7. Neighborhood metadata fields to propagate
##
## Structural fields are handled explicitly below.
## These fields describe functional evidence or prior family mapping.
## ================================================================== ##

evidence_columns = [
    "is_integrated_candidate",

    "candidate_source",

    "fegenie_positive",
    "fegenie_HMMs",
    "fegenie_categories",

    "findmehemes_positive",
    "number_of_hemes",

    "signalp_prediction",
    "signalp_cs_position",

    "deeptmhmm_class",
    "deeptmhmm_n_tm_helices",

    "export_evidence",
    "localization_class",

    "fegenie_unannotated_heme_candidate",
    "fegenie_unannotated_exported_heme_candidate",
    "strong_fegenie_unannotated_soluble_periplasmic_like_candidate",

    "globdb_assignment_status",
    "globdb_n_passing_hits",
    "globdb_n_clusters_hit",

    "globdb_top_scoring_cluster",
    "globdb_top_module",

    "globdb_top_target_protein",
    "globdb_top_cluster_representative",

    "globdb_top_pident",
    "globdb_top_query_coverage",
    "globdb_top_target_coverage",
    "globdb_top_evalue",
    "globdb_top_bits",

    "globdb_all_clusters_hit",

    "priority_cyc2",
    "priority_ge5_hemes",
    "priority_reason",

    "partial",
    "start_type",
]


evidence_columns = [
    column

    for column
    in evidence_columns

    if column
    in catalog.columns
]


## ================================================================== ##
## 8. Extract neighborhoods
## ================================================================== ##

print()
print("Extracting observed focal neighborhoods...")


neighborhood_rows = []

focal_summary_rows = []

self_rows = 0


grouped = {
    key:
        group.sort_values(
            [
                "contig_gene_rank",
                "start",
                "end",
                "protein_id",
            ],
            kind="stable",
        ).reset_index(
            drop=True
        )

    for key, group
    in context.groupby(
        [
            "genome",
            "contig",
        ],
        sort=False,
    )
}


for focal_number, focal in enumerate(
    focals.to_dict(
        orient="records"
    ),
    start=1,
):

    genome = focal[
        "genome"
    ]

    focal_protein = focal[
        "protein_id"
    ]

    contig = focal[
        "contig"
    ]

    focal_id = focal[
        "focal_id"
    ]


    key = (
        genome,
        contig,
    )


    if key not in grouped:

        fail(
            f"Focal contig absent from working context:\n"
            f"{genome}\t{contig}"
        )


    group = grouped[
        key
    ]


    matches = np.flatnonzero(
        (
            group[
                "protein_id"
            ]
            ==
            focal_protein
        )
        .to_numpy()
    )


    if len(
        matches
    ) != 1:

        fail(
            f"Expected exactly one focal match for "
            f"{genome} / {focal_protein}; "
            f"found {len(matches)}."
        )


    focal_index = int(
        matches[0]
    )


    focal_record = group.iloc[
        focal_index
    ]


    focal_rank = int(
        focal_record[
            "contig_gene_rank"
        ]
    )

    focal_gene_count = int(
        focal_record[
            "contig_gene_count"
        ]
    )

    focal_start = int(
        focal_record[
            "start"
        ]
    )

    focal_end = int(
        focal_record[
            "end"
        ]
    )

    focal_strand = focal_record[
        "strand"
    ]


    if focal_strand not in {
        "+",
        "-",
    }:

        fail(
            f"Unexpected focal strand "
            f"{focal_strand!r} for "
            f"{genome} / {focal_protein}"
        )


    ranks = group[
        "contig_gene_rank"
    ].to_numpy(
        dtype=np.int64
    )


    starts = group[
        "start"
    ].to_numpy(
        dtype=np.int64
    )


    ends = group[
        "end"
    ].to_numpy(
        dtype=np.int64
    )


    ## -------------------------------------------------------------- ##
    ## Extraction masks
    ## -------------------------------------------------------------- ##

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


    ## -------------------------------------------------------------- ##
    ## Gene-space edge availability
    ## -------------------------------------------------------------- ##

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

        n_upstream = n_left

        n_downstream = n_right

    else:

        n_upstream = n_right

        n_downstream = n_left


    ## -------------------------------------------------------------- ##
    ## Counts for focal summary
    ## -------------------------------------------------------------- ##

    n_integrated_neighbors = 0

    n_fegenie_neighbors = 0

    n_heme_neighbors = 0

    n_globdb_clustered_neighbors = 0

    n_other_priority_focals = 0


    for neighbor_index in neighbor_indices:

        if neighbor_index == focal_index:

            continue


        neighbor = group.iloc[
            neighbor_index
        ]


        if (
            clean(
                neighbor.get(
                    "is_integrated_candidate",
                    ""
                )
            )
            ==
            "1"
        ):

            n_integrated_neighbors += 1


        if (
            clean(
                neighbor.get(
                    "fegenie_positive",
                    ""
                )
            )
            ==
            "1"
        ):

            n_fegenie_neighbors += 1


        try:

            n_hemes = int(
                float(
                    clean(
                        neighbor.get(
                            "number_of_hemes",
                            "0"
                        )
                    )
                    or
                    0
                )
            )

        except Exception:

            n_hemes = 0


        if n_hemes > 0:

            n_heme_neighbors += 1


        if clean(
            neighbor.get(
                "globdb_top_scoring_cluster",
                ""
            )
        ):

            n_globdb_clustered_neighbors += 1


        if clean(
            neighbor.get(
                "priority_reason",
                ""
            )
        ):

            n_other_priority_focals += 1


    focal_summary = dict(
        focal
    )


    focal_summary.update(
        {
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

            "n_other_integrated_candidates_in_envelope":
                n_integrated_neighbors,

            "n_other_fegenie_positive_in_envelope":
                n_fegenie_neighbors,

            "n_other_heme_positive_in_envelope":
                n_heme_neighbors,

            "n_other_globdb_clustered_in_envelope":
                n_globdb_clustered_neighbors,

            "n_other_priority_focals_in_envelope":
                n_other_priority_focals,

            "n_genes_left_available":
                n_left,

            "n_genes_right_available":
                n_right,

            "n_oriented_upstream_genes_available":
                n_upstream,

            "n_oriented_downstream_genes_available":
                n_downstream,
        }
    )


    focal_summary_rows.append(
        focal_summary
    )


    ## -------------------------------------------------------------- ##
    ## Orientation normalization
    ## -------------------------------------------------------------- ##

    direction = (
        1
        if focal_strand == "+"
        else -1
    )


    ## -------------------------------------------------------------- ##
    ## Write every observed gene
    ## -------------------------------------------------------------- ##

    for neighbor_index in neighbor_indices:

        neighbor = group.iloc[
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

        neighbor_strand = neighbor[
            "strand"
        ]


        genomic_gene_offset = (
            neighbor_rank
            -
            focal_rank
        )


        oriented_gene_offset = (
            genomic_gene_offset
            *
            direction
        )


        genomic_midpoint_offset_bp = (
            midpoint_offset(
                neighbor_start,
                neighbor_end,
                focal_start,
                focal_end,
            )
        )


        oriented_midpoint_offset_bp = (
            genomic_midpoint_offset_bp
            *
            direction
        )


        (
            intergenic_gap_bp,
            cds_overlap_nt,
        ) = cds_relationship(
            focal_start,
            focal_end,
            neighbor_start,
            neighbor_end,
        )


        in_gene_window = int(
            abs(
                genomic_gene_offset
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

            extraction_basis = "both"

        elif in_gene_window:

            extraction_basis = "gene_window"

        elif in_bp_window:

            extraction_basis = "bp_window"

        else:

            fail(
                "Internal extraction-mask error."
            )


        is_focal = int(
            neighbor_index
            ==
            focal_index
        )


        if is_focal:

            self_rows += 1


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


        same_strand = int(
            neighbor_strand
            ==
            focal_strand
        )


        row = {
            ## ------------------------------------------------------ ##
            ## Focal identity
            ## ------------------------------------------------------ ##

            "focal_id":
                focal_id,

            "genome":
                genome,

            "focal_protein_id":
                focal_protein,

            "focal_priority_reason":
                focal[
                    "priority_reason"
                ],

            "focal_priority_cyc2":
                focal[
                    "priority_cyc2"
                ],

            "focal_priority_ge5_hemes":
                focal[
                    "priority_ge5_hemes"
                ],

            "focal_number_of_hemes":
                focal[
                    "number_of_hemes"
                ],

            "focal_fegenie_HMMs":
                focal.get(
                    "fegenie_HMMs",
                    ""
                ),

            "focal_globdb_top_cluster":
                focal.get(
                    "globdb_top_scoring_cluster",
                    ""
                ),

            "focal_globdb_top_module":
                focal.get(
                    "globdb_top_module",
                    ""
                ),

            ## ------------------------------------------------------ ##
            ## Focal coordinates
            ## ------------------------------------------------------ ##

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

            ## ------------------------------------------------------ ##
            ## Focal observability
            ## ------------------------------------------------------ ##

            "focal_full_5kb_context_observable":
                focal.get(
                    "full_5kb_context_observable",
                    ""
                ),

            "focal_full_10kb_context_observable":
                focal.get(
                    "full_10kb_context_observable",
                    ""
                ),

            "focal_full_20kb_context_observable":
                focal.get(
                    "full_20kb_context_observable",
                    ""
                ),

            "focal_full_5genes_context_observable":
                focal.get(
                    "full_5genes_context_observable",
                    ""
                ),

            "focal_full_10genes_context_observable":
                focal.get(
                    "full_10genes_context_observable",
                    ""
                ),

            "focal_full_20genes_context_observable":
                focal.get(
                    "full_20genes_context_observable",
                    ""
                ),

            "focal_full_extraction_context_observable":
                focal.get(
                    "full_extraction_context_observable",
                    ""
                ),

            ## ------------------------------------------------------ ##
            ## Neighbor identity + coordinates
            ## ------------------------------------------------------ ##

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

            ## ------------------------------------------------------ ##
            ## Relative position
            ## ------------------------------------------------------ ##

            "genomic_gene_offset":
                genomic_gene_offset,

            "oriented_gene_offset":
                oriented_gene_offset,

            "genomic_midpoint_offset_bp":
                genomic_midpoint_offset_bp,

            "oriented_midpoint_offset_bp":
                oriented_midpoint_offset_bp,

            "intergenic_gap_bp":
                intergenic_gap_bp,

            "cds_overlap_nt":
                cds_overlap_nt,

            "within_20_genes":
                in_gene_window,

            "within_20kb":
                in_bp_window,

            "extraction_basis":
                extraction_basis,

            "neighbor_same_strand_as_focal":
                same_strand,

            "oriented_neighbor_strand":
                oriented_neighbor_strand,

            "is_focal":
                is_focal,
        }


        ## ---------------------------------------------------------- ##
        ## Propagate neighbor evidence
        ## ---------------------------------------------------------- ##

        for column in evidence_columns:

            row[
                "neighbor_"
                +
                column
            ] = neighbor.get(
                column,
                ""
            )


        neighborhood_rows.append(
            row
        )


    print(
        f"  {focal_number:>2}/"
        f"{len(focals)}  "
        f"{genome} / "
        f"{focal_protein}: "
        f"{int(keep_mask.sum())} genes"
    )


## ================================================================== ##
## 9. Convert outputs to data frames
## ================================================================== ##

neighborhoods = pd.DataFrame(
    neighborhood_rows
)


focal_summary = pd.DataFrame(
    focal_summary_rows
)


## ================================================================== ##
## 10. Strict neighborhood QC
## ================================================================== ##

if len(
    focal_summary
) != EXPECTED_FOCALS:

    fail(
        f"Expected {EXPECTED_FOCALS} focal summary rows "
        f"but created {len(focal_summary)}."
    )


if self_rows != EXPECTED_FOCALS:

    fail(
        f"Expected exactly one focal-self row per "
        f"focal ({EXPECTED_FOCALS}) but found "
        f"{self_rows}."
    )


self_counts = (
    neighborhoods[
        neighborhoods[
            "is_focal"
        ]
        ==
        1
    ]
    .groupby(
        "focal_id"
    )
    .size()
)


if len(
    self_counts
) != EXPECTED_FOCALS:

    fail(
        "Not every focal has a self row."
    )


if (
    self_counts
    !=
    1
).any():

    fail(
        "One or more focals have !=1 self row."
    )


## Every neighborhood gene must come from the original catalogue.

neighbor_keys = set(
    zip(
        neighborhoods[
            "genome"
        ],
        neighborhoods[
            "neighbor_protein_id"
        ],
    )
)


if not neighbor_keys.issubset(
    catalog_keys
):

    fail(
        "Neighborhood output contains genes absent "
        "from the authoritative gene catalogue."
    )


## ================================================================== ##
## 11. Unique actual genes represented anywhere in the neighborhoods
## ================================================================== ##

unique_keys = (
    neighborhoods[
        [
            "genome",
            "neighbor_protein_id",
        ]
    ]
    .drop_duplicates()
    .rename(
        columns={
            "neighbor_protein_id":
                "protein_id"
        }
    )
)


unique_genes = unique_keys.merge(
    catalog,
    on=[
        "genome",
        "protein_id",
    ],
    how="left",
    validate="one_to_one",
)


if unique_genes[
    "contig"
].isna().any():

    fail(
        "Failed to recover one or more unique "
        "neighborhood genes from catalogue."
    )


## Number of focal neighborhoods in which each gene occurs.

gene_occurrence = (
    neighborhoods[
        [
            "genome",
            "neighbor_protein_id",
            "focal_id",
        ]
    ]
    .drop_duplicates()
    .groupby(
        [
            "genome",
            "neighbor_protein_id",
        ],
        as_index=False,
    )
    .agg(
        n_priority_neighborhoods=(
            "focal_id",
            "nunique",
        )
    )
    .rename(
        columns={
            "neighbor_protein_id":
                "protein_id"
        }
    )
)


unique_genes = unique_genes.merge(
    gene_occurrence,
    on=[
        "genome",
        "protein_id",
    ],
    how="left",
    validate="one_to_one",
)


unique_genes = (
    unique_genes
    .sort_values(
        [
            "genome",
            "contig",
            "contig_gene_rank",
            "protein_id",
        ],
        kind="stable",
    )
    .reset_index(
        drop=True
    )
)


## ================================================================== ##
## 12. Sort and write outputs
## ================================================================== ##

focal_summary = (
    focal_summary
    .sort_values(
        [
            "focal_id",
        ],
        kind="stable",
    )
    .reset_index(
        drop=True
    )
)


neighborhoods = (
    neighborhoods
    .sort_values(
        [
            "focal_id",
            "oriented_gene_offset",
            "oriented_midpoint_offset_bp",
            "neighbor_protein_id",
        ],
        kind="stable",
    )
    .reset_index(
        drop=True
    )
)


focal_summary.to_csv(
    OUT_FOCALS,
    sep="\t",
    index=False,
    na_rep="",
)


neighborhoods.to_csv(
    OUT_NEIGHBORHOODS,
    sep="\t",
    index=False,
    na_rep="",
    float_format="%.1f",
)


unique_genes.to_csv(
    OUT_UNIQUE_GENES,
    sep="\t",
    index=False,
    na_rep="",
)


## ================================================================== ##
## 13. QC summary
## ================================================================== ##

n_unique_genes = len(
    unique_genes
)


n_neighborhood_rows = len(
    neighborhoods
)


n_integrated_rows = int(
    (
        neighborhoods[
            "neighbor_is_integrated_candidate"
        ]
        .astype(str)
        ==
        "1"
    ).sum()
) if (
    "neighbor_is_integrated_candidate"
    in neighborhoods.columns
) else 0


n_heme_rows = 0


if "neighbor_number_of_hemes" in neighborhoods.columns:

    n_heme_rows = int(
        (
            pd.to_numeric(
                neighborhoods[
                    "neighbor_number_of_hemes"
                ],
                errors="coerce",
            )
            .fillna(0)
            >
            0
        ).sum()
    )


n_fegenie_rows = 0


if "neighbor_fegenie_positive" in neighborhoods.columns:

    n_fegenie_rows = int(
        (
            neighborhoods[
                "neighbor_fegenie_positive"
            ]
            .astype(str)
            ==
            "1"
        ).sum()
    )


n_clustered_rows = 0


if "neighbor_globdb_top_scoring_cluster" in neighborhoods.columns:

    n_clustered_rows = int(
        (
            neighborhoods[
                "neighbor_globdb_top_scoring_cluster"
            ]
            .astype(str)
            .str.strip()
            !=
            ""
        ).sum()
    )


qc_rows = [
    (
        "gene_radius",
        GENE_RADIUS,
    ),

    (
        "bp_radius",
        BP_RADIUS,
    ),

    (
        "priority_focals",
        len(
            focal_summary
        ),
    ),

    (
        "focal_contigs",
        len(
            target_contigs
        ),
    ),

    (
        "neighborhood_rows",
        n_neighborhood_rows,
    ),

    (
        "unique_actual_neighborhood_genes",
        n_unique_genes,
    ),

    (
        "focal_self_rows",
        self_rows,
    ),

    (
        "neighborhood_rows_integrated_candidate",
        n_integrated_rows,
    ),

    (
        "neighborhood_rows_fegenie_positive",
        n_fegenie_rows,
    ),

    (
        "neighborhood_rows_heme_positive",
        n_heme_rows,
    ),

    (
        "neighborhood_rows_with_globdb_cluster",
        n_clustered_rows,
    ),

    (
        "focals_full_extraction_context",
        int(
            pd.to_numeric(
                focal_summary[
                    "full_extraction_context_observable"
                ],
                errors="coerce",
            )
            .fillna(0)
            .sum()
        ),
    ),

    (
        "focals_censored_extraction_context",
        int(
            (
                pd.to_numeric(
                    focal_summary[
                        "full_extraction_context_observable"
                    ],
                    errors="coerce",
                )
                .fillna(0)
                ==
                0
            ).sum()
        ),
    ),
]


qc = pd.DataFrame(
    qc_rows,
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
## 14. Terminal summary
## ================================================================== ##

print()
print("Neighborhood extraction summary")

print(
    f"  Priority focal genes:          "
    f"{len(focal_summary):,}"
)

print(
    f"  Focal-containing contigs:      "
    f"{len(target_contigs):,}"
)

print(
    f"  Neighborhood observation rows: "
    f"{n_neighborhood_rows:,}"
)

print(
    f"  Unique actual genes:           "
    f"{n_unique_genes:,}"
)

print(
    f"  Focal self rows:               "
    f"{self_rows:,}"
)

print(
    f"  Integrated-candidate rows:     "
    f"{n_integrated_rows:,}"
)

print(
    f"  FeGenie-positive rows:         "
    f"{n_fegenie_rows:,}"
)

print(
    f"  Heme-positive rows:            "
    f"{n_heme_rows:,}"
)

print(
    f"  GlobDB-clustered rows:         "
    f"{n_clustered_rows:,}"
)

print(
    f"  Fully observable envelopes:    "
    f"{int(pd.to_numeric(focal_summary['full_extraction_context_observable'], errors='coerce').fillna(0).sum()):,}"
    f"/{len(focal_summary)}"
)

print()
print("=" * 80)
print("STAGE 52 COMPLETE")
print("=" * 80)

print(
    f"Focal occurrences:\n  "
    f"{OUT_FOCALS}"
)

print(
    f"Observed neighborhoods:\n  "
    f"{OUT_NEIGHBORHOODS}"
)

print(
    f"Unique neighborhood genes:\n  "
    f"{OUT_UNIQUE_GENES}"
)

print(
    f"QC:\n  "
    f"{OUT_QC}"
)
