#!/usr/bin/env python3

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd


## ================================================================== ##
## STAGE 63
##
## BUILD COMBINED GLOBDB + METAGENOME DATASET FOR CLUSTER_00121
##
## Biological design
## -----------------
##
## GlobDB:
##     19 proteins belonging to old Cluster_00121
##
## Metagenome:
##     4 new 9-heme proteins assigned uniquely to Cluster_00121
##
## Total expected focal proteins:
##     23
##
##
## For GlobDB, neighborhood extraction reproduces the validated
## workflow:
##
##     +/- 20 genes
##
##          OR
##
##     CDS intersects:
##
##         focal_start - 20 kb
##         ...
##         focal_end   + 20 kb
##
##
## The output deliberately retains ACTUAL GENES.
##
## A second shared sequence-family clustering will subsequently be
## performed across all genes from all 23 neighborhoods so that
## GlobDB and metagenome loci use one common local-family vocabulary.
## ================================================================== ##


## ================================================================== ##
## Paths
## ================================================================== ##

HOME = Path.home()

ROOT = (
    HOME
    / "methanotrophs"
    / "methanotroph_project"
    / "jeppe"
)

OLD_WF = (
    ROOT
    / "genome_analysis_workflow"
)

MG_WF = (
    ROOT
    / "metagenome"
)


OLD_CLUSTER_MEMBERSHIP = (
    OLD_WF
    / "06_clustering"
    / "cluster_membership.tsv"
)


OLD_GENE_CATALOG = (
    OLD_WF
    / "11_genome_annotations"
    / "annotated_gene_catalog.tsv"
)


OLD_ORF_DIR = (
    OLD_WF
    / "FeGenie_conda_full_run_20260803"
    / "fegenie_results"
    / "ORF_calls"
)


MG_GENE_TABLE = (
    MG_WF
    / "functional_analysis"
    / "methylococcales_neighborhoods"
    / "priority_observed_neighborhood_genes_local_families_COG20.tsv"
)


MG_PROTEINS = (
    MG_WF
    / "functional_analysis"
    / "methylococcales_neighborhoods"
    / "priority_neighborhood_proteins.faa"
)


OUT_DIR = (
    MG_WF
    / "comparative_analysis"
    / "Cluster_00121"
)


OUT_DIR.mkdir(
    parents=True,
    exist_ok=True,
)


OUT_OLD_FOCALS = (
    OUT_DIR
    / "Cluster_00121_GlobDB_focal_proteins.tsv"
)


OUT_COMBINED_GENES = (
    OUT_DIR
    / "Cluster_00121_combined_neighborhood_genes.tsv"
)


OUT_PROTEIN_METADATA = (
    OUT_DIR
    / "Cluster_00121_combined_protein_metadata.tsv"
)


OUT_PROTEINS = (
    OUT_DIR
    / "Cluster_00121_combined_neighborhood_proteins.faa"
)


OUT_FOCAL_METADATA = (
    OUT_DIR
    / "Cluster_00121_focal_23_metadata.tsv"
)


OUT_FOCAL_FASTA = (
    OUT_DIR
    / "Cluster_00121_focal_23.faa"
)


OUT_QC = (
    OUT_DIR
    / "Cluster_00121_build_qc.tsv"
)


## ================================================================== ##
## Constants
## ================================================================== ##

TARGET_CLUSTER = "Cluster_00121"

TARGET_MG_GROUP = "Multiheme__Cluster_00121"

GENE_RADIUS = 20

BP_RADIUS = 20_000

EXPECTED_OLD_FOCALS = 19

EXPECTED_MG_FOCALS = 4

EXPECTED_TOTAL_FOCALS = 23


## ================================================================== ##
## Helpers
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
            f"Required file not found:\n{path}"
        )


def clean(value):

    if value is None:

        return ""

    if pd.isna(value):

        return ""

    return str(
        value
    ).strip()


def flag(value):

    value = clean(
        value
    ).lower()

    return int(
        value
        in {
            "1",
            "true",
            "yes",
        }
    )


def flip_strand(strand):

    strand = clean(
        strand
    )

    if strand == "+":
        return "-"

    if strand == "-":
        return "+"

    return strand


def resolve_column(
    columns,
    candidates,
    description,
    required=True,
):

    for candidate in candidates:

        if candidate in columns:

            return candidate


    if required:

        fail(
            f"Could not resolve {description}. "
            f"Tried: {', '.join(candidates)}"
        )


    return None


def fasta_iter(path):

    header = None

    seq = []


    with path.open(
        "r",
        encoding="utf-8",
    ) as handle:

        for raw in handle:

            line = raw.strip()

            if not line:

                continue


            if line.startswith(">"):

                if header is not None:

                    yield (
                        header,
                        "".join(
                            seq
                        ),
                    )


                header = line[
                    1:
                ]

                seq = []


            else:

                seq.append(
                    line
                )


    if header is not None:

        yield (
            header,
            "".join(
                seq
            ),
        )


def fasta_dict(path):

    result = {}


    for header, sequence in fasta_iter(
        path
    ):

        protein_id = (
            header
            .split()[0]
        )


        if protein_id in result:

            fail(
                f"Duplicate FASTA ID {protein_id} in {path}"
            )


        result[
            protein_id
        ] = sequence


    return result


def resolve_old_proteome(genome):

    exact_candidates = [
        OLD_ORF_DIR
        / f"{genome}.fa-proteins.faa",

        OLD_ORF_DIR
        / f"{genome}.fna-proteins.faa",

        OLD_ORF_DIR
        / f"{genome}.fasta-proteins.faa",
    ]


    exact = [
        path

        for path
        in exact_candidates

        if path.is_file()
    ]


    if len(
        exact
    ) == 1:

        return exact[
            0
        ]


    matches = sorted(
        OLD_ORF_DIR.glob(
            f"{genome}*-proteins.faa"
        )
    )


    if len(
        matches
    ) != 1:

        fail(
            f"Expected exactly one old ORF FASTA for {genome}; "
            f"found {len(matches)}:\n"
            +
            "\n".join(
                str(x)
                for x
                in matches
            )
        )


    return matches[
        0
    ]


## ================================================================== ##
## Input validation
## ================================================================== ##

print("=" * 80)

print(
    "STAGE 63 - BUILD COMBINED CLUSTER_00121 DATASET"
)

print("=" * 80)


for path in [
    OLD_CLUSTER_MEMBERSHIP,
    OLD_GENE_CATALOG,
    MG_GENE_TABLE,
    MG_PROTEINS,
]:

    require_file(
        path
    )


if not OLD_ORF_DIR.is_dir():

    fail(
        f"Old ORF directory does not exist:\n{OLD_ORF_DIR}"
    )


## ================================================================== ##
## 1. Identify the 19 GlobDB Cluster_00121 focal proteins
## ================================================================== ##

print()
print(
    "Identifying GlobDB Cluster_00121 proteins..."
)


membership = pd.read_csv(
    OLD_CLUSTER_MEMBERSHIP,
    sep="\t",
    dtype=str,
    keep_default_na=False,
)


required_membership = {
    "cluster",
    "genome",
    "protein_id",
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
        "Old cluster membership is missing columns: "
        +
        ", ".join(
            sorted(
                missing
            )
        )
    )


old_focals = (
    membership[
        membership[
            "cluster"
        ]
        ==
        TARGET_CLUSTER
    ]
    .copy()
    .reset_index(
        drop=True
    )
)


if len(
    old_focals
) != EXPECTED_OLD_FOCALS:

    fail(
        f"Expected {EXPECTED_OLD_FOCALS} old "
        f"{TARGET_CLUSTER} proteins; found {len(old_focals)}."
    )


if old_focals[
    [
        "genome",
        "protein_id",
    ]
].duplicated().any():

    fail(
        "Duplicate GlobDB focal protein keys."
    )


print(
    f"  GlobDB focal proteins: "
    f"{len(old_focals)}"
)

print(
    f"  GlobDB genomes:        "
    f"{old_focals['genome'].nunique()}"
)


old_focals.to_csv(
    OUT_OLD_FOCALS,
    sep="\t",
    index=False,
)


old_genomes = set(
    old_focals[
        "genome"
    ]
)


## ================================================================== ##
## 2. Resolve Stage-11 gene catalogue schema
## ================================================================== ##

print()
print(
    "Resolving GlobDB Stage-11 gene catalogue..."
)


header = pd.read_csv(
    OLD_GENE_CATALOG,
    sep="\t",
    nrows=0,
)


columns = header.columns


COL_GENOME = resolve_column(
    columns,
    [
        "genome",
    ],
    "genome",
)


COL_PROTEIN = resolve_column(
    columns,
    [
        "protein_id",
        "prodigal_protein_id",
    ],
    "protein ID",
)


COL_CONTIG = resolve_column(
    columns,
    [
        "contig",
    ],
    "contig",
)


COL_START = resolve_column(
    columns,
    [
        "start",
    ],
    "start",
)


COL_END = resolve_column(
    columns,
    [
        "end",
    ],
    "end",
)


COL_STRAND = resolve_column(
    columns,
    [
        "strand",
    ],
    "strand",
)


COL_RANK = resolve_column(
    columns,
    [
        "gene_rank",
        "contig_gene_rank",
    ],
    "gene rank",
)


optional_candidates = {
    "globdb_cog":
        [
            "globdb_cog",
            "cog",
        ],

    "globdb_gene":
        [
            "globdb_gene",
            "gene",
        ],

    "globdb_product":
        [
            "globdb_product",
            "product",
        ],

    "annotation_accepted":
        [
            "annotation_accepted",
        ],

    "old_mmseq_cluster":
        [
            "mmseq_cluster",
            "cluster",
        ],

    "old_mcl_module":
        [
            "mcl_module",
            "module",
        ],

    "fegenie_positive":
        [
            "fegenie_positive",
        ],

    "fegenie_HMMs":
        [
            "fegenie_HMMs",
            "fegenie_hmms",
        ],

    "findmehemes_positive":
        [
            "findmehemes_positive",
        ],

    "number_of_hemes":
        [
            "number_of_hemes",
            "heme_count",
            "CXXCH_count",
        ],
}


optional_columns = {}


for output_name, candidates in optional_candidates.items():

    optional_columns[
        output_name
    ] = resolve_column(
        columns,
        candidates,
        output_name,
        required=False,
    )


selected_columns = [
    COL_GENOME,
    COL_PROTEIN,
    COL_CONTIG,
    COL_START,
    COL_END,
    COL_STRAND,
    COL_RANK,
]


selected_columns += [
    value

    for value
    in optional_columns.values()

    if value is not None
]


selected_columns = list(
    dict.fromkeys(
        selected_columns
    )
)


## ================================================================== ##
## 3. Read only the 19 relevant GlobDB genomes
## ================================================================== ##

print()
print(
    "Reading the 19 focal GlobDB genomes from "
    "the complete 2,002,656-gene catalogue..."
)


catalog_parts = []


for chunk in pd.read_csv(
    OLD_GENE_CATALOG,
    sep="\t",
    dtype=str,
    keep_default_na=False,
    usecols=selected_columns,
    chunksize=200_000,
):

    keep = chunk[
        chunk[
            COL_GENOME
        ].isin(
            old_genomes
        )
    ]


    if len(
        keep
    ):

        catalog_parts.append(
            keep.copy()
        )


if not catalog_parts:

    fail(
        "No Stage-11 genes recovered for the "
        "Cluster_00121 genomes."
    )


catalog = pd.concat(
    catalog_parts,
    ignore_index=True,
)


rename = {
    COL_GENOME:
        "genome",

    COL_PROTEIN:
        "protein_id",

    COL_CONTIG:
        "contig",

    COL_START:
        "start",

    COL_END:
        "end",

    COL_STRAND:
        "strand",

    COL_RANK:
        "gene_rank",
}


for output_name, old_name in optional_columns.items():

    if old_name is not None:

        rename[
            old_name
        ] = output_name


catalog = catalog.rename(
    columns=rename
)


for field in optional_candidates:

    if field not in catalog.columns:

        catalog[
            field
        ] = ""


for field in [
    "start",
    "end",
    "gene_rank",
]:

    catalog[
        field
    ] = pd.to_numeric(
        catalog[
            field
        ],
        errors="raise",
    ).astype(int)


if catalog[
    [
        "genome",
        "protein_id",
    ]
].duplicated().any():

    fail(
        "Duplicate genome + protein_id rows "
        "in selected old gene catalogue."
    )


print(
    f"  Genes recovered: "
    f"{len(catalog):,}"
)

print(
    f"  Genomes recovered: "
    f"{catalog['genome'].nunique()}"
)


## ================================================================== ##
## 4. Validate all 19 focal proteins against Stage-11 coordinates
## ================================================================== ##

focal_keys = old_focals[
    [
        "genome",
        "protein_id",
    ]
].copy()


focal_catalog = focal_keys.merge(
    catalog,
    on=[
        "genome",
        "protein_id",
    ],
    how="left",
    validate="one_to_one",
)


if focal_catalog[
    "contig"
].isna().any():

    missing_focals = focal_catalog[
        focal_catalog[
            "contig"
        ].isna()
    ]


    fail(
        "One or more old focal proteins were not "
        "found in Stage-11:\n"
        +
        missing_focals.to_string(
            index=False
        )
    )


## ================================================================== ##
## 5. Extract the 19 GlobDB neighborhoods
## ================================================================== ##

print()
print(
    "Extracting GlobDB +/-20-gene OR +/-20-kb neighborhoods..."
)


old_neighborhood_rows = []


for focal_number, focal in enumerate(
    focal_catalog.itertuples(
        index=False
    ),
    start=1,
):

    genome = focal.genome

    focal_protein = focal.protein_id

    focal_contig = focal.contig

    focal_start = int(
        focal.start
    )

    focal_end = int(
        focal.end
    )

    focal_rank = int(
        focal.gene_rank
    )

    focal_strand = clean(
        focal.strand
    )


    if focal_strand not in {
        "+",
        "-",
    }:

        fail(
            f"Unexpected focal strand for "
            f"{focal_protein}: {focal_strand}"
        )


    focal_id = (
        f"GlobDB_C121_{focal_number:02d}"
    )


    contig_genes = (
        catalog[
            (
                catalog[
                    "genome"
                ]
                ==
                genome
            )
            &
            (
                catalog[
                    "contig"
                ]
                ==
                focal_contig
            )
        ]
        .sort_values(
            [
                "gene_rank",
                "start",
                "protein_id",
            ],
            kind="stable",
        )
        .copy()
    )


    ranks = contig_genes[
        "gene_rank"
    ].to_numpy(
        dtype=int
    )


    starts = contig_genes[
        "start"
    ].to_numpy(
        dtype=int
    )


    ends = contig_genes[
        "end"
    ].to_numpy(
        dtype=int
    )


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


    keep = (
        within_gene_window
        |
        within_bp_window
    )


    selected = contig_genes.loc[
        keep
    ].copy()


    focal_midpoint = (
        focal_start
        +
        focal_end
    ) / 2.0


    for row in selected.itertuples(
        index=False
    ):

        neighbor_start = int(
            row.start
        )

        neighbor_end = int(
            row.end
        )

        neighbor_rank = int(
            row.gene_rank
        )

        neighbor_strand = clean(
            row.strand
        )


        genomic_gene_offset = (
            neighbor_rank
            -
            focal_rank
        )


        neighbor_midpoint = (
            neighbor_start
            +
            neighbor_end
        ) / 2.0


        midpoint_delta = (
            neighbor_midpoint
            -
            focal_midpoint
        )


        if focal_strand == "+":

            oriented_gene_offset = (
                genomic_gene_offset
            )

            oriented_midpoint = (
                midpoint_delta
            )

            oriented_neighbor_strand = (
                neighbor_strand
            )


        else:

            oriented_gene_offset = (
                -genomic_gene_offset
            )

            oriented_midpoint = (
                -midpoint_delta
            )

            oriented_neighbor_strand = (
                flip_strand(
                    neighbor_strand
                )
            )


        in_gene = int(
            abs(
                genomic_gene_offset
            )
            <=
            GENE_RADIUS
        )


        in_bp = int(
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
            in_gene
            and
            in_bp
        ):

            extraction_basis = (
                "both"
            )


        elif in_gene:

            extraction_basis = (
                "gene_window"
            )


        elif in_bp:

            extraction_basis = (
                "bp_window"
            )


        else:

            fail(
                "Internal extraction-mask error."
            )


        accepted = flag(
            getattr(
                row,
                "annotation_accepted",
                "",
            )
        )


        raw_cog = clean(
            getattr(
                row,
                "globdb_cog",
                "",
            )
        )


        raw_gene = clean(
            getattr(
                row,
                "globdb_gene",
                "",
            )
        )


        raw_product = clean(
            getattr(
                row,
                "globdb_product",
                "",
            )
        )


        ## If annotation_accepted is unavailable entirely,
        ## accept existing annotation fields as legacy output.
        if (
            optional_columns[
                "annotation_accepted"
            ]
            is None
        ):

            use_annotation = (
                bool(
                    raw_cog
                    or
                    raw_gene
                    or
                    raw_product
                )
            )


        else:

            use_annotation = (
                accepted
                ==
                1
            )


        comparison_cog = (
            raw_cog
            if use_annotation
            else ""
        )


        comparison_product = (
            raw_product
            if use_annotation
            else ""
        )


        old_neighborhood_rows.append(
            {
                "source_dataset":
                    "GlobDB",

                "focal_id":
                    focal_id,

                "focal_cluster":
                    TARGET_CLUSTER,

                "genome":
                    genome,

                "focal_protein_id":
                    focal_protein,

                "focal_contig":
                    focal_contig,

                "focal_start":
                    focal_start,

                "focal_end":
                    focal_end,

                "focal_strand":
                    focal_strand,

                "focal_gene_rank":
                    focal_rank,

                "neighbor_protein_id":
                    row.protein_id,

                "neighbor_contig":
                    row.contig,

                "neighbor_start":
                    neighbor_start,

                "neighbor_end":
                    neighbor_end,

                "neighbor_strand":
                    neighbor_strand,

                "neighbor_gene_rank":
                    neighbor_rank,

                "genomic_gene_offset":
                    genomic_gene_offset,

                "oriented_gene_offset":
                    oriented_gene_offset,

                "genomic_midpoint_offset_bp":
                    midpoint_delta,

                "oriented_midpoint_offset_bp":
                    oriented_midpoint,

                "oriented_neighbor_strand":
                    oriented_neighbor_strand,

                "within_20_genes":
                    in_gene,

                "within_20kb":
                    in_bp,

                "extraction_basis":
                    extraction_basis,

                "is_focal":
                    int(
                        row.protein_id
                        ==
                        focal_protein
                    ),

                "globdb_cog":
                    raw_cog,

                "globdb_gene":
                    raw_gene,

                "globdb_product":
                    raw_product,

                "globdb_annotation_accepted":
                    accepted,

                "direct_COG20_cog":
                    "",

                "direct_COG20_function":
                    "",

                "comparison_cog":
                    comparison_cog,

                "comparison_product":
                    comparison_product,

                "comparison_annotation_source":
                    (
                        "GlobDB_COG20"
                        if comparison_cog
                        or comparison_product
                        else ""
                    ),

                "existing_local_family":
                    "",

                "old_global_mmseq_cluster":
                    clean(
                        getattr(
                            row,
                            "old_mmseq_cluster",
                            "",
                        )
                    ),

                "old_mcl_module":
                    clean(
                        getattr(
                            row,
                            "old_mcl_module",
                            "",
                        )
                    ),

                "fegenie_positive":
                    clean(
                        getattr(
                            row,
                            "fegenie_positive",
                            "",
                        )
                    ),

                "fegenie_HMMs":
                    clean(
                        getattr(
                            row,
                            "fegenie_HMMs",
                            "",
                        )
                    ),

                "findmehemes_positive":
                    clean(
                        getattr(
                            row,
                            "findmehemes_positive",
                            "",
                        )
                    ),

                "number_of_hemes":
                    clean(
                        getattr(
                            row,
                            "number_of_hemes",
                            "",
                        )
                    ),

                "neighborhood_protein_id":
                    "",
            }
        )


old_neighborhoods = pd.DataFrame(
    old_neighborhood_rows
)


if (
    old_neighborhoods[
        "is_focal"
    ].sum()
    !=
    EXPECTED_OLD_FOCALS
):

    fail(
        "Expected exactly 19 GlobDB focal-self rows."
    )


print(
    f"  GlobDB neighborhood rows: "
    f"{len(old_neighborhoods):,}"
)


## ================================================================== ##
## 6. Recover the 4 metagenome Cluster_00121 neighborhoods
## ================================================================== ##

print()
print(
    "Reading metagenome Cluster_00121 neighborhoods..."
)


mg = pd.read_csv(
    MG_GENE_TABLE,
    sep="\t",
    dtype=str,
    keep_default_na=False,
)


required_mg = {
    "focal_group",
    "focal_id",

    "genome",

    "focal_protein_id",
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
    "oriented_neighbor_strand",

    "within_20_genes",
    "within_20kb",
    "extraction_basis",

    "is_focal",

    "neighborhood_protein_id",
    "local_family_id",

    "direct_COG20_primary_cog",
    "direct_COG20_primary_function",
}


missing = (
    required_mg
    -
    set(
        mg.columns
    )
)


if missing:

    fail(
        "Metagenome gene table is missing columns:\n"
        +
        "\n".join(
            sorted(
                missing
            )
        )
    )


mg = (
    mg[
        mg[
            "focal_group"
        ]
        ==
        TARGET_MG_GROUP
    ]
    .copy()
    .reset_index(
        drop=True
    )
)


if mg[
    "focal_id"
].nunique() != EXPECTED_MG_FOCALS:

    fail(
        f"Expected {EXPECTED_MG_FOCALS} metagenome "
        f"Cluster_00121 focal regions; found "
        f"{mg['focal_id'].nunique()}."
    )


if (
    pd.to_numeric(
        mg[
            "is_focal"
        ],
        errors="raise",
    ).sum()
    !=
    EXPECTED_MG_FOCALS
):

    fail(
        "Expected exactly four metagenome focal-self rows."
    )


def mg_column(
    row,
    name,
    default="",
):

    if name not in mg.columns:

        return default

    return clean(
        row[
            name
        ]
    )


mg_rows = []


for _, row in mg.iterrows():

    direct_cog = mg_column(
        row,
        "direct_COG20_primary_cog",
    )


    direct_function = mg_column(
        row,
        "direct_COG20_primary_function",
    )


    transferred_cog = mg_column(
        row,
        "neighbor_transferred_cog",
    )


    transferred_product = mg_column(
        row,
        "neighbor_transferred_product",
    )


    if direct_cog:

        comparison_cog = (
            direct_cog
        )

        comparison_product = (
            direct_function
        )

        annotation_source = (
            "direct_COG20"
        )


    elif transferred_cog:

        comparison_cog = (
            transferred_cog
        )

        comparison_product = (
            transferred_product
        )

        annotation_source = (
            "GlobDB_homology_transfer"
        )


    else:

        comparison_cog = ""

        comparison_product = ""

        annotation_source = ""


    mg_rows.append(
        {
            "source_dataset":
                "Metagenome",

            "focal_id":
                clean(
                    row[
                        "focal_id"
                    ]
                ),

            "focal_cluster":
                TARGET_CLUSTER,

            "genome":
                clean(
                    row[
                        "genome"
                    ]
                ),

            "focal_protein_id":
                clean(
                    row[
                        "focal_protein_id"
                    ]
                ),

            "focal_contig":
                clean(
                    row[
                        "focal_contig"
                    ]
                ),

            "focal_start":
                row[
                    "focal_start"
                ],

            "focal_end":
                row[
                    "focal_end"
                ],

            "focal_strand":
                clean(
                    row[
                        "focal_strand"
                    ]
                ),

            "focal_gene_rank":
                row[
                    "focal_gene_rank"
                ],

            "neighbor_protein_id":
                clean(
                    row[
                        "neighbor_protein_id"
                    ]
                ),

            "neighbor_contig":
                mg_column(
                    row,
                    "neighbor_contig",
                    clean(
                        row[
                            "focal_contig"
                        ]
                    ),
                ),

            "neighbor_start":
                row[
                    "neighbor_start"
                ],

            "neighbor_end":
                row[
                    "neighbor_end"
                ],

            "neighbor_strand":
                clean(
                    row[
                        "neighbor_strand"
                    ]
                ),

            "neighbor_gene_rank":
                row[
                    "neighbor_gene_rank"
                ],

            "genomic_gene_offset":
                row[
                    "genomic_gene_offset"
                ],

            "oriented_gene_offset":
                row[
                    "oriented_gene_offset"
                ],

            "genomic_midpoint_offset_bp":
                row[
                    "genomic_midpoint_offset_bp"
                ],

            "oriented_midpoint_offset_bp":
                row[
                    "oriented_midpoint_offset_bp"
                ],

            "oriented_neighbor_strand":
                clean(
                    row[
                        "oriented_neighbor_strand"
                    ]
                ),

            "within_20_genes":
                row[
                    "within_20_genes"
                ],

            "within_20kb":
                row[
                    "within_20kb"
                ],

            "extraction_basis":
                clean(
                    row[
                        "extraction_basis"
                    ]
                ),

            "is_focal":
                row[
                    "is_focal"
                ],

            "globdb_cog":
                transferred_cog,

            "globdb_gene":
                "",

            "globdb_product":
                transferred_product,

            "globdb_annotation_accepted":
                "",

            "direct_COG20_cog":
                direct_cog,

            "direct_COG20_function":
                direct_function,

            "comparison_cog":
                comparison_cog,

            "comparison_product":
                comparison_product,

            "comparison_annotation_source":
                annotation_source,

            "existing_local_family":
                clean(
                    row[
                        "local_family_id"
                    ]
                ),

            "old_global_mmseq_cluster":
                mg_column(
                    row,
                    "globdb_top_scoring_cluster",
                ),

            "old_mcl_module":
                "",

            "fegenie_positive":
                mg_column(
                    row,
                    "neighbor_fegenie_positive",
                ),

            "fegenie_HMMs":
                mg_column(
                    row,
                    "neighbor_fegenie_HMMs",
                ),

            "findmehemes_positive":
                mg_column(
                    row,
                    "neighbor_findmehemes_positive",
                ),

            "number_of_hemes":
                mg_column(
                    row,
                    "neighbor_number_of_hemes",
                ),

            "neighborhood_protein_id":
                clean(
                    row[
                        "neighborhood_protein_id"
                    ]
                ),
        }
    )


mg_neighborhoods = pd.DataFrame(
    mg_rows
)


print(
    f"  Metagenome neighborhood rows: "
    f"{len(mg_neighborhoods):,}"
)


## ================================================================== ##
## 7. Combine neighborhood observations
## ================================================================== ##

combined = pd.concat(
    [
        old_neighborhoods,
        mg_neighborhoods,
    ],
    ignore_index=True,
)


if combined[
    "focal_id"
].nunique() != EXPECTED_TOTAL_FOCALS:

    fail(
        f"Expected {EXPECTED_TOTAL_FOCALS} total focal regions; "
        f"found {combined['focal_id'].nunique()}."
    )


if (
    pd.to_numeric(
        combined[
            "is_focal"
        ],
        errors="coerce",
    )
    .fillna(0)
    .astype(int)
    .sum()
    !=
    EXPECTED_TOTAL_FOCALS
):

    fail(
        "Expected exactly 23 focal-self rows in "
        "combined neighborhoods."
    )


## ================================================================== ##
## 8. Assign stable IDs to UNIQUE ACTUAL neighborhood proteins
##
## Same physical gene receives the same combined ID even if it were
## ever observed in >1 focal neighborhood.
## ================================================================== ##

protein_keys = (
    combined[
        [
            "source_dataset",
            "genome",
            "neighbor_protein_id",
            "neighborhood_protein_id",
        ]
    ]
    .drop_duplicates(
        [
            "source_dataset",
            "genome",
            "neighbor_protein_id",
        ]
    )
    .sort_values(
        [
            "source_dataset",
            "genome",
            "neighbor_protein_id",
        ],
        kind="stable",
    )
    .reset_index(
        drop=True
    )
)


protein_keys[
    "combined_protein_id"
] = [
    f"C121P{index:06d}"

    for index
    in range(
        1,
        len(
            protein_keys
        )
        + 1,
    )
]


combined = combined.merge(
    protein_keys,
    on=[
        "source_dataset",
        "genome",
        "neighbor_protein_id",
        "neighborhood_protein_id",
    ],
    how="left",
    validate="many_to_one",
)


if combined[
    "combined_protein_id"
].isna().any():

    fail(
        "Failed assigning combined protein IDs."
    )


## ================================================================== ##
## 9. Recover protein sequences
## ================================================================== ##

print()
print(
    "Recovering all combined neighborhood protein sequences..."
)


mg_sequences = fasta_dict(
    MG_PROTEINS
)


old_sequence_cache = {}


protein_metadata_rows = []

sequence_by_combined_id = {}


for row in protein_keys.itertuples(
    index=False
):

    source = row.source_dataset

    genome = row.genome

    protein_id = row.neighbor_protein_id

    combined_id = row.combined_protein_id


    if source == "Metagenome":

        nbp = row.neighborhood_protein_id


        if not nbp:

            fail(
                f"Metagenome protein {genome}/{protein_id} "
                f"has no NBP identifier."
            )


        sequence = mg_sequences.get(
            nbp
        )


        if sequence is None:

            fail(
                f"Could not find {nbp} in "
                f"{MG_PROTEINS}"
            )


        source_fasta = str(
            MG_PROTEINS
        )


    elif source == "GlobDB":

        if genome not in old_sequence_cache:

            proteome_path = resolve_old_proteome(
                genome
            )


            old_sequence_cache[
                genome
            ] = (
                proteome_path,
                fasta_dict(
                    proteome_path
                ),
            )


        proteome_path, sequences = (
            old_sequence_cache[
                genome
            ]
        )


        sequence = sequences.get(
            protein_id
        )


        if sequence is None:

            fail(
                f"Could not recover old protein "
                f"{genome}/{protein_id} from "
                f"{proteome_path}"
            )


        source_fasta = str(
            proteome_path
        )


    else:

        fail(
            f"Unknown source dataset: {source}"
        )


    sequence = sequence.rstrip(
        "*"
    )


    if not sequence:

        fail(
            f"Empty sequence for {combined_id}"
        )


    if "*" in sequence:

        fail(
            f"Internal stop in {combined_id}"
        )


    sequence_by_combined_id[
        combined_id
    ] = sequence


    protein_metadata_rows.append(
        {
            "combined_protein_id":
                combined_id,

            "source_dataset":
                source,

            "genome":
                genome,

            "protein_id":
                protein_id,

            "neighborhood_protein_id":
                row.neighborhood_protein_id,

            "sequence_length":
                len(
                    sequence
                ),

            "source_fasta":
                source_fasta,
        }
    )


protein_metadata = pd.DataFrame(
    protein_metadata_rows
)


if len(
    protein_metadata
) != len(
    protein_keys
):

    fail(
        "Protein sequence recovery changed protein count."
    )


## ================================================================== ##
## 10. Write all-neighborhood FASTA
## ================================================================== ##

with OUT_PROTEINS.open(
    "w",
    encoding="utf-8",
) as handle:

    for row in protein_metadata.itertuples(
        index=False
    ):

        sequence = sequence_by_combined_id[
            row.combined_protein_id
        ]


        handle.write(
            f">{row.combined_protein_id}\n"
        )


        for start in range(
            0,
            len(
                sequence
            ),
            80,
        ):

            handle.write(
                sequence[
                    start:
                    start
                    +
                    80
                ]
                +
                "\n"
            )


## ================================================================== ##
## 11. Write protein metadata and neighborhood table
## ================================================================== ##

protein_metadata.to_csv(
    OUT_PROTEIN_METADATA,
    sep="\t",
    index=False,
)


combined = (
    combined
    .sort_values(
        [
            "source_dataset",
            "genome",
            "focal_id",
            "oriented_gene_offset",
            "oriented_midpoint_offset_bp",
        ],
        kind="stable",
    )
    .reset_index(
        drop=True
    )
)


combined.to_csv(
    OUT_COMBINED_GENES,
    sep="\t",
    index=False,
)


## ================================================================== ##
## 12. Build the 23-protein focal FASTA for later phylogeny
## ================================================================== ##

focal_rows = (
    combined[
        pd.to_numeric(
            combined[
                "is_focal"
            ],
            errors="coerce",
        )
        .fillna(0)
        .astype(int)
        ==
        1
    ]
    .copy()
)


focal_metadata = (
    focal_rows[
        [
            "combined_protein_id",
            "source_dataset",
            "focal_id",
            "genome",
            "focal_protein_id",
            "focal_contig",
            "number_of_hemes",
            "comparison_cog",
            "comparison_product",
        ]
    ]
    .drop_duplicates(
        "combined_protein_id"
    )
    .reset_index(
        drop=True
    )
)


if len(
    focal_metadata
) != EXPECTED_TOTAL_FOCALS:

    fail(
        f"Expected {EXPECTED_TOTAL_FOCALS} unique focal proteins; "
        f"found {len(focal_metadata)}."
    )


focal_metadata.to_csv(
    OUT_FOCAL_METADATA,
    sep="\t",
    index=False,
)


with OUT_FOCAL_FASTA.open(
    "w",
    encoding="utf-8",
) as handle:

    for row in focal_metadata.itertuples(
        index=False
    ):

        sequence = sequence_by_combined_id[
            row.combined_protein_id
        ]


        handle.write(
            f">{row.combined_protein_id}\n"
        )


        for start in range(
            0,
            len(
                sequence
            ),
            80,
        ):

            handle.write(
                sequence[
                    start:
                    start
                    +
                    80
                ]
                +
                "\n"
            )


## ================================================================== ##
## 13. QC
## ================================================================== ##

qc = pd.DataFrame(
    [
        (
            "GlobDB_focal_proteins",
            old_neighborhoods[
                "focal_id"
            ].nunique(),
        ),

        (
            "metagenome_focal_proteins",
            mg_neighborhoods[
                "focal_id"
            ].nunique(),
        ),

        (
            "total_focal_proteins",
            combined[
                "focal_id"
            ].nunique(),
        ),

        (
            "GlobDB_neighborhood_rows",
            len(
                old_neighborhoods
            ),
        ),

        (
            "metagenome_neighborhood_rows",
            len(
                mg_neighborhoods
            ),
        ),

        (
            "combined_neighborhood_rows",
            len(
                combined
            ),
        ),

        (
            "unique_combined_proteins",
            len(
                protein_metadata
            ),
        ),

        (
            "focal_FASTA_sequences",
            len(
                focal_metadata
            ),
        ),

        (
            "rows_with_comparison_COG",
            int(
                (
                    combined[
                        "comparison_cog"
                    ]
                    .astype(str)
                    .str.strip()
                    !=
                    ""
                ).sum()
            ),
        ),
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
## Final report
## ================================================================== ##

print()
print("=" * 80)

print(
    "STAGE 63 COMPLETE"
)

print("=" * 80)

print()
print(
    f"GlobDB focal regions:          "
    f"{old_neighborhoods['focal_id'].nunique()}"
)

print(
    f"Metagenome focal regions:      "
    f"{mg_neighborhoods['focal_id'].nunique()}"
)

print(
    f"Total focal regions:           "
    f"{combined['focal_id'].nunique()}"
)

print()
print(
    f"GlobDB neighborhood rows:      "
    f"{len(old_neighborhoods):,}"
)

print(
    f"Metagenome neighborhood rows:  "
    f"{len(mg_neighborhoods):,}"
)

print(
    f"Combined neighborhood rows:    "
    f"{len(combined):,}"
)

print(
    f"Unique neighborhood proteins:  "
    f"{len(protein_metadata):,}"
)

print()
print(
    f"Focal sequences for phylogeny: "
    f"{len(focal_metadata)}"
)

print()
print("Combined neighborhood table:")

print(
    f"  {OUT_COMBINED_GENES}"
)

print()
print("Combined protein FASTA:")

print(
    f"  {OUT_PROTEINS}"
)

print()
print("23 focal proteins:")

print(
    f"  {OUT_FOCAL_FASTA}"
)

print()
print("QC:")

print(
    f"  {OUT_QC}"
)
