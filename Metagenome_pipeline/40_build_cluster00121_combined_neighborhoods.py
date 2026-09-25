#!/usr/bin/env python3

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd


## ================================================================== ##
## STAGE 63
##
## BUILD COMBINED GLOBDB + METAGENOME NEIGHBORHOODS FOR CLUSTER_00121
##
## Goal
## ----
##
## Compare genomic architecture around:
##
##     19 GlobDB Cluster_00121 proteins
##      4 metagenome Cluster_00121 proteins
##
## Total:
##
##     23 focal loci
##
##
## Neighborhood definition
## -----------------------
##
## Retain every CDS satisfying either:
##
##     abs(gene-rank offset) <= 20
##
## OR
##
##     CDS intersects:
##
##         focal_start - 20,000 bp
##         ...
##         focal_end   + 20,000 bp
##
##
## The output contains ACTUAL GENES.
##
## No representative selection is performed here.
## No phylogeny is performed here.
##
## All unique neighborhood proteins are exported for a subsequent
## shared 40% identity / 80% coverage MMseqs clustering.
## ================================================================== ##


## ================================================================== ##
## Paths
## ================================================================== ##

ROOT = (
    Path.home()
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


MG_PROTEIN_FASTA = (
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


OUT_GLOBDB_FOCALS = (
    OUT_DIR
    / "Cluster_00121_GlobDB_focals.tsv"
)


OUT_COMBINED = (
    OUT_DIR
    / "Cluster_00121_combined_neighborhood_genes.tsv"
)


OUT_PROTEIN_METADATA = (
    OUT_DIR
    / "Cluster_00121_combined_neighborhood_protein_metadata.tsv"
)


OUT_FASTA = (
    OUT_DIR
    / "Cluster_00121_combined_neighborhood_proteins.faa"
)


OUT_FOCAL_SUMMARY = (
    OUT_DIR
    / "Cluster_00121_focal_regions.tsv"
)


OUT_QC = (
    OUT_DIR
    / "Cluster_00121_neighborhood_build_qc.tsv"
)


## ================================================================== ##
## Constants
## ================================================================== ##

TARGET_CLUSTER = "Cluster_00121"

TARGET_MG_GROUP = "Multiheme__Cluster_00121"

GENE_RADIUS = 20

BP_RADIUS = 20_000

EXPECTED_GLOBDB = 19

EXPECTED_METAGENOME = 4

EXPECTED_TOTAL = 23


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


def as_flag(value):

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


def flip_strand(value):

    value = clean(
        value
    )

    if value == "+":
        return "-"

    if value == "-":
        return "+"

    return value


def resolve_column(
    columns,
    alternatives,
    description,
    required=True,
):

    for name in alternatives:

        if name in columns:

            return name


    if required:

        fail(
            f"Could not resolve column for {description}. "
            f"Tried: {', '.join(alternatives)}"
        )


    return None


def fasta_iter(path):

    header = None

    sequence = []


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
                            sequence
                        ),
                    )


                header = line[
                    1:
                ]

                sequence = []


            else:

                sequence.append(
                    line
                )


    if header is not None:

        yield (
            header,
            "".join(
                sequence
            ),
        )


def read_fasta(path):

    result = {}


    for header, sequence in fasta_iter(
        path
    ):

        identifier = (
            header
            .split()[0]
        )


        if identifier in result:

            fail(
                f"Duplicate FASTA identifier "
                f"{identifier} in {path}"
            )


        result[
            identifier
        ] = sequence.rstrip(
            "*"
        )


    return result


def find_old_proteome(genome):

    matches = sorted(
        OLD_ORF_DIR.glob(
            f"{genome}*-proteins.faa"
        )
    )


    if len(
        matches
    ) != 1:

        fail(
            f"Expected exactly one protein FASTA for "
            f"{genome}; found {len(matches)}:\n"
            +
            "\n".join(
                str(path)
                for path
                in matches
            )
        )


    return matches[
        0
    ]


## ================================================================== ##
## Start
## ================================================================== ##

print("=" * 80)

print(
    "STAGE 63 - BUILD COMBINED CLUSTER_00121 NEIGHBORHOODS"
)

print("=" * 80)


for path in [
    OLD_CLUSTER_MEMBERSHIP,
    OLD_GENE_CATALOG,
    MG_GENE_TABLE,
    MG_PROTEIN_FASTA,
]:

    require_file(
        path
    )


if not OLD_ORF_DIR.is_dir():

    fail(
        f"Old ORF directory missing:\n{OLD_ORF_DIR}"
    )


## ================================================================== ##
## 1. Identify the 19 old Cluster_00121 proteins
## ================================================================== ##

print()
print(
    "Identifying GlobDB Cluster_00121 focal proteins..."
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
        "cluster_membership.tsv is missing: "
        +
        ", ".join(
            sorted(
                missing
            )
        )
    )


globdb_focals = (
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
    globdb_focals
) != EXPECTED_GLOBDB:

    fail(
        f"Expected {EXPECTED_GLOBDB} GlobDB "
        f"{TARGET_CLUSTER} proteins; found "
        f"{len(globdb_focals)}."
    )


if globdb_focals[
    [
        "genome",
        "protein_id",
    ]
].duplicated().any():

    fail(
        "Duplicate GlobDB Cluster_00121 focal proteins."
    )


globdb_focals.to_csv(
    OUT_GLOBDB_FOCALS,
    sep="\t",
    index=False,
)


globdb_genomes = set(
    globdb_focals[
        "genome"
    ]
)


print(
    f"  Focal proteins: "
    f"{len(globdb_focals)}"
)

print(
    f"  Genomes:        "
    f"{globdb_focals['genome'].nunique()}"
)


## ================================================================== ##
## 2. Resolve old complete-gene catalogue schema
## ================================================================== ##

print()
print(
    "Resolving old complete-gene catalogue schema..."
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


optional = {
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
        ],

    "findmehemes_positive":
        [
            "findmehemes_positive",
        ],

    "number_of_hemes":
        [
            "number_of_hemes",
        ],
}


resolved_optional = {}


for output_name, alternatives in optional.items():

    resolved_optional[
        output_name
    ] = resolve_column(
        columns,
        alternatives,
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
    column

    for column
    in resolved_optional.values()

    if column is not None
]


selected_columns = list(
    dict.fromkeys(
        selected_columns
    )
)


## ================================================================== ##
## 3. Read only the 19 relevant genomes from Stage 11
## ================================================================== ##

print()
print(
    "Reading relevant GlobDB genomes from Stage-11 catalogue..."
)


catalog_chunks = []


for chunk in pd.read_csv(
    OLD_GENE_CATALOG,
    sep="\t",
    dtype=str,
    keep_default_na=False,
    usecols=selected_columns,
    chunksize=200_000,
):

    selected = chunk[
        chunk[
            COL_GENOME
        ].isin(
            globdb_genomes
        )
    ]


    if len(
        selected
    ) > 0:

        catalog_chunks.append(
            selected.copy()
        )


if not catalog_chunks:

    fail(
        "No genes recovered for Cluster_00121 genomes."
    )


catalog = pd.concat(
    catalog_chunks,
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


for output_name, old_column in resolved_optional.items():

    if old_column is not None:

        rename[
            old_column
        ] = output_name


catalog = catalog.rename(
    columns=rename
)


for output_name in optional:

    if output_name not in catalog.columns:

        catalog[
            output_name
        ] = ""


for column in [
    "start",
    "end",
    "gene_rank",
]:

    catalog[
        column
    ] = pd.to_numeric(
        catalog[
            column
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
        "Duplicate genome + protein_id rows in "
        "selected GlobDB catalogue."
    )


print(
    f"  Genes:   "
    f"{len(catalog):,}"
)

print(
    f"  Genomes: "
    f"{catalog['genome'].nunique()}"
)


## ================================================================== ##
## 4. Recover focal coordinates
## ================================================================== ##

focal_coordinates = (
    globdb_focals[
        [
            "genome",
            "protein_id",
        ]
    ]
    .merge(
        catalog,
        on=[
            "genome",
            "protein_id",
        ],
        how="left",
        validate="one_to_one",
    )
)


if focal_coordinates[
    "contig"
].isna().any():

    fail(
        "One or more GlobDB focal proteins could not "
        "be located in the Stage-11 catalogue."
    )


## ================================================================== ##
## 5. Extract GlobDB neighborhoods
## ================================================================== ##

print()
print(
    "Extracting GlobDB neighborhoods..."
)


globdb_rows = []


for index, focal in enumerate(
    focal_coordinates.itertuples(
        index=False
    ),
    start=1,
):

    genome = focal.genome

    focal_protein = focal.protein_id

    contig = focal.contig

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
        f"GlobDB_C121_{index:02d}"
    )


    genes = (
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
                contig
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


    rank_offsets = (
        genes[
            "gene_rank"
        ].to_numpy(
            dtype=int
        )
        -
        focal_rank
    )


    within_gene_window = (
        np.abs(
            rank_offsets
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
            genes[
                "end"
            ].to_numpy(
                dtype=int
            )
            >=
            window_left
        )
        &
        (
            genes[
                "start"
            ].to_numpy(
                dtype=int
            )
            <=
            window_right
        )
    )


    genes = genes.loc[
        (
            within_gene_window
            |
            within_bp_window
        )
    ].copy()


    focal_midpoint = (
        focal_start
        +
        focal_end
    ) / 2.0


    for gene in genes.itertuples(
        index=False
    ):

        start = int(
            gene.start
        )

        end = int(
            gene.end
        )

        rank = int(
            gene.gene_rank
        )

        strand = clean(
            gene.strand
        )


        genomic_gene_offset = (
            rank
            -
            focal_rank
        )


        midpoint = (
            start
            +
            end
        ) / 2.0


        genomic_midpoint_offset = (
            midpoint
            -
            focal_midpoint
        )


        if focal_strand == "+":

            oriented_gene_offset = (
                genomic_gene_offset
            )

            oriented_midpoint_offset = (
                genomic_midpoint_offset
            )

            oriented_strand = (
                strand
            )


        else:

            oriented_gene_offset = (
                -genomic_gene_offset
            )

            oriented_midpoint_offset = (
                -genomic_midpoint_offset
            )

            oriented_strand = (
                flip_strand(
                    strand
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
                end
                >=
                window_left
            )
            and
            (
                start
                <=
                window_right
            )
        )


        if (
            in_gene
            and
            in_bp
        ):

            basis = "both"


        elif in_gene:

            basis = "gene_window"


        else:

            basis = "bp_window"


        annotation_accepted = as_flag(
            getattr(
                gene,
                "annotation_accepted",
                "",
            )
        )


        cog = clean(
            getattr(
                gene,
                "globdb_cog",
                "",
            )
        )


        product = clean(
            getattr(
                gene,
                "globdb_product",
                "",
            )
        )


        if (
            resolved_optional[
                "annotation_accepted"
            ]
            is None
        ):

            comparison_cog = cog

            comparison_product = product


        elif annotation_accepted == 1:

            comparison_cog = cog

            comparison_product = product


        else:

            comparison_cog = ""

            comparison_product = ""


        globdb_rows.append(
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
                    gene.protein_id,

                "neighbor_contig":
                    gene.contig,

                "neighbor_start":
                    start,

                "neighbor_end":
                    end,

                "neighbor_strand":
                    strand,

                "neighbor_gene_rank":
                    rank,

                "genomic_gene_offset":
                    genomic_gene_offset,

                "oriented_gene_offset":
                    oriented_gene_offset,

                "genomic_midpoint_offset_bp":
                    genomic_midpoint_offset,

                "oriented_midpoint_offset_bp":
                    oriented_midpoint_offset,

                "oriented_neighbor_strand":
                    oriented_strand,

                "within_20_genes":
                    in_gene,

                "within_20kb":
                    in_bp,

                "extraction_basis":
                    basis,

                "is_focal":
                    int(
                        gene.protein_id
                        ==
                        focal_protein
                    ),

                "comparison_cog":
                    comparison_cog,

                "comparison_product":
                    comparison_product,

                "comparison_annotation_source":
                    (
                        "GlobDB"
                        if (
                            comparison_cog
                            or
                            comparison_product
                        )
                        else
                        ""
                    ),

                "globdb_cog":
                    cog,

                "globdb_gene":
                    clean(
                        getattr(
                            gene,
                            "globdb_gene",
                            "",
                        )
                    ),

                "globdb_product":
                    product,

                "direct_COG20_cog":
                    "",

                "direct_COG20_function":
                    "",

                "existing_local_family":
                    "",

                "old_global_mmseq_cluster":
                    clean(
                        getattr(
                            gene,
                            "old_mmseq_cluster",
                            "",
                        )
                    ),

                "old_mcl_module":
                    clean(
                        getattr(
                            gene,
                            "old_mcl_module",
                            "",
                        )
                    ),

                "fegenie_positive":
                    clean(
                        getattr(
                            gene,
                            "fegenie_positive",
                            "",
                        )
                    ),

                "fegenie_HMMs":
                    clean(
                        getattr(
                            gene,
                            "fegenie_HMMs",
                            "",
                        )
                    ),

                "findmehemes_positive":
                    clean(
                        getattr(
                            gene,
                            "findmehemes_positive",
                            "",
                        )
                    ),

                "number_of_hemes":
                    clean(
                        getattr(
                            gene,
                            "number_of_hemes",
                            "",
                        )
                    ),

                "neighborhood_protein_id":
                    "",
            }
        )


globdb = pd.DataFrame(
    globdb_rows
)


if (
    globdb[
        "is_focal"
    ].sum()
    !=
    EXPECTED_GLOBDB
):

    fail(
        "GlobDB neighborhood extraction did not "
        "produce exactly 19 focal-self rows."
    )


print(
    f"  GlobDB neighborhood rows: "
    f"{len(globdb):,}"
)


## ================================================================== ##
## 6. Extract the four metagenome Cluster_00121 regions
## ================================================================== ##

print()
print(
    "Extracting metagenome Cluster_00121 regions..."
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
        "Metagenome neighborhood table is missing:\n"
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
].nunique() != EXPECTED_METAGENOME:

    fail(
        f"Expected {EXPECTED_METAGENOME} metagenome "
        f"Cluster_00121 regions; found "
        f"{mg['focal_id'].nunique()}."
    )


def value(
    row,
    column,
):

    if column not in mg.columns:

        return ""

    return clean(
        row[
            column
        ]
    )


mg_rows = []


for _, row in mg.iterrows():

    direct_cog = value(
        row,
        "direct_COG20_primary_cog",
    )


    direct_function = value(
        row,
        "direct_COG20_primary_function",
    )


    transferred_cog = value(
        row,
        "neighbor_transferred_cog",
    )


    transferred_product = value(
        row,
        "neighbor_transferred_product",
    )


    ## Prefer directly predicted COG20 annotation.
    ## Homology transfer remains visible separately.
    if direct_cog:

        comparison_cog = direct_cog

        comparison_product = direct_function

        annotation_source = (
            "direct_COG20"
        )


    elif transferred_cog:

        comparison_cog = transferred_cog

        comparison_product = transferred_product

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
                value(
                    row,
                    "neighbor_contig",
                )
                or
                clean(
                    row[
                        "focal_contig"
                    ]
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

            "comparison_cog":
                comparison_cog,

            "comparison_product":
                comparison_product,

            "comparison_annotation_source":
                annotation_source,

            "globdb_cog":
                transferred_cog,

            "globdb_gene":
                "",

            "globdb_product":
                transferred_product,

            "direct_COG20_cog":
                direct_cog,

            "direct_COG20_function":
                direct_function,

            "existing_local_family":
                clean(
                    row[
                        "local_family_id"
                    ]
                ),

            "old_global_mmseq_cluster":
                value(
                    row,
                    "globdb_top_scoring_cluster",
                ),

            "old_mcl_module":
                "",

            "fegenie_positive":
                value(
                    row,
                    "neighbor_fegenie_positive",
                ),

            "fegenie_HMMs":
                value(
                    row,
                    "neighbor_fegenie_HMMs",
                ),

            "findmehemes_positive":
                value(
                    row,
                    "neighbor_findmehemes_positive",
                ),

            "number_of_hemes":
                value(
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


metagenome = pd.DataFrame(
    mg_rows
)


if (
    pd.to_numeric(
        metagenome[
            "is_focal"
        ],
        errors="coerce",
    )
    .fillna(0)
    .astype(int)
    .sum()
    !=
    EXPECTED_METAGENOME
):

    fail(
        "Metagenome subset does not contain "
        "exactly four focal-self rows."
    )


print(
    f"  Metagenome neighborhood rows: "
    f"{len(metagenome):,}"
)


## ================================================================== ##
## 7. Combine both datasets
## ================================================================== ##

combined = pd.concat(
    [
        globdb,
        metagenome,
    ],
    ignore_index=True,
)


if combined[
    "focal_id"
].nunique() != EXPECTED_TOTAL:

    fail(
        f"Expected {EXPECTED_TOTAL} focal regions; "
        f"found {combined['focal_id'].nunique()}."
    )


## ================================================================== ##
## 8. Assign stable combined protein IDs
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
    f"C121P{i:06d}"

    for i in range(
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
        "Failed to assign combined protein IDs."
    )


## ================================================================== ##
## 9. Recover protein sequences
## ================================================================== ##

print()
print(
    "Recovering neighborhood protein sequences..."
)


mg_sequences = read_fasta(
    MG_PROTEIN_FASTA
)


old_cache = {}


metadata_rows = []

sequences = {}


for row in protein_keys.itertuples(
    index=False
):

    source = row.source_dataset

    genome = row.genome

    protein_id = row.neighbor_protein_id

    combined_id = row.combined_protein_id


    if source == "Metagenome":

        nbp = clean(
            row.neighborhood_protein_id
        )


        if not nbp:

            fail(
                f"Metagenome protein {protein_id} "
                f"has no neighborhood protein ID."
            )


        if nbp not in mg_sequences:

            fail(
                f"{nbp} not found in metagenome FASTA."
            )


        sequence = mg_sequences[
            nbp
        ]

        source_fasta = (
            MG_PROTEIN_FASTA
        )


    else:

        if genome not in old_cache:

            proteome = find_old_proteome(
                genome
            )


            old_cache[
                genome
            ] = (
                proteome,
                read_fasta(
                    proteome
                ),
            )


        proteome, proteome_sequences = (
            old_cache[
                genome
            ]
        )


        if protein_id not in proteome_sequences:

            fail(
                f"{protein_id} not found in "
                f"{proteome}"
            )


        sequence = proteome_sequences[
            protein_id
        ]

        source_fasta = proteome


    if not sequence:

        fail(
            f"Empty protein sequence for {combined_id}"
        )


    if "*" in sequence:

        fail(
            f"Internal stop in {combined_id}"
        )


    sequences[
        combined_id
    ] = sequence


    metadata_rows.append(
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
                str(
                    source_fasta
                ),
        }
    )


protein_metadata = pd.DataFrame(
    metadata_rows
)


## ================================================================== ##
## 10. Write shared neighborhood protein FASTA
## ================================================================== ##

with OUT_FASTA.open(
    "w",
    encoding="utf-8",
) as handle:

    for row in protein_metadata.itertuples(
        index=False
    ):

        sequence = sequences[
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


protein_metadata.to_csv(
    OUT_PROTEIN_METADATA,
    sep="\t",
    index=False,
)


## ================================================================== ##
## 11. Write combined neighborhood table
## ================================================================== ##

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
    OUT_COMBINED,
    sep="\t",
    index=False,
)


## ================================================================== ##
## 12. Focal-region summary
## ================================================================== ##

focal_summary = (
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
    [
        [
            "source_dataset",
            "focal_id",
            "genome",
            "focal_protein_id",
            "focal_contig",
            "focal_start",
            "focal_end",
            "focal_strand",
            "number_of_hemes",
            "combined_protein_id",
        ]
    ]
    .drop_duplicates()
    .sort_values(
        [
            "source_dataset",
            "genome",
            "focal_protein_id",
        ],
        kind="stable",
    )
    .reset_index(
        drop=True
    )
)


if len(
    focal_summary
) != EXPECTED_TOTAL:

    fail(
        f"Expected {EXPECTED_TOTAL} focal-summary rows; "
        f"found {len(focal_summary)}."
    )


focal_summary.to_csv(
    OUT_FOCAL_SUMMARY,
    sep="\t",
    index=False,
)


## ================================================================== ##
## 13. QC
## ================================================================== ##

qc = pd.DataFrame(
    [
        (
            "GlobDB_focal_regions",
            globdb[
                "focal_id"
            ].nunique(),
        ),

        (
            "metagenome_focal_regions",
            metagenome[
                "focal_id"
            ].nunique(),
        ),

        (
            "total_focal_regions",
            combined[
                "focal_id"
            ].nunique(),
        ),

        (
            "GlobDB_neighborhood_rows",
            len(
                globdb
            ),
        ),

        (
            "metagenome_neighborhood_rows",
            len(
                metagenome
            ),
        ),

        (
            "combined_neighborhood_rows",
            len(
                combined
            ),
        ),

        (
            "unique_neighborhood_proteins",
            len(
                protein_metadata
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
## Report
## ================================================================== ##

print()
print("=" * 80)

print(
    "STAGE 63 COMPLETE"
)

print("=" * 80)

print()
print(
    f"GlobDB focal regions:         "
    f"{globdb['focal_id'].nunique()}"
)

print(
    f"Metagenome focal regions:     "
    f"{metagenome['focal_id'].nunique()}"
)

print(
    f"Combined focal regions:       "
    f"{combined['focal_id'].nunique()}"
)

print()
print(
    f"GlobDB neighborhood rows:     "
    f"{len(globdb):,}"
)

print(
    f"Metagenome neighborhood rows: "
    f"{len(metagenome):,}"
)

print(
    f"Combined neighborhood rows:   "
    f"{len(combined):,}"
)

print(
    f"Unique neighborhood proteins: "
    f"{len(protein_metadata):,}"
)

print()
print("Combined neighborhood table:")

print(
    f"  {OUT_COMBINED}"
)

print()
print("Combined neighborhood proteins:")

print(
    f"  {OUT_FASTA}"
)

print()
print("Focal-region table:")

print(
    f"  {OUT_FOCAL_SUMMARY}"
)

print()
print("QC:")

print(
    f"  {OUT_QC}"
)
