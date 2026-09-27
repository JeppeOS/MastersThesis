#!/usr/bin/env python3

from __future__ import annotations

import gzip
import sys
from pathlib import Path

import pandas as pd


## ================================================================== ##
## STAGE 51
##
## PRIORITY FOCAL CONTEXT OBSERVABILITY
##
## Purpose
## -------
##
## For the 21 priority Methylococcales focal proteins:
##
##   - recover nucleotide contig lengths from the authoritative
##     SemiBin2 MAG FASTAs;
##
##   - validate all 59,786 Stage-50 Prodigal gene coordinates against
##     those FASTAs;
##
##   - calculate physical bp availability to each contig end;
##
##   - calculate gene-count availability to each contig end;
##
##   - orient both quantities relative to focal transcription;
##
##   - explicitly flag complete 5/10/20 kb context;
##
##   - explicitly flag complete 5/10/20-gene context.
##
## This stage does NOT extract neighborhoods and does NOT interpret
## missing neighboring genes biologically.
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


MAG_MANIFEST = (
    PROJECT
    / "functional_analysis"
    / "SemiBin2_187_MAG_manifest.tsv"
)


OUT_CONTIG_LENGTHS = (
    WORK_DIR
    / "methylococcales_contig_lengths.tsv"
)


OUT_OBSERVABILITY = (
    WORK_DIR
    / "priority_context_observability.tsv"
)


OUT_QC = (
    WORK_DIR
    / "context_observability_qc.tsv"
)


OUT_MISSING_CONTIGS = (
    WORK_DIR
    / "missing_gene_contigs.tsv"
)


OUT_BAD_COORDINATES = (
    WORK_DIR
    / "genes_exceeding_contig_length.tsv"
)


## ================================================================== ##
## Locked values established by Stage 50
## ================================================================== ##

EXPECTED_MAGS = 16

EXPECTED_GENES = 59_786

EXPECTED_GENE_CONTIGS = 714

EXPECTED_FOCALS = 21


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


def resolve_path(value):

    path = Path(
        str(value).strip()
    ).expanduser()

    if not path.is_absolute():

        path = (
            PROJECT
            / path
        )

    return path.resolve()


def open_maybe_gzip(path):

    ## -------------------------------------------------------------- ##
    ## SemiBin-derived FASTAs have previously been encountered where
    ## compression did not agree with the filename extension.
    ##
    ## Detect gzip from the actual magic bytes:
    ##
    ##     1f 8b
    ## -------------------------------------------------------------- ##

    with path.open(
        "rb"
    ) as handle:

        magic = handle.read(
            2
        )


    if magic == b"\x1f\x8b":

        return gzip.open(
            path,
            "rt",
        )


    return path.open(
        "r",
        encoding="utf-8",
    )


def parse_fasta_lengths(
    genome,
    path,
):

    rows = []

    seen = set()

    current_contig = None

    current_length = 0


    def flush():

        nonlocal current_contig
        nonlocal current_length


        if current_contig is None:

            return


        if current_contig in seen:

            fail(
                f"Duplicate contig ID {current_contig} "
                f"in FASTA:\n{path}"
            )


        seen.add(
            current_contig
        )


        rows.append(
            {
                "genome":
                    genome,

                "contig":
                    current_contig,

                "contig_length_nt":
                    current_length,

                "source_nucleotide_fasta":
                    str(
                        path
                    ),
            }
        )


    with open_maybe_gzip(
        path
    ) as handle:

        for raw in handle:

            if raw.startswith(
                ">"
            ):

                flush()


                header = (
                    raw[1:]
                    .strip()
                )


                if not header:

                    fail(
                        f"Empty FASTA header in:\n"
                        f"{path}"
                    )


                current_contig = (
                    header.split()[0]
                )


                current_length = 0

                continue


            sequence = "".join(
                raw.split()
            )


            if not sequence:

                continue


            if current_contig is None:

                fail(
                    "Sequence encountered before "
                    f"first FASTA header:\n{path}"
                )


            current_length += len(
                sequence
            )


    flush()


    if not rows:

        fail(
            f"No sequences parsed from:\n"
            f"{path}"
        )


    return rows


## ================================================================== ##
## Start
## ================================================================== ##

print("=" * 80)
print("STAGE 51 - BUILD PRIORITY CONTEXT OBSERVABILITY")
print("=" * 80)


for path in (
    GENE_CATALOG,
    FOCALS,
    MAG_MANIFEST,
):

    require_file(
        path
    )


## ================================================================== ##
## 1. Read Stage-50 catalogue and focal genes
## ================================================================== ##

print()
print("Reading Stage-50 gene catalogue...")


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
    "start",
    "end",
    "strand",
    "contig_gene_rank",
    "contig_gene_count",
}


missing = (
    required_catalog
    - set(
        catalog.columns
    )
)


if missing:

    fail(
        "Gene catalogue is missing columns:\n"
        + "\n".join(
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
        "in gene catalogue."
    )


for column in [
    "start",
    "end",
    "contig_gene_rank",
    "contig_gene_count",
]:

    catalog[
        column
    ] = pd.to_numeric(
        catalog[
            column
        ],
        errors="raise",
    ).astype(int)


gene_contigs = (
    catalog[
        [
            "genome",
            "contig",
        ]
    ]
    .drop_duplicates()
)


if len(
    gene_contigs
) != EXPECTED_GENE_CONTIGS:

    fail(
        f"Expected {EXPECTED_GENE_CONTIGS:,} "
        f"gene-containing contigs but found "
        f"{len(gene_contigs):,}."
    )


genomes = sorted(
    catalog[
        "genome"
    ].unique()
)


if len(
    genomes
) != EXPECTED_MAGS:

    fail(
        f"Expected {EXPECTED_MAGS} MAGs but "
        f"found {len(genomes)}."
    )


print(
    f"  Genes:                  "
    f"{len(catalog):,}"
)

print(
    f"  Gene-containing contigs:"
    f" {len(gene_contigs):,}"
)

print(
    f"  MAGs:                   "
    f"{len(genomes):,}"
)


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
    "start",
    "end",
    "strand",
    "contig_gene_rank",
    "contig_gene_count",
    "priority_reason",
}


missing = (
    required_focals
    - set(
        focals.columns
    )
)


if missing:

    fail(
        "Priority focal table is missing columns:\n"
        + "\n".join(
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
    "start",
    "end",
    "contig_gene_rank",
    "contig_gene_count",
]:

    focals[
        column
    ] = pd.to_numeric(
        focals[
            column
        ],
        errors="raise",
    ).astype(int)


## ================================================================== ##
## 2. Locate nucleotide FASTAs from authoritative SemiBin2 manifest
## ================================================================== ##

print()
print("Locating authoritative SemiBin2 nucleotide FASTAs...")


manifest = pd.read_csv(
    MAG_MANIFEST,
    sep="\t",
    dtype=str,
    keep_default_na=False,
)


required_manifest = {
    "Genome_ID",
    "Functional_input_FASTA",
}


missing = (
    required_manifest
    - set(
        manifest.columns
    )
)


if missing:

    fail(
        "SemiBin2 manifest is missing columns:\n"
        + "\n".join(
            sorted(
                missing
            )
        )
    )


manifest = manifest[
    manifest[
        "Genome_ID"
    ].isin(
        genomes
    )
].copy()


if len(
    manifest
) != EXPECTED_MAGS:

    fail(
        f"Expected {EXPECTED_MAGS} selected MAG rows "
        f"in SemiBin2 manifest but found "
        f"{len(manifest)}."
    )


if manifest[
    "Genome_ID"
].duplicated().any():

    fail(
        "Duplicate Genome_ID in selected "
        "SemiBin2 manifest rows."
    )


genome_fastas = {}


for row in manifest.itertuples(
    index=False
):

    genome = str(
        row.Genome_ID
    ).strip()


    fasta = resolve_path(
        row.Functional_input_FASTA
    )


    if not fasta.is_file():

        fail(
            f"Nucleotide FASTA does not exist for "
            f"{genome}:\n{fasta}"
        )


    genome_fastas[
        genome
    ] = fasta


missing_fastas = sorted(
    set(
        genomes
    )
    - set(
        genome_fastas
    )
)


if missing_fastas:

    fail(
        "MAGs lacking nucleotide FASTAs:\n"
        + "\n".join(
            missing_fastas
        )
    )


print(
    f"  Nucleotide FASTAs: "
    f"{len(genome_fastas)}"
)


## ================================================================== ##
## 3. Parse all contig lengths
## ================================================================== ##

print()
print("Parsing nucleotide contig lengths...")


contig_rows = []


for number, genome in enumerate(
    genomes,
    start=1,
):

    rows = parse_fasta_lengths(
        genome,
        genome_fastas[
            genome
        ],
    )


    contig_rows.extend(
        rows
    )


    print(
        f"  Parsed {number:>2}/"
        f"{len(genomes)}: "
        f"{genome} "
        f"({len(rows):,} contigs)"
    )


contig_lengths = pd.DataFrame(
    contig_rows
)


if contig_lengths[
    [
        "genome",
        "contig",
    ]
].duplicated().any():

    fail(
        "Duplicate genome + contig keys "
        "in FASTA length table."
    )


if (
    contig_lengths[
        "contig_length_nt"
    ]
    <= 0
).any():

    fail(
        "One or more FASTA contigs have "
        "non-positive length."
    )


contig_lengths = (
    contig_lengths
    .sort_values(
        [
            "genome",
            "contig",
        ],
        kind="stable",
    )
    .reset_index(
        drop=True
    )
)


contig_lengths.to_csv(
    OUT_CONTIG_LENGTHS,
    sep="\t",
    index=False,
)


print(
    f"  Total FASTA contigs: "
    f"{len(contig_lengths):,}"
)


## ================================================================== ##
## 4. Validate every Prodigal gene against nucleotide contig lengths
## ================================================================== ##

print()
print("Validating all Stage-50 gene coordinates against FASTAs...")


catalog_contigs = (
    catalog
    .groupby(
        [
            "genome",
            "contig",
        ],
        as_index=False,
    )
    .agg(
        n_genes=(
            "protein_id",
            "size",
        ),

        minimum_gene_start=(
            "start",
            "min",
        ),

        maximum_gene_end=(
            "end",
            "max",
        ),
    )
)


if len(
    catalog_contigs
) != EXPECTED_GENE_CONTIGS:

    fail(
        f"Expected {EXPECTED_GENE_CONTIGS:,} "
        f"gene-containing contigs but found "
        f"{len(catalog_contigs):,}."
    )


catalog_check = catalog_contigs.merge(
    contig_lengths[
        [
            "genome",
            "contig",
            "contig_length_nt",
        ]
    ],
    on=[
        "genome",
        "contig",
    ],
    how="left",
    validate="one_to_one",
)


missing_contigs = catalog_check[
    catalog_check[
        "contig_length_nt"
    ].isna()
].copy()


missing_contigs.to_csv(
    OUT_MISSING_CONTIGS,
    sep="\t",
    index=False,
)


if len(
    missing_contigs
) > 0:

    fail(
        f"{len(missing_contigs):,} Prodigal "
        f"gene-containing contigs are absent "
        f"from their nucleotide FASTAs.\n"
        f"See:\n{OUT_MISSING_CONTIGS}"
    )


catalog_check[
    "contig_length_nt"
] = (
    catalog_check[
        "contig_length_nt"
    ]
    .astype(int)
)


bad_coordinates = catalog_check[
    (
        catalog_check[
            "minimum_gene_start"
        ]
        < 1
    )
    |
    (
        catalog_check[
            "maximum_gene_end"
        ]
        >
        catalog_check[
            "contig_length_nt"
        ]
    )
].copy()


bad_coordinates.to_csv(
    OUT_BAD_COORDINATES,
    sep="\t",
    index=False,
)


if len(
    bad_coordinates
) > 0:

    fail(
        f"{len(bad_coordinates):,} contigs contain "
        f"Prodigal coordinates outside the "
        f"nucleotide FASTA bounds.\n"
        f"See:\n{OUT_BAD_COORDINATES}"
    )


print(
    f"  Prodigal genes validated: "
    f"{len(catalog):,}"
)

print(
    f"  Missing gene contigs:     "
    f"{len(missing_contigs):,}"
)

print(
    f"  Coordinate violations:    "
    f"{len(bad_coordinates):,}"
)


## ================================================================== ##
## 5. Join contig lengths onto priority focal genes
## ================================================================== ##

print()
print("Calculating focal context observability...")


observability = focals.merge(
    contig_lengths[
        [
            "genome",
            "contig",
            "contig_length_nt",
            "source_nucleotide_fasta",
        ]
    ],
    on=[
        "genome",
        "contig",
    ],
    how="left",
    validate="many_to_one",
)


if len(
    observability
) != EXPECTED_FOCALS:

    fail(
        "Contig-length join changed focal row count."
    )


if observability[
    "contig_length_nt"
].isna().any():

    fail(
        "One or more focal contigs could not "
        "be assigned a nucleotide length."
    )


observability[
    "contig_length_nt"
] = (
    observability[
        "contig_length_nt"
    ]
    .astype(int)
)


## ================================================================== ##
## 6. Physical bp availability to contig boundaries
##
## Coordinates are 1-based and inclusive.
## ================================================================== ##

observability[
    "genomic_left_bp_available"
] = (
    observability[
        "start"
    ]
    - 1
)


observability[
    "genomic_right_bp_available"
] = (
    observability[
        "contig_length_nt"
    ]
    -
    observability[
        "end"
    ]
)


if (
    observability[
        "genomic_left_bp_available"
    ]
    < 0
).any():

    fail(
        "Negative left-side bp availability detected."
    )


if (
    observability[
        "genomic_right_bp_available"
    ]
    < 0
).any():

    fail(
        "A focal gene extends beyond its contig."
    )


## ================================================================== ##
## 7. Gene-count availability to contig boundaries
## ================================================================== ##

observability[
    "genomic_left_genes_available"
] = (
    observability[
        "contig_gene_rank"
    ]
    - 1
)


observability[
    "genomic_right_genes_available"
] = (
    observability[
        "contig_gene_count"
    ]
    -
    observability[
        "contig_gene_rank"
    ]
)


if (
    observability[
        "genomic_left_genes_available"
    ]
    < 0
).any():

    fail(
        "Negative left-side gene availability detected."
    )


if (
    observability[
        "genomic_right_genes_available"
    ]
    < 0
).any():

    fail(
        "Negative right-side gene availability detected."
    )


## ================================================================== ##
## 8. Orient upstream/downstream relative to focal strand
##
## '+' :
##     upstream   = genomic left
##     downstream = genomic right
##
## '-' :
##     upstream   = genomic right
##     downstream = genomic left
## ================================================================== ##

if not observability[
    "strand"
].isin(
    [
        "+",
        "-",
    ]
).all():

    fail(
        "Unexpected focal strand value."
    )


plus = (
    observability[
        "strand"
    ]
    == "+"
)


observability[
    "oriented_upstream_bp_available"
] = observability[
    "genomic_left_bp_available"
]


observability[
    "oriented_downstream_bp_available"
] = observability[
    "genomic_right_bp_available"
]


observability[
    "oriented_upstream_genes_available"
] = observability[
    "genomic_left_genes_available"
]


observability[
    "oriented_downstream_genes_available"
] = observability[
    "genomic_right_genes_available"
]


observability.loc[
    ~plus,
    "oriented_upstream_bp_available"
] = observability.loc[
    ~plus,
    "genomic_right_bp_available"
]


observability.loc[
    ~plus,
    "oriented_downstream_bp_available"
] = observability.loc[
    ~plus,
    "genomic_left_bp_available"
]


observability.loc[
    ~plus,
    "oriented_upstream_genes_available"
] = observability.loc[
    ~plus,
    "genomic_right_genes_available"
]


observability.loc[
    ~plus,
    "oriented_downstream_genes_available"
] = observability.loc[
    ~plus,
    "genomic_left_genes_available"
]


## ================================================================== ##
## 9. Explicit bp-window observability
## ================================================================== ##

for radius in [
    5_000,
    10_000,
    20_000,
]:

    label = (
        f"{radius // 1000}kb"
    )


    up_col = (
        f"upstream_{label}_observable"
    )


    down_col = (
        f"downstream_{label}_observable"
    )


    full_col = (
        f"full_{label}_context_observable"
    )


    observability[
        up_col
    ] = (
        observability[
            "oriented_upstream_bp_available"
        ]
        >= radius
    ).astype(int)


    observability[
        down_col
    ] = (
        observability[
            "oriented_downstream_bp_available"
        ]
        >= radius
    ).astype(int)


    observability[
        full_col
    ] = (
        (
            observability[
                up_col
            ]
            == 1
        )
        &
        (
            observability[
                down_col
            ]
            == 1
        )
    ).astype(int)


    ## Physical genomic left/right retained separately.

    observability[
        f"left_{label}_observable"
    ] = (
        observability[
            "genomic_left_bp_available"
        ]
        >= radius
    ).astype(int)


    observability[
        f"right_{label}_observable"
    ] = (
        observability[
            "genomic_right_bp_available"
        ]
        >= radius
    ).astype(int)


## ================================================================== ##
## 10. Explicit gene-window observability
## ================================================================== ##

for radius in [
    5,
    10,
    20,
]:

    label = (
        f"{radius}genes"
    )


    up_col = (
        f"upstream_{label}_observable"
    )


    down_col = (
        f"downstream_{label}_observable"
    )


    full_col = (
        f"full_{label}_context_observable"
    )


    observability[
        up_col
    ] = (
        observability[
            "oriented_upstream_genes_available"
        ]
        >= radius
    ).astype(int)


    observability[
        down_col
    ] = (
        observability[
            "oriented_downstream_genes_available"
        ]
        >= radius
    ).astype(int)


    observability[
        full_col
    ] = (
        (
            observability[
                up_col
            ]
            == 1
        )
        &
        (
            observability[
                down_col
            ]
            == 1
        )
    ).astype(int)


## ================================================================== ##
## 11. Useful combined observability for our extraction rule
##
## Later neighborhood extraction retains a gene when:
##
##     abs(gene offset) <= 20 genes
##
## OR
##
##     its CDS intersects focal CDS +/- 20 kb
##
## A focal has fully observable extraction context only when BOTH:
##
##     +/-20 genes
##     +/-20 kb
##
## are observable.
##
## This is deliberately conservative.
## ================================================================== ##

observability[
    "full_extraction_context_observable"
] = (
    (
        observability[
            "full_20genes_context_observable"
        ]
        == 1
    )
    &
    (
        observability[
            "full_20kb_context_observable"
        ]
        == 1
    )
).astype(int)


## ================================================================== ##
## 12. Final column ordering
## ================================================================== ##

first_columns = [
    "genome",
    "protein_id",
    "priority_reason",
    "priority_cyc2",
    "priority_ge5_hemes",

    "number_of_hemes",
    "fegenie_HMMs",

    "globdb_assignment_status",
    "globdb_top_scoring_cluster",
    "globdb_top_module",

    "contig",
    "contig_length_nt",

    "start",
    "end",
    "strand",

    "contig_gene_rank",
    "contig_gene_count",

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

    "source_nucleotide_fasta",
]


first_columns = [
    column

    for column
    in first_columns

    if column
    in observability.columns
]


remaining_columns = [
    column

    for column
    in observability.columns

    if column
    not in first_columns
]


observability = observability[
    first_columns
    +
    remaining_columns
]


observability = (
    observability
    .sort_values(
        [
            "priority_reason",
            "number_of_hemes",
            "genome",
            "contig",
            "contig_gene_rank",
            "protein_id",
        ],
        ascending=[
            True,
            False,
            True,
            True,
            True,
            True,
        ],
        kind="stable",
    )
    .reset_index(
        drop=True
    )
)


observability.to_csv(
    OUT_OBSERVABILITY,
    sep="\t",
    index=False,
    na_rep="",
)


## ================================================================== ##
## 13. QC summary
## ================================================================== ##

n_fasta_contigs = len(
    contig_lengths
)


n_fasta_contigs_without_genes = (
    contig_lengths[
        [
            "genome",
            "contig",
        ]
    ]
    .merge(
        gene_contigs,
        on=[
            "genome",
            "contig",
        ],
        how="left",
        indicator=True,
    )
    .query(
        "_merge == 'left_only'"
    )
    .shape[0]
)


qc_rows = [
    (
        "MAGs",
        len(
            genomes
        ),
    ),

    (
        "nucleotide_FASTAs",
        len(
            genome_fastas
        ),
    ),

    (
        "FASTA_contigs_total",
        n_fasta_contigs,
    ),

    (
        "contigs_with_Prodigal_genes",
        len(
            gene_contigs
        ),
    ),

    (
        "FASTA_contigs_without_Prodigal_genes",
        n_fasta_contigs_without_genes,
    ),

    (
        "Prodigal_genes_validated",
        len(
            catalog
        ),
    ),

    (
        "missing_gene_contigs",
        len(
            missing_contigs
        ),
    ),

    (
        "coordinate_violations",
        len(
            bad_coordinates
        ),
    ),

    (
        "priority_focals",
        len(
            observability
        ),
    ),
]


for label in [
    "5kb",
    "10kb",
    "20kb",
]:

    for direction in [
        "upstream",
        "downstream",
    ]:

        column = (
            f"{direction}_{label}_observable"
        )


        qc_rows.append(
            (
                f"focals_{direction}_{label}_observable",
                int(
                    observability[
                        column
                    ].sum()
                ),
            )
        )


    qc_rows.append(
        (
            f"focals_full_{label}_context_observable",
            int(
                observability[
                    f"full_{label}_context_observable"
                ].sum()
            ),
        )
    )


for label in [
    "5genes",
    "10genes",
    "20genes",
]:

    for direction in [
        "upstream",
        "downstream",
    ]:

        column = (
            f"{direction}_{label}_observable"
        )


        qc_rows.append(
            (
                f"focals_{direction}_{label}_observable",
                int(
                    observability[
                        column
                    ].sum()
                ),
            )
        )


    qc_rows.append(
        (
            f"focals_full_{label}_context_observable",
            int(
                observability[
                    f"full_{label}_context_observable"
                ].sum()
            ),
        )
    )


qc_rows.append(
    (
        "focals_full_extraction_context_observable",
        int(
            observability[
                "full_extraction_context_observable"
            ].sum()
        ),
    )
)


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
print("FASTA / coordinate validation")

print(
    f"  MAGs:                           "
    f"{len(genomes):,}"
)

print(
    f"  FASTA contigs:                  "
    f"{n_fasta_contigs:,}"
)

print(
    f"  Contigs with Prodigal genes:    "
    f"{len(gene_contigs):,}"
)

print(
    f"  FASTA contigs without genes:    "
    f"{n_fasta_contigs_without_genes:,}"
)

print(
    f"  Prodigal genes validated:       "
    f"{len(catalog):,}"
)

print(
    f"  Missing gene contigs:           "
    f"{len(missing_contigs):,}"
)

print(
    f"  Coordinate violations:          "
    f"{len(bad_coordinates):,}"
)


print()
print("Priority focal observability")


for label in [
    "5kb",
    "10kb",
    "20kb",
]:

    print(
        f"  Full +/-{label:<7} context:     "
        f"{int(observability[f'full_{label}_context_observable'].sum()):>2}"
        f"/{len(observability)}"
    )


for label in [
    "5genes",
    "10genes",
    "20genes",
]:

    pretty = label.replace(
        "genes",
        " genes"
    )

    print(
        f"  Full +/-{pretty:<8} context:    "
        f"{int(observability[f'full_{label}_context_observable'].sum()):>2}"
        f"/{len(observability)}"
    )


print(
    f"  Full extraction context:        "
    f"{int(observability['full_extraction_context_observable'].sum()):>2}"
    f"/{len(observability)}"
)


print()
print("=" * 80)
print("STAGE 51 COMPLETE")
print("=" * 80)

print(
    f"Contig lengths:\n  "
    f"{OUT_CONTIG_LENGTHS}"
)

print(
    f"Focal observability:\n  "
    f"{OUT_OBSERVABILITY}"
)

print(
    f"QC:\n  "
    f"{OUT_QC}"
)
