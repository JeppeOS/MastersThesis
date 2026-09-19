#!/usr/bin/env python3

from pathlib import Path
import gzip
import sys

import pandas as pd


## ================================================================== ##
## STAGE 12B
##
## Determine nucleotide contig lengths and quantify how much genomic
## sequence is actually observable around every focal module protein.
##
## This stage does NOT alter the Stage-12A neighborhood observations.
## It adds censoring information that will later be used when calculating
## neighborhood conservation.
##
## Example:
##
##   If a focal gene is only 3 kb from a contig end, absence of a gene
##   10 kb in that direction cannot be treated as biological absence.
## ================================================================== ##


## ================================================================== ##
## Paths
## ================================================================== ##

HERE = Path(__file__).resolve().parent
WORKFLOW = HERE.parent


FASTA_DIR = Path(
    "~/methanotrophs/methanotroph_project/jeppe/"
    "methylococcales_fastas"
).expanduser()


GENE_CATALOG = (
    WORKFLOW
    / "11_genome_annotations"
    / "annotated_gene_catalog.tsv"
)


FOCAL_OCCURRENCES = (
    HERE
    / "focal_module_occurrences.tsv"
)


OUT_CONTIG_LENGTHS = (
    HERE
    / "contig_lengths.tsv"
)


OUT_OBSERVABILITY = (
    HERE
    / "focal_context_observability.tsv"
)


OUT_QC = (
    HERE
    / "context_observability_qc.tsv"
)


OUT_MISSING_CONTIGS = (
    HERE
    / "missing_catalog_contigs.tsv"
)


OUT_BAD_COORDINATES = (
    HERE
    / "genes_exceeding_contig_length.tsv"
)


## ================================================================== ##
## Expected upstream values
## ================================================================== ##

EXPECTED_GENOMES = 631
EXPECTED_GENES = 2_002_656
EXPECTED_GENE_CONTIGS = 176_128

EXPECTED_MODULE_PROTEINS = 10_537
EXPECTED_FOCAL_CONTIGS = 6_941

CHUNK_SIZE = 250_000


## ================================================================== ##
## Utility functions
## ================================================================== ##

def fail(message):

    print(
        f"\nERROR: {message}",
        file=sys.stderr
    )

    sys.exit(1)


def require_columns(
    df,
    required,
    source,
):

    missing = (
        set(required)
        -
        set(df.columns)
    )

    if missing:

        fail(
            f"{source} is missing required column(s):\n"
            +
            ", ".join(
                sorted(missing)
            )
        )


def find_genome_fasta(
    genome,
):

    candidates = [
        FASTA_DIR / f"{genome}.fa.gz",
        FASTA_DIR / f"{genome}.fna.gz",
        FASTA_DIR / f"{genome}.fa",
        FASTA_DIR / f"{genome}.fna",
    ]


    found = [
        path
        for path
        in candidates
        if path.exists()
    ]


    if len(found) == 0:

        fail(
            f"No FASTA found for genome "
            f"{genome}."
        )


    if len(found) > 1:

        fail(
            f"Multiple possible FASTAs found "
            f"for genome {genome}:\n"
            +
            "\n".join(
                str(x)
                for x
                in found
            )
        )


    return found[0]


def open_fasta(
    path,
):

    if path.suffix == ".gz":

        return gzip.open(
            path,
            "rt"
        )

    return open(
        path,
        "r"
    )


def parse_fasta_lengths(
    genome,
    path,
):

    """
    Stream one FASTA and return:

        contig
        contig_length_nt

    Headers such as:

        > GCA_002134785_000000000001

    are handled deliberately. The contig ID is the first
    whitespace-delimited token after removing '>'.
    """

    rows = []

    seen = set()

    current_contig = None
    current_length = 0


    def flush_current():

        nonlocal current_contig
        nonlocal current_length


        if current_contig is None:
            return


        if current_contig in seen:

            fail(
                f"Duplicate FASTA contig ID "
                f"{current_contig} in:\n"
                f"{path}"
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

                "source_fasta":
                    str(
                        path
                    ),
            }
        )


    with open_fasta(
        path
    ) as handle:

        for raw_line in handle:

            if raw_line.startswith(">"):

                flush_current()


                header = (
                    raw_line[1:]
                    .strip()
                )


                if not header:

                    fail(
                        f"Empty FASTA header in:\n"
                        f"{path}"
                    )


                current_contig = (
                    header
                    .split()[0]
                )


                current_length = 0

                continue


            sequence = "".join(
                raw_line.split()
            )


            if not sequence:
                continue


            if current_contig is None:

                fail(
                    f"Sequence encountered before "
                    f"first FASTA header in:\n"
                    f"{path}"
                )


            current_length += len(
                sequence
            )


    flush_current()


    if len(rows) == 0:

        fail(
            f"No sequences parsed from:\n"
            f"{path}"
        )


    return rows


## ================================================================== ##
## Start
## ================================================================== ##

print("=" * 80)
print("STAGE 12B - BUILD CONTIG OBSERVABILITY")
print("=" * 80)


## ================================================================== ##
## 1. Read focal module occurrences
## ================================================================== ##

print()
print(
    "Reading focal module occurrences..."
)


focals = pd.read_csv(
    FOCAL_OCCURRENCES,
    sep="\t",
    dtype=str,
    keep_default_na=False
)


require_columns(
    focals,
    [
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
    ],
    FOCAL_OCCURRENCES.name
)


if len(focals) != EXPECTED_MODULE_PROTEINS:

    fail(
        f"Expected "
        f"{EXPECTED_MODULE_PROTEINS:,} "
        f"focal proteins but found "
        f"{len(focals):,}."
    )


if (
    focals[
        [
            "genome",
            "focal_protein_id",
        ]
    ]
    .duplicated()
    .any()
):

    fail(
        "Duplicate focal-protein keys."
    )


genomes = sorted(
    focals[
        "genome"
    ]
    .unique()
)


if len(genomes) != EXPECTED_GENOMES:

    fail(
        f"Expected {EXPECTED_GENOMES} "
        f"genomes but found "
        f"{len(genomes)}."
    )


n_focal_contigs = (
    focals[
        [
            "genome",
            "contig",
        ]
    ]
    .drop_duplicates()
    .shape[0]
)


if n_focal_contigs != EXPECTED_FOCAL_CONTIGS:

    fail(
        f"Expected "
        f"{EXPECTED_FOCAL_CONTIGS:,} "
        f"focal contigs but found "
        f"{n_focal_contigs:,}."
    )


print(
    f"  Focal proteins: "
    f"{len(focals):,}"
)

print(
    f"  Genomes:        "
    f"{len(genomes):,}"
)

print(
    f"  Focal contigs:  "
    f"{n_focal_contigs:,}"
)


## ================================================================== ##
## 2. Locate exactly one original nucleotide FASTA per genome
## ================================================================== ##

print()
print(
    "Locating original GlobDB genome FASTAs..."
)


genome_fastas = {
    genome:
    find_genome_fasta(
        genome
    )

    for genome
    in genomes
}


print(
    f"  FASTAs found: "
    f"{len(genome_fastas):,}"
)


## ================================================================== ##
## 3. Parse nucleotide length of every contig
## ================================================================== ##

print()
print(
    "Parsing contig lengths..."
)


contig_rows = []


for number, genome in enumerate(
    genomes,
    start=1
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


    if (
        number % 50 == 0
        or
        number == len(
            genomes
        )
    ):

        print(
            f"  Parsed "
            f"{number:>3}/"
            f"{len(genomes)} genomes"
        )


contig_lengths = pd.DataFrame(
    contig_rows
)


if (
    contig_lengths[
        [
            "genome",
            "contig",
        ]
    ]
    .duplicated()
    .any()
):

    fail(
        "Duplicate genome + contig keys "
        "in FASTA-derived length table."
    )


if (
    contig_lengths[
        "contig_length_nt"
    ]
    <=
    0
).any():

    fail(
        "One or more FASTA contigs have "
        "non-positive sequence length."
    )


contig_lengths = (
    contig_lengths
    .sort_values(
        [
            "genome",
            "contig",
        ],
        kind="stable"
    )
    .reset_index(
        drop=True
    )
)


contig_lengths.to_csv(
    OUT_CONTIG_LENGTHS,
    sep="\t",
    index=False
)


print(
    f"  FASTA contigs:  "
    f"{len(contig_lengths):,}"
)


## ================================================================== ##
## 4. Independently summarize all contigs used by Prodigal genes
##
## We verify:
##
##   - all 2,002,656 genes are seen
##   - all 176,128 gene-containing contigs exist in the nucleotide FASTAs
##   - no Prodigal gene extends beyond the nucleotide contig
## ================================================================== ##

print()
print(
    "Validating all Prodigal gene coordinates "
    "against FASTA lengths..."
)


catalog_header = pd.read_csv(
    GENE_CATALOG,
    sep="\t",
    nrows=0
)


require_columns(
    catalog_header,
    [
        "genome",
        "contig",
        "start",
        "end",
    ],
    GENE_CATALOG.name
)


catalog_parts = []

total_gene_rows = 0


for chunk in pd.read_csv(
    GENE_CATALOG,
    sep="\t",
    usecols=[
        "genome",
        "contig",
        "start",
        "end",
    ],
    dtype={
        "genome":
            str,

        "contig":
            str,

        "start":
            int,

        "end":
            int,
    },
    chunksize=CHUNK_SIZE,
):

    total_gene_rows += len(
        chunk
    )


    summary = (
        chunk
        .groupby(
            [
                "genome",
                "contig",
            ],
            as_index=False
        )
        .agg(
            n_genes=(
                "start",
                "size"
            ),

            minimum_gene_start=(
                "start",
                "min"
            ),

            maximum_gene_end=(
                "end",
                "max"
            ),
        )
    )


    catalog_parts.append(
        summary
    )


if total_gene_rows != EXPECTED_GENES:

    fail(
        f"Expected "
        f"{EXPECTED_GENES:,} "
        f"catalog genes but read "
        f"{total_gene_rows:,}."
    )


catalog_contigs = pd.concat(
    catalog_parts,
    ignore_index=True
)


## Groups may cross chunk boundaries, so aggregate once more. ##

catalog_contigs = (
    catalog_contigs
    .groupby(
        [
            "genome",
            "contig",
        ],
        as_index=False
    )
    .agg(
        n_genes=(
            "n_genes",
            "sum"
        ),

        minimum_gene_start=(
            "minimum_gene_start",
            "min"
        ),

        maximum_gene_end=(
            "maximum_gene_end",
            "max"
        ),
    )
)


if len(catalog_contigs) != EXPECTED_GENE_CONTIGS:

    fail(
        f"Expected "
        f"{EXPECTED_GENE_CONTIGS:,} "
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
    validate="one_to_one"
)


missing_contigs = catalog_check[
    catalog_check[
        "contig_length_nt"
    ]
    .isna()
].copy()


missing_contigs.to_csv(
    OUT_MISSING_CONTIGS,
    sep="\t",
    index=False
)


if len(missing_contigs) > 0:

    fail(
        f"{len(missing_contigs):,} "
        f"Prodigal gene-containing contigs "
        f"were not found in the original "
        f"FASTA files."
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
        <
        1
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
    index=False
)


if len(bad_coordinates) > 0:

    fail(
        f"{len(bad_coordinates):,} "
        f"contigs contain Prodigal coordinates "
        f"outside the nucleotide FASTA bounds."
    )


print(
    f"  Prodigal genes validated: "
    f"{total_gene_rows:,}"
)

print(
    f"  Gene-containing contigs:   "
    f"{len(catalog_check):,}"
)

print(
    f"  Missing contigs:            "
    f"{len(missing_contigs):,}"
)

print(
    f"  Coordinate violations:      "
    f"{len(bad_coordinates):,}"
)


## ================================================================== ##
## 5. Join nucleotide contig lengths onto focal occurrences
## ================================================================== ##

print()
print(
    "Calculating focal context observability..."
)


for column in [
    "focal_start",
    "focal_end",
    "contig_gene_rank",
    "contig_gene_count",
]:

    focals[
        column
    ] = pd.to_numeric(
        focals[
            column
        ],
        errors="raise"
    ).astype(int)


observability = focals.merge(
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
    validate="many_to_one"
)


if len(observability) != EXPECTED_MODULE_PROTEINS:

    fail(
        "Contig-length join changed focal "
        "row count."
    )


if (
    observability[
        "contig_length_nt"
    ]
    .isna()
    .any()
):

    fail(
        "One or more focal contigs could not "
        "be assigned a FASTA length."
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
## 6. Genomic distance to contig boundaries
##
## Coordinates are 1-based and inclusive.
##
## For a CDS:
##
##     start = 100
##
## there are exactly 99 nucleotides before the CDS.
## ================================================================== ##

observability[
    "genomic_left_bp_available"
] = (
    observability[
        "focal_start"
    ]
    -
    1
)


observability[
    "genomic_right_bp_available"
] = (
    observability[
        "contig_length_nt"
    ]
    -
    observability[
        "focal_end"
    ]
)


if (
    observability[
        "genomic_left_bp_available"
    ]
    <
    0
).any():

    fail(
        "Negative left-side sequence availability "
        "detected."
    )


if (
    observability[
        "genomic_right_bp_available"
    ]
    <
    0
).any():

    fail(
        "A focal CDS extends beyond its "
        "nucleotide contig."
    )


## ================================================================== ##
## 7. Normalize upstream/downstream by focal strand
##
## '+' focal:
##
##     upstream   = genomic left
##     downstream = genomic right
##
## '-' focal:
##
##     upstream   = genomic right
##     downstream = genomic left
## ================================================================== ##

valid_strands = (
    observability[
        "focal_strand"
    ]
    .isin(
        [
            "+",
            "-",
        ]
    )
)


if not valid_strands.all():

    fail(
        "Unexpected focal strand value."
    )


plus = (
    observability[
        "focal_strand"
    ]
    ==
    "+"
)


observability[
    "oriented_upstream_bp_available"
] = (
    observability[
        "genomic_left_bp_available"
    ]
)


observability[
    "oriented_downstream_bp_available"
] = (
    observability[
        "genomic_right_bp_available"
    ]
)


observability.loc[
    ~plus,
    "oriented_upstream_bp_available"
] = (
    observability.loc[
        ~plus,
        "genomic_right_bp_available"
    ]
)


observability.loc[
    ~plus,
    "oriented_downstream_bp_available"
] = (
    observability.loc[
        ~plus,
        "genomic_left_bp_available"
    ]
)


## ================================================================== ##
## 8. Explicit observability at useful physical windows
## ================================================================== ##

for radius in [
    5_000,
    10_000,
    20_000,
]:

    label = (
        f"{radius // 1000}kb"
    )


    upstream_column = (
        f"upstream_{label}_observable"
    )


    downstream_column = (
        f"downstream_{label}_observable"
    )


    full_column = (
        f"full_{label}_context_observable"
    )


    observability[
        upstream_column
    ] = (
        observability[
            "oriented_upstream_bp_available"
        ]
        >=
        radius
    ).astype(int)


    observability[
        downstream_column
    ] = (
        observability[
            "oriented_downstream_bp_available"
        ]
        >=
        radius
    ).astype(int)


    observability[
        full_column
    ] = (
        (
            observability[
                upstream_column
            ]
            ==
            1
        )
        &
        (
            observability[
                downstream_column
            ]
            ==
            1
        )
    ).astype(int)


## ================================================================== ##
## 9. Useful genomic, non-oriented observability flags
##
## These are retained because some later analyses may care about
## physical left/right rather than transcription-oriented
## upstream/downstream.
## ================================================================== ##

for radius in [
    5_000,
    10_000,
    20_000,
]:

    label = (
        f"{radius // 1000}kb"
    )


    observability[
        f"left_{label}_observable"
    ] = (
        observability[
            "genomic_left_bp_available"
        ]
        >=
        radius
    ).astype(int)


    observability[
        f"right_{label}_observable"
    ] = (
        observability[
            "genomic_right_bp_available"
        ]
        >=
        radius
    ).astype(int)


## ================================================================== ##
## 10. Final ordering and output
## ================================================================== ##

first_columns = [
    "genome",
    "focal_protein_id",
    "focal_cluster",
    "focal_module",

    "contig",
    "contig_length_nt",

    "focal_start",
    "focal_end",
    "focal_strand",

    "contig_gene_rank",
    "contig_gene_count",

    "genomic_left_bp_available",
    "genomic_right_bp_available",

    "oriented_upstream_bp_available",
    "oriented_downstream_bp_available",

    "upstream_5kb_observable",
    "downstream_5kb_observable",
    "full_5kb_context_observable",

    "upstream_10kb_observable",
    "downstream_10kb_observable",
    "full_10kb_context_observable",

    "upstream_20kb_observable",
    "downstream_20kb_observable",
    "full_20kb_context_observable",

    "left_5kb_observable",
    "right_5kb_observable",

    "left_10kb_observable",
    "right_10kb_observable",

    "left_20kb_observable",
    "right_20kb_observable",
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


observability.to_csv(
    OUT_OBSERVABILITY,
    sep="\t",
    index=False
)


## ================================================================== ##
## 11. QC statistics
## ================================================================== ##

n_fasta_contigs = len(
    contig_lengths
)


n_fasta_contigs_without_prodigal_genes = (
    contig_lengths[
        [
            "genome",
            "contig",
        ]
    ]
    .merge(
        catalog_contigs[
            [
                "genome",
                "contig",
            ]
        ],
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
    [
        "genomes",
        len(
            genomes
        ),
    ],

    [
        "original_genome_fastas",
        len(
            genome_fastas
        ),
    ],

    [
        "fasta_contigs_total",
        n_fasta_contigs,
    ],

    [
        "contigs_with_prodigal_genes",
        len(
            catalog_contigs
        ),
    ],

    [
        "fasta_contigs_without_prodigal_genes",
        n_fasta_contigs_without_prodigal_genes,
    ],

    [
        "prodigal_genes_validated",
        total_gene_rows,
    ],

    [
        "missing_gene_contigs",
        len(
            missing_contigs
        ),
    ],

    [
        "gene_contigs_with_coordinate_violations",
        len(
            bad_coordinates
        ),
    ],

    [
        "focal_module_proteins",
        len(
            observability
        ),
    ],

    [
        "focal_contigs",
        n_focal_contigs,
    ],
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


        count = int(
            observability[
                column
            ]
            .sum()
        )


        qc_rows.append(
            [
                f"focals_{direction}_{label}_observable",
                count,
            ]
        )


    full_column = (
        f"full_{label}_context_observable"
    )


    qc_rows.append(
        [
            f"focals_full_{label}_context_observable",
            int(
                observability[
                    full_column
                ]
                .sum()
            ),
        ]
    )


qc = pd.DataFrame(
    qc_rows,
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
## 12. Terminal summary
## ================================================================== ##

print()
print(
    "FASTA / coordinate validation"
)

print(
    f"  Genomes:                         "
    f"{len(genomes):,}"
)

print(
    f"  FASTA contigs:                   "
    f"{n_fasta_contigs:,}"
)

print(
    f"  Contigs with Prodigal genes:     "
    f"{len(catalog_contigs):,}"
)

print(
    f"  FASTA contigs without genes:     "
    f"{n_fasta_contigs_without_prodigal_genes:,}"
)

print(
    f"  Prodigal genes validated:        "
    f"{total_gene_rows:,}"
)

print(
    f"  Missing gene contigs:            "
    f"{len(missing_contigs):,}"
)

print(
    f"  Coordinate violations:           "
    f"{len(bad_coordinates):,}"
)


print()
print(
    "Focal physical-context observability"
)


for label in [
    "5kb",
    "10kb",
    "20kb",
]:

    upstream = int(
        observability[
            f"upstream_{label}_observable"
        ]
        .sum()
    )

    downstream = int(
        observability[
            f"downstream_{label}_observable"
        ]
        .sum()
    )

    full = int(
        observability[
            f"full_{label}_context_observable"
        ]
        .sum()
    )


    print()
    print(
        f"  {label}:"
    )

    print(
        f"    upstream observable:           "
        f"{upstream:,}/"
        f"{len(observability):,} "
        f"({100 * upstream / len(observability):.2f}%)"
    )

    print(
        f"    downstream observable:         "
        f"{downstream:,}/"
        f"{len(observability):,} "
        f"({100 * downstream / len(observability):.2f}%)"
    )

    print(
        f"    complete both sides:           "
        f"{full:,}/"
        f"{len(observability):,} "
        f"({100 * full / len(observability):.2f}%)"
    )


print()
print("=" * 80)
print("SUCCESS")
print("=" * 80)


print(
    f"Contig lengths:      "
    f"{OUT_CONTIG_LENGTHS}"
)

print(
    f"Focal observability: "
    f"{OUT_OBSERVABILITY}"
)

print(
    f"QC:                  "
    f"{OUT_QC}"
)

print(
    f"Missing contigs QC:  "
    f"{OUT_MISSING_CONTIGS}"
)

print(
    f"Coordinate QC:       "
    f"{OUT_BAD_COORDINATES}"
)
