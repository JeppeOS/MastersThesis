#!/usr/bin/env python3

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd


## ================================================================== ##
## STAGE 53
##
## PREPARE UNIQUE PRIORITY-NEIGHBORHOOD PROTEINS FOR ANNOTATION
##
## Purpose
## -------
##
## Recover the amino-acid sequence of every unique actual gene
## represented in the 21 Stage-52 priority neighborhoods.
##
## Outputs:
##
##   priority_neighborhood_proteins.faa
##   priority_neighborhood_protein_metadata.tsv
##
## Protein FASTA identifiers are replaced with stable unique IDs:
##
##     NBP000001
##     NBP000002
##     ...
##
## genome + original protein_id remain in the metadata sidecar.
##
## This avoids relying on contig-derived protein IDs being globally
## unique across MAGs.
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


UNIQUE_GENES = (
    WORK_DIR
    / "priority_neighborhood_unique_genes.tsv"
)


ORF_DIR = (
    PROJECT
    / "functional_analysis"
    / "FeGenie_SemiBin2_full_run_20260918"
    / "fegenie_results"
    / "ORF_calls"
)


OUT_FASTA = (
    WORK_DIR
    / "priority_neighborhood_proteins.faa"
)


OUT_METADATA = (
    WORK_DIR
    / "priority_neighborhood_protein_metadata.tsv"
)


OUT_QC = (
    WORK_DIR
    / "priority_neighborhood_protein_export_qc.tsv"
)


EXPECTED_UNIQUE_GENES = 845


## ================================================================== ##
## Utilities
## ================================================================== ##

def fail(message):

    print(
        f"\nERROR: {message}",
        file=sys.stderr,
    )

    sys.exit(1)


def normalise_genome_filename(path):

    name = path.name

    suffixes = (
        "-proteins.faa",
        ".fasta",
        ".fna",
        ".faa",
        ".fa",
    )


    changed = True


    while changed:

        changed = False

        for suffix in suffixes:

            if name.endswith(
                suffix
            ):

                name = name[
                    :-len(suffix)
                ]

                changed = True

                break


    return name


def fasta_records(path):

    header = None

    sequence_parts = []


    with path.open(
        "r",
        encoding="utf-8",
    ) as handle:

        for raw in handle:

            line = raw.rstrip(
                "\r\n"
            )


            if not line:

                continue


            if line.startswith(
                ">"
            ):

                if header is not None:

                    yield (
                        header,
                        "".join(
                            sequence_parts
                        ),
                    )


                header = line[
                    1:
                ]


                sequence_parts = []


            else:

                if header is None:

                    fail(
                        "Sequence encountered before "
                        f"FASTA header in:\n{path}"
                    )


                sequence_parts.append(
                    line.strip()
                )


    if header is not None:

        yield (
            header,
            "".join(
                sequence_parts
            ),
        )


## ================================================================== ##
## Start
## ================================================================== ##

print("=" * 80)
print("STAGE 53 - PREPARE PRIORITY NEIGHBORHOOD PROTEINS")
print("=" * 80)


if not UNIQUE_GENES.is_file():

    fail(
        f"Missing Stage-52 unique-gene table:\n"
        f"{UNIQUE_GENES}"
    )


if not ORF_DIR.is_dir():

    fail(
        f"Missing FeGenie ORF directory:\n"
        f"{ORF_DIR}"
    )


## ================================================================== ##
## 1. Read neighborhood gene catalogue
## ================================================================== ##

genes = pd.read_csv(
    UNIQUE_GENES,
    sep="\t",
    dtype=str,
    keep_default_na=False,
)


required = {
    "genome",
    "protein_id",
}


missing = (
    required
    -
    set(
        genes.columns
    )
)


if missing:

    fail(
        "Unique neighborhood-gene table is missing:\n"
        +
        "\n".join(
            sorted(
                missing
            )
        )
    )


if len(
    genes
) != EXPECTED_UNIQUE_GENES:

    fail(
        f"Expected {EXPECTED_UNIQUE_GENES} unique "
        f"neighborhood genes but found "
        f"{len(genes)}."
    )


if genes[
    [
        "genome",
        "protein_id",
    ]
].duplicated().any():

    fail(
        "Duplicate genome + protein_id keys "
        "in unique neighborhood-gene table."
    )


target_keys = set(
    zip(
        genes[
            "genome"
        ],
        genes[
            "protein_id"
        ],
    )
)


target_genomes = set(
    genes[
        "genome"
    ]
)


print(
    f"Unique neighborhood genes: "
    f"{len(genes):,}"
)

print(
    f"MAGs represented:          "
    f"{len(target_genomes):,}"
)


## ================================================================== ##
## 2. Locate authoritative proteomes
## ================================================================== ##

genome_to_faa = {}


for path in sorted(
    ORF_DIR.glob(
        "*-proteins.faa"
    )
):

    if (
        not path.is_file()
        or
        path.stat().st_size == 0
    ):

        continue


    genome = normalise_genome_filename(
        path
    )


    if genome not in target_genomes:

        continue


    if genome in genome_to_faa:

        fail(
            f"More than one protein FASTA found "
            f"for genome:\n{genome}"
        )


    genome_to_faa[
        genome
    ] = path


missing_genomes = sorted(
    target_genomes
    -
    set(
        genome_to_faa
    )
)


if missing_genomes:

    fail(
        "Neighborhood genomes lacking proteomes:\n"
        +
        "\n".join(
            missing_genomes
        )
    )


print(
    f"Authoritative proteomes:    "
    f"{len(genome_to_faa):,}"
)


## ================================================================== ##
## 3. Recover target protein sequences
## ================================================================== ##

sequence_lookup = {}


for number, genome in enumerate(
    sorted(
        target_genomes
    ),
    start=1,
):

    faa = genome_to_faa[
        genome
    ]


    genome_target_ids = set(
        genes.loc[
            genes[
                "genome"
            ]
            ==
            genome,
            "protein_id",
        ]
    )


    recovered = set()


    for full_header, sequence in fasta_records(
        faa
    ):

        protein_id = (
            full_header
            .split()[0]
        )


        if protein_id not in genome_target_ids:

            continue


        key = (
            genome,
            protein_id,
        )


        if key in sequence_lookup:

            fail(
                "Protein recovered more than once:\n"
                f"{genome}\t{protein_id}"
            )


        if sequence.endswith(
            "*"
        ):

            sequence = sequence[
                :-1
            ]


        if "*" in sequence:

            fail(
                f"Internal stop symbol in "
                f"{genome} / {protein_id}"
            )


        if not sequence:

            fail(
                f"Empty amino-acid sequence for "
                f"{genome} / {protein_id}"
            )


        sequence_lookup[
            key
        ] = sequence


        recovered.add(
            protein_id
        )


    missing_ids = sorted(
        genome_target_ids
        -
        recovered
    )


    if missing_ids:

        fail(
            f"{len(missing_ids)} neighborhood proteins "
            f"were not recovered from {genome}.\n"
            f"First missing IDs:\n"
            +
            "\n".join(
                missing_ids[
                    :20
                ]
            )
        )


    print(
        f"  {number:>2}/"
        f"{len(target_genomes)} "
        f"{genome}: "
        f"{len(recovered):,} proteins"
    )


## ================================================================== ##
## 4. Global recovery QC
## ================================================================== ##

missing_keys = sorted(
    target_keys
    -
    set(
        sequence_lookup
    )
)


unexpected_keys = sorted(
    set(
        sequence_lookup
    )
    -
    target_keys
)


if missing_keys:

    fail(
        f"{len(missing_keys)} target proteins "
        f"were not recovered."
    )


if unexpected_keys:

    fail(
        f"{len(unexpected_keys)} unexpected proteins "
        f"were recovered."
    )


if len(
    sequence_lookup
) != EXPECTED_UNIQUE_GENES:

    fail(
        f"Recovered {len(sequence_lookup)} sequences; "
        f"expected {EXPECTED_UNIQUE_GENES}."
    )


## ================================================================== ##
## 5. Assign stable annotation IDs
## ================================================================== ##

genes = (
    genes
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


genes[
    "neighborhood_protein_id"
] = [
    f"NBP{index:06d}"

    for index
    in range(
        1,
        len(
            genes
        )
        + 1,
    )
]


genes[
    "annotation_sequence_length"
] = [
    len(
        sequence_lookup[
            (
                row.genome,
                row.protein_id,
            )
        ]
    )

    for row
    in genes[
        [
            "genome",
            "protein_id",
        ]
    ].itertuples(
        index=False
    )
]


## ================================================================== ##
## 6. Write FASTA
## ================================================================== ##

with OUT_FASTA.open(
    "w",
    encoding="utf-8",
) as handle:

    for row in genes[
        [
            "neighborhood_protein_id",
            "genome",
            "protein_id",
        ]
    ].itertuples(
        index=False
    ):

        sequence = sequence_lookup[
            (
                row.genome,
                row.protein_id,
            )
        ]


        handle.write(
            f">{row.neighborhood_protein_id}\n"
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
                    start + 80
                ]
                +
                "\n"
            )


## ================================================================== ##
## 7. Write metadata
## ================================================================== ##

first_columns = [
    "neighborhood_protein_id",
    "genome",
    "protein_id",
    "annotation_sequence_length",
    "contig",
    "contig_gene_rank",
    "start",
    "end",
    "strand",
    "n_priority_neighborhoods",

    "is_integrated_candidate",

    "fegenie_positive",
    "fegenie_HMMs",

    "findmehemes_positive",
    "number_of_hemes",

    "signalp_prediction",
    "deeptmhmm_class",
    "localization_class",

    "globdb_top_scoring_cluster",
    "globdb_top_module",

    "priority_reason",
]


first_columns = [
    column

    for column
    in first_columns

    if column
    in genes.columns
]


remaining_columns = [
    column

    for column
    in genes.columns

    if column
    not in first_columns
]


genes = genes[
    first_columns
    +
    remaining_columns
]


genes.to_csv(
    OUT_METADATA,
    sep="\t",
    index=False,
    na_rep="",
)


## ================================================================== ##
## 8. QC
## ================================================================== ##

qc = pd.DataFrame(
    [
        (
            "unique_neighborhood_genes",
            len(
                genes
            ),
        ),

        (
            "proteins_recovered",
            len(
                sequence_lookup
            ),
        ),

        (
            "MAGs_represented",
            len(
                target_genomes
            ),
        ),

        (
            "stable_FASTA_ids",
            genes[
                "neighborhood_protein_id"
            ].nunique(),
        ),

        (
            "minimum_protein_length",
            genes[
                "annotation_sequence_length"
            ].min(),
        ),

        (
            "maximum_protein_length",
            genes[
                "annotation_sequence_length"
            ].max(),
        ),

        (
            "mean_protein_length",
            f"{genes['annotation_sequence_length'].mean():.2f}",
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
print("STAGE 53 COMPLETE")
print("=" * 80)

print(
    f"Neighborhood proteins:  "
    f"{len(genes):,}"
)

print(
    f"Sequences recovered:    "
    f"{len(sequence_lookup):,}"
)

print(
    f"MAGs represented:       "
    f"{len(target_genomes):,}"
)

print(
    f"Minimum protein length: "
    f"{genes['annotation_sequence_length'].min():,} aa"
)

print(
    f"Maximum protein length: "
    f"{genes['annotation_sequence_length'].max():,} aa"
)

print()
print(
    f"FASTA:\n  "
    f"{OUT_FASTA}"
)

print(
    f"Metadata:\n  "
    f"{OUT_METADATA}"
)

print(
    f"QC:\n  "
    f"{OUT_QC}"
)
