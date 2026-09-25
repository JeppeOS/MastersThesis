#!/usr/bin/env python3

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd


## ================================================================== ##
## STAGE 54
##
## BUILD FIXED GLOBDB METHYLOCOCCALES ANNOTATION REFERENCE
##
## Purpose
## -------
##
## Build a protein reference from the same 631 Methylococcales genomes
## used in the original GlobDB workflow.
##
## The corresponding Stage-11 GlobDB annotation fields are retained:
##
##     globdb_cog
##     globdb_gene
##     globdb_product
##
## No new annotation is generated.
## ================================================================== ##


BASE = (
    Path.home()
    / "methanotrophs"
    / "methanotroph_project"
    / "jeppe"
)


OLD_WORKFLOW = (
    BASE
    / "genome_analysis_workflow"
)


OLD_ORF_DIR = (
    OLD_WORKFLOW
    / "FeGenie_conda_full_run_20260803"
    / "fegenie_results"
    / "ORF_calls"
)


OLD_CATALOG = (
    OLD_WORKFLOW
    / "11_genome_annotations"
    / "annotated_gene_catalog.tsv"
)


PROJECT = (
    BASE
    / "metagenome"
)


OUT_DIR = (
    PROJECT
    / "functional_analysis"
    / "methylococcales_neighborhoods"
    / "globdb_annotation_transfer"
)


REFERENCE_FASTA = (
    OUT_DIR
    / "GlobDB631_all_proteins.faa"
)


REFERENCE_METADATA = (
    OUT_DIR
    / "GlobDB631_annotation_metadata.tsv"
)


QC_OUT = (
    OUT_DIR
    / "GlobDB631_reference_qc.tsv"
)


EXPECTED_GENOMES = 631

EXPECTED_PROTEINS = 2_002_656


def fail(message):

    print(
        f"\nERROR: {message}",
        file=sys.stderr,
    )

    sys.exit(1)


def fasta_records(path):

    header = None
    parts = []

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

            if line.startswith(">"):

                if header is not None:

                    yield (
                        header,
                        "".join(parts),
                    )

                header = line[1:]

                parts = []

            else:

                if header is None:

                    fail(
                        f"Sequence before header in:\n{path}"
                    )

                parts.append(
                    line.strip()
                )

    if header is not None:

        yield (
            header,
            "".join(parts),
        )


def normalise_genome(path):

    name = path.name

    suffixes = (
        "-proteins.faa",
        ".faa",
        ".fa",
        ".fna",
        ".fasta",
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


print("=" * 80)
print("STAGE 54 - BUILD GLOBDB ANNOTATION REFERENCE")
print("=" * 80)


if not OLD_ORF_DIR.is_dir():

    fail(
        f"Missing old ORF directory:\n{OLD_ORF_DIR}"
    )


if not OLD_CATALOG.is_file():

    fail(
        f"Missing Stage-11 gene catalogue:\n{OLD_CATALOG}"
    )


OUT_DIR.mkdir(
    parents=True,
    exist_ok=True,
)


## ================================================================== ##
## 1. Read existing Stage-11 annotation metadata
## ================================================================== ##

print()
print("Reading original Stage-11 gene catalogue...")


catalog = pd.read_csv(
    OLD_CATALOG,
    sep="\t",
    dtype=str,
    keep_default_na=False,
)


required = {
    "genome",
    "protein_id",
    "globdb_cog",
    "globdb_gene",
    "globdb_product",
}


missing = (
    required
    -
    set(
        catalog.columns
    )
)


if missing:

    fail(
        "Old gene catalogue is missing:\n"
        +
        "\n".join(
            sorted(
                missing
            )
        )
    )


if len(
    catalog
) != EXPECTED_PROTEINS:

    fail(
        f"Expected {EXPECTED_PROTEINS:,} old catalogue genes "
        f"but found {len(catalog):,}."
    )


if catalog[
    [
        "genome",
        "protein_id",
    ]
].duplicated().any():

    fail(
        "Duplicate genome + protein_id keys "
        "in old Stage-11 catalogue."
    )


metadata_columns = [
    "genome",
    "protein_id",
    "globdb_cog",
    "globdb_gene",
    "globdb_product",

    "annotation_match_type",
    "annotation_accepted",
    "reciprocal_overlap",
]


metadata_columns = [
    column

    for column
    in metadata_columns

    if column
    in catalog.columns
]


metadata = catalog[
    metadata_columns
].copy()


metadata[
    "has_globdb_annotation"
] = (
    (
        metadata[
            "globdb_cog"
        ].str.strip()
        != ""
    )
    |
    (
        metadata[
            "globdb_gene"
        ].str.strip()
        != ""
    )
    |
    (
        metadata[
            "globdb_product"
        ].str.strip()
        != ""
    )
).astype(int)


metadata.to_csv(
    REFERENCE_METADATA,
    sep="\t",
    index=False,
)


print(
    f"  Stage-11 proteins:       "
    f"{len(metadata):,}"
)

print(
    f"  With COG/gene/product:   "
    f"{int(metadata['has_globdb_annotation'].sum()):,}"
)


## ================================================================== ##
## 2. Discover authoritative old proteomes
## ================================================================== ##

proteomes = sorted(
    path

    for path
    in OLD_ORF_DIR.glob(
        "*-proteins.faa"
    )

    if (
        path.is_file()
        and
        path.stat().st_size > 0
    )
)


if len(
    proteomes
) != EXPECTED_GENOMES:

    fail(
        f"Expected {EXPECTED_GENOMES} old proteomes "
        f"but found {len(proteomes)}."
    )


print(
    f"  Old proteomes:           "
    f"{len(proteomes):,}"
)


## ================================================================== ##
## 3. Build combined reference FASTA
##
## Protein IDs are retained unchanged because the original GlobDB /
## Prodigal IDs are globally genome-prefixed.
## ================================================================== ##

print()
print("Building combined GlobDB protein reference...")


n_proteins = 0

seen_ids = set()


with REFERENCE_FASTA.open(
    "w",
    encoding="utf-8",
) as output:

    for number, proteome in enumerate(
        proteomes,
        start=1,
    ):

        genome = normalise_genome(
            proteome
        )


        genome_count = 0


        for full_header, sequence in fasta_records(
            proteome
        ):

            protein_id = (
                full_header
                .split()[0]
            )


            if not protein_id:

                fail(
                    f"Empty protein ID in:\n{proteome}"
                )


            if protein_id in seen_ids:

                fail(
                    f"Protein ID is not globally unique:\n"
                    f"{protein_id}"
                )


            seen_ids.add(
                protein_id
            )


            if sequence.endswith(
                "*"
            ):

                sequence = sequence[
                    :-1
                ]


            if "*" in sequence:

                fail(
                    f"Internal stop in old reference protein:\n"
                    f"{protein_id}"
                )


            if not sequence:

                fail(
                    f"Empty sequence:\n"
                    f"{protein_id}"
                )


            output.write(
                f">{protein_id}\n"
            )


            for start in range(
                0,
                len(sequence),
                80,
            ):

                output.write(
                    sequence[
                        start:
                        start + 80
                    ]
                    +
                    "\n"
                )


            n_proteins += 1

            genome_count += 1


        if (
            number % 50 == 0
            or
            number == len(
                proteomes
            )
        ):

            print(
                f"  Processed {number:>3}/"
                f"{len(proteomes)} proteomes"
            )


if n_proteins != EXPECTED_PROTEINS:

    fail(
        f"Combined FASTA contains {n_proteins:,} proteins; "
        f"expected {EXPECTED_PROTEINS:,}."
    )


## ================================================================== ##
## 4. Validate FASTA IDs against Stage-11 annotation catalogue
## ================================================================== ##

catalog_ids = set(
    metadata[
        "protein_id"
    ]
)


missing_metadata = (
    seen_ids
    -
    catalog_ids
)


missing_fasta = (
    catalog_ids
    -
    seen_ids
)


if missing_metadata:

    fail(
        f"{len(missing_metadata):,} FASTA proteins lack "
        f"Stage-11 metadata."
    )


if missing_fasta:

    fail(
        f"{len(missing_fasta):,} Stage-11 proteins are "
        f"missing from the combined FASTA."
    )


## ================================================================== ##
## 5. QC
## ================================================================== ##

qc = pd.DataFrame(
    [
        (
            "old_proteomes",
            len(
                proteomes
            ),
        ),

        (
            "reference_proteins",
            n_proteins,
        ),

        (
            "Stage11_metadata_rows",
            len(
                metadata
            ),
        ),

        (
            "proteins_with_GlobDB_annotation",
            int(
                metadata[
                    "has_globdb_annotation"
                ].sum()
            ),
        ),

        (
            "FASTA_proteins_missing_metadata",
            len(
                missing_metadata
            ),
        ),

        (
            "metadata_proteins_missing_FASTA",
            len(
                missing_fasta
            ),
        ),
    ],
    columns=[
        "metric",
        "value",
    ],
)


qc.to_csv(
    QC_OUT,
    sep="\t",
    index=False,
)


print()
print("=" * 80)
print("STAGE 54 COMPLETE")
print("=" * 80)

print(
    f"Reference proteins:       "
    f"{n_proteins:,}"
)

print(
    f"Proteins with annotation: "
    f"{int(metadata['has_globdb_annotation'].sum()):,}"
)

print()
print(
    f"Reference FASTA:\n  "
    f"{REFERENCE_FASTA}"
)

print(
    f"Reference metadata:\n  "
    f"{REFERENCE_METADATA}"
)

print(
    f"QC:\n  "
    f"{QC_OUT}"
)
