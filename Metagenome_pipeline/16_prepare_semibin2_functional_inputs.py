#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import gzip
import shutil


## ------------------------------------------------------------
## Paths
## ------------------------------------------------------------

PROJECT = (
    Path.home()
    / "methanotrophs"
    / "methanotroph_project"
    / "jeppe"
    / "metagenome"
)

MASTER = (
    PROJECT
    / "semibin2_analysis"
    / "MAG_master_table_SemiBin2.tsv"
)

COVERAGE = (
    PROJECT
    / "semibin2_analysis"
    / "SemiBin2_MAG_coverage_per_bin.tsv"
)

OUTROOT = (
    PROJECT
    / "functional_analysis"
)

INPUT_ROOT = (
    OUTROOT
    / "input_bins"
)

MANIFEST_OUT = (
    OUTROOT
    / "SemiBin2_187_MAG_manifest.tsv"
)


## ------------------------------------------------------------
## Read inputs
## ------------------------------------------------------------

master = pd.read_csv(
    MASTER,
    sep="\t"
)

coverage = pd.read_csv(
    COVERAGE,
    sep="\t"
)


if len(master) != 187:

    raise RuntimeError(
        f"Expected 187 Flye + SemiBin2 MAGs, "
        f"found {len(master)}."
    )


## ------------------------------------------------------------
## Merge genome properties
## ------------------------------------------------------------

cov_cols = [
    "Genome_ID",
    "Genome_length_bp",
    "N_contigs",
    "Mean_coverage_x"
]

coverage_small = coverage[
    cov_cols
].copy()


x = master.merge(
    coverage_small,
    on="Genome_ID",
    how="left",
    validate="one_to_one"
)


if x["Genome_length_bp"].isna().any():

    raise RuntimeError(
        "Some MAGs are missing coverage/genome statistics."
    )


## ------------------------------------------------------------
## Recreate staged FASTA directory only
## ------------------------------------------------------------

if INPUT_ROOT.exists():

    shutil.rmtree(
        INPUT_ROOT
    )


INPUT_ROOT.mkdir(
    parents=True,
    exist_ok=True
)


## ------------------------------------------------------------
## gzip detection
## ------------------------------------------------------------

def is_gzip(path):

    with open(
        path,
        "rb"
    ) as handle:

        return (
            handle.read(2)
            == b"\x1f\x8b"
        )


## ------------------------------------------------------------
## Copy/decompress FASTA
## ------------------------------------------------------------

def stage_fasta(
    source,
    destination
):

    source = Path(
        source
    )


    if not source.exists():

        raise RuntimeError(
            f"Missing source FASTA:\n{source}"
        )


    if is_gzip(
        source
    ):

        with gzip.open(
            source,
            "rb"
        ) as src:

            with open(
                destination,
                "wb"
            ) as dst:

                shutil.copyfileobj(
                    src,
                    dst
                )

    else:

        shutil.copy2(
            source,
            destination
        )


## ------------------------------------------------------------
## Stage all MAGs
## ------------------------------------------------------------

staged_paths = []


for _, row in x.iterrows():

    genome_id = str(
        row["Genome_ID"]
    )

    sample = str(
        row["Sample"]
    )

    source = Path(
        str(
            row["Original_path"]
        )
    )


    sample_dir = (
        INPUT_ROOT
        / sample
    )

    sample_dir.mkdir(
        parents=True,
        exist_ok=True
    )


    destination = (
        sample_dir
        / f"{genome_id}.fa"
    )


    stage_fasta(
        source,
        destination
    )


    staged_paths.append(
        str(
            destination
        )
    )


x[
    "Functional_input_FASTA"
] = staged_paths


## ------------------------------------------------------------
## Flag Methylococcales
## ------------------------------------------------------------

x[
    "Is_Methylococcales"
] = (
    x["Order"]
    .fillna("")
    .eq(
        "Methylococcales"
    )
)


## ------------------------------------------------------------
## Save manifest
## ------------------------------------------------------------

x = x.sort_values(
    [
        "Sample",
        "Genome_ID"
    ]
)


x.to_csv(
    MANIFEST_OUT,
    sep="\t",
    index=False,
    na_rep=""
)


## ------------------------------------------------------------
## Validation
## ------------------------------------------------------------

fastas = list(
    INPUT_ROOT.glob(
        "*/*.fa"
    )
)


if len(fastas) != 187:

    raise RuntimeError(
        f"Expected 187 staged FASTAs, "
        f"found {len(fastas)}."
    )


n_meth = int(
    x[
        "Is_Methylococcales"
    ].sum()
)


if n_meth != 16:

    raise RuntimeError(
        f"Expected 16 Methylococcales MAGs, "
        f"found {n_meth}."
    )


## ------------------------------------------------------------
## Print summary
## ------------------------------------------------------------

print()

print(
    "Functional-analysis input set"
)

print(
    "=" * 70
)

print(
    x.groupby(
        "Sample"
    )
    .size()
    .to_string()
)

print()

print(
    f"Total MAGs:            {len(x)}"
)

print(
    f"Methylococcales MAGs:  {n_meth}"
)

print()

print(
    f"Input directory:\n"
    f"{INPUT_ROOT}"
)

print()

print(
    f"Manifest:\n"
    f"{MANIFEST_OUT}"
)
