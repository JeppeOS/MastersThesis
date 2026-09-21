#!/usr/bin/env python3


from pathlib import Path
import pandas as pd
import gzip
import hashlib
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

OUTDIR = (
    PROJECT
    / "semibin2_analysis"
    / "final_Methylococcales"
)

FASTA_DIR = (
    OUTDIR
    / "fasta"
)

FINAL_TABLE = (
    OUTDIR
    / "SemiBin2_Methylococcales_FINAL.tsv"
)

CHECKSUM_TABLE = (
    OUTDIR
    / "SemiBin2_Methylococcales_FINAL.sha256.tsv"
)


## ------------------------------------------------------------
## Recreate frozen output directory
## ------------------------------------------------------------

if OUTDIR.exists():

    shutil.rmtree(
        OUTDIR
    )


FASTA_DIR.mkdir(
    parents=True,
    exist_ok=True
)


## ------------------------------------------------------------
## Read master table
## ------------------------------------------------------------

master = pd.read_csv(
    MASTER,
    sep="\t"
)


coverage = pd.read_csv(
    COVERAGE,
    sep="\t"
)


## ------------------------------------------------------------
## Methylococcales only
## ------------------------------------------------------------

meth = master[
    master["Order"]
    .fillna("")
    .eq(
        "Methylococcales"
    )
].copy()


if len(meth) != 16:

    raise RuntimeError(
        f"Expected 16 Methylococcales MAGs, "
        f"found {len(meth)}."
    )


## ------------------------------------------------------------
## Coverage columns to merge
## ------------------------------------------------------------

coverage_cols = [
    "Genome_ID",
    "Genome_length_bp",
    "N_contigs",
    "Mean_coverage_x"
]


missing_cov = set(
    coverage_cols
) - set(
    coverage.columns
)


if missing_cov:

    raise RuntimeError(
        "Missing coverage columns: "
        + ", ".join(
            sorted(
                missing_cov
            )
        )
    )


coverage_small = coverage[
    coverage_cols
].copy()


## ------------------------------------------------------------
## One-to-one merge
## ------------------------------------------------------------

meth = meth.merge(
    coverage_small,
    on="Genome_ID",
    how="left",
    validate="one_to_one"
)


if meth[
    "Mean_coverage_x"
].isna().any():

    missing = meth.loc[
        meth["Mean_coverage_x"].isna(),
        "Genome_ID"
    ].tolist()

    raise RuntimeError(
        "Missing coverage for:\n"
        + "\n".join(
            missing
        )
    )


## ------------------------------------------------------------
## Detect gzip compression by magic bytes
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
## Materialize one frozen plain-text FASTA
## ------------------------------------------------------------

def freeze_fasta(source, destination):

    source = Path(
        source
    )

    if not source.exists():

        raise RuntimeError(
            f"Missing source FASTA:\n"
            f"{source}"
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
## SHA256
## ------------------------------------------------------------

def sha256_file(path):

    digest = hashlib.sha256()

    with open(
        path,
        "rb"
    ) as handle:

        while True:

            block = handle.read(
                1024 * 1024
            )

            if not block:
                break

            digest.update(
                block
            )

    return digest.hexdigest()


## ------------------------------------------------------------
## Freeze FASTAs
## ------------------------------------------------------------

frozen_paths = []
checksums = []


for _, row in meth.iterrows():

    genome_id = str(
        row["Genome_ID"]
    )

    source = Path(
        str(
            row["Original_path"]
        )
    )


    destination = (
        FASTA_DIR
        / f"{genome_id}.fa"
    )


    freeze_fasta(
        source,
        destination
    )


    checksum = sha256_file(
        destination
    )


    frozen_paths.append(
        str(
            destination
        )
    )


    checksums.append({

        "Genome_ID":
            genome_id,

        "Frozen_FASTA":
            str(
                destination
            ),

        "SHA256":
            checksum
    })


meth[
    "Frozen_FASTA"
] = frozen_paths


## ------------------------------------------------------------
## Authoritative column order
## ------------------------------------------------------------

preferred_columns = [

    ## Identity
    "Genome_ID",
    "Sample",
    "Sample_short",

    ## MAG quality
    "Quality",
    "Completeness",
    "Contamination",

    ## Genome properties
    "Genome_length_bp",
    "N_contigs",
    "Mean_coverage_x",

    ## GTDB taxonomy
    "Domain",
    "Phylum",
    "Class",
    "Order",
    "Family",
    "Genus",
    "Species",
    "classification",

    ## GTDB classification details
    "classification_method",
    "fastani_reference",
    "fastani_reference_radius",
    "fastani_taxonomy",
    "fastani_ani",
    "fastani_af",
    "closest_placement_reference",
    "closest_placement_taxonomy",
    "closest_placement_ani",
    "closest_placement_af",
    "msa_percent",
    "red_value",
    "warnings",

    ## Provenance
    "Binner",
    "Assembly",
    "Original_path",
    "Frozen_FASTA"
]


preferred_columns = [
    col
    for col in preferred_columns
    if col in meth.columns
]


extra_columns = [
    col
    for col in meth.columns
    if col not in preferred_columns
]


meth = meth[
    preferred_columns
    + extra_columns
]


## ------------------------------------------------------------
## Sort authoritative table
## ------------------------------------------------------------

meth = meth.sort_values(
    [
        "Sample",
        "Genus",
        "Mean_coverage_x"
    ],
    ascending=[
        True,
        True,
        False
    ]
).reset_index(
    drop=True
)


## ------------------------------------------------------------
## Save authoritative table
## ------------------------------------------------------------

meth.to_csv(
    FINAL_TABLE,
    sep="\t",
    index=False,
    na_rep=""
)


## ------------------------------------------------------------
## Save checksums
## ------------------------------------------------------------

checksum_df = pd.DataFrame(
    checksums
)

checksum_df = checksum_df[
    checksum_df["Genome_ID"]
    .isin(
        meth["Genome_ID"]
    )
].copy()


checksum_df = (
    meth[
        [
            "Genome_ID"
        ]
    ]
    .merge(
        checksum_df,
        on="Genome_ID",
        how="left",
        validate="one_to_one"
    )
)


checksum_df.to_csv(
    CHECKSUM_TABLE,
    sep="\t",
    index=False
)


## ------------------------------------------------------------
## Validation
## ------------------------------------------------------------

n_fasta = len(
    list(
        FASTA_DIR.glob(
            "*.fa"
        )
    )
)


if n_fasta != 16:

    raise RuntimeError(
        f"Expected 16 frozen FASTAs, "
        f"found {n_fasta}."
    )


if meth["Genome_ID"].duplicated().any():

    raise RuntimeError(
        "Duplicate Genome_ID values "
        "in final table."
    )


## ------------------------------------------------------------
## Print concise final table
## ------------------------------------------------------------

print()

print(
    "FINAL SemiBin2 Methylococcales dataset"
)

print(
    "=" * 110
)


display_columns = [
    "Sample_short",
    "Genome_ID",
    "Genus",
    "Species",
    "Completeness",
    "Contamination",
    "Genome_length_bp",
    "N_contigs",
    "Mean_coverage_x"
]


print(

    meth[
        display_columns
    ].to_string(
        index=False
    )
)


print()

print(
    f"MAGs: {len(meth)}"
)

print(
    f"Frozen FASTAs: {n_fasta}"
)

print()

print(
    f"Authoritative table:\n"
    f"{FINAL_TABLE}"
)

print()

print(
    f"Checksums:\n"
    f"{CHECKSUM_TABLE}"
)

print()

print(
    f"Frozen FASTAs:\n"
    f"{FASTA_DIR}"
)
