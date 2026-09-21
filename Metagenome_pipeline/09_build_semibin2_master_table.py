#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import os


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

CHECKM_MASTER = (
    PROJECT
    / "binner_comparison"
    / "CheckM2_binner_master.tsv"
)

GTDB_DIR = (
    PROJECT
    / "binner_comparison"
    / "gtdbtk"
)

BAC = (
    GTDB_DIR
    / "gtdbtk.bac120.summary.tsv"
)

ARC = (
    GTDB_DIR
    / "gtdbtk.ar53.summary.tsv"
)

OUTDIR = (
    PROJECT
    / "semibin2_analysis"
)

MASTER_OUT = (
    OUTDIR
    / "MAG_master_table_SemiBin2.tsv"
)

METH_OUT = (
    OUTDIR
    / "Methylococcales_MAGs_SemiBin2.tsv"
)

METH_BINDIR = (
    OUTDIR
    / "Methylococcales_bins"
)

OUTDIR.mkdir(
    parents=True,
    exist_ok=True
)

METH_BINDIR.mkdir(
    parents=True,
    exist_ok=True
)


## ------------------------------------------------------------
## Read CheckM2 master table
## ------------------------------------------------------------

checkm = pd.read_csv(
    CHECKM_MASTER,
    sep="\t"
)


required = {
    "Genome_ID",
    "Binner",
    "Assembly",
    "Sample",
    "Original_path",
    "Completeness",
    "Contamination"
}

missing = required - set(
    checkm.columns
)

if missing:

    raise RuntimeError(
        "Missing columns from CheckM2 master table: "
        + ", ".join(
            sorted(missing)
        )
    )


## ------------------------------------------------------------
## Keep Flye + SemiBin2 only
## ------------------------------------------------------------

df = checkm[
    (
        checkm["Assembly"]
        .astype(str)
        .str.lower()
        == "flye"
    )
    &
    (
        checkm["Binner"]
        .astype(str)
        .str.lower()
        == "semibin2"
    )
].copy()


print(
    f"Flye + SemiBin2 CheckM2 bins: {len(df)}"
)


if len(df) != 187:

    raise RuntimeError(
        f"Expected 187 Flye + SemiBin2 bins, "
        f"found {len(df)}."
    )


## ------------------------------------------------------------
## Read GTDB-Tk summaries
## ------------------------------------------------------------

gtdb_tables = []


if BAC.exists():

    bac = pd.read_csv(
        BAC,
        sep="\t"
    )

    bac["GTDB_domain_source"] = "Bacteria"

    gtdb_tables.append(
        bac
    )


if ARC.exists():

    arc = pd.read_csv(
        ARC,
        sep="\t"
    )

    arc["GTDB_domain_source"] = "Archaea"

    gtdb_tables.append(
        arc
    )


if not gtdb_tables:

    raise RuntimeError(
        "No GTDB-Tk summary files found."
    )


gtdb = pd.concat(
    gtdb_tables,
    ignore_index=True
)


if "user_genome" not in gtdb.columns:

    raise RuntimeError(
        "GTDB-Tk table does not contain user_genome."
    )


gtdb = gtdb.rename(
    columns={
        "user_genome":
            "Genome_ID"
    }
)


## ------------------------------------------------------------
## Keep GTDB columns useful for downstream work
## ------------------------------------------------------------

wanted_gtdb = [
    "Genome_ID",
    "classification",
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
    "GTDB_domain_source"
]

wanted_gtdb = [
    col
    for col in wanted_gtdb
    if col in gtdb.columns
]

gtdb = gtdb[
    wanted_gtdb
].copy()


## ------------------------------------------------------------
## Merge CheckM2 + GTDB-Tk
## ------------------------------------------------------------

df = df.merge(
    gtdb,
    on="Genome_ID",
    how="left",
    validate="one_to_one"
)


missing_gtdb = (
    df["classification"]
    .isna()
    .sum()
)

print(
    f"Bins without GTDB classification row: "
    f"{missing_gtdb}"
)


## ------------------------------------------------------------
## Parse GTDB taxonomy
## ------------------------------------------------------------

taxonomy_ranks = {
    "d__": "Domain",
    "p__": "Phylum",
    "c__": "Class",
    "o__": "Order",
    "f__": "Family",
    "g__": "Genus",
    "s__": "Species"
}


for column in taxonomy_ranks.values():

    df[column] = pd.NA


for idx, classification in df["classification"].items():

    if pd.isna(classification):
        continue

    parts = str(
        classification
    ).split(";")

    for part in parts:

        for prefix, column in taxonomy_ranks.items():

            if part.startswith(prefix):

                value = part[
                    len(prefix):
                ]

                if value != "":
                    df.at[
                        idx,
                        column
                    ] = value

                break


## ------------------------------------------------------------
## Quality categories
##
## Same definitions as previous analyses.
## ------------------------------------------------------------

def quality_class(row):

    comp = row["Completeness"]
    cont = row["Contamination"]

    if comp >= 90 and cont <= 5:

        return "High quality"

    if comp >= 50 and cont <= 10:

        return "Medium quality"

    if cont > 10:

        return "High contamination"

    return "Low completeness"


df["Quality"] = df.apply(
    quality_class,
    axis=1
)


## ------------------------------------------------------------
## Convenient sample labels
## ------------------------------------------------------------

df["Sample_short"] = (
    df["Sample"]
    .astype(str)
    .str.replace(
        "_seqs$",
        "",
        regex=True
    )
)


## ------------------------------------------------------------
## Make sure original FASTA paths still exist
## ------------------------------------------------------------

df["FASTA_exists"] = (
    df["Original_path"]
    .apply(
        lambda x:
            Path(str(x)).exists()
    )
)


missing_fasta = (
    (~df["FASTA_exists"])
    .sum()
)

print(
    f"Missing FASTA files: {missing_fasta}"
)


## ------------------------------------------------------------
## Reorder useful columns first
## ------------------------------------------------------------

front = [
    "Genome_ID",
    "Sample",
    "Sample_short",
    "Binner",
    "Assembly",
    "Quality",
    "Completeness",
    "Contamination",
    "Domain",
    "Phylum",
    "Class",
    "Order",
    "Family",
    "Genus",
    "Species",
    "classification",
    "classification_method",
    "Original_path",
    "FASTA_exists"
]

front = [
    col
    for col in front
    if col in df.columns
]

rest = [
    col
    for col in df.columns
    if col not in front
]

df = df[
    front
    + rest
]


## ------------------------------------------------------------
## Sort
## ------------------------------------------------------------

df = df.sort_values(
    [
        "Sample",
        "Order",
        "Genus",
        "Completeness"
    ],
    ascending=[
        True,
        True,
        True,
        False
    ],
    na_position="last"
)


## ------------------------------------------------------------
## Write complete SemiBin2 master table
## ------------------------------------------------------------

df.to_csv(
    MASTER_OUT,
    sep="\t",
    index=False
)


## ------------------------------------------------------------
## Methylococcales subset
## ------------------------------------------------------------

meth = df[
    df["Order"]
    .fillna("")
    .eq(
        "Methylococcales"
    )
].copy()


meth = meth.sort_values(
    [
        "Sample",
        "Genus",
        "Completeness"
    ],
    ascending=[
        True,
        True,
        False
    ]
)


meth.to_csv(
    METH_OUT,
    sep="\t",
    index=False
)


## ------------------------------------------------------------
## Create links to Methylococcales FASTAs
##
## Useful later for FeGenie and other downstream analyses.
## ------------------------------------------------------------

for _, row in meth.iterrows():

    source = Path(
        row["Original_path"]
    )

    destination = (
        METH_BINDIR
        / f"{row['Genome_ID']}.fa"
    )

    if destination.exists() or destination.is_symlink():

        destination.unlink()

    if source.exists():

        destination.symlink_to(
            source.resolve()
        )


## ------------------------------------------------------------
## Quality summary
## ------------------------------------------------------------

print()
print(
    "SemiBin2 quality summary"
)

print(
    "=" * 60
)

quality_counts = (
    df["Quality"]
    .value_counts()
)

for quality in [
    "High quality",
    "Medium quality",
    "Low completeness",
    "High contamination"
]:

    print(
        f"{quality:20s} "
        f"{quality_counts.get(quality, 0)}"
    )


## ------------------------------------------------------------
## Methylococcales summary
## ------------------------------------------------------------

print()
print(
    f"Methylococcales MAGs: "
    f"{len(meth)}"
)

print()


if len(meth) > 0:

    display_cols = [
        "Sample_short",
        "Genome_ID",
        "Genus",
        "Species",
        "Completeness",
        "Contamination",
        "Quality"
    ]

    print(
        meth[
            display_cols
        ].to_string(
            index=False
        )
    )


print()
print(
    f"Master table:\n{MASTER_OUT}"
)

print()
print(
    f"Methylococcales table:\n{METH_OUT}"
)

print()
print(
    f"Methylococcales FASTAs:\n{METH_BINDIR}"
)
