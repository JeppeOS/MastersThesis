#!/usr/bin/env python3

from pathlib import Path
import pandas as pd


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

MANIFEST = (
    PROJECT
    / "binner_comparison"
    / "bin_manifest.tsv"
)

CHECKM = (
    PROJECT
    / "binner_comparison"
    / "checkm2"
    / "quality_report.tsv"
)

OUTDIR = (
    PROJECT
    / "binner_comparison"
)

SUMMARY_OUT = (
    OUTDIR
    / "CheckM2_binner_summary.tsv"
)

SAMPLE_OUT = (
    OUTDIR
    / "CheckM2_binner_summary_by_sample.tsv"
)

MASTER_OUT = (
    OUTDIR
    / "CheckM2_binner_master.tsv"
)


## ------------------------------------------------------------
## Read files
## ------------------------------------------------------------

manifest = pd.read_csv(
    MANIFEST,
    sep="\t"
)

checkm = pd.read_csv(
    CHECKM,
    sep="\t"
)


## ------------------------------------------------------------
## Standardize CheckM2 genome ID
## ------------------------------------------------------------

if "Name" not in checkm.columns:
    raise RuntimeError(
        f"Could not find 'Name' in CheckM2 columns:\n"
        f"{list(checkm.columns)}"
    )


checkm["Genome_ID"] = (
    checkm["Name"]
    .astype(str)
    .str.replace(r"\.(fa|fasta|fna)$", "", regex=True)
)


## ------------------------------------------------------------
## Merge
## ------------------------------------------------------------

df = manifest.merge(
    checkm,
    on="Genome_ID",
    how="left",
    validate="one_to_one"
)


missing = df["Completeness"].isna().sum()

if missing:
    raise RuntimeError(
        f"{missing} manifest bins had no CheckM2 result."
    )


## ------------------------------------------------------------
## Quality categories
##
## Same definitions used previously:
##
## High:
##   completeness >= 90
##   contamination <= 5
##
## Medium:
##   completeness >= 50
##   contamination <= 10
##   excluding High
##
## High contamination:
##   contamination > 10
##
## Low completeness:
##   everything else
## ------------------------------------------------------------

def classify(row):

    comp = row["Completeness"]
    cont = row["Contamination"]

    if comp >= 90 and cont <= 5:
        return "High quality"

    if comp >= 50 and cont <= 10:
        return "Medium quality"

    if cont > 10:
        return "High contamination"

    return "Low completeness"


df["Quality_class"] = df.apply(
    classify,
    axis=1
)


## ------------------------------------------------------------
## Save full merged table
## ------------------------------------------------------------

df.to_csv(
    MASTER_OUT,
    sep="\t",
    index=False
)


## ------------------------------------------------------------
## Group summary
## ------------------------------------------------------------

quality_order = [
    "High quality",
    "Medium quality",
    "Low completeness",
    "High contamination",
]


summary = (
    df
    .groupby(
        [
            "Assembly",
            "Binner",
            "Quality_class"
        ]
    )
    .size()
    .unstack(
        fill_value=0
    )
    .reindex(
        columns=quality_order,
        fill_value=0
    )
    .reset_index()
)


summary["Total"] = (
    summary[quality_order]
    .sum(axis=1)
)


summary["High_plus_medium"] = (
    summary["High quality"]
    + summary["Medium quality"]
)


summary["High_plus_medium_pct"] = (
    100
    * summary["High_plus_medium"]
    / summary["Total"]
)


summary["High_quality_pct"] = (
    100
    * summary["High quality"]
    / summary["Total"]
)


summary.to_csv(
    SUMMARY_OUT,
    sep="\t",
    index=False
)


## ------------------------------------------------------------
## Per-sample summary
## ------------------------------------------------------------

sample_summary = (
    df
    .groupby(
        [
            "Assembly",
            "Binner",
            "Sample",
            "Quality_class"
        ]
    )
    .size()
    .unstack(
        fill_value=0
    )
    .reindex(
        columns=quality_order,
        fill_value=0
    )
    .reset_index()
)


sample_summary["Total"] = (
    sample_summary[quality_order]
    .sum(axis=1)
)


sample_summary["High_plus_medium"] = (
    sample_summary["High quality"]
    + sample_summary["Medium quality"]
)


sample_summary.to_csv(
    SAMPLE_OUT,
    sep="\t",
    index=False
)


## ------------------------------------------------------------
## Print concise summary
## ------------------------------------------------------------

print()
print("CheckM2 summary")
print("=" * 85)

for _, row in summary.iterrows():

    print(
        f"{row['Assembly']:12s} "
        f"{row['Binner']:10s} | "
        f"HQ {int(row['High quality']):3d} | "
        f"MQ {int(row['Medium quality']):3d} | "
        f"Low {int(row['Low completeness']):3d} | "
        f">10% contam {int(row['High contamination']):3d} | "
        f"HQ+MQ {int(row['High_plus_medium']):3d} "
        f"({row['High_plus_medium_pct']:.1f}%)"
    )

print()
print(f"Written: {SUMMARY_OUT}")
print(f"Written: {SAMPLE_OUT}")
print(f"Written: {MASTER_OUT}")
