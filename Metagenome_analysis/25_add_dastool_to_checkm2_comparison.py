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

DAS_CHECKM = (
    PROJECT
    / "dastool_evaluation"
    / "checkm2"
    / "quality_report.tsv"
)

BASE_TABLE = (
    PROJECT
    / "binner_comparison"
    / "CheckM2_all_combinations.tsv"
)

OUTPUT = (
    PROJECT
    / "binner_comparison"
    / "CheckM2_all_combinations_with_DASTool.tsv"
)


## ------------------------------------------------------------
## Read DAS Tool CheckM2 results
## ------------------------------------------------------------

df = pd.read_csv(
    DAS_CHECKM,
    sep="\t"
)


if len(df) != 80:
    raise RuntimeError(
        f"Expected 80 DAS Tool bins, found {len(df)}."
    )


## ------------------------------------------------------------
## Quality classification
##
## Same definitions used throughout:
##
## HQ:
##   completeness >= 90
##   contamination <= 5
##
## MQ:
##   completeness >= 50
##   contamination <= 10
##   excluding HQ
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

counts = df["Quality_class"].value_counts()


high = int(
    counts.get(
        "High quality",
        0
    )
)

medium = int(
    counts.get(
        "Medium quality",
        0
    )
)

low = int(
    counts.get(
        "Low completeness",
        0
    )
)

high_contam = int(
    counts.get(
        "High contamination",
        0
    )
)

total = len(df)

high_medium = (
    high
    + medium
)


## ------------------------------------------------------------
## Create DAS Tool row
## ------------------------------------------------------------

das_row = pd.DataFrame([{
    "Assembly":
        "Flye",

    "Binner":
        "DAS_Tool",

    "High_quality":
        high,

    "Medium_quality":
        medium,

    "High_plus_medium":
        high_medium,

    "Low_completeness":
        low,

    "High_contamination":
        high_contam,

    "Total_bins":
        total,

    "High_quality_pct":
        round(
            100 * high / total,
            2
        ),

    "Medium_quality_pct":
        round(
            100 * medium / total,
            2
        ),

    "High_plus_medium_pct":
        round(
            100 * high_medium / total,
            2
        ),

    "Low_completeness_pct":
        round(
            100 * low / total,
            2
        ),

    "High_contamination_pct":
        round(
            100 * high_contam / total,
            2
        ),
}])


## ------------------------------------------------------------
## Read existing six-way comparison
## ------------------------------------------------------------

comparison = pd.read_csv(
    BASE_TABLE,
    sep="\t"
)


## Remove old DAS Tool row if script is rerun
comparison = comparison[
    comparison["Binner"] != "DAS_Tool"
].copy()


## ------------------------------------------------------------
## Add DAS Tool
## ------------------------------------------------------------

comparison = pd.concat(
    [
        comparison,
        das_row
    ],
    ignore_index=True
)


## ------------------------------------------------------------
## Sort
## ------------------------------------------------------------

assembly_order = {
    "Flye": 0,
    "Strainberry": 1
}

binner_order = {
    "MetaBAT2": 0,
    "SemiBin2": 1,
    "MaxBin2": 2,
    "DAS_Tool": 3
}


comparison["_assembly"] = (
    comparison["Assembly"]
    .map(assembly_order)
)

comparison["_binner"] = (
    comparison["Binner"]
    .map(binner_order)
)


comparison = (
    comparison
    .sort_values(
        [
            "_assembly",
            "_binner"
        ]
    )
    .drop(
        columns=[
            "_assembly",
            "_binner"
        ]
    )
    .reset_index(
        drop=True
    )
)


## ------------------------------------------------------------
## Save
## ------------------------------------------------------------

comparison.to_csv(
    OUTPUT,
    sep="\t",
    index=False
)


## ------------------------------------------------------------
## Print DAS Tool result
## ------------------------------------------------------------

print()
print("DAS Tool CheckM2")
print("=" * 70)

print(
    f"HQ:             {high}"
)

print(
    f"MQ:             {medium}"
)

print(
    f"HQ + MQ:        {high_medium} "
    f"({100 * high_medium / total:.1f}%)"
)

print(
    f"Low:            {low}"
)

print(
    f">10% contam:    {high_contam}"
)

print(
    f"Total:          {total}"
)

print()
print(
    f"Written: {OUTPUT}"
)
