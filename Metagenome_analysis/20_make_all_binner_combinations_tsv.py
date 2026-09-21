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

## Original Flye + MetaBAT2 master table
ORIGINAL_MASTER = (
    PROJECT
    / "MAG_master_table.tsv"
)

## Strainberry + MetaBAT2 CheckM2
STRAINBERRY_METABAT = (
    PROJECT
    / "strainberry_checkm2_results"
    / "checkm2"
    / "quality_report.tsv"
)

## Flye/Strainberry + SemiBin2/MaxBin2 summary
OTHER_BINNERS = (
    PROJECT
    / "binner_comparison"
    / "CheckM2_binner_summary.tsv"
)

OUTPUT = (
    PROJECT
    / "binner_comparison"
    / "CheckM2_all_combinations.tsv"
)


## ------------------------------------------------------------
## Quality classification
##
## Same definitions used throughout the analysis.
## ------------------------------------------------------------

def quality_class(comp, cont):

    if comp >= 90 and cont <= 5:
        return "High quality"

    if comp >= 50 and cont <= 10:
        return "Medium quality"

    if cont > 10:
        return "High contamination"

    return "Low completeness"


## ------------------------------------------------------------
## Summarize one dataframe
## ------------------------------------------------------------

def summarize(df, assembly, binner):

    classes = df.apply(
        lambda row: quality_class(
            row["Completeness"],
            row["Contamination"]
        ),
        axis=1
    )

    counts = classes.value_counts()

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

    return {
        "Assembly": assembly,
        "Binner": binner,

        "High_quality": high,
        "Medium_quality": medium,
        "High_plus_medium": high_medium,
        "Low_completeness": low,
        "High_contamination": high_contam,
        "Total_bins": total,

        "High_quality_pct":
            100 * high / total,

        "Medium_quality_pct":
            100 * medium / total,

        "High_plus_medium_pct":
            100 * high_medium / total,

        "Low_completeness_pct":
            100 * low / total,

        "High_contamination_pct":
            100 * high_contam / total,
    }


rows = []


## ------------------------------------------------------------
## 1. Flye + MetaBAT2
##
## Restrict original MAG master table to barcode10-16.
## barcode09 is deliberately excluded.
## ------------------------------------------------------------

original = pd.read_csv(
    ORIGINAL_MASTER,
    sep="\t"
)


required = {
    "Sample",
    "Completeness",
    "Contamination"
}

missing = required - set(
    original.columns
)

if missing:

    raise RuntimeError(
        "Missing columns in MAG_master_table.tsv: "
        + ", ".join(
            sorted(missing)
        )
    )


original = original[
    original["Sample"]
    .astype(str)
    .str.match(
        r"^barcode(10|11|12|13|14|15|16)(_seqs)?$"
    )
].copy()


rows.append(
    summarize(
        original,
        assembly="Flye",
        binner="MetaBAT2"
    )
)


## ------------------------------------------------------------
## 2. Strainberry + MetaBAT2
## ------------------------------------------------------------

sb_metabat = pd.read_csv(
    STRAINBERRY_METABAT,
    sep="\t"
)


rows.append(
    summarize(
        sb_metabat,
        assembly="Strainberry",
        binner="MetaBAT2"
    )
)


## ------------------------------------------------------------
## 3-6. SemiBin2 + MaxBin2
##
## These were already summarized by script 20.
## ------------------------------------------------------------

others = pd.read_csv(
    OTHER_BINNERS,
    sep="\t"
)


## Standardize capitalization
others["Assembly"] = (
    others["Assembly"]
    .replace({
        "flye": "Flye",
        "strainberry": "Strainberry"
    })
)


## Convert each existing summary row into the same format
for _, row in others.iterrows():

    total = int(
        row["Total"]
    )

    high = int(
        row["High quality"]
    )

    medium = int(
        row["Medium quality"]
    )

    low = int(
        row["Low completeness"]
    )

    high_contam = int(
        row["High contamination"]
    )

    high_medium = (
        high
        + medium
    )

    rows.append({
        "Assembly":
            row["Assembly"],

        "Binner":
            row["Binner"],

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
            100 * high / total,

        "Medium_quality_pct":
            100 * medium / total,

        "High_plus_medium_pct":
            100 * high_medium / total,

        "Low_completeness_pct":
            100 * low / total,

        "High_contamination_pct":
            100 * high_contam / total,
    })


## ------------------------------------------------------------
## Final table
## ------------------------------------------------------------

result = pd.DataFrame(
    rows
)


## Desired ordering
assembly_order = {
    "Flye": 0,
    "Strainberry": 1
}

binner_order = {
    "MetaBAT2": 0,
    "SemiBin2": 1,
    "MaxBin2": 2
}


result["_assembly_order"] = (
    result["Assembly"]
    .map(assembly_order)
)

result["_binner_order"] = (
    result["Binner"]
    .map(binner_order)
)


result = (
    result
    .sort_values(
        [
            "_assembly_order",
            "_binner_order"
        ]
    )
    .drop(
        columns=[
            "_assembly_order",
            "_binner_order"
        ]
    )
    .reset_index(
        drop=True
    )
)


## ------------------------------------------------------------
## Round percentages for readability
## ------------------------------------------------------------

pct_columns = [
    "High_quality_pct",
    "Medium_quality_pct",
    "High_plus_medium_pct",
    "Low_completeness_pct",
    "High_contamination_pct",
]

result[pct_columns] = (
    result[pct_columns]
    .round(2)
)


## ------------------------------------------------------------
## Sanity checks
## ------------------------------------------------------------

if len(result) != 6:

    raise RuntimeError(
        f"Expected 6 assembly/binner combinations, "
        f"but found {len(result)}."
    )


for _, row in result.iterrows():

    calculated_total = (
        row["High_quality"]
        + row["Medium_quality"]
        + row["Low_completeness"]
        + row["High_contamination"]
    )

    if calculated_total != row["Total_bins"]:

        raise RuntimeError(
            f"Quality classes do not sum correctly for "
            f"{row['Assembly']} + {row['Binner']}."
        )


## ------------------------------------------------------------
## Save
## ------------------------------------------------------------

result.to_csv(
    OUTPUT,
    sep="\t",
    index=False
)


print()
print("All assembly × binner combinations")
print("=" * 100)

print(
    result[
        [
            "Assembly",
            "Binner",
            "High_quality",
            "Medium_quality",
            "High_plus_medium",
            "Low_completeness",
            "High_contamination",
            "Total_bins",
            "High_plus_medium_pct"
        ]
    ].to_string(
        index=False
    )
)

print()
print(
    f"Written: {OUTPUT}"
)
