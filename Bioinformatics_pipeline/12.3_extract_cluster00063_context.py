#!/usr/bin/env python3

from pathlib import Path
import pandas as pd


###############################################################################
## Configuration
###############################################################################

WF = (
    Path.home()
    / "methanotrophs"
    / "methanotroph_project"
    / "jeppe"
    / "genome_analysis_workflow"
)

TARGET = "Cluster_00063"

STAGE12 = (
    WF
    / "12_neighborhood_extraction"
)

OUTDIR = (
    STAGE12
    / "Cluster_00063_MCA0421"
)

OUTDIR.mkdir(
    parents=True,
    exist_ok=True
)


###############################################################################
## Input tables
###############################################################################

FILES = {
    "occurrences": (
        STAGE12
        / "focal_module_occurrences.tsv"
    ),

    "genes": (
        STAGE12
        / "observed_neighborhood_genes.tsv"
    ),

    "observability": (
        STAGE12
        / "focal_context_observability.tsv"
    ),

    "coordinates": (
        WF
        / "10_genomic_coordinates"
        / "module_protein_coordinates.tsv"
    ),
}


###############################################################################
## Output tables
###############################################################################

OUTPUTS = {
    "occurrences": (
        OUTDIR
        / "occurrences_direct_cluster00063.tsv"
    ),

    "genes": (
        OUTDIR
        / "genes_direct_cluster00063.tsv"
    ),

    "observability": (
        OUTDIR
        / "observability_direct_cluster00063.tsv"
    ),

    "coordinates": (
        OUTDIR
        / "coordinates_direct_cluster00063.tsv"
    ),
}


###############################################################################
## Extract exact Cluster_00063 rows
###############################################################################

results = {}

for name, infile in FILES.items():

    if not infile.is_file():

        raise FileNotFoundError(
            f"Required input not found: {infile}"
        )

    print()
    print("=" * 100)
    print(name.upper())
    print(infile)
    print("=" * 100)

    df = pd.read_csv(
        infile,
        sep="\t",
        dtype=str,
        low_memory=False
    )

    print(
        f"Rows: {len(df)}"
    )

    print(
        "Columns:"
    )

    print(
        df.columns.tolist()
    )


    ###########################################################################
    ## Find columns containing an exact Cluster_00063 value
    ###########################################################################

    matching_columns = []

    for column in df.columns:

        if (
            df[column]
            .fillna("")
            .eq(TARGET)
            .any()
        ):

            matching_columns.append(
                column
            )

    print(
        f"Columns containing exact "
        f"{TARGET}: "
        f"{matching_columns}"
    )

    if not matching_columns:

        raise RuntimeError(
            f"No column in {infile} contains "
            f"an exact {TARGET} value."
        )


    ###########################################################################
    ## Retain rows with Cluster_00063 in ANY matching column
    ###########################################################################

    mask = pd.Series(
        False,
        index=df.index
    )

    for column in matching_columns:

        mask = (
            mask
            |
            df[column]
            .fillna("")
            .eq(TARGET)
        )

    direct = (
        df.loc[mask]
        .copy()
    )

    results[name] = direct


    ###########################################################################
    ## Write table
    ###########################################################################

    outfile = OUTPUTS[name]

    direct.to_csv(
        outfile,
        sep="\t",
        index=False
    )


    ###########################################################################
    ## Report
    ###########################################################################

    print()
    print("=" * 100)
    print(
        f"{name}: DIRECT "
        f"{TARGET} ROWS"
    )
    print("=" * 100)

    print(
        f"Matching columns: "
        f"{matching_columns}"
    )

    print(
        f"Rows: {len(direct)}"
    )

    print()

    print(
        direct.head(10)
        .to_string(index=False)
    )


###############################################################################
## Cross-check the focal Cluster_00063 tables
###############################################################################

for name in [
    "occurrences",
    "observability",
    "coordinates",
]:

    direct = results[name]

    if len(direct) != 46:

        print(
            f"WARNING: {name} contains "
            f"{len(direct)} direct rows; "
            "the original analysis contained 46."
        )


###############################################################################
## Final output report
###############################################################################

print()
print("=" * 100)
print("OUTPUT DIRECTORY")
print("=" * 100)

print(
    OUTDIR
)

print()

for name, path in OUTPUTS.items():

    print(
        f"{name:15s}  "
        f"{len(results[name]):7d} rows  "
        f"{path.name}"
    )
