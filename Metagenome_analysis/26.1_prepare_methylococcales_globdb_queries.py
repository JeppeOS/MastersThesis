#!/usr/bin/env python3

from __future__ import annotations

import csv
from pathlib import Path


###############################################################################
## Paths
###############################################################################

PROJECT = (
    Path.home()
    / "methanotrophs"
    / "methanotroph_project"
    / "jeppe"
    / "metagenome"
)

MASTER = (
    PROJECT
    / "functional_analysis"
    / "integrated_annotations"
    / "Methylococcales"
    / "Methylococcales_MASTER_PROTEIN_ANNOTATIONS.tsv"
)

OUT_DIR = (
    PROJECT
    / "functional_analysis"
    / "globdb_cluster_mapping"
)

FASTA_OUT = (
    OUT_DIR
    / "Methylococcales_heme_queries.faa"
)

METADATA_OUT = (
    OUT_DIR
    / "Methylococcales_heme_query_metadata.tsv"
)

SUMMARY_OUT = (
    OUT_DIR
    / "Methylococcales_heme_query_summary.tsv"
)


###############################################################################
## Helpers
###############################################################################

def as_int(value) -> int:

    value = str(value).strip()

    if value in {
        "1",
        "TRUE",
        "True",
        "true"
    }:

        return 1

    if value in {
        "0",
        "FALSE",
        "False",
        "false",
        "",
        "nan",
        "NaN"
    }:

        return 0

    return int(
        float(
            value
        )
    )


###############################################################################
## Validate input
###############################################################################

if not MASTER.is_file():

    raise SystemExit(
        f"ERROR: Methylococcales master table not found:\n{MASTER}"
    )


OUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)


###############################################################################
## Read master table
###############################################################################

with MASTER.open(
    "r",
    encoding="utf-8",
    newline=""
) as handle:

    reader = csv.DictReader(
        handle,
        delimiter="\t"
    )

    rows = list(
        reader
    )

    fields = (
        reader.fieldnames
        or []
    )


required = {
    "genome",
    "protein_id",
    "findmehemes_positive",
    "number_of_hemes",
    "export_evidence",
    "fegenie_positive",
    "fegenie_HMMs",
    "signalp_prediction",
    "deeptmhmm_class",
    "deeptmhmm_sequence",
    "localization_class"
}


missing = required.difference(
    fields
)


if missing:

    raise SystemExit(

        "ERROR: Methylococcales master table "
        f"is missing columns: {sorted(missing)}"
    )


###############################################################################
## Select all heme-positive Methylococcales proteins
###############################################################################

selected = []


for row in rows:

    if as_int(
        row[
            "findmehemes_positive"
        ]
    ) != 1:

        continue


    genome = row[
        "genome"
    ].strip()


    protein_id = row[
        "protein_id"
    ].strip()


    if not genome or not protein_id:

        raise SystemExit(
            "ERROR: selected row has empty genome or protein_id."
        )


    sequence = row[
        "deeptmhmm_sequence"
    ].strip()


    ## Prodigal translations retain a terminal stop symbol.
    if sequence.endswith(
        "*"
    ):

        sequence = sequence[
            :-1
        ]


    if "*" in sequence:

        raise SystemExit(

            f"ERROR: internal stop character found in "
            f"{genome} / {protein_id}"
        )


    if not sequence:

        raise SystemExit(

            f"ERROR: empty sequence for "
            f"{genome} / {protein_id}"
        )


    hemes = as_int(
        row[
            "number_of_hemes"
        ]
    )


    hmms = {

        x.strip()

        for x in row[
            "fegenie_HMMs"
        ].split(
            ";"
        )

        if x.strip()
    }


    is_cyc2 = int(
        "Cyc2_repCluster2"
        in hmms
    )


    is_high_heme = int(
        hemes >= 5
    )


    reasons = []


    if is_cyc2:

        reasons.append(
            "Cyc2"
        )


    if is_high_heme:

        reasons.append(
            ">=5_hemes"
        )


    copied = dict(
        row
    )


    copied[
        "clustering_sequence"
    ] = sequence


    copied[
        "clustering_sequence_length"
    ] = str(
        len(
            sequence
        )
    )


    copied[
        "priority_cyc2"
    ] = str(
        is_cyc2
    )


    copied[
        "priority_ge5_hemes"
    ] = str(
        is_high_heme
    )


    copied[
        "priority_reason"
    ] = ";".join(
        reasons
    )


    selected.append(
        copied
    )


###############################################################################
## Stable ordering and unique query IDs
###############################################################################

selected.sort(

    key=lambda row: (
        row[
            "genome"
        ],
        row[
            "protein_id"
        ]
    )
)


for index, row in enumerate(
    selected,
    start=1
):

    row[
        "query_id"
    ] = f"MGQ{index:06d}"


###############################################################################
## QC
###############################################################################

if not selected:

    raise SystemExit(
        "ERROR: no Methylococcales heme-positive proteins selected."
    )


if len(
    selected
) != 712:

    print(
        f"WARNING: expected 712 heme-positive proteins, "
        f"selected {len(selected)}."
    )


n_cyc2 = sum(

    as_int(
        row[
            "priority_cyc2"
        ]
    )

    for row in selected
)


n_high_heme = sum(

    as_int(
        row[
            "priority_ge5_hemes"
        ]
    )

    for row in selected
)


if n_cyc2 != 6:

    print(
        f"WARNING: expected 6 Cyc2 candidates, found {n_cyc2}."
    )


if n_high_heme != 15:

    print(
        f"WARNING: expected 15 >=5-heme candidates, "
        f"found {n_high_heme}."
    )


###############################################################################
## FASTA
###############################################################################

with FASTA_OUT.open(
    "w",
    encoding="utf-8"
) as handle:

    for row in selected:

        handle.write(
            f">{row['query_id']}\n"
        )


        sequence = row[
            "clustering_sequence"
        ]


        for start in range(
            0,
            len(
                sequence
            ),
            80
        ):

            handle.write(

                sequence[
                    start:start + 80
                ]
                + "\n"
            )


###############################################################################
## Metadata
###############################################################################

preferred_fields = [

    "query_id",
    "genome",
    "protein_id",
    "length",
    "clustering_sequence_length",
    "number_of_hemes",

    "priority_cyc2",
    "priority_ge5_hemes",
    "priority_reason",

    "candidate_source",

    "fegenie_positive",
    "fegenie_HMMs",
    "fegenie_categories",

    "findmehemes_positive",

    "signalp_prediction",
    "signalp_cs_position",

    "deeptmhmm_class",
    "deeptmhmm_n_tm_helices",

    "export_evidence",
    "localization_class",

    "source_proteome"
]


metadata_fields = [

    field

    for field
    in preferred_fields

    if (
        field in selected[0]
        or field == "query_id"
        or field.startswith(
            "priority_"
        )
        or field == "clustering_sequence_length"
    )
]


with METADATA_OUT.open(
    "w",
    encoding="utf-8",
    newline=""
) as handle:

    writer = csv.DictWriter(

        handle,

        fieldnames=
            metadata_fields,

        delimiter="\t",

        lineterminator="\n",

        extrasaction="ignore"
    )


    writer.writeheader()

    writer.writerows(
        selected
    )


###############################################################################
## Summary
###############################################################################

n_exported = sum(

    as_int(
        row[
            "export_evidence"
        ]
    )

    for row in selected
)


n_nonexported = (
    len(
        selected
    )
    - n_exported
)


summary_rows = [

    (
        "heme_positive_queries",
        len(
            selected
        )
    ),

    (
        "export_evidence_positive",
        n_exported
    ),

    (
        "export_evidence_negative",
        n_nonexported
    ),

    (
        "priority_Cyc2",
        n_cyc2
    ),

    (
        "priority_ge5_hemes",
        n_high_heme
    )
]


with SUMMARY_OUT.open(
    "w",
    encoding="utf-8",
    newline=""
) as handle:

    writer = csv.writer(
        handle,
        delimiter="\t",
        lineterminator="\n"
    )


    writer.writerow(
        [
            "metric",
            "value"
        ]
    )


    writer.writerows(
        summary_rows
    )


###############################################################################
## Report
###############################################################################

print(
    "Methylococcales GlobDB-query preparation complete."
)

print(
    f"Heme-positive queries:       {len(selected)}"
)

print(
    f"With export evidence:        {n_exported}"
)

print(
    f"Without export evidence:     {n_nonexported}"
)

print(
    f"Cyc2 priority proteins:      {n_cyc2}"
)

print(
    f">=5-heme priority proteins:  {n_high_heme}"
)

print()

print(
    f"FASTA:    {FASTA_OUT}"
)

print(
    f"Metadata: {METADATA_OUT}"
)

print(
    f"Summary:  {SUMMARY_OUT}"
)
