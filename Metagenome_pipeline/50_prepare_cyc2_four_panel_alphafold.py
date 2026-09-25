#!/usr/bin/env python3

from pathlib import Path
import hashlib
import json
import re
import sys

import pandas as pd


## ================================================================== ##
## STAGE 73
##
## PREPARE FOUR HEME-AWARE CYC2 STRUCTURES
##
## Panel set:
##
##   1. MAG Cluster_00024
##      Methylobacter_A
##      barcode14 SemiBin_0
##      contig_38_1172
##
##   2. MAG Cluster_00050
##      barcode10 SemiBin_1
##      contig_25_317
##
##   3. GlobDB Cluster_00050 representative
##
##   4. GlobDB Cluster_00206 representative
##
## Each Cyc2 candidate is expected to contain one canonical CXXCH
## motif and is submitted with one HEC ligand.
##
## Cluster_00206 is explicitly inspected here rather than assumed.
## ================================================================== ##


HOME = Path.home()

META_PROJECT = (
    HOME
    / "methanotrophs"
    / "methanotroph_project"
    / "jeppe"
    / "metagenome"
)

OLD_PROJECT = (
    HOME
    / "methanotrophs"
    / "methanotroph_project"
    / "jeppe"
    / "genome_analysis_workflow"
)


## ================================================================== ##
## Inputs
## ================================================================== ##

MAG_METADATA = (
    META_PROJECT
    / "functional_analysis"
    / "methylococcales_neighborhoods"
    / "priority_neighborhood_protein_metadata.tsv"
)

MAG_FASTA = (
    META_PROJECT
    / "functional_analysis"
    / "methylococcales_neighborhoods"
    / "priority_neighborhood_proteins.faa"
)

OLD_MEMBERSHIP = (
    OLD_PROJECT
    / "06_clustering"
    / "cluster_membership.tsv"
)

OLD_FASTA = (
    OLD_PROJECT
    / "06_clustering"
    / "exported_heme_candidates.faa"
)


## ================================================================== ##
## Outputs
## ================================================================== ##

OUTDIR = (
    META_PROJECT
    / "comparative_analysis"
    / "Cyc2_four_panel"
    / "structure_analysis"
    / "alphafold"
)

OUTDIR.mkdir(
    parents=True,
    exist_ok=True,
)

OUT_TABLE = (
    OUTDIR
    / "Cyc2_four_panel_selected_proteins.tsv"
)

OUT_FASTA = (
    OUTDIR
    / "Cyc2_four_panel_selected_proteins.faa"
)

OUT_JSON = (
    OUTDIR
    / "Cyc2_four_panel_hemeC1_alphafoldserver.json"
)

OUT_CLUSTER_CHECK = (
    OUTDIR
    / "Cyc2_GlobDB_cluster_representative_check.tsv"
)


## ================================================================== ##
## Helpers
## ================================================================== ##

def fail(message):

    print(
        f"\nERROR: {message}",
        file=sys.stderr,
    )

    sys.exit(1)


def read_fasta(path):

    records = {}

    header = None
    sequence = []


    with open(path) as handle:

        for line in handle:

            line = line.strip()

            if not line:
                continue


            if line.startswith(">"):

                if header is not None:

                    records[
                        header
                    ] = "".join(
                        sequence
                    )


                header = (
                    line[1:]
                    .split()[0]
                )

                sequence = []


            else:

                sequence.append(
                    line
                )


    if header is not None:

        records[
            header
        ] = "".join(
            sequence
        )


    return records


def deterministic_seed(name):

    digest = hashlib.sha256(
        name.encode()
    ).hexdigest()

    return str(
        (
            int(
                digest[:8],
                16,
            )
            %
            2147483646
        )
        +
        1
    )


def resolve_column(
    columns,
    candidates,
    required=False,
):

    lower = {
        str(c).lower():
            c
        for c
        in columns
    }


    for candidate in candidates:

        if candidate in columns:
            return candidate

        if candidate.lower() in lower:
            return lower[
                candidate.lower()
            ]


    if required:

        fail(
            "Could not resolve required column from: "
            +
            ", ".join(
                candidates
            )
        )


    return None


def find_nbp_id(row):

    values = [
        str(value)

        for value
        in row.values
    ]


    hits = sorted(
        {
            value

            for value
            in values

            if re.fullmatch(
                r"NBP\d+",
                value,
            )
        }
    )


    if len(hits) != 1:

        fail(
            f"Could not uniquely identify NBP ID: {hits}"
        )


    return hits[0]


def get_mag_sequence(
    metadata,
    sequences,
    protein_id,
):

    mask = pd.Series(
        False,
        index=metadata.index,
    )


    for column in metadata.columns:

        mask |= (
            metadata[column]
            .astype(str)
            ==
            protein_id
        )


    hits = metadata[
        mask
    ]


    if len(hits) != 1:

        fail(
            f"{protein_id}: expected exactly one "
            f"Stage-53 metadata row; found {len(hits)}."
        )


    nbp_id = find_nbp_id(
        hits.iloc[0]
    )


    if nbp_id not in sequences:

        fail(
            f"{protein_id}: {nbp_id} absent from MAG FASTA."
        )


    return (
        nbp_id,
        sequences[
            nbp_id
        ].replace(
            "*",
            "",
        ).upper(),
    )


def get_cluster_representative(
    membership,
    cluster_id,
):

    cluster_col = resolve_column(
        membership.columns,
        [
            "cluster",
            "mmseq_cluster",
        ],
        required=True,
    )


    rep_col = resolve_column(
        membership.columns,
        [
            "representative_protein",
            "representative",
        ],
        required=True,
    )


    protein_col = resolve_column(
        membership.columns,
        [
            "protein_id",
            "protein",
        ],
        required=True,
    )


    subset = membership[
        membership[
            cluster_col
        ]
        ==
        cluster_id
    ].copy()


    if len(subset) == 0:

        fail(
            f"{cluster_id}: no membership rows."
        )


    representatives = sorted(
        set(
            subset[
                rep_col
            ]
        )
        -
        {""}
    )


    if len(
        representatives
    ) != 1:

        fail(
            f"{cluster_id}: expected one representative ID; "
            f"found {representatives}."
        )


    representative = representatives[
        0
    ]


    rep_rows = subset[
        subset[
            protein_col
        ]
        ==
        representative
    ].copy()


    if len(
        rep_rows
    ) != 1:

        fail(
            f"{cluster_id}: representative {representative} "
            f"is not represented by exactly one membership row."
        )


    return (
        representative,
        rep_rows.iloc[0],
        subset,
    )


def motif_starts(sequence):

    return [
        m.start() + 1

        for m
        in re.finditer(
            r"C..CH",
            sequence,
        )
    ]


## ================================================================== ##
## Load inputs
## ================================================================== ##

for path in [
    MAG_METADATA,
    MAG_FASTA,
    OLD_MEMBERSHIP,
    OLD_FASTA,
]:

    if not path.is_file():

        fail(
            f"Missing input:\n{path}"
        )


mag_metadata = pd.read_csv(
    MAG_METADATA,
    sep="\t",
    dtype=str,
).fillna("")


mag_sequences = read_fasta(
    MAG_FASTA
)


old_membership = pd.read_csv(
    OLD_MEMBERSHIP,
    sep="\t",
    dtype=str,
).fillna("")


old_sequences = read_fasta(
    OLD_FASTA
)


print("=" * 80)
print("STAGE 73 - CYC2 FOUR-PANEL ALPHAFOLD PREPARATION")
print("=" * 80)


## ================================================================== ##
## 1. MAG proteins
## ================================================================== ##

selected = []


MAG_TARGETS = [

    {
        "panel":
            "MAG_Cluster_00024_Methylobacter",

        "cluster":
            "Cluster_00024",

        "protein_id":
            "contig_38_1172",

        "source":
            "MAG",

        "taxon":
            "Methylobacter_A",

        "genome":
            "barcode14 SemiBin_0",
    },

    {
        "panel":
            "MAG_Cluster_00050",

        "cluster":
            "Cluster_00050",

        "protein_id":
            "contig_25_317",

        "source":
            "MAG",

        "taxon":
            "Methylumidiphilus",

        "genome":
            "barcode10 SemiBin_1",
    },
]


for target in MAG_TARGETS:

    nbp_id, sequence = get_mag_sequence(
        mag_metadata,
        mag_sequences,
        target[
            "protein_id"
        ],
    )


    selected.append(
        {
            **target,

            "sequence_id":
                nbp_id,

            "sequence":
                sequence,

            "n_CXXCH":
                len(
                    motif_starts(
                        sequence
                    )
                ),

            "CXXCH_start_positions":
                ";".join(
                    map(
                        str,
                        motif_starts(
                            sequence
                        ),
                    )
                ),
        }
    )


## ================================================================== ##
## 2. GlobDB representatives
## ================================================================== ##

cluster_check_rows = []


for cluster_id in [
    "Cluster_00050",
    "Cluster_00206",
]:

    representative, rep_row, cluster_subset = (
        get_cluster_representative(
            old_membership,
            cluster_id,
        )
    )


    if representative not in old_sequences:

        fail(
            f"{cluster_id}: representative "
            f"{representative} absent from:\n{OLD_FASTA}"
        )


    sequence = (
        old_sequences[
            representative
        ]
        .replace(
            "*",
            "",
        )
        .upper()
    )


    heme_col = resolve_column(
        cluster_subset.columns,
        [
            "number_of_hemes",
            "heme_count",
            "n_hemes",
        ],
    )


    fe_col = resolve_column(
        cluster_subset.columns,
        [
            "FeGenie_HMMs",
            "fegenie_HMMs",
            "fegenie_hmms",
        ],
    )


    signalp_col = resolve_column(
        cluster_subset.columns,
        [
            "signalp_prediction",
            "SignalP_prediction",
        ],
    )


    localization_col = resolve_column(
        cluster_subset.columns,
        [
            "localization_class",
            "upstream_localization_class",
        ],
    )


    cluster_check_rows.append(
        {
            "cluster":
                cluster_id,

            "n_members":
                len(
                    cluster_subset
                ),

            "representative_protein":
                representative,

            "representative_length":
                len(
                    sequence
                ),

            "representative_CXXCH":
                len(
                    motif_starts(
                        sequence
                    )
                ),

            "representative_heme_annotation":
                (
                    rep_row[
                        heme_col
                    ]
                    if heme_col
                    else
                    ""
                ),

            "representative_FeGenie_HMMs":
                (
                    rep_row[
                        fe_col
                    ]
                    if fe_col
                    else
                    ""
                ),

            "representative_SignalP":
                (
                    rep_row[
                        signalp_col
                    ]
                    if signalp_col
                    else
                    ""
                ),

            "representative_localization":
                (
                    rep_row[
                        localization_col
                    ]
                    if localization_col
                    else
                    ""
                ),
        }
    )


    selected.append(
        {
            "panel":
                (
                    "GlobDB_"
                    +
                    cluster_id
                ),

            "cluster":
                cluster_id,

            "protein_id":
                representative,

            "source":
                "GlobDB",

            "taxon":
                "",

            "genome":
                (
                    str(
                        rep_row.get(
                            "genome",
                            "",
                        )
                    )
                ),

            "sequence_id":
                representative,

            "sequence":
                sequence,

            "n_CXXCH":
                len(
                    motif_starts(
                        sequence
                    )
                ),

            "CXXCH_start_positions":
                ";".join(
                    map(
                        str,
                        motif_starts(
                            sequence
                        ),
                    )
                ),
        }
    )


cluster_check = pd.DataFrame(
    cluster_check_rows
)


cluster_check.to_csv(
    OUT_CLUSTER_CHECK,
    sep="\t",
    index=False,
)


## ================================================================== ##
## 3. Strict heme-motif QC
## ================================================================== ##

if len(
    selected
) != 4:

    fail(
        f"Expected four selected proteins; "
        f"found {len(selected)}."
    )


for row in selected:

    if row[
        "n_CXXCH"
    ] != 1:

        fail(
            f"{row['panel']} / {row['protein_id']}: "
            f"found {row['n_CXXCH']} canonical CXXCH motifs. "
            "Expected exactly one for this Cyc2 comparison.\n"
            "Inspect this protein before proceeding."
        )


## ================================================================== ##
## 4. Build AlphaFold jobs
## ================================================================== ##

jobs = []


for row in selected:

    job_name = (
        "Cyc2_"
        +
        row[
            "panel"
        ]
        +
        "_hemeC1"
    )


    seed = deterministic_seed(
        job_name
    )


    row[
        "job_name"
    ] = job_name

    row[
        "alphafold_seed"
    ] = seed


    jobs.append(
        {
            "name":
                job_name,

            "modelSeeds": [
                seed
            ],

            "sequences": [

                {
                    "proteinChain": {

                        "sequence":
                            row[
                                "sequence"
                            ],

                        "count":
                            1,

                        "useStructureTemplate":
                            False,
                    }
                },

                {
                    "ligand": {

                        "ligand":
                            "CCD_HEC",

                        "count":
                            1,
                    }
                },
            ],

            "dialect":
                "alphafoldserver",

            "version":
                3,
        }
    )


## ================================================================== ##
## 5. Write outputs
## ================================================================== ##

table = pd.DataFrame(
    selected
)


table[
    [
        "panel",
        "source",
        "cluster",
        "genome",
        "taxon",
        "protein_id",
        "sequence_id",
        "n_CXXCH",
        "CXXCH_start_positions",
        "job_name",
        "alphafold_seed",
    ]
].to_csv(
    OUT_TABLE,
    sep="\t",
    index=False,
)


with open(
    OUT_FASTA,
    "w",
) as handle:

    for row in selected:

        handle.write(
            f">{row['job_name']}\n"
        )


        sequence = row[
            "sequence"
        ]


        for i in range(
            0,
            len(
                sequence
            ),
            80,
        ):

            handle.write(
                sequence[
                    i:
                    i + 80
                ]
                +
                "\n"
            )


with open(
    OUT_JSON,
    "w",
) as handle:

    json.dump(
        jobs,
        handle,
        indent=2,
    )


## ================================================================== ##
## Report
## ================================================================== ##

print()
print("GlobDB representative check:")
print()

print(
    cluster_check
    .to_string(
        index=False
    )
)


print()
print("Selected four-panel set:")
print()

print(
    table[
        [
            "panel",
            "source",
            "cluster",
            "genome",
            "taxon",
            "protein_id",
            "sequence",
            "n_CXXCH",
        ]
    ]
    .assign(
        sequence=lambda x:
            x[
                "sequence"
            ].str.len()
    )
    .rename(
        columns={
            "sequence":
                "length"
        }
    )
    .to_string(
        index=False
    )
)


print()
print("=" * 80)
print("STAGE 73 COMPLETE")
print("=" * 80)

print()
print("AlphaFold jobs: 4")
print("HEC per job:    1")
print("Apo jobs:       0")

print()
print(
    f"Representative check:\n  {OUT_CLUSTER_CHECK}"
)

print()
print(
    f"Selected proteins:\n  {OUT_TABLE}"
)

print()
print(
    f"FASTA:\n  {OUT_FASTA}"
)

print()
print(
    f"AlphaFold Server JSON:\n  {OUT_JSON}"
)
