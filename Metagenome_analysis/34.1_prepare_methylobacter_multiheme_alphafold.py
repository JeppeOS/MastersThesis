#!/usr/bin/env python3

from pathlib import Path
import hashlib
import json
import re
import sys

import pandas as pd


## ================================================================== ##
## STAGE 70
##
## PREPARE HEME-AWARE ALPHAFOLD SERVER INPUTS FOR THE THREE
## METHYLOBACTER MULTIHEME CANDIDATES
##
## Targets:
##
## barcode14 SemiBin_1
##     contig_75_3433
##     5-heme
##     GlobDB Cluster_00042
##
## barcode15 SemiBin_3
##     contig_171_1
##     7-heme
##     no GlobDB family match
##
## barcode16 SemiBin_1
##     contig_483_99
##     8-heme
##     top Cluster_00288
##     secondary Cluster_00218
##
## Heme-aware only. No apo controls.
##
## The script recovers the proteins from the Stage 53 neighborhood
## protein collection, verifies the expected canonical CXXCH motif
## count, and generates AlphaFold Server v3 JSON.
## ================================================================== ##


PROJECT = (
    Path.home()
    / "methanotrophs"
    / "methanotroph_project"
    / "jeppe"
    / "metagenome"
)


NEIGHBOR_ROOT = (
    PROJECT
    / "functional_analysis"
    / "methylococcales_neighborhoods"
)


METADATA_FILE = (
    NEIGHBOR_ROOT
    / "priority_neighborhood_protein_metadata.tsv"
)


FASTA_FILE = (
    NEIGHBOR_ROOT
    / "priority_neighborhood_proteins.faa"
)


OUTDIR = (
    PROJECT
    / "comparative_analysis"
    / "Methylobacter_multiheme"
    / "structure_analysis"
    / "alphafold"
)


OUTDIR.mkdir(
    parents=True,
    exist_ok=True,
)


OUT_TABLE = (
    OUTDIR
    / "Methylobacter_multiheme_AF_selected_proteins.tsv"
)


OUT_FASTA = (
    OUTDIR
    / "Methylobacter_multiheme_AF_selected_proteins.faa"
)


OUT_JSON = (
    OUTDIR
    / "Methylobacter_multiheme_heme_alphafoldserver.json"
)


## ================================================================== ##
## Candidate definitions
## ================================================================== ##

TARGETS = [

    {
        "candidate":
            "MB_C00042_b14SB1",

        "protein_id":
            "contig_75_3433",

        "taxon":
            "Methylobacter_A",

        "mag":
            "barcode14 SemiBin_1",

        "expected_hemes":
            5,

        "globdb_mapping":
            "Cluster_00042",

        "caution":
            "",
    },

    {
        "candidate":
            "MB_nomatch_b15SB3",

        "protein_id":
            "contig_171_1",

        "taxon":
            "Methylobacter_C sp002256465",

        "mag":
            "barcode15 SemiBin_3",

        "expected_hemes":
            7,

        "globdb_mapping":
            "no_match",

        "caution":
            "",
    },

    {
        "candidate":
            "MB_C00288_b16SB1",

        "protein_id":
            "contig_483_99",

        "taxon":
            "Methylobacter_C sp002256465",

        "mag":
            "barcode16 SemiBin_1",

        "expected_hemes":
            8,

        "globdb_mapping":
            "Cluster_00288; secondary Cluster_00218",

        "caution":
            "MAG CheckM2 contamination 17.61%; interpret locus cautiously",
    },
]


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


    with open(
        path,
        "r",
    ) as handle:

        for line in handle:

            line = line.strip()


            if not line:
                continue


            if line.startswith(
                ">"
            ):

                if header is not None:

                    records[
                        header
                    ] = "".join(
                        sequence
                    )


                header = (
                    line[
                        1:
                    ]
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
        name.encode(
            "utf-8"
        )
    ).hexdigest()


    value = int(
        digest[
            :8
        ],
        16,
    )


    return str(
        (
            value
            %
            2147483646
        )
        +
        1
    )


def find_exact_row(
    table,
    protein_id,
):

    mask = pd.Series(
        False,
        index=table.index,
    )


    for column in table.columns:

        mask = (
            mask
            |
            (
                table[
                    column
                ]
                .astype(str)
                ==
                protein_id
            )
        )


    hits = table[
        mask
    ].copy()


    if len(
        hits
    ) != 1:

        fail(
            f"{protein_id}: expected exactly one row "
            f"in {METADATA_FILE}, found {len(hits)}."
        )


    return hits.iloc[
        0
    ]


def find_nbp_id(row):

    matches = []


    for value in row.values:

        value = str(
            value
        )


        if re.fullmatch(
            r"NBP\d+",
            value,
        ):

            matches.append(
                value
            )


    matches = sorted(
        set(
            matches
        )
    )


    if len(
        matches
    ) != 1:

        fail(
            "Could not uniquely determine NBP identifier "
            f"from metadata row. Found: {matches}"
        )


    return matches[
        0
    ]


## ================================================================== ##
## Validate inputs
## ================================================================== ##

for path in [
    METADATA_FILE,
    FASTA_FILE,
]:

    if not path.is_file():

        fail(
            f"Missing input file:\n{path}"
        )


metadata = pd.read_csv(
    METADATA_FILE,
    sep="\t",
    dtype=str,
).fillna(
    ""
)


sequences = read_fasta(
    FASTA_FILE
)


print("=" * 80)

print(
    "STAGE 70 - METHYLOBACTER MULTIHEME ALPHAFOLD PREPARATION"
)

print("=" * 80)

print()

print(
    f"Neighborhood metadata rows: {len(metadata):,}"
)

print(
    f"Neighborhood FASTA entries: {len(sequences):,}"
)

print()


## ================================================================== ##
## Recover candidates and validate CXXCH motifs
## ================================================================== ##

output_rows = []

fasta_records = []

jobs = []


for target in TARGETS:

    row = find_exact_row(
        metadata,
        target[
            "protein_id"
        ],
    )


    nbp_id = find_nbp_id(
        row
    )


    if nbp_id not in sequences:

        fail(
            f"{target['protein_id']}: "
            f"{nbp_id} not present in FASTA."
        )


    sequence = (
        sequences[
            nbp_id
        ]
        .replace(
            "*",
            "",
        )
        .upper()
    )


    motifs = list(
        re.finditer(
            r"C..CH",
            sequence,
        )
    )


    observed_hemes = len(
        motifs
    )


    if (
        observed_hemes
        !=
        target[
            "expected_hemes"
        ]
    ):

        fail(
            f"{target['protein_id']}: expected "
            f"{target['expected_hemes']} canonical CXXCH motifs "
            f"but found {observed_hemes}.\n"
            "Inspect the sequence rather than automatically changing "
            "the heme count."
        )


    motif_positions = ";".join(
        str(
            match.start()
            +
            1
        )

        for match
        in motifs
    )


    job_name = (
        target[
            "candidate"
        ]
        +
        f"_hemeC{observed_hemes}"
    )


    seed = deterministic_seed(
        job_name
    )


    output_rows.append(
        {
            "candidate":
                target[
                    "candidate"
                ],

            "job_name":
                job_name,

            "mag":
                target[
                    "mag"
                ],

            "taxon":
                target[
                    "taxon"
                ],

            "protein_id":
                target[
                    "protein_id"
                ],

            "nbp_id":
                nbp_id,

            "sequence_length":
                len(
                    sequence
                ),

            "canonical_CXXCH":
                observed_hemes,

            "CXXCH_start_positions":
                motif_positions,

            "globdb_mapping":
                target[
                    "globdb_mapping"
                ],

            "alphafold_seed":
                seed,

            "caution":
                target[
                    "caution"
                ],
        }
    )


    fasta_records.append(
        (
            job_name,
            sequence,
        )
    )


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
                            sequence,

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
                            observed_hemes,
                    }
                },
            ],

            "dialect":
                "alphafoldserver",

            "version":
                3,
        }
    )


    print(
        f"{target['protein_id']}: "
        f"{len(sequence)} aa, "
        f"{observed_hemes} CXXCH motifs, "
        f"{nbp_id}"
    )


## ================================================================== ##
## Write outputs
## ================================================================== ##

out = pd.DataFrame(
    output_rows
)


out.to_csv(
    OUT_TABLE,
    sep="\t",
    index=False,
)


with open(
    OUT_FASTA,
    "w",
) as handle:

    for name, sequence in fasta_records:

        handle.write(
            f">{name}\n"
        )


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
## Final QC
## ================================================================== ##

if len(
    jobs
) != 3:

    fail(
        "Expected exactly three AlphaFold jobs."
    )


if sum(
    row[
        "canonical_CXXCH"
    ]

    for row
    in output_rows
) != (
    5
    +
    7
    +
    8
):

    fail(
        "Unexpected total heme count."
    )


print()

print("=" * 80)

print(
    "STAGE 70 COMPLETE"
)

print("=" * 80)

print()

print(
    out[
        [
            "candidate",
            "taxon",
            "protein_id",
            "sequence_length",
            "canonical_CXXCH",
            "globdb_mapping",
        ]
    ]
    .to_string(
        index=False
    )
)

print()

print(
    "AlphaFold jobs: 3"
)

print(
    "HEC counts:     5, 7, 8"
)

print(
    "Apo jobs:       0"
)

print()

print(
    f"Metadata:\n  {OUT_TABLE}"
)

print()

print(
    f"FASTA:\n  {OUT_FASTA}"
)

print()

print(
    f"AlphaFold Server JSON:\n  {OUT_JSON}"
)
