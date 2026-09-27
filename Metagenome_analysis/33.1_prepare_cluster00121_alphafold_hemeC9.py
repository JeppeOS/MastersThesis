#!/usr/bin/env python3

from pathlib import Path
import hashlib
import json
import re
import sys

import pandas as pd


## ================================================================== ##
## STAGE 67
##
## PREPARE ALPHAFOLD SERVER JOBS FOR CLUSTER_00121
##
## Selected proteins:
##
##     4 metagenome Cluster_00121 focal proteins
##     1 GlobDB Cluster_00121 MMseqs representative
##
## Prediction:
##
##     holo only
##     9 x HEC per protein
##
## No apo models are generated.
## ================================================================== ##


ROOT = (
    Path.home()
    / "methanotrophs"
    / "methanotroph_project"
    / "jeppe"
)


MG = ROOT / "metagenome"

OLD = ROOT / "genome_analysis_workflow"


CROSSWALK = (
    MG
    / "comparative_analysis"
    / "Cluster_00121"
    / "local_clustering"
    / "C121_focal_family_crosswalk.tsv"
)


COMBINED_FASTA = (
    MG
    / "comparative_analysis"
    / "Cluster_00121"
    / "Cluster_00121_combined_neighborhood_proteins.faa"
)


OLD_MEMBERSHIP = (
    OLD
    / "06_clustering"
    / "cluster_membership.tsv"
)


OUTDIR = (
    MG
    / "comparative_analysis"
    / "Cluster_00121"
    / "structure_analysis"
    / "alphafold"
)


OUTDIR.mkdir(
    parents=True,
    exist_ok=True,
)


OUT_SELECTED = (
    OUTDIR
    / "Cluster_00121_AF_selected_proteins.tsv"
)


OUT_FASTA = (
    OUTDIR
    / "Cluster_00121_AF_selected_proteins.faa"
)


OUT_JSON = (
    OUTDIR
    / "Cluster_00121_hemeC9_alphafoldserver.json"
)


TARGET_CLUSTER = "Cluster_00121"

EXPECTED_MAG = 4

EXPECTED_TOTAL = 5

EXPECTED_HEMES = 9


## ================================================================== ##
## Helpers
## ================================================================== ##

def fail(message):

    print(
        f"ERROR: {message}",
        file=sys.stderr,
    )

    sys.exit(1)


def read_fasta(path):

    seqs = {}

    header = None
    seq = []


    with path.open() as handle:

        for raw in handle:

            line = raw.strip()


            if not line:
                continue


            if line.startswith(">"):

                if header is not None:

                    seqs[
                        header.split()[0]
                    ] = "".join(seq).rstrip("*")


                header = line[1:]

                seq = []


            else:

                seq.append(line)


    if header is not None:

        seqs[
            header.split()[0]
        ] = "".join(seq).rstrip("*")


    return seqs


def safe_name(text):

    return re.sub(
        r"[^A-Za-z0-9_.-]+",
        "_",
        str(text),
    )


def deterministic_seed(text):

    digest = hashlib.sha256(
        text.encode("utf-8")
    ).hexdigest()

    ## Keep comfortably inside signed 32-bit integer range.
    return str(
        int(
            digest[:8],
            16,
        )
        %
        2_000_000_000
        +
        1
    )


## ================================================================== ##
## Validate inputs
## ================================================================== ##

for path in [
    CROSSWALK,
    COMBINED_FASTA,
    OLD_MEMBERSHIP,
]:

    if not path.is_file():

        fail(
            f"Missing input:\n{path}"
        )


## ================================================================== ##
## 1. Read focal crosswalk
## ================================================================== ##

cross = pd.read_csv(
    CROSSWALK,
    sep="\t",
    dtype=str,
    keep_default_na=False,
)


required = {
    "source_dataset",
    "genome",
    "focal_protein_id",
    "combined_protein_id",
    "local_family_id",
    "number_of_hemes",
}


missing = required - set(
    cross.columns
)


if missing:

    fail(
        "Crosswalk missing columns: "
        +
        ", ".join(
            sorted(missing)
        )
    )


if set(
    cross["local_family_id"]
) != {
    "C121_local_001"
}:

    fail(
        "Not all focal proteins belong to C121_local_001."
    )


## ================================================================== ##
## 2. Select four MAG proteins
## ================================================================== ##

mag = (
    cross[
        cross[
            "source_dataset"
        ]
        ==
        "Metagenome"
    ]
    .copy()
)


if len(mag) != EXPECTED_MAG:

    fail(
        f"Expected {EXPECTED_MAG} MAG focal proteins; "
        f"found {len(mag)}."
    )


## ================================================================== ##
## 3. Identify actual GlobDB MMseqs representative
## ================================================================== ##

old = pd.read_csv(
    OLD_MEMBERSHIP,
    sep="\t",
    dtype=str,
    keep_default_na=False,
)


required_old = {
    "cluster",
    "representative_protein",
    "protein_id",
}


missing = required_old - set(
    old.columns
)


if missing:

    fail(
        "Old membership missing columns: "
        +
        ", ".join(
            sorted(missing)
        )
    )


cluster = old[
    old[
        "cluster"
    ]
    ==
    TARGET_CLUSTER
].copy()


representatives = sorted(
    set(
        x.strip()

        for x
        in cluster[
            "representative_protein"
        ]

        if x.strip()
    )
)


if len(representatives) != 1:

    fail(
        f"Expected one GlobDB representative for "
        f"{TARGET_CLUSTER}; found {representatives}"
    )


globdb_rep = representatives[0]


globdb = cross[
    (
        cross[
            "source_dataset"
        ]
        ==
        "GlobDB"
    )
    &
    (
        cross[
            "focal_protein_id"
        ]
        ==
        globdb_rep
    )
].copy()


if len(globdb) != 1:

    fail(
        f"Could not map GlobDB representative "
        f"{globdb_rep} uniquely into focal crosswalk."
    )


## ================================================================== ##
## 4. Combine selected proteins
## ================================================================== ##

selected = pd.concat(
    [
        globdb,
        mag,
    ],
    ignore_index=True,
)


if len(selected) != EXPECTED_TOTAL:

    fail(
        f"Expected {EXPECTED_TOTAL} selected proteins; "
        f"found {len(selected)}."
    )


## ================================================================== ##
## 5. Recover sequences and count CXXCH motifs
## ================================================================== ##

seqs = read_fasta(
    COMBINED_FASTA
)


records = []

jobs = []


for _, row in selected.iterrows():

    combined_id = row[
        "combined_protein_id"
    ]


    if combined_id not in seqs:

        fail(
            f"{combined_id} not found in combined FASTA."
        )


    sequence = seqs[
        combined_id
    ]


    motifs = [
        match.start() + 1

        for match
        in re.finditer(
            r"C..CH",
            sequence,
        )
    ]


    n_motifs = len(
        motifs
    )


    if n_motifs != EXPECTED_HEMES:

        fail(
            f"{row['focal_protein_id']} has "
            f"{n_motifs} canonical CXXCH motifs, "
            f"expected {EXPECTED_HEMES}."
        )


    if row[
        "source_dataset"
    ] == "Metagenome":

        source_label = (
            "MAG_"
            +
            safe_name(
                row[
                    "genome"
                ]
            )
        )


    else:

        source_label = (
            "GlobDB_rep_"
            +
            safe_name(
                row[
                    "genome"
                ]
            )
        )


    job_name = (
        "C121_"
        +
        source_label
        +
        "_hemeC9"
    )


    seed = deterministic_seed(
        row[
            "focal_protein_id"
        ]
    )


    jobs.append(
        {
            "name":
                job_name,

            "modelSeeds":
                [
                    seed
                ],

            "sequences":
                [
                    {
                        "proteinChain":
                            {
                                "sequence":
                                    sequence,

                                "count":
                                    1,

                                "useStructureTemplate":
                                    False,
                            }
                    },

                    {
                        "ligand":
                            {
                                "ligand":
                                    "CCD_HEC",

                                "count":
                                    EXPECTED_HEMES,
                            }
                    },
                ],

            "dialect":
                "alphafoldserver",

            "version":
                3,
        }
    )


    records.append(
        {
            "job_name":
                job_name,

            "source_dataset":
                row[
                    "source_dataset"
                ],

            "genome":
                row[
                    "genome"
                ],

            "protein_id":
                row[
                    "focal_protein_id"
                ],

            "combined_protein_id":
                combined_id,

            "sequence_length":
                len(
                    sequence
                ),

            "canonical_CXXCH_count":
                n_motifs,

            "CXXCH_start_positions":
                ";".join(
                    map(
                        str,
                        motifs,
                    )
                ),

            "HEC_count":
                EXPECTED_HEMES,

            "model_seed":
                seed,
        }
    )


selected_out = pd.DataFrame(
    records
)


## ================================================================== ##
## 6. Write metadata
## ================================================================== ##

selected_out.to_csv(
    OUT_SELECTED,
    sep="\t",
    index=False,
)


## ================================================================== ##
## 7. Write selected FASTA
## ================================================================== ##

with OUT_FASTA.open(
    "w"
) as handle:

    for record in records:

        sequence = seqs[
            record[
                "combined_protein_id"
            ]
        ]


        handle.write(
            f">{record['job_name']} "
            f"protein_id={record['protein_id']} "
            f"source={record['source_dataset']}\n"
        )


        for start in range(
            0,
            len(sequence),
            80,
        ):

            handle.write(
                sequence[
                    start:
                    start + 80
                ]
                +
                "\n"
            )


## ================================================================== ##
## 8. Write one AlphaFold Server JSON containing all five jobs
## ================================================================== ##

with OUT_JSON.open(
    "w"
) as handle:

    json.dump(
        jobs,
        handle,
        indent=2,
    )

    handle.write("\n")


## ================================================================== ##
## Report
## ================================================================== ##

print("=" * 80)

print(
    "STAGE 67 - CLUSTER_00121 HEME-AWARE ALPHAFOLD JOBS"
)

print("=" * 80)

print()

print(
    selected_out[
        [
            "source_dataset",
            "genome",
            "protein_id",
            "sequence_length",
            "canonical_CXXCH_count",
            "HEC_count",
        ]
    ].to_string(
        index=False
    )
)

print()

print(
    f"Selected proteins: {len(selected_out)}"
)

print(
    f"Heme molecules/job: {EXPECTED_HEMES}"
)

print(
    "Apo jobs: 0"
)

print()

print(
    f"AlphaFold Server JSON:\n  {OUT_JSON}"
)

print()

print(
    f"Metadata:\n  {OUT_SELECTED}"
)

print()

print(
    f"FASTA:\n  {OUT_FASTA}"
)
