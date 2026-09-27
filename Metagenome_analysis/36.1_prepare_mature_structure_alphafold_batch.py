#!/usr/bin/env python3

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import pandas as pd


## ================================================================== ##
## STAGE 75
##
## PREPARE SIGNAL-PEPTIDE-TRIMMED ALPHAFOLD BATCH
##
## Final structure set:
##
##   Cyc2:
##       4 proteins
##
##   Cluster_00121:
##       5 proteins
##
##   Methylobacter multiheme:
##       3 proteins
##
## Total:
##       12 proteins
##
##
## Strategy:
##
##   1. Reuse the exact protein sequences and HEC counts from the
##      already validated AlphaFold Server JSON files.
##
##   2. Retrieve SignalP prediction + cleavage site from the
##      authoritative integrated annotation tables.
##
##   3. If SignalP predicts a cleavable signal peptide:
##
##          SP / LIPO
##
##      remove residues up to the cleavage site.
##
##      Example:
##
##          25-26
##
##      means residues 1-25 are removed and the mature protein begins
##      at residue 26.
##
##   4. Proteins with OTHER / no SignalP cleavage remain unchanged.
##
##   5. Verify that the mature sequence still contains exactly the
##      expected number of canonical CXXCH motifs.
##
##   6. Write:
##
##      - one combined AlphaFold Server JSON containing all 12 jobs
##      - mature FASTA
##      - complete trimming/QC table
##
##
## No apo jobs.
## HEC counts are retained exactly as in the validated original jobs.
## ================================================================== ##


HOME = Path.home()


PROJECT = (
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
## Original validated AlphaFold job files
## ================================================================== ##

C121_JSON = (
    PROJECT
    / "comparative_analysis"
    / "Cluster_00121"
    / "structure_analysis"
    / "alphafold"
    / "Cluster_00121_hemeC9_alphafoldserver.json"
)


MB_JSON = (
    PROJECT
    / "comparative_analysis"
    / "Methylobacter_multiheme"
    / "structure_analysis"
    / "alphafold"
    / "Methylobacter_multiheme_heme_alphafoldserver.json"
)


CYC2_JSON = (
    PROJECT
    / "comparative_analysis"
    / "Cyc2_four_panel"
    / "structure_analysis"
    / "alphafold"
    / "Cyc2_four_panel_hemeC1_alphafoldserver.json"
)


## ================================================================== ##
## Authoritative SignalP annotation tables
## ================================================================== ##

MAG_MASTER = (
    PROJECT
    / "functional_analysis"
    / "integrated_annotations"
    / "MASTER_PROTEIN_ANNOTATIONS.tsv"
)


GLOBDB_MASTER = (
    OLD_PROJECT
    / "05_integrated_annotations"
    / "MASTER_PROTEIN_ANNOTATIONS.tsv"
)


## ================================================================== ##
## Output
## ================================================================== ##

OUTDIR = (
    PROJECT
    / "comparative_analysis"
    / "final_structure_figures"
    / "00_mature_alphafold"
)


OUTDIR.mkdir(
    parents=True,
    exist_ok=True,
)


OUT_JSON = (
    OUTDIR
    / "final_12_mature_holo_alphafoldserver.json"
)


OUT_FASTA = (
    OUTDIR
    / "final_12_mature_holo_proteins.faa"
)


OUT_TABLE = (
    OUTDIR
    / "final_12_mature_holo_proteins.tsv"
)


## ================================================================== ##
## Final protein set
##
## json_group tells us which validated original JSON contains the
## sequence.
##
## job_pattern identifies the original AlphaFold job.
##
## protein_id is used only for annotation lookup.
## ================================================================== ##

TARGETS = [

    ## -------------------------------------------------------------- ##
    ## Cyc2
    ## -------------------------------------------------------------- ##

    {
        "panel_set": "cyc2",
        "protein_id": "contig_38_1172",
        "genome_hint": "barcode14",
        "annotation_source": "MAG",
        "json_group": "cyc2",
        "job_pattern": r"cyc2.*mag.*cluster_00024",
        "expected_hemes": 1,
    },

    {
        "panel_set": "cyc2",
        "protein_id": "contig_25_317",
        "genome_hint": "barcode10",
        "annotation_source": "MAG",
        "json_group": "cyc2",
        "job_pattern": r"cyc2.*mag.*cluster_00050",
        "expected_hemes": 1,
    },

    {
        "panel_set": "cyc2",
        "protein_id": "GCA_023229605_000000000176_34",
        "genome_hint": "GCA_023229605",
        "annotation_source": "GlobDB",
        "json_group": "cyc2",
        "job_pattern": r"cyc2.*globdb.*cluster_00050",
        "expected_hemes": 1,
    },

    {
        "panel_set": "cyc2",
        "protein_id": "GCA_002843155_000000000119_52",
        "genome_hint": "GCA_002843155",
        "annotation_source": "GlobDB",
        "json_group": "cyc2",
        "job_pattern": r"cyc2.*globdb.*cluster_00206",
        "expected_hemes": 1,
    },


    ## -------------------------------------------------------------- ##
    ## Cluster_00121
    ## -------------------------------------------------------------- ##

    {
        "panel_set": "cluster00121",
        "protein_id": "BCRBG_22465_000000000369_1",
        "genome_hint": "BCRBG_22465",
        "annotation_source": "GlobDB",
        "json_group": "cluster00121",
        "job_pattern": r"c121.*globdb|globdb.*c121|globdb_rep",
        "expected_hemes": 9,
    },

    {
        "panel_set": "cluster00121",
        "protein_id": "contig_25_1899",
        "genome_hint": "barcode10",
        "annotation_source": "MAG",
        "json_group": "cluster00121",
        "job_pattern": r"c121.*barcode10.*semibin_?1(?!\d)",
        "expected_hemes": 9,
    },

    {
        "panel_set": "cluster00121",
        "protein_id": "contig_1431_201",
        "genome_hint": "barcode10",
        "annotation_source": "MAG",
        "json_group": "cluster00121",
        "job_pattern": r"c121.*barcode10.*semibin_?8(?!\d)",
        "expected_hemes": 9,
    },

    {
        "panel_set": "cluster00121",
        "protein_id": "contig_122_3517",
        "genome_hint": "barcode12",
        "annotation_source": "MAG",
        "json_group": "cluster00121",
        "job_pattern": r"c121.*barcode12.*semibin_?0(?!\d)",
        "expected_hemes": 9,
    },

    {
        "panel_set": "cluster00121",
        "protein_id": "contig_871_20",
        "genome_hint": "barcode15",
        "annotation_source": "MAG",
        "json_group": "cluster00121",
        "job_pattern": r"c121.*barcode15.*semibin_?10(?!\d)",
        "expected_hemes": 9,
    },


    ## -------------------------------------------------------------- ##
    ## Methylobacter multihemes
    ## -------------------------------------------------------------- ##

    {
        "panel_set": "methylobacter",
        "protein_id": "contig_75_3433",
        "genome_hint": "barcode14",
        "annotation_source": "MAG",
        "json_group": "methylobacter",
        "job_pattern": r"mb.*c00042",
        "expected_hemes": 5,
    },

    {
        "panel_set": "methylobacter",
        "protein_id": "contig_483_99",
        "genome_hint": "barcode16",
        "annotation_source": "MAG",
        "json_group": "methylobacter",
        "job_pattern": r"mb.*c00288",
        "expected_hemes": 8,
    },

    {
        "panel_set": "methylobacter",
        "protein_id": "contig_171_1",
        "genome_hint": "barcode15",
        "annotation_source": "MAG",
        "json_group": "methylobacter",
        "job_pattern": r"mb.*nomatch",
        "expected_hemes": 7,
    },
]


## ================================================================== ##
## Utilities
## ================================================================== ##

def fail(message):

    print(
        f"\nERROR: {message}",
        file=sys.stderr,
    )

    sys.exit(1)


def require_file(path):

    if not path.is_file():

        fail(
            f"Required file missing:\n{path}"
        )


def resolve_column(
    columns,
    aliases,
    label,
):

    columns = list(
        columns
    )


    lower = {
        str(column).lower():
            column

        for column
        in columns
    }


    for alias in aliases:

        if alias in columns:

            return alias


        if alias.lower() in lower:

            return lower[
                alias.lower()
            ]


    fail(
        f"Could not resolve column '{label}'.\n"
        f"Tried: {aliases}\n"
        f"Available columns:\n"
        +
        "\n".join(
            columns
        )
    )


def load_json(path):

    require_file(
        path
    )


    data = json.loads(
        path.read_text()
    )


    if isinstance(
        data,
        dict,
    ):

        data = [
            data
        ]


    if not isinstance(
        data,
        list,
    ):

        fail(
            f"Unexpected AlphaFold JSON format:\n{path}"
        )


    return data


def protein_sequence_from_job(job):

    hits = []


    for sequence_entry in job.get(
        "sequences",
        [],
    ):

        if "proteinChain" in sequence_entry:

            protein = sequence_entry[
                "proteinChain"
            ]


            sequence = str(
                protein.get(
                    "sequence",
                    "",
                )
            ).strip().upper()


            if sequence:

                hits.append(
                    sequence
                )


    if len(
        hits
    ) != 1:

        fail(
            f"Expected exactly one protein chain in job "
            f"{job.get('name', '')}; found {len(hits)}."
        )


    return hits[
        0
    ]


def hec_count_from_job(job):

    counts = []


    for sequence_entry in job.get(
        "sequences",
        [],
    ):

        if "ligand" not in sequence_entry:
            continue


        ligand = sequence_entry[
            "ligand"
        ]


        if (
            str(
                ligand.get(
                    "ligand",
                    "",
                )
            ).upper()
            ==
            "CCD_HEC"
        ):

            counts.append(
                int(
                    ligand.get(
                        "count",
                        1,
                    )
                )
            )


    if len(
        counts
    ) != 1:

        fail(
            f"Expected exactly one HEC ligand entry in "
            f"{job.get('name', '')}; found {counts}."
        )


    return counts[
        0
    ]


def match_job(
    jobs,
    regex,
):

    pattern = re.compile(
        regex,
        flags=re.IGNORECASE,
    )


    hits = [
        job

        for job
        in jobs

        if pattern.search(
            str(
                job.get(
                    "name",
                    "",
                )
            )
        )
    ]


    if len(
        hits
    ) != 1:

        names = [
            str(
                job.get(
                    "name",
                    "",
                )
            )

            for job
            in jobs
        ]


        fail(
            f"Pattern did not identify exactly one AlphaFold job:\n"
            f"  {regex}\n\n"
            f"Matches: {len(hits)}\n\n"
            f"Available jobs:\n"
            +
            "\n".join(
                names
            )
        )


    return hits[
        0
    ]


def parse_cleavage_position(value):

    value = str(
        value
    ).strip()


    if not value:

        return None


    ## SignalP examples:
    ##
    ##   25-26.
    ##   25-26
    ##   24-25. Pr: 0.999
    ##

    match = re.search(
        r"(\d+)\s*-\s*(\d+)",
        value,
    )


    if match is None:

        fail(
            f"Could not parse SignalP cleavage site: '{value}'"
        )


    left = int(
        match.group(
            1
        )
    )


    right = int(
        match.group(
            2
        )
    )


    if right != left + 1:

        fail(
            f"Unexpected cleavage-site coordinates: '{value}'"
        )


    ## Remove residues 1..left.
    return left


def lookup_annotation(
    df,
    protein_id,
    genome_hint,
    protein_col,
    genome_col,
):

    hits = df[
        df[
            protein_col
        ].astype(str)
        ==
        protein_id
    ].copy()


    if len(
        hits
    ) > 1:

        hits2 = hits[
            hits[
                genome_col
            ]
            .astype(str)
            .str.contains(
                genome_hint,
                regex=False,
            )
        ].copy()


        if len(
            hits2
        ) == 1:

            hits = hits2


    if len(
        hits
    ) != 1:

        fail(
            f"{protein_id}: expected one annotation row; "
            f"found {len(hits)}.\n"
            f"Genome hint: {genome_hint}"
        )


    return hits.iloc[
        0
    ]


## ================================================================== ##
## Load inputs
## ================================================================== ##

for path in [
    C121_JSON,
    MB_JSON,
    CYC2_JSON,
    MAG_MASTER,
    GLOBDB_MASTER,
]:

    require_file(
        path
    )


json_groups = {

    "cluster00121":
        load_json(
            C121_JSON
        ),

    "methylobacter":
        load_json(
            MB_JSON
        ),

    "cyc2":
        load_json(
            CYC2_JSON
        ),
}


mag = pd.read_csv(
    MAG_MASTER,
    sep="\t",
    dtype=str,
    keep_default_na=False,
)


old = pd.read_csv(
    GLOBDB_MASTER,
    sep="\t",
    dtype=str,
    keep_default_na=False,
)


## ================================================================== ##
## Resolve annotation columns independently
## ================================================================== ##

def annotation_columns(df):

    return {

        "protein":
            resolve_column(
                df.columns,
                [
                    "protein_id",
                    "protein",
                ],
                "protein ID",
            ),

        "genome":
            resolve_column(
                df.columns,
                [
                    "genome",
                    "Genome_ID",
                ],
                "genome",
            ),

        "prediction":
            resolve_column(
                df.columns,
                [
                    "signalp_prediction",
                    "SignalP_prediction",
                    "signalp_class",
                ],
                "SignalP prediction",
            ),

        "cleavage":
            resolve_column(
                df.columns,
                [
                    "signalp_cs_position",
                    "SignalP_cs_position",
                    "signalp_cleavage_site",
                    "cleavage_site",
                ],
                "SignalP cleavage site",
            ),
    }


mag_cols = annotation_columns(
    mag
)


old_cols = annotation_columns(
    old
)


## ================================================================== ##
## Build mature jobs
## ================================================================== ##

new_jobs = []

table_rows = []

fasta_records = []


for target in TARGETS:

    jobs = json_groups[
        target[
            "json_group"
        ]
    ]


    original_job = match_job(
        jobs,
        target[
            "job_pattern"
        ],
    )


    original_name = str(
        original_job[
            "name"
        ]
    )


    full_sequence = protein_sequence_from_job(
        original_job
    )


    hec_count = hec_count_from_job(
        original_job
    )


    if (
        hec_count
        !=
        target[
            "expected_hemes"
        ]
    ):

        fail(
            f"{target['protein_id']}: original AlphaFold HEC count "
            f"{hec_count} != expected "
            f"{target['expected_hemes']}."
        )


    if (
        target[
            "annotation_source"
        ]
        ==
        "MAG"
    ):

        annotation_df = mag

        cols = mag_cols


    else:

        annotation_df = old

        cols = old_cols


    annotation = lookup_annotation(

        annotation_df,

        target[
            "protein_id"
        ],

        target[
            "genome_hint"
        ],

        cols[
            "protein"
        ],

        cols[
            "genome"
        ],
    )


    signalp_prediction = str(
        annotation[
            cols[
                "prediction"
            ]
        ]
    ).strip()


    signalp_cs = str(
        annotation[
            cols[
                "cleavage"
            ]
        ]
    ).strip()


    prediction_upper = signalp_prediction.upper()


    cleavable = (
        prediction_upper
        in {
            "SP",
            "LIPO",
        }
    )


    if cleavable:

        trim_n = parse_cleavage_position(
            signalp_cs
        )


        if trim_n is None:

            fail(
                f"{target['protein_id']}: "
                f"SignalP={signalp_prediction} "
                "but no cleavage position was available."
            )


    else:

        trim_n = 0


    if trim_n >= len(
        full_sequence
    ):

        fail(
            f"{target['protein_id']}: invalid trim length "
            f"{trim_n} for sequence length "
            f"{len(full_sequence)}."
        )


    mature_sequence = full_sequence[
        trim_n:
    ]


    motif_count = len(
        re.findall(
            r"C..CH",
            mature_sequence,
        )
    )


    if (
        motif_count
        !=
        target[
            "expected_hemes"
        ]
    ):

        fail(
            f"{target['protein_id']}: after trimming, "
            f"found {motif_count} canonical CXXCH motifs; "
            f"expected {target['expected_hemes']}.\n"
            f"SignalP: {signalp_prediction}\n"
            f"Cleavage: {signalp_cs}\n"
            f"Removed: {trim_n} aa"
        )


    new_name = (
        original_name
        +
        "_mature"
    )


    ## -------------------------------------------------------------- ##
    ## Preserve original seeds / dialect / version.
    ##
    ## Only protein sequence and job name change.
    ## -------------------------------------------------------------- ##

    new_job = json.loads(
        json.dumps(
            original_job
        )
    )


    new_job[
        "name"
    ] = new_name


    replaced = 0


    for sequence_entry in new_job[
        "sequences"
    ]:

        if "proteinChain" in sequence_entry:

            sequence_entry[
                "proteinChain"
            ][
                "sequence"
            ] = mature_sequence

            replaced += 1


    if replaced != 1:

        fail(
            f"{original_name}: could not replace exactly one "
            "protein sequence."
        )


    new_jobs.append(
        new_job
    )


    fasta_records.append(
        (
            new_name,
            mature_sequence,
        )
    )


    table_rows.append(
        {
            "panel_set":
                target[
                    "panel_set"
                ],

            "protein_id":
                target[
                    "protein_id"
                ],

            "annotation_source":
                target[
                    "annotation_source"
                ],

            "original_job":
                original_name,

            "mature_job":
                new_name,

            "signalp_prediction":
                signalp_prediction,

            "signalp_cs_position":
                signalp_cs,

            "removed_N_terminal_aa":
                trim_n,

            "original_length":
                len(
                    full_sequence
                ),

            "mature_length":
                len(
                    mature_sequence
                ),

            "expected_hemes":
                target[
                    "expected_hemes"
                ],

            "mature_CXXCH_count":
                motif_count,

            "sequence_changed":
                int(
                    trim_n
                    >
                    0
                ),
        }
    )


## ================================================================== ##
## Strict final QC
## ================================================================== ##

if len(
    new_jobs
) != 12:

    fail(
        f"Expected 12 mature AlphaFold jobs; "
        f"created {len(new_jobs)}."
    )


names = [
    job[
        "name"
    ]

    for job
    in new_jobs
]


if len(
    names
) != len(
    set(
        names
    )
):

    fail(
        "Duplicate mature AlphaFold job names."
    )


## ================================================================== ##
## Write outputs
## ================================================================== ##

OUT_JSON.write_text(
    json.dumps(
        new_jobs,
        indent=2,
    )
    +
    "\n"
)


with OUT_FASTA.open(
    "w",
    encoding="utf-8",
) as handle:

    for name, sequence in fasta_records:

        handle.write(
            f">{name}\n"
        )


        for start in range(
            0,
            len(
                sequence
            ),
            80,
        ):

            handle.write(
                sequence[
                    start:
                    start
                    +
                    80
                ]
                +
                "\n"
            )


table = pd.DataFrame(
    table_rows
)


table.to_csv(
    OUT_TABLE,
    sep="\t",
    index=False,
)


## ================================================================== ##
## Report
## ================================================================== ##

print("=" * 100)

print(
    "STAGE 75 COMPLETE - MATURE STRUCTURE ALPHAFOLD BATCH"
)

print("=" * 100)

print()


print(
    table[
        [
            "panel_set",
            "protein_id",
            "signalp_prediction",
            "signalp_cs_position",
            "removed_N_terminal_aa",
            "original_length",
            "mature_length",
            "expected_hemes",
            "mature_CXXCH_count",
        ]
    ]
    .to_string(
        index=False
    )
)


print()

print(
    f"Total jobs:             {len(new_jobs)}"
)

print(
    "Signal peptides removed: "
    f"{int((table['removed_N_terminal_aa'] > 0).sum())}"
)

print(
    "Sequences unchanged:      "
    f"{int((table['removed_N_terminal_aa'] == 0).sum())}"
)

print(
    "Apo jobs:                 0"
)

print()

print(
    f"AlphaFold Server JSON:\n  {OUT_JSON}"
)

print()

print(
    f"Mature FASTA:\n  {OUT_FASTA}"
)

print()

print(
    f"Trimming/QC table:\n  {OUT_TABLE}"
)
