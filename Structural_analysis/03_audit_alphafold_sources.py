#!/usr/bin/env python3

from pathlib import Path
import json
import re
import pandas as pd


## ================================================================== ##
## PURPOSE
##
## Audit the AlphaFold3 structures used in the four thesis figures:
##
##   1. Cluster_00035 MtoA / FeGenie-negative
##   2. Cyc2
##   3. MCA0421 / Cluster_00063
##   4. Cluster_00121
##
## This script:
##
##   - keeps original workflow provenance explicit
##   - does NOT modify any original files
##   - identifies the selected model CIF
##   - locates the corresponding AlphaFold summary_confidences JSON
##   - extracts the actual AF3 pTM value
##   - compares it with the pTM currently reported in the thesis
##
## The result becomes the source audit used to build the unified
## structure manifest in the next stage.
## ================================================================== ##


BASE = (
    Path.home()
    / "methanotrophs"
    / "methanotroph_project"
    / "jeppe"
)


GENOME_WF = (
    BASE
    / "genome_analysis_workflow"
    / "18_structure_prediction"
)


META_WF = (
    BASE
    / "metagenome"
    / "comparative_analysis"
    / "final_structure_figures"
)


OUTROOT = (
    BASE
    / "structure_standardization"
)


OUTDIR = (
    OUTROOT
    / "00_audit"
)


OUTDIR.mkdir(
    parents=True,
    exist_ok=True
)


OUTFILE = (
    OUTDIR
    / "alphafold_structure_source_audit.tsv"
)


## ================================================================== ##
## Source paths
## ================================================================== ##


MTOA_ROOT = (
    GENOME_WF
    / "09_mature_heme_structure_figure"
)


MTOA_INPUT_MANIFEST = (
    MTOA_ROOT
    / "01_inputs"
    / "mature_heme_manifest.tsv"
)


MTOA_SELECTED_MANIFEST = (
    MTOA_ROOT
    / "02_results"
    / "selected_model_manifest.tsv"
)


MCA_ROOT = (
    GENOME_WF
    / "10_MCA0421_Cluster00063"
)


MCA_UNPACKED = (
    MCA_ROOT
    / "04_alphafold_results"
    / "unpacked"
)


MCA_CONFIDENCE = (
    MCA_ROOT
    / "05_model_qc"
    / "alphafold3_model_confidence_summary.tsv"
)


META_SELECTED = (
    META_WF
    / "00_mature_alphafold"
    / "qc"
    / "final_12_mature_AF_selected_models.tsv"
)


## ================================================================== ##
## Expected thesis pTM values
##
## These are used ONLY as QC.
##
## Actual pTM values are always read from AlphaFold output JSON.
## ================================================================== ##


EXPECTED_PTM = {

    ## Cluster_00035
    ("cluster00035", "A"):
        0.88,

    ("cluster00035", "B"):
        0.88,


    ## Cyc2
    ("cyc2", "A"):
        0.91,

    ("cyc2", "B"):
        0.88,

    ("cyc2", "C"):
        0.89,

    ("cyc2", "D"):
        0.93,


    ## MCA0421
    ("mca0421", "A"):
        0.92,

    ("mca0421", "B"):
        0.92,

    ("mca0421", "C"):
        0.92,

    ("mca0421", "D"):
        0.95,


    ## Cluster_00121
    ("cluster00121", "A"):
        0.91,

    ("cluster00121", "B"):
        0.90,

    ("cluster00121", "C"):
        0.90,

    ("cluster00121", "D"):
        0.95,

    ("cluster00121", "E"):
        0.88,
}


## ================================================================== ##
## Expected panel definitions
##
## These preserve the current biological panel ordering.
## ================================================================== ##


META_PANELS = {

    "cyc2": [

        {
            "panel": "A",
            "protein_id":
                "contig_38_1172",
        },

        {
            "panel": "B",
            "protein_id":
                "contig_25_317",
        },

        {
            "panel": "C",
            "protein_id":
                "GCA_023229605_000000000176_34",
        },

        {
            "panel": "D",
            "protein_id":
                "GCA_002843155_000000000119_52",
        },
    ],


    "cluster00121": [

        {
            "panel": "A",
            "protein_id":
                "BCRBG_22465_000000000369_1",
        },

        {
            "panel": "B",
            "protein_id":
                "contig_25_1899",
        },

        {
            "panel": "C",
            "protein_id":
                "contig_1431_201",
        },

        {
            "panel": "D",
            "protein_id":
                "contig_122_3517",
        },

        {
            "panel": "E",
            "protein_id":
                "contig_871_20",
        },
    ],
}


MCA_PANELS = [

    {
        "panel": "A",
        "job":
            "fold_bath_mc_9hec",
        "protein_id":
            "MCA0421",
        "genome_label":
            "Methylococcus capsulatus Bath",
    },

    {
        "panel": "B",
        "job":
            "fold_gcf_000424685_000000000006_23_9hec",
        "protein_id":
            "GCF_000424685_000000000006_23",
        "genome_label":
            "Methylococcus capsulatus",
    },

    {
        "panel": "C",
        "job":
            "fold_gcf_036864975_000000000001_380_9hec",
        "protein_id":
            "GCF_036864975_000000000001_380",
        "genome_label":
            "Methylococcus sp. 016106025",
    },

    {
        "panel": "D",
        "job":
            "fold_gca_026399025_000000000114_35_9hec",
        "protein_id":
            "GCA_026399025_000000000114_35",
        "genome_label":
            "Methylobacter_A sp. 026399025",
    },
]


## ================================================================== ##
## Helpers
## ================================================================== ##


def fail(message):

    raise RuntimeError(
        message
    )


def normalize_job(value):

    value = str(
        value
    ).strip().lower()

    value = re.sub(
        r"[^a-z0-9]+",
        "_",
        value
    )

    return value.strip(
        "_"
    )


def read_ptm(json_file):

    with open(
        json_file,
        "r"
    ) as handle:

        data = json.load(
            handle
        )


    ## AlphaFold3 summary_confidences JSON
    ## normally contains a top-level "ptm".

    for key in (
        "ptm",
        "pTM",
        "PTM",
    ):

        if key in data:

            return float(
                data[key]
            )


    fail(
        "No pTM field found in:\n"
        f"{json_file}"
    )


def find_summary_json(
    root,
    job,
    model_number
):

    job_norm = normalize_job(
        job
    )


    model_number = int(
        model_number
    )


    candidates = []


    for path in root.rglob(
        "*summary_confidences*.json"
    ):

        name_norm = normalize_job(
            path.name
        )


        parent_norm = normalize_job(
            path.parent.name
        )


        if (
            job_norm in name_norm
            or job_norm == parent_norm
        ):

            model_pattern = re.compile(
                rf"summary_confidences[_-]{model_number}(?:\D|$)"
            )


            if model_pattern.search(
                path.name
            ):

                candidates.append(
                    path
                )


    candidates = sorted(
        set(
            candidates
        )
    )


    if len(candidates) == 1:

        return candidates[0]


    if len(candidates) == 0:

        fail(
            "Could not locate summary_confidences JSON for:\n"
            f"  job: {job}\n"
            f"  model: {model_number}\n"
            f"  root: {root}"
        )


    fail(
        "Multiple matching summary JSON files found for:\n"
        f"  job: {job}\n"
        f"  model: {model_number}\n\n"
        +
        "\n".join(
            str(x)
            for x in candidates
        )
    )


def check_file(
    path,
    label
):

    if not Path(
        path
    ).is_file():

        fail(
            f"Missing {label}:\n"
            f"{path}"
        )


## ================================================================== ##
## Input validation
## ================================================================== ##


check_file(
    MTOA_INPUT_MANIFEST,
    "Cluster_00035 mature input manifest"
)


check_file(
    MTOA_SELECTED_MANIFEST,
    "Cluster_00035 selected-model manifest"
)


check_file(
    META_SELECTED,
    "metagenome Stage-76 selected-model table"
)


check_file(
    MCA_CONFIDENCE,
    "MCA0421 AlphaFold confidence table"
)


## ================================================================== ##
## Audit rows
## ================================================================== ##


rows = []


## ================================================================== ##
## 1. Cluster_00035
##
## Original genome_analysis_workflow
## ================================================================== ##


mtoa_inputs = pd.read_csv(
    MTOA_INPUT_MANIFEST,
    sep="\t",
    dtype=str,
    keep_default_na=False
)


mtoa_selected = pd.read_csv(
    MTOA_SELECTED_MANIFEST,
    sep="\t",
    dtype=str,
    keep_default_na=False
)


for panel, label in [
    (
        "A",
        "MtoA"
    ),
    (
        "B",
        "FeGenie-negative"
    ),
]:

    input_hit = mtoa_inputs[
        mtoa_inputs["old_label"]
        == label
    ].copy()


    if len(
        input_hit
    ) != 1:

        fail(
            f"Expected one Cluster_00035 "
            f"{label} input row; "
            f"found {len(input_hit)}."
        )


    selected_hit = mtoa_selected[
        (
            mtoa_selected["family"]
            == "C00035"
        )
        &
        (
            mtoa_selected["label"]
            == label
        )
    ].copy()


    if len(
        selected_hit
    ) != 1:

        fail(
            f"Expected one selected Cluster_00035 "
            f"{label} model; "
            f"found {len(selected_hit)}."
        )


    inp = input_hit.iloc[
        0
    ]


    sel = selected_hit.iloc[
        0
    ]


    job = str(
        sel["job_name"]
    )


    model_number = int(
        sel[
            "selected_model_index"
        ]
    )


    cif = Path(
        sel[
            "best_model"
        ]
    )


    check_file(
        cif,
        f"Cluster_00035 {label} CIF"
    )


    summary_json = find_summary_json(
        MTOA_ROOT,
        job,
        model_number
    )


    ptm = read_ptm(
        summary_json
    )


    expected = EXPECTED_PTM[
        (
            "cluster00035",
            panel
        )
    ]


    rows.append({

        "figure_set":
            "cluster00035",

        "panel":
            panel,

        "source_workflow":
            "genome_analysis_workflow",

        "source_root":
            str(
                MTOA_ROOT
            ),

        "protein_id":
            inp[
                "protein_id"
            ],

        "genome_id":
            inp[
                "genome"
            ],

        "genome_label":
            "",

        "annotation":
            (
                "MtoA"
                if label == "MtoA"
                else
                "FeGenie-negative homologue"
            ),

        "job":
            job,

        "model_number":
            model_number,

        "structure_cif":
            str(
                cif
            ),

        "summary_json":
            str(
                summary_json
            ),

        "ptm":
            ptm,

        "expected_thesis_ptm":
            expected,

        "ptm_matches_thesis":
            abs(
                ptm - expected
            ) < 0.005,

        "sequence_length":
            inp.get(
                "mature_length",
                ""
            ),

        "expected_hemes":
            inp.get(
                "cxxch_count",
                ""
            ),
    })


## ================================================================== ##
## 2. Cyc2 and Cluster_00121
##
## Metagenome Stage-76 authoritative mature structures
## ================================================================== ##


meta = pd.read_csv(
    META_SELECTED,
    sep="\t",
    dtype=str,
    keep_default_na=False
)


for figure_set in (
    "cyc2",
    "cluster00121"
):

    for panel_def in META_PANELS[
        figure_set
    ]:

        panel = panel_def[
            "panel"
        ]


        protein_id = panel_def[
            "protein_id"
        ]


        hit = meta[
            (
                meta["panel_set"]
                == figure_set
            )
            &
            (
                meta["protein_id"]
                == protein_id
            )
        ].copy()


        if len(
            hit
        ) != 1:

            fail(
                f"{figure_set} panel {panel}: "
                f"expected one Stage-76 row for "
                f"{protein_id}; found {len(hit)}."
            )


        row = hit.iloc[
            0
        ]


        job = row[
            "job"
        ]


        model_number = int(
            row[
                "model_number"
            ]
        )


        cif = Path(
            row[
                "selected_holo_cif"
            ]
        )


        check_file(
            cif,
            f"{figure_set} {panel} CIF"
        )


        summary_json = find_summary_json(
            META_WF,
            job,
            model_number
        )


        ptm = read_ptm(
            summary_json
        )


        expected = EXPECTED_PTM[
            (
                figure_set,
                panel
            )
        ]


        rows.append({

            "figure_set":
                figure_set,

            "panel":
                panel,

            "source_workflow":
                "metagenome",

            "source_root":
                str(
                    META_WF
                ),

            "protein_id":
                protein_id,

            "genome_id":
                "",

            "genome_label":
                "",

            "annotation":
                (
                    "FeGenie: Cyc2_repCluster2"
                    if figure_set == "cyc2"
                    else
                    "Cluster_00121 9-heme cytochrome"
                ),

            "job":
                job,

            "model_number":
                model_number,

            "structure_cif":
                str(
                    cif
                ),

            "summary_json":
                str(
                    summary_json
                ),

            "ptm":
                ptm,

            "expected_thesis_ptm":
                expected,

            "ptm_matches_thesis":
                abs(
                    ptm - expected
                ) < 0.005,

            "sequence_length":
                row.get(
                    "sequence_length",
                    ""
                ),

            "expected_hemes":
                row.get(
                    "expected_hemes",
                    ""
                ),
        })


## ================================================================== ##
## 3. MCA0421 / Cluster_00063
##
## Original genome_analysis_workflow
##
## The existing figure used model 0 for these four AF3 jobs.
## pTM is taken directly from summary_confidences_0.json.
## ================================================================== ##


mca_conf = pd.read_csv(
    MCA_CONFIDENCE,
    sep="\t",
    dtype=str,
    keep_default_na=False
)


for panel_def in MCA_PANELS:

    panel = panel_def[
        "panel"
    ]


    job = panel_def[
        "job"
    ]


    model_number = 0


    job_dir = (
        MCA_UNPACKED
        / job
    )


    if not job_dir.is_dir():

        fail(
            f"Missing MCA0421 AlphaFold job directory:\n"
            f"{job_dir}"
        )


    cif = (
        job_dir
        / f"{job}_model_{model_number}.cif"
    )


    check_file(
        cif,
        f"MCA0421 panel {panel} CIF"
    )


    summary_json = (
        job_dir
        / (
            f"{job}_summary_confidences_"
            f"{model_number}.json"
        )
    )


    check_file(
        summary_json,
        f"MCA0421 panel {panel} confidence JSON"
    )


    ptm = read_ptm(
        summary_json
    )


    expected = EXPECTED_PTM[
        (
            "mca0421",
            panel
        )
    ]


    rows.append({

        "figure_set":
            "mca0421",

        "panel":
            panel,

        "source_workflow":
            "genome_analysis_workflow",

        "source_root":
            str(
                MCA_ROOT
            ),

        "protein_id":
            panel_def[
                "protein_id"
            ],

        "genome_id":
            "",

        "genome_label":
            panel_def[
                "genome_label"
            ],

        "annotation":
            "MCA0421 / Cluster_00063 9-heme cytochrome",

        "job":
            job,

        "model_number":
            model_number,

        "structure_cif":
            str(
                cif
            ),

        "summary_json":
            str(
                summary_json
            ),

        "ptm":
            ptm,

        "expected_thesis_ptm":
            expected,

        "ptm_matches_thesis":
            abs(
                ptm - expected
            ) < 0.005,

        "sequence_length":
            "",

        "expected_hemes":
            9,
    })


## ================================================================== ##
## Final audit table
## ================================================================== ##


audit = pd.DataFrame(
    rows
)


figure_order = {
    "cluster00035":
        1,
    "cyc2":
        2,
    "mca0421":
        3,
    "cluster00121":
        4,
}


audit[
    "_figure_order"
] = audit[
    "figure_set"
].map(
    figure_order
)


audit = audit.sort_values(
    [
        "_figure_order",
        "panel"
    ]
).drop(
    columns=[
        "_figure_order"
    ]
).reset_index(
    drop=True
)


## ================================================================== ##
## Sanity checks
## ================================================================== ##


if len(
    audit
) != 15:

    fail(
        f"Expected 15 structures across the four figures; "
        f"found {len(audit)}."
    )


if audit[
    "structure_cif"
].duplicated().any():

    duplicated = audit.loc[
        audit[
            "structure_cif"
        ].duplicated(
            keep=False
        ),
        [
            "figure_set",
            "panel",
            "structure_cif"
        ]
    ]

    fail(
        "Duplicate selected CIFs detected:\n"
        +
        duplicated.to_string(
            index=False
        )
    )


audit.to_csv(
    OUTFILE,
    sep="\t",
    index=False
)


## ================================================================== ##
## Console report
## ================================================================== ##


print()

print(
    "=" * 120
)

print(
    "ALPHAFOLD STRUCTURE SOURCE AUDIT"
)

print(
    "=" * 120
)

print()


display = audit[
    [
        "figure_set",
        "panel",
        "source_workflow",
        "protein_id",
        "genome_label",
        "model_number",
        "ptm",
        "expected_thesis_ptm",
        "ptm_matches_thesis"
    ]
].copy()


print(
    display.to_string(
        index=False
    )
)


print()

print(
    "=" * 120
)


n_match = int(
    audit[
        "ptm_matches_thesis"
    ].sum()
)


print(
    f"Structures audited:        {len(audit)}"
)

print(
    f"pTM matching thesis:       {n_match}/{len(audit)}"
)

print()


if n_match != len(
    audit
):

    print(
        "WARNING: Some AlphaFold pTM values differ from "
        "the values currently reported in the thesis."
    )

    print()

    print(
        audit.loc[
            ~audit[
                "ptm_matches_thesis"
            ],
            [
                "figure_set",
                "panel",
                "protein_id",
                "ptm",
                "expected_thesis_ptm",
                "summary_json"
            ]
        ].to_string(
            index=False
        )
    )

    print()


print(
    f"Audit table:\n{OUTFILE}"
)

print()

print(
    "No source files were modified."
)
