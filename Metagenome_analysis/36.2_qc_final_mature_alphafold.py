#!/usr/bin/env python3

from __future__ import annotations

import json
import math
import re
import shutil
import sys
import zipfile
from functools import lru_cache
from pathlib import Path

import pandas as pd

try:
    from Bio.PDB import MMCIFParser
    from Bio.SeqUtils import seq1

except ImportError:

    sys.exit(
        "\nERROR: Biopython is required.\n"
        "Run this script from the genome_viz environment.\n"
    )


## ================================================================== ##
## STAGE 76
##
## QC FINAL SIGNAL-PEPTIDE-TRIMMED ALPHAFOLD STRUCTURES
##
## Final figure set:
##
##   Cyc2                         4 proteins x 1 HEC
##   Cluster_00121               5 proteins x 9 HEC
##   Methylobacter multiheme     5 + 8 + 7 HEC
##
## Total:
##
##   12 proteins
##   69 expected CXXCH / HEC sites
##
##
## For every AlphaFold model:
##
##   - recover the mature protein chain
##   - identify canonical CXXCH motifs
##   - identify HEC ligands
##   - assign CXXCH motifs to HEC ligands one-to-one
##   - measure:
##
##         Cys1-S -> nearest HEC atom
##         Cys2-S -> nearest HEC atom
##         His-NE2 -> HEC Fe
##
## A heme site passes when all three distances are <= 3.0 Å.
##
##
## Model selection:
##
##   1. retain models in which ALL expected heme sites pass
##   2. choose highest AlphaFold ranking score
##   3. model number breaks ties
##
##
## Selected mature holo CIFs are copied unchanged into:
##
##   00_mature_alphafold/selected_models/
##
## These become the authoritative structures for the three final
## structure figures.
## ================================================================== ##


PROJECT = (
    Path.home()
    / "methanotrophs"
    / "methanotroph_project"
    / "jeppe"
    / "metagenome"
)


ROOT = (
    PROJECT
    / "comparative_analysis"
    / "final_structure_figures"
    / "00_mature_alphafold"
)


RAW_DIR = (
    ROOT
    / "raw_data"
)


EXTRACTED_DIR = (
    ROOT
    / "extracted"
)


QC_DIR = (
    ROOT
    / "qc"
)


SELECTED_DIR = (
    ROOT
    / "selected_models"
)


TRIM_TABLE = (
    ROOT
    / "final_12_mature_holo_proteins.tsv"
)


OUT_ALL = (
    QC_DIR
    / "final_12_mature_AF_heme_geometry_all_models.tsv"
)


OUT_SELECTED = (
    QC_DIR
    / "final_12_mature_AF_selected_models.tsv"
)


OUT_SELECTED_SITES = (
    QC_DIR
    / "final_12_mature_AF_selected_model_heme_sites.tsv"
)


OUT_ALL_SITES = (
    QC_DIR
    / "final_12_mature_AF_heme_sites_all_models.tsv"
)


EXPECTED_JOBS = 12

EXPECTED_MODELS_PER_JOB = 5

EXPECTED_TOTAL_HEMES = 69

DISTANCE_THRESHOLD = 3.0


## ================================================================== ##
## Helpers
## ================================================================== ##

def fail(message):

    print(
        f"\nERROR: {message}",
        file=sys.stderr,
    )

    sys.exit(1)


def normalize_name(value):

    value = str(
        value
    ).strip().lower()


    value = re.sub(
        r"[^a-z0-9]+",
        "_",
        value,
    )


    value = value.strip(
        "_"
    )


    ## AlphaFold Server prefixes downloaded result files and CIFs with
    ## "fold_", whereas the original submitted job name does not
    ## necessarily contain this prefix. Treat both forms as the same
    ## biological job.
    value = re.sub(
        r"^fold_",
        "",
        value,
    )


    return value


def distance(
    atom1,
    atom2,
):

    return float(
        atom1
        -
        atom2
    )


def min_distance_to_residue(
    atom,
    residue,
):

    values = [

        distance(
            atom,
            other,
        )

        for other
        in residue.get_atoms()
    ]


    if not values:

        return math.inf


    return min(
        values
    )


def get_protein_chain(
    structure,
):

    model = next(
        structure.get_models()
    )


    ## Prefer chain A, as used by AlphaFold Server.

    if "A" in model:

        chain = model[
            "A"
        ]


        residues = [

            residue

            for residue
            in chain.get_residues()

            if residue.id[
                0
            ]
            ==
            " "
        ]


        if residues:

            return (
                chain,
                residues,
            )


    ## Fallback if AlphaFold used another chain ID.

    for chain in model:

        residues = [

            residue

            for residue
            in chain.get_residues()

            if residue.id[
                0
            ]
            ==
            " "
        ]


        if residues:

            return (
                chain,
                residues,
            )


    fail(
        "Could not identify protein polymer chain."
    )


def residues_to_sequence(
    residues,
):

    sequence = ""


    for residue in residues:

        try:

            aa = seq1(
                residue.resname
            )

        except Exception:

            aa = "X"


        sequence += aa


    return sequence


def get_hemes(
    structure,
):

    return [

        residue

        for residue
        in structure.get_residues()

        if (
            residue.resname
            .strip()
            .upper()
            ==
            "HEC"
        )
    ]


def get_motifs(
    sequence,
    residues,
):

    matches = list(
        re.finditer(
            r"C..CH",
            sequence,
        )
    )


    motifs = []


    for number, match in enumerate(
        matches,
        start=1,
    ):

        start = match.start()


        if (
            start
            +
            4
            >=
            len(
                residues
            )
        ):

            fail(
                "CXXCH motif exceeds protein residue list."
            )


        cys1 = residues[
            start
        ]


        cys2 = residues[
            start
            +
            3
        ]


        his = residues[
            start
            +
            4
        ]


        if (
            cys1.resname.upper()
            !=
            "CYS"
            or
            cys2.resname.upper()
            !=
            "CYS"
            or
            his.resname.upper()
            !=
            "HIS"
        ):

            fail(
                f"Sequence/residue mismatch at CXXCH motif "
                f"{number}."
            )


        motifs.append(
            {
                "motif_number":
                    number,

                "sequence_start_1based":
                    start
                    +
                    1,

                "cys1":
                    cys1,

                "cys2":
                    cys2,

                "his":
                    his,
            }
        )


    return motifs


def pair_geometry(
    motif,
    heme,
):

    cys1 = motif[
        "cys1"
    ]


    cys2 = motif[
        "cys2"
    ]


    his = motif[
        "his"
    ]


    if "SG" not in cys1:

        return None


    if "SG" not in cys2:

        return None


    if "NE2" not in his:

        return None


    if "FE" not in heme:

        return None


    d_cys1 = min_distance_to_residue(
        cys1[
            "SG"
        ],
        heme,
    )


    d_cys2 = min_distance_to_residue(
        cys2[
            "SG"
        ],
        heme,
    )


    d_his_fe = distance(
        his[
            "NE2"
        ],
        heme[
            "FE"
        ],
    )


    return {
        "cys1_heme_distance":
            d_cys1,

        "cys2_heme_distance":
            d_cys2,

        "his_fe_distance":
            d_his_fe,

        "assignment_score":
            (
                d_cys1
                +
                d_cys2
                +
                d_his_fe
            ),
    }


## ================================================================== ##
## Optimal one-to-one motif -> HEC assignment
##
## Maximum n = 9, so exact bit-mask dynamic programming is trivial.
## ================================================================== ##

def optimal_assignment(
    cost_matrix,
):

    n = len(
        cost_matrix
    )


    @lru_cache(
        maxsize=None
    )
    def solve(
        motif_index,
        used_mask,
    ):

        if motif_index == n:

            return (
                0.0,
                (),
            )


        best_cost = math.inf

        best_assignment = None


        for heme_index in range(
            n
        ):

            bit = (
                1
                <<
                heme_index
            )


            if (
                used_mask
                &
                bit
            ):

                continue


            remaining_cost, remaining_assignment = solve(

                motif_index
                +
                1,

                used_mask
                |
                bit,
            )


            total = (
                cost_matrix[
                    motif_index
                ][
                    heme_index
                ]
                +
                remaining_cost
            )


            if total < best_cost:

                best_cost = total

                best_assignment = (
                    heme_index,
                ) + remaining_assignment


        return (
            best_cost,
            best_assignment,
        )


    return solve(
        0,
        0,
    )


def read_ranking_score(
    cif_path,
    model_number,
):

    parent = cif_path.parent


    prefix = re.sub(
        rf"_model_{model_number}\.cif$",
        "",
        cif_path.name,
    )


    candidates = [

        parent
        / f"{prefix}_summary_confidences_{model_number}.json",

        parent
        / f"{prefix}_summary_confidence_{model_number}.json",
    ]


    candidates += sorted(
        parent.glob(
            f"*summary*{model_number}.json"
        )
    )


    seen = set()


    for path in candidates:

        if path in seen:

            continue


        seen.add(
            path
        )


        if not path.is_file():

            continue


        try:

            data = json.loads(
                path.read_text()
            )

        except Exception:

            continue


        for key in [

            "ranking_score",

            "ranking_confidence",
        ]:

            if key not in data:

                continue


            try:

                return (
                    float(
                        data[
                            key
                        ]
                    ),
                    str(
                        path
                    ),
                )

            except Exception:

                pass


    return (
        math.nan,
        "",
    )


## ================================================================== ##
## Input QC
## ================================================================== ##

if not RAW_DIR.is_dir():

    fail(
        f"Raw AlphaFold directory does not exist:\n"
        f"{RAW_DIR}"
    )


if not TRIM_TABLE.is_file():

    fail(
        f"Stage-75 trimming table does not exist:\n"
        f"{TRIM_TABLE}"
    )


trim = pd.read_csv(
    TRIM_TABLE,
    sep="\t",
    dtype=str,
    keep_default_na=False,
)


required_columns = [

    "panel_set",

    "protein_id",

    "mature_job",

    "mature_length",

    "expected_hemes",

    "mature_CXXCH_count",
]


missing_columns = [

    column

    for column
    in required_columns

    if column not in trim.columns
]


if missing_columns:

    fail(
        "Stage-75 table is missing column(s):\n"
        +
        "\n".join(
            missing_columns
        )
    )


if len(
    trim
) != EXPECTED_JOBS:

    fail(
        f"Expected {EXPECTED_JOBS} Stage-75 proteins; "
        f"found {len(trim)}."
    )


trim[
    "_normalized_job"
] = trim[
    "mature_job"
].apply(
    normalize_name
)


if trim[
    "_normalized_job"
].duplicated().any():

    fail(
        "Duplicate normalized mature-job names in Stage-75 table."
    )


metadata = {

    row[
        "_normalized_job"
    ]:
        row

    for _,
    row
    in trim.iterrows()
}


## ================================================================== ##
## Prepare clean output directories
## ================================================================== ##

QC_DIR.mkdir(
    parents=True,
    exist_ok=True,
)


SELECTED_DIR.mkdir(
    parents=True,
    exist_ok=True,
)


EXTRACTED_DIR.mkdir(
    parents=True,
    exist_ok=True,
)


## Remove stale selected CIFs from previous attempts.

for path in SELECTED_DIR.glob(
    "*.cif"
):

    path.unlink()


## ================================================================== ##
## Extract all 12 AlphaFold archives
## ================================================================== ##

archives = sorted(
    RAW_DIR.glob(
        "*.zip"
    )
)


if len(
    archives
) != EXPECTED_JOBS:

    fail(
        f"Expected {EXPECTED_JOBS} AlphaFold ZIP archives; "
        f"found {len(archives)}."
    )


print("=" * 100)

print(
    "STAGE 76 - QC FINAL MATURE ALPHAFOLD STRUCTURES"
)

print("=" * 100)

print()

print(
    f"AlphaFold archives: {len(archives)}"
)


for archive in archives:

    destination = (
        EXTRACTED_DIR
        / archive.stem
    )


    ## Clean the per-job extraction directory on reruns.

    if destination.exists():

        shutil.rmtree(
            destination
        )


    destination.mkdir(
        parents=True,
        exist_ok=True,
    )


    with zipfile.ZipFile(
        archive,
        "r",
    ) as handle:

        handle.extractall(
            destination
        )


    print(
        f"  extracted: {archive.name}"
    )


## ================================================================== ##
## Discover CIF models
## ================================================================== ##

pattern = re.compile(
    r"^(?P<job>.+)_model_(?P<model>\d+)\.cif$"
)


discovered = []


for cif in sorted(
    EXTRACTED_DIR.rglob(
        "*_model_*.cif"
    )
):

    match = pattern.match(
        cif.name
    )


    if match is None:

        continue


    job = match.group(
        "job"
    )


    normalized_job = normalize_name(
        job
    )


    if normalized_job not in metadata:

        fail(
            f"AlphaFold job not found in Stage-75 metadata:\n"
            f"  CIF job: {job}\n"
            f"  normalized: {normalized_job}"
        )


    meta = metadata[
        normalized_job
    ]


    discovered.append(
        {
            "job":
                job,

            "normalized_job":
                normalized_job,

            "model_number":
                int(
                    match.group(
                        "model"
                    )
                ),

            "expected_hemes":
                int(
                    meta[
                        "expected_hemes"
                    ]
                ),

            "expected_length":
                int(
                    meta[
                        "mature_length"
                    ]
                ),

            "panel_set":
                meta[
                    "panel_set"
                ],

            "protein_id":
                meta[
                    "protein_id"
                ],

            "cif":
                cif,
        }
    )


jobs = sorted(
    set(
        row[
            "normalized_job"
        ]

        for row
        in discovered
    )
)


if len(
    jobs
) != EXPECTED_JOBS:

    fail(
        f"Expected {EXPECTED_JOBS} mature AlphaFold jobs; "
        f"found {len(jobs)}."
    )


print()

print(
    f"Discovered mature jobs: {len(jobs)}"
)


## ================================================================== ##
## QC every model
## ================================================================== ##

parser = MMCIFParser(
    QUIET=True
)


model_rows = []

site_rows = []


for normalized_job in jobs:

    job_models = sorted(

        [

            row

            for row
            in discovered

            if row[
                "normalized_job"
            ]
            ==
            normalized_job
        ],

        key=lambda row:
            row[
                "model_number"
            ],
    )


    if len(
        job_models
    ) != EXPECTED_MODELS_PER_JOB:

        fail(
            f"{normalized_job}: expected "
            f"{EXPECTED_MODELS_PER_JOB} models; "
            f"found {len(job_models)}."
        )


    first = job_models[
        0
    ]


    expected_hemes = first[
        "expected_hemes"
    ]


    expected_length = first[
        "expected_length"
    ]


    print()

    print(
        f"{first['job']} "
        f"[{first['panel_set']} | "
        f"{first['protein_id']} | "
        f"{expected_hemes} HEC]"
    )


    for item in job_models:

        model_number = item[
            "model_number"
        ]


        cif = item[
            "cif"
        ]


        structure = parser.get_structure(
            (
                item[
                    "job"
                ]
                +
                f"_{model_number}"
            ),
            str(
                cif
            ),
        )


        chain, residues = get_protein_chain(
            structure
        )


        sequence = residues_to_sequence(
            residues
        )


        if len(
            sequence
        ) != expected_length:

            fail(
                f"{item['job']} model {model_number}: "
                f"structure length {len(sequence)} != "
                f"Stage-75 mature length {expected_length}."
            )


        motifs = get_motifs(
            sequence,
            residues,
        )


        hemes = get_hemes(
            structure
        )


        if len(
            motifs
        ) != expected_hemes:

            fail(
                f"{item['job']} model {model_number}: "
                f"found {len(motifs)} CXXCH motifs; "
                f"expected {expected_hemes}."
            )


        if len(
            hemes
        ) != expected_hemes:

            fail(
                f"{item['job']} model {model_number}: "
                f"found {len(hemes)} HEC residues; "
                f"expected {expected_hemes}."
            )


        geometry_matrix = []

        geometry_details = []


        for motif in motifs:

            cost_row = []

            detail_row = []


            for heme in hemes:

                geometry = pair_geometry(
                    motif,
                    heme,
                )


                if geometry is None:

                    fail(
                        f"{item['job']} model {model_number}: "
                        "required SG / NE2 / FE atom missing."
                    )


                cost_row.append(
                    geometry[
                        "assignment_score"
                    ]
                )


                detail_row.append(
                    geometry
                )


            geometry_matrix.append(
                cost_row
            )


            geometry_details.append(
                detail_row
            )


        total_cost, assignment = optimal_assignment(
            geometry_matrix
        )


        current_sites = []

        n_pass = 0


        for motif_index, heme_index in enumerate(
            assignment
        ):

            motif = motifs[
                motif_index
            ]


            heme = hemes[
                heme_index
            ]


            geometry = geometry_details[
                motif_index
            ][
                heme_index
            ]


            passed = int(

                geometry[
                    "cys1_heme_distance"
                ]
                <=
                DISTANCE_THRESHOLD

                and

                geometry[
                    "cys2_heme_distance"
                ]
                <=
                DISTANCE_THRESHOLD

                and

                geometry[
                    "his_fe_distance"
                ]
                <=
                DISTANCE_THRESHOLD
            )


            n_pass += passed


            site = {

                "panel_set":
                    item[
                        "panel_set"
                    ],

                "protein_id":
                    item[
                        "protein_id"
                    ],

                "job":
                    item[
                        "job"
                    ],

                "model_number":
                    model_number,

                "expected_hemes":
                    expected_hemes,

                "motif_number":
                    motif[
                        "motif_number"
                    ],

                "motif_start_1based":
                    motif[
                        "sequence_start_1based"
                    ],

                "heme_assignment":
                    heme_index
                    +
                    1,

                "heme_chain":
                    heme.get_parent().id,

                "heme_residue_number":
                    heme.id[
                        1
                    ],

                "cys1_heme_distance":
                    geometry[
                        "cys1_heme_distance"
                    ],

                "cys2_heme_distance":
                    geometry[
                        "cys2_heme_distance"
                    ],

                "his_fe_distance":
                    geometry[
                        "his_fe_distance"
                    ],

                "site_pass_3A":
                    passed,
            }


            current_sites.append(
                site
            )


            site_rows.append(
                site
            )


        ranking_score, summary_file = read_ranking_score(
            cif,
            model_number,
        )


        model_rows.append(
            {
                "panel_set":
                    item[
                        "panel_set"
                    ],

                "protein_id":
                    item[
                        "protein_id"
                    ],

                "job":
                    item[
                        "job"
                    ],

                "normalized_job":
                    normalized_job,

                "model_number":
                    model_number,

                "sequence_length":
                    len(
                        sequence
                    ),

                "expected_hemes":
                    expected_hemes,

                "n_CXXCH":
                    len(
                        motifs
                    ),

                "n_HEC":
                    len(
                        hemes
                    ),

                "sites_passing_3A":
                    n_pass,

                "all_sites_pass":
                    int(
                        n_pass
                        ==
                        expected_hemes
                    ),

                "ranking_score":
                    ranking_score,

                "geometry_assignment_score":
                    total_cost,

                "max_cys1_heme_distance":
                    max(
                        row[
                            "cys1_heme_distance"
                        ]

                        for row
                        in current_sites
                    ),

                "max_cys2_heme_distance":
                    max(
                        row[
                            "cys2_heme_distance"
                        ]

                        for row
                        in current_sites
                    ),

                "max_his_fe_distance":
                    max(
                        row[
                            "his_fe_distance"
                        ]

                        for row
                        in current_sites
                    ),

                "protein_chain":
                    chain.id,

                "cif":
                    str(
                        cif
                    ),

                "summary_file":
                    summary_file,
            }
        )


        ranking_text = (

            f"{ranking_score:.4f}"

            if not math.isnan(
                ranking_score
            )

            else
            "NA"
        )


        print(
            f"  model_{model_number}: "
            f"{n_pass}/{expected_hemes} sites <=3 Å, "
            f"ranking={ranking_text}"
        )


models = pd.DataFrame(
    model_rows
)


sites = pd.DataFrame(
    site_rows
)


models.to_csv(
    OUT_ALL,
    sep="\t",
    index=False,
    na_rep="",
    float_format="%.6f",
)


sites.to_csv(
    OUT_ALL_SITES,
    sep="\t",
    index=False,
    na_rep="",
    float_format="%.6f",
)


## ================================================================== ##
## Select best fully passing model
## ================================================================== ##

selected_rows = []

selected_site_frames = []


for normalized_job, group in models.groupby(
    "normalized_job",
    sort=True,
):

    group = group.copy()


    passing = group[
        group[
            "all_sites_pass"
        ]
        ==
        1
    ].copy()


    if len(
        passing
    ) == 0:

        fail(
            f"{normalized_job}: no model passed all "
            "heme-geometry sites."
        )


    ## Highest ranking score.
    ##
    ## If ranking score is missing, model number becomes the fallback.
    ## Stable sorting makes model_0 win exact ties.

    passing[
        "_ranking_sort"
    ] = pd.to_numeric(
        passing[
            "ranking_score"
        ],
        errors="coerce",
    ).fillna(
        -math.inf
    )


    passing = passing.sort_values(
        [
            "_ranking_sort",
            "model_number",
        ],
        ascending=[
            False,
            True,
        ],
        kind="stable",
    )


    best = passing.iloc[
        0
    ].copy()


    source_cif = Path(
        best[
            "cif"
        ]
    )


    destination_cif = (
        SELECTED_DIR
        /
        (
            str(
                best[
                    "job"
                ]
            )
            +
            "_selected_model_"
            +
            str(
                int(
                    best[
                        "model_number"
                    ]
                )
            )
            +
            ".cif"
        )
    )


    shutil.copy2(
        source_cif,
        destination_cif,
    )


    best[
        "selected_holo_cif"
    ] = str(
        destination_cif
    )


    ## Add Stage-75 trimming metadata.

    stage75 = metadata[
        normalized_job
    ]


    for column in [

        "signalp_prediction",

        "signalp_cs_position",

        "removed_N_terminal_aa",

        "original_length",

        "mature_length",

        "sequence_changed",
    ]:

        if column in stage75.index:

            best[
                column
            ] = stage75[
                column
            ]


    selected_rows.append(
        best
    )


    selected_site_frames.append(

        sites[
            (
                sites[
                    "normalized_job"
                ]
                if "normalized_job"
                in sites.columns
                else
                sites[
                    "job"
                ].apply(
                    normalize_name
                )
            )
            ==
            normalized_job
        ][
            sites[
                "model_number"
            ]
            ==
            int(
                best[
                    "model_number"
                ]
            )
        ].copy()
    )


selected = pd.DataFrame(
    selected_rows
)


## Rebuild selected site table explicitly to avoid chained-indexing
## ambiguity.

selected_site_parts = []


for _, best in selected.iterrows():

    block = sites[
        (
            sites[
                "job"
            ]
            ==
            best[
                "job"
            ]
        )
        &
        (
            sites[
                "model_number"
            ]
            ==
            int(
                best[
                    "model_number"
                ]
            )
        )
    ].copy()


    selected_site_parts.append(
        block
    )


selected_sites = pd.concat(
    selected_site_parts,
    ignore_index=True,
)


selected.to_csv(
    OUT_SELECTED,
    sep="\t",
    index=False,
    na_rep="",
    float_format="%.6f",
)


selected_sites.to_csv(
    OUT_SELECTED_SITES,
    sep="\t",
    index=False,
    na_rep="",
    float_format="%.6f",
)


## ================================================================== ##
## Final locked QC
## ================================================================== ##

if len(
    selected
) != EXPECTED_JOBS:

    fail(
        f"Expected {EXPECTED_JOBS} selected structures; "
        f"found {len(selected)}."
    )


observed_total_hemes = int(

    pd.to_numeric(
        selected[
            "expected_hemes"
        ]
    ).sum()
)


if (
    observed_total_hemes
    !=
    EXPECTED_TOTAL_HEMES
):

    fail(
        f"Selected structures contain "
        f"{observed_total_hemes} expected hemes; "
        f"expected {EXPECTED_TOTAL_HEMES}."
    )


observed_passing_sites = int(

    pd.to_numeric(
        selected[
            "sites_passing_3A"
        ]
    ).sum()
)


if (
    observed_passing_sites
    !=
    EXPECTED_TOTAL_HEMES
):

    fail(
        f"Selected structures have "
        f"{observed_passing_sites}/"
        f"{EXPECTED_TOTAL_HEMES} passing sites."
    )


selected_cifs = sorted(
    SELECTED_DIR.glob(
        "*.cif"
    )
)


if len(
    selected_cifs
) != EXPECTED_JOBS:

    fail(
        f"Expected {EXPECTED_JOBS} selected CIFs; "
        f"found {len(selected_cifs)}."
    )


## Figure-set counts.

set_counts = (
    selected[
        "panel_set"
    ]
    .value_counts()
    .to_dict()
)


expected_set_counts = {

    "cyc2":
        4,

    "cluster00121":
        5,

    "methylobacter":
        3,
}


if set_counts != expected_set_counts:

    fail(
        f"Unexpected figure-set composition:\n"
        f"Observed: {set_counts}\n"
        f"Expected: {expected_set_counts}"
    )


## ================================================================== ##
## Report
## ================================================================== ##

print()

print("=" * 100)

print(
    "STAGE 76 COMPLETE"
)

print("=" * 100)

print()

print(
    "Selected mature models:"
)

print()


display = selected[
    [
        "panel_set",
        "protein_id",
        "sequence_length",
        "expected_hemes",
        "model_number",
        "ranking_score",
        "sites_passing_3A",
        "max_cys1_heme_distance",
        "max_cys2_heme_distance",
        "max_his_fe_distance",
    ]
].copy()


print(
    display.to_string(
        index=False
    )
)


print()

print(
    f"Selected structures:          {len(selected)}"
)

print(
    f"Expected heme sites:          {EXPECTED_TOTAL_HEMES}"
)

print(
    f"Selected sites passing QC:    "
    f"{observed_passing_sites}/"
    f"{EXPECTED_TOTAL_HEMES}"
)

print()

print(
    "Figure sets:"
)

for figure_set in [

    "cyc2",

    "cluster00121",

    "methylobacter",
]:

    print(
        f"  {figure_set:<16} "
        f"{set_counts[figure_set]}"
    )


print()

print(
    "Selected mature holo CIFs:"
)


for path in selected_cifs:

    print(
        f"  {path}"
    )


print()

print(
    "All-model QC:"
)

print(
    f"  {OUT_ALL}"
)

print()

print(
    "All heme-site QC:"
)

print(
    f"  {OUT_ALL_SITES}"
)

print()

print(
    "Selected-model QC:"
)

print(
    f"  {OUT_SELECTED}"
)

print()

print(
    "Selected heme-site QC:"
)

print(
    f"  {OUT_SELECTED_SITES}"
)
