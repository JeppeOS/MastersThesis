#!/usr/bin/env python3

from __future__ import annotations

import json
import math
import re
import shutil
import sys
import zipfile
from pathlib import Path

import pandas as pd

try:
    from Bio.PDB import MMCIFParser
    from Bio.SeqUtils import seq1
except ImportError:
    sys.exit(
        "\nERROR: Biopython is required.\n"
        "Run this script from an environment containing Biopython.\n"
    )


## ================================================================== ##
## STAGE 74
##
## QC HEME-AWARE ALPHAFOLD MODELS FOR THE FOUR CYC2 PROTEINS
##
## Each structure is expected to contain:
##
##     1 canonical CXXCH motif
##     1 HEC ligand
##
## Geometry checked:
##
##     Cys1-S -> nearest HEC atom
##     Cys2-S -> nearest HEC atom
##     His-NE2 -> HEC Fe
##
## Pass criterion:
##
##     all three distances <= 3.0 Å
##
## Model selection:
##
##     highest AlphaFold ranking score
##     model number used as tie-breaker
##
## Selected holo CIF files are copied unchanged for later visualization.
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
    / "Cyc2_four_panel"
    / "structure_analysis"
    / "alphafold"
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


OUT_ALL = (
    QC_DIR
    / "Cyc2_four_panel_AF_heme_geometry_all_models.tsv"
)


OUT_SELECTED = (
    QC_DIR
    / "Cyc2_four_panel_AF_selected_models.tsv"
)


OUT_SITES = (
    QC_DIR
    / "Cyc2_four_panel_AF_selected_model_heme_sites.tsv"
)


EXPECTED_JOBS = 4
EXPECTED_MODELS_PER_JOB = 5
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


def distance(atom1, atom2):

    return float(
        atom1 - atom2
    )


def min_distance_to_residue(atom, residue):

    distances = [
        distance(
            atom,
            other,
        )

        for other
        in residue.get_atoms()
    ]

    if not distances:
        return math.inf

    return min(
        distances
    )


def read_summary_score(
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


    candidates += list(
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

            if key in data:

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


def get_protein_chain(structure):

    model = next(
        structure.get_models()
    )


    if "A" in model:

        chain = model[
            "A"
        ]


        polymer = [
            residue

            for residue
            in chain.get_residues()

            if residue.id[
                0
            ] == " "
        ]


        if polymer:

            return (
                chain,
                polymer,
            )


    for chain in model:

        polymer = [
            residue

            for residue
            in chain.get_residues()

            if residue.id[
                0
            ] == " "
        ]


        if polymer:

            return (
                chain,
                polymer,
            )


    fail(
        "Could not identify protein polymer chain."
    )


def residues_to_sequence(residues):

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


def get_hemes(structure):

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


## ================================================================== ##
## Prepare directories
## ================================================================== ##

for directory in [
    EXTRACTED_DIR,
    QC_DIR,
    SELECTED_DIR,
]:

    directory.mkdir(
        parents=True,
        exist_ok=True,
    )


## ================================================================== ##
## Extract AlphaFold archives
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
        f"Expected {EXPECTED_JOBS} ZIP archives; "
        f"found {len(archives)} in:\n{RAW_DIR}"
    )


print(
    f"Found {len(archives)} AlphaFold archives."
)


for archive in archives:

    destination = (
        EXTRACTED_DIR
        / archive.stem
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
## Discover AlphaFold CIF models
## ================================================================== ##

cifs = sorted(
    EXTRACTED_DIR.rglob(
        "*_model_*.cif"
    )
)


if not cifs:

    fail(
        "No AlphaFold CIF files found."
    )


pattern = re.compile(
    r"^(?P<job>.+)_model_(?P<model>\d+)\.cif$"
)


discovered = []


for cif in cifs:

    match = pattern.match(
        cif.name
    )


    if match is None:
        continue


    discovered.append(
        {
            "job":
                match.group(
                    "job"
                ),

            "model_number":
                int(
                    match.group(
                        "model"
                    )
                ),

            "cif":
                cif,
        }
    )


jobs = sorted(
    set(
        row[
            "job"
        ]

        for row
        in discovered
    )
)


if len(
    jobs
) != EXPECTED_JOBS:

    fail(
        f"Expected {EXPECTED_JOBS} jobs; "
        f"found {len(jobs)}:\n"
        +
        "\n".join(
            jobs
        )
    )


print()

print(
    f"Discovered {len(jobs)} Cyc2 jobs."
)


## ================================================================== ##
## QC all models
## ================================================================== ##

parser = MMCIFParser(
    QUIET=True
)


model_rows = []
site_rows = []


for job in jobs:

    models = sorted(
        [
            row

            for row
            in discovered

            if row[
                "job"
            ]
            ==
            job
        ],
        key=lambda x:
            x[
                "model_number"
            ],
    )


    if len(
        models
    ) != EXPECTED_MODELS_PER_JOB:

        fail(
            f"{job}: expected "
            f"{EXPECTED_MODELS_PER_JOB} models; "
            f"found {len(models)}."
        )


    print()

    print(
        job
    )


    for row in models:

        model_number = row[
            "model_number"
        ]

        cif = row[
            "cif"
        ]


        structure = parser.get_structure(
            f"{job}_{model_number}",
            str(
                cif
            ),
        )


        chain, polymer = get_protein_chain(
            structure
        )


        sequence = residues_to_sequence(
            polymer
        )


        motifs = list(
            re.finditer(
                r"C..CH",
                sequence,
            )
        )


        if len(
            motifs
        ) != 1:

            fail(
                f"{job} model {model_number}: "
                f"expected exactly one CXXCH motif; "
                f"found {len(motifs)}."
            )


        motif = motifs[
            0
        ]


        start = motif.start()


        cys1 = polymer[
            start
        ]

        cys2 = polymer[
            start + 3
        ]

        his = polymer[
            start + 4
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
                f"{job} model {model_number}: "
                "CXXCH sequence/residue mapping failed."
            )


        if "SG" not in cys1:

            fail(
                f"{job} model {model_number}: "
                "Cys1 SG atom missing."
            )


        if "SG" not in cys2:

            fail(
                f"{job} model {model_number}: "
                "Cys2 SG atom missing."
            )


        if "NE2" not in his:

            fail(
                f"{job} model {model_number}: "
                "His NE2 atom missing."
            )


        hemes = get_hemes(
            structure
        )


        if len(
            hemes
        ) != 1:

            fail(
                f"{job} model {model_number}: "
                f"expected exactly one HEC; "
                f"found {len(hemes)}."
            )


        heme = hemes[
            0
        ]


        if "FE" not in heme:

            fail(
                f"{job} model {model_number}: "
                "HEC FE atom missing."
            )


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


        passed = int(
            d_cys1
            <=
            DISTANCE_THRESHOLD
            and
            d_cys2
            <=
            DISTANCE_THRESHOLD
            and
            d_his_fe
            <=
            DISTANCE_THRESHOLD
        )


        ranking_score, summary_file = read_summary_score(
            cif,
            model_number,
        )


        model_rows.append(
            {
                "job":
                    job,

                "model_number":
                    model_number,

                "sequence_length":
                    len(
                        sequence
                    ),

                "protein_chain":
                    chain.id,

                "CXXCH_start_1based":
                    start + 1,

                "n_CXXCH":
                    1,

                "n_HEC":
                    1,

                "cys1_heme_distance":
                    d_cys1,

                "cys2_heme_distance":
                    d_cys2,

                "his_fe_distance":
                    d_his_fe,

                "site_pass_3A":
                    passed,

                "ranking_score":
                    ranking_score,

                "cif":
                    str(
                        cif
                    ),

                "summary_file":
                    summary_file,
            }
        )


        site_rows.append(
            {
                "job":
                    job,

                "model_number":
                    model_number,

                "CXXCH_start_1based":
                    start + 1,

                "heme_chain":
                    heme.get_parent().id,

                "heme_residue_number":
                    heme.id[
                        1
                    ],

                "cys1_heme_distance":
                    d_cys1,

                "cys2_heme_distance":
                    d_cys2,

                "his_fe_distance":
                    d_his_fe,

                "site_pass_3A":
                    passed,
            }
        )


        ranking_text = (
            f", ranking={ranking_score:.4f}"
            if not math.isnan(
                ranking_score
            )
            else
            ""
        )


        print(
            f"  model_{model_number}: "
            f"{passed}/1 site <=3 Å"
            f"{ranking_text}"
        )


models_df = pd.DataFrame(
    model_rows
)


sites_df = pd.DataFrame(
    site_rows
)


models_df.to_csv(
    OUT_ALL,
    sep="\t",
    index=False,
    na_rep="",
    float_format="%.6f",
)


## ================================================================== ##
## Select highest-ranked model per Cyc2 protein
## ================================================================== ##

selected_rows = []
selected_sites = []


for job, group in models_df.groupby(
    "job",
    sort=True,
):

    group = group.copy()


    if group[
        "ranking_score"
    ].notna().any():

        group = group.sort_values(
            [
                "ranking_score",
                "site_pass_3A",
                "model_number",
            ],
            ascending=[
                False,
                False,
                True,
            ],
            kind="stable",
        )


    else:

        group = group.sort_values(
            [
                "site_pass_3A",
                "model_number",
            ],
            ascending=[
                False,
                True,
            ],
            kind="stable",
        )


    best = group.iloc[
        0
    ].copy()


    if int(
        best[
            "site_pass_3A"
        ]
    ) != 1:

        fail(
            f"Selected model for {job} "
            "does not pass heme geometry QC."
        )


    source_cif = Path(
        best[
            "cif"
        ]
    )


    selected_cif = (
        SELECTED_DIR
        /
        f"{job}_selected_model_{int(best['model_number'])}.cif"
    )


    shutil.copy2(
        source_cif,
        selected_cif,
    )


    best[
        "selected_holo_cif"
    ] = str(
        selected_cif
    )


    selected_rows.append(
        best
    )


    site = sites_df[
        (
            sites_df[
                "job"
            ]
            ==
            job
        )
        &
        (
            sites_df[
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


    selected_sites.append(
        site
    )


selected_df = pd.DataFrame(
    selected_rows
)


selected_sites_df = pd.concat(
    selected_sites,
    ignore_index=True,
)


selected_df.to_csv(
    OUT_SELECTED,
    sep="\t",
    index=False,
    na_rep="",
    float_format="%.6f",
)


selected_sites_df.to_csv(
    OUT_SITES,
    sep="\t",
    index=False,
    na_rep="",
    float_format="%.6f",
)


## ================================================================== ##
## Final QC
## ================================================================== ##

if len(
    selected_df
) != EXPECTED_JOBS:

    fail(
        f"Expected {EXPECTED_JOBS} selected models; "
        f"found {len(selected_df)}."
    )


if int(
    selected_df[
        "site_pass_3A"
    ].sum()
) != EXPECTED_JOBS:

    fail(
        "Not all four selected Cyc2 models passed geometry QC."
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
        f"Expected four selected CIF files; "
        f"found {len(selected_cifs)}."
    )


## ================================================================== ##
## Report
## ================================================================== ##

print()

print("=" * 80)

print(
    "STAGE 74 COMPLETE"
)

print("=" * 80)

print()

print(
    "Selected models:"
)

print(
    selected_df[
        [
            "job",
            "sequence_length",
            "model_number",
            "ranking_score",
            "cys1_heme_distance",
            "cys2_heme_distance",
            "his_fe_distance",
        ]
    ]
    .to_string(
        index=False
    )
)


print()

print(
    f"Selected Cyc2 structures: {len(selected_df)}"
)

print(
    "Selected heme sites passing QC: "
    f"{int(selected_df['site_pass_3A'].sum())}/4"
)

print()

print(
    "Selected holo CIF files:"
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
    f"  {OUT_SITES}"
)
