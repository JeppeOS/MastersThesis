#!/usr/bin/env python3

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd


## ================================================================== ##
## STAGE 77
##
## BUILD FINAL STRUCTURE-FIGURE MANIFESTS
##
## Uses ONLY the authoritative mature-protein structures selected by
## Stage 76.
##
## Final figures:
##
##   1. Cyc2 four-panel
##   2. Cluster_00121 five-panel
##   3. Methylobacter multiheme three-panel
##
## Output manifests are written directly to the existing final figure
## workflow so that 02_render_structures.py and 03_assemble_figure.py
## can use them.
## ================================================================== ##


PROJECT = (
    Path.home()
    / "methanotrophs"
    / "methanotroph_project"
    / "jeppe"
    / "metagenome"
)


FIGROOT = (
    PROJECT
    / "comparative_analysis"
    / "final_structure_figures"
)


SELECTED_TSV = (
    FIGROOT
    / "00_mature_alphafold"
    / "qc"
    / "final_12_mature_AF_selected_models.tsv"
)


## ================================================================== ##
## Figure definitions
##
## Foldseek labels:
##
## Cluster_00121:
##   top AFDB50 annotations are predominantly IPT/TIG-like.
##   The recurrent MtrC/MtrF-like structural hit is retained as the
##   biologically relevant structural context.
##
## Methylobacter:
##   Cluster_00042 -> penta-heme cytochrome c552 PDB100 hit
##   Cluster_00288 -> OcwA PDB100 hit
##   no-match      -> OcwA PDB100 hit
##
## Cyc2:
##   no Foldseek analysis was performed; use the FeGenie Cyc2 call.
## ================================================================== ##


FIGURES = {

    "cyc2": {

        "ncols": 2,

        "align_to_first": 1,

        "panels": [

            {
                "panel_id": "A",
                "protein_id": "contig_38_1172",
                "title": "MAG Cluster_00024",
                "subtitle": "Methylobacter_A",
                "evidence": "FeGenie: Cyc2_repCluster2",
            },

            {
                "panel_id": "B",
                "protein_id": "contig_25_317",
                "title": "MAG Cluster_00050",
                "subtitle": "Methylumidiphilus",
                "evidence": "FeGenie: Cyc2_repCluster2",
            },

            {
                "panel_id": "C",
                "protein_id": "GCA_023229605_000000000176_34",
                "title": "GlobDB Cluster_00050",
                "subtitle": "GlobDB representative",
                "evidence": "FeGenie: Cyc2_repCluster2",
            },

            {
                "panel_id": "D",
                "protein_id": "GCA_002843155_000000000119_52",
                "title": "GlobDB Cluster_00206",
                "subtitle": "GlobDB representative",
                "evidence": "FeGenie: Cyc2_repCluster2",
            },
        ],
    },


    "cluster00121": {

        "ncols": 3,

        "align_to_first": 1,

        "panels": [

            {
                "panel_id": "A",
                "protein_id": "BCRBG_22465_000000000369_1",
                "title": "Cluster_00121",
                "subtitle": "GlobDB representative",
                "evidence": "Foldseek: IPT/TIG-like; MtrC/MtrF-like hit",
            },

            {
                "panel_id": "B",
                "protein_id": "contig_25_1899",
                "title": "Cluster_00121",
                "subtitle": "Methylumidiphilus",
                "evidence": "Foldseek: IPT/TIG-like; MtrC/MtrF-like hit",
            },

            {
                "panel_id": "C",
                "protein_id": "contig_1431_201",
                "title": "Cluster_00121",
                "subtitle": "Methylumidiphilus",
                "evidence": "Foldseek: IPT/TIG-like; MtrC/MtrF-like hit",
            },

            {
                "panel_id": "D",
                "protein_id": "contig_122_3517",
                "title": "Cluster_00121",
                "subtitle": "Methylovulum",
                "evidence": "Foldseek: cell-surface IPT/TIG-like; MtrC/MtrF-like hit",
            },

            {
                "panel_id": "E",
                "protein_id": "contig_871_20",
                "title": "Cluster_00121",
                "subtitle": "Methylomonas",
                "evidence": "Foldseek: IPT/TIG-like; MtrC/MtrF-like hit",
            },
        ],
    },


    "methylobacter": {

        "ncols": 3,

        ## These three proteins are not one homologous protein family.
        ## Do not force them into one structural alignment.
        "align_to_first": 0,

        "panels": [

            {
                "panel_id": "A",
                "protein_id": "contig_75_3433",
                "title": "Cluster_00042",
                "subtitle": "Methylobacter_A",
                "evidence": "Foldseek PDB100: penta-heme cytochrome c552",
            },

            {
                "panel_id": "B",
                "protein_id": "contig_483_99",
                "title": "Cluster_00288",
                "subtitle": "Methylobacter_C sp002256465",
                "evidence": "Foldseek: TPR protein; PDB100: OcwA",
            },

            {
                "panel_id": "C",
                "protein_id": "contig_171_1",
                "title": "No GlobDB cluster",
                "subtitle": "Methylobacter_C sp002256465",
                "evidence": "Foldseek: TPR protein; PDB100: OcwA",
            },
        ],
    },
}


def fail(message):

    print(
        f"\nERROR: {message}",
        file=sys.stderr,
    )

    sys.exit(1)


## ================================================================== ##
## Load Stage-76 authoritative selections
## ================================================================== ##

if not SELECTED_TSV.is_file():

    fail(
        f"Stage-76 selected-model table not found:\n"
        f"{SELECTED_TSV}"
    )


selected = pd.read_csv(
    SELECTED_TSV,
    sep="\t",
    dtype=str,
    keep_default_na=False,
)


required = [

    "panel_set",
    "protein_id",
    "job",
    "model_number",
    "ranking_score",
    "sequence_length",
    "expected_hemes",
    "selected_holo_cif",
]


missing = [

    column

    for column
    in required

    if column not in selected.columns
]


if missing:

    fail(
        "Stage-76 table is missing column(s):\n"
        +
        "\n".join(
            missing
        )
    )


if len(selected) != 12:

    fail(
        f"Expected 12 selected mature structures; "
        f"found {len(selected)}."
    )


## ================================================================== ##
## Build one manifest per figure
## ================================================================== ##

for figure_name, config in FIGURES.items():

    work_dir = (
        FIGROOT
        / "02_work"
        / figure_name
    )


    manifest_dir = (
        work_dir
        / "01_manifest"
    )


    manifest_dir.mkdir(
        parents=True,
        exist_ok=True,
    )


    output_rows = []


    for panel_order, panel in enumerate(
        config["panels"],
        start=1,
    ):

        hits = selected[
            (
                selected["panel_set"]
                ==
                figure_name
            )
            &
            (
                selected["protein_id"]
                ==
                panel["protein_id"]
            )
        ].copy()


        if len(hits) != 1:

            fail(
                f"{figure_name} / {panel['protein_id']}: "
                f"expected one Stage-76 structure; "
                f"found {len(hits)}."
            )


        row = hits.iloc[0]


        cif = Path(
            row["selected_holo_cif"]
        )


        if not cif.is_file():

            fail(
                f"Selected mature CIF missing:\n"
                f"{cif}"
            )


        n_hemes = int(
            row["expected_hemes"]
        )


        heme_word = (
            "heme"
            if n_hemes == 1
            else
            "hemes"
        )


        detail = (
            f"{panel['protein_id']} | "
            f"{row['sequence_length']} aa | "
            f"{n_hemes} {heme_word}"
        )


        output_rows.append(
            {
                "panel_order": panel_order,
                "panel_id": panel["panel_id"],
                "job": row["job"],
                "structure_path": str(cif),
                "model_number": row["model_number"],
                "ranking_score": row["ranking_score"],
                "title": panel["title"],
                "subtitle": panel["subtitle"],
                "detail": detail,
                "evidence": panel["evidence"],
                "ncols": config["ncols"],
                "align_to_first": config["align_to_first"],
            }
        )


    manifest = pd.DataFrame(
        output_rows
    )


    out = (
        manifest_dir
        / f"{figure_name}_manifest.tsv"
    )


    manifest.to_csv(
        out,
        sep="\t",
        index=False,
    )


    print("=" * 80)

    print(
        f"FINAL MATURE MANIFEST: {figure_name}"
    )

    print("=" * 80)

    print()

    print(
        manifest[
            [
                "panel_id",
                "title",
                "subtitle",
                "detail",
                "evidence",
                "model_number",
            ]
        ].to_string(
            index=False
        )
    )

    print()

    print(
        f"Manifest: {out}"
    )

    print()


print("=" * 80)

print(
    "STAGE 77 COMPLETE"
)

print("=" * 80)
