#!/usr/bin/env python3

from pathlib import Path
import shutil
import subprocess
import pandas as pd


## ================================================================== ##
## PURPOSE
##
## Render all 15 standardized AlphaFold3 holo structures.
##
## IMPORTANT SEPARATION:
##
##   US-align:
##       protein-only PDBs
##
##   Final figures:
##       ORIGINAL mature holo AlphaFold3 CIFs
##
## Therefore Heme C ligands and Fe atoms are retained in the figures.
##
##
## Standardized rendering:
##
##   protein cartoon colored by pLDDT:
##
##       >90       #0053D6
##       70-90     #65CBF3
##       50-70     #FFDB13
##       <50       #FF7D45
##
##   HEC:
##       dark red sticks
##
##   Fe:
##       gold/orange spheres
##
## Within each figure set:
##
##   - panel A is the visual reference
##   - B/C/D/E are superposed onto A using protein coordinates
##   - all panels use exactly the same camera/view and scale
##
## PyMOL superposition is used ONLY to standardize orientation.
## Reported TM-score/RMSD values remain the independent US-align values.
## ================================================================== ##


BASE = (
    Path.home()
    / "methanotrophs"
    / "methanotroph_project"
    / "jeppe"
)


MANIFEST = (
    BASE
    / "structure_standardization"
    / "03_plot_manifest"
    / "final_structure_plot_manifest.tsv"
)


ROOT = (
    BASE
    / "structure_standardization"
    / "04_rendered"
)


HELPERS = (
    ROOT
    / "pymol_helpers"
)


PANELS = (
    ROOT
    / "panels"
)


LOGS = (
    ROOT
    / "logs"
)


TABLES = (
    ROOT
    / "tables"
)


for directory in (
    HELPERS,
    PANELS,
    LOGS,
    TABLES,
):

    directory.mkdir(
        parents=True,
        exist_ok=True
    )


QC_OUT = (
    TABLES
    / "rendering_qc.tsv"
)


## ================================================================== ##
## Input
## ================================================================== ##


df = pd.read_csv(
    MANIFEST,
    sep="\t"
)


if len(df) != 15:

    raise RuntimeError(
        f"Expected 15 structures; found {len(df)}."
    )


## ================================================================== ##
## PyMOL
## ================================================================== ##


PYMOL = shutil.which(
    "pymol"
)


if PYMOL is None:

    raise RuntimeError(
        "Could not find pymol in PATH. "
        "Activate genome_viz first."
    )


print()
print(
    f"PyMOL: {PYMOL}"
)


## ================================================================== ##
## Expected figure sizes
## ================================================================== ##


EXPECTED = {

    "cluster00035":
        2,

    "cyc2":
        4,

    "mca0421":
        4,

    "cluster00121":
        5,
}


## ================================================================== ##
## Build and run one PyMOL helper per figure set
## ================================================================== ##


qc_rows = []


for figure_set, group in df.groupby(
    "figure_set",
    sort=False
):

    group = group.sort_values(
        "panel"
    )


    if len(group) != EXPECTED[
        figure_set
    ]:

        raise RuntimeError(
            f"{figure_set}: expected "
            f"{EXPECTED[figure_set]} panels; "
            f"found {len(group)}."
        )


    # ------------------------------------------------------------- #
    # Panel A is the structural reference for each figure set.
    #
    # Stage 43 changes which biological structure is assigned panel A
    # for Cyc2 and Cluster_00121. Therefore no separate reference flag
    # is required.
    # ------------------------------------------------------------- #

    reference = group.loc[
        group["panel"]
        .astype(str)
        .str.strip()
        .eq("A")
    ]

    if len(reference) != 1:

        raise RuntimeError(
            f"{figure_set}: expected exactly "
            f"one reference structure."
        )


    if (
        reference.iloc[
            0
        ][
            "panel"
        ]
        != "A"
    ):

        raise RuntimeError(
            f"{figure_set}: panel A must be reference."
        )


    figure_panel_dir = (
        PANELS
        / figure_set
    )


    figure_panel_dir.mkdir(
        parents=True,
        exist_ok=True
    )


    helper = (
        HELPERS
        / f"{figure_set}.py"
    )


    log_file = (
        LOGS
        / f"{figure_set}.log"
    )


    records = []


    for _, row in group.iterrows():

        cif = Path(
            row[
                "render_structure"
            ]
        )


        if not cif.is_file():

            raise RuntimeError(
                f"Missing holo CIF:\n{cif}"
            )


        records.append(
            {
                "panel":
                    str(
                        row[
                            "panel"
                        ]
                    ),

                "cif":
                    str(
                        cif
                    ),

                "expected_hemes":
                    int(
                        row[
                            "expected_hemes"
                        ]
                    ),

                "output_png":
                    str(
                        figure_panel_dir
                        / (
                            f"{row['panel']}.png"
                        )
                    ),
            }
        )


    ## ============================================================= ##
    ## PyMOL helper
    ##
    ## Note:
    ##
    ## PyMOL's parser in this installation did not accept >= in the
    ## previous workflow. We therefore use strict inequalities.
    ## ============================================================= ##


    helper_code = f'''
from pymol import cmd


FIGURE_SET = {figure_set!r}

PANELS = {records!r}

QC_FILE = {str(TABLES / (figure_set + "_rendering_qc.tsv"))!r}


## ================================================================ ##
## Global rendering configuration
## ================================================================ ##


cmd.reinitialize()

cmd.bg_color(
    "white"
)

cmd.set(
    "ray_opaque_background",
    1
)

cmd.set(
    "orthoscopic",
    1
)

cmd.set(
    "depth_cue",
    0
)

cmd.set(
    "ray_shadows",
    0
)

cmd.set(
    "antialias",
    2
)

cmd.set(
    "cartoon_fancy_helices",
    1
)

cmd.set(
    "cartoon_smooth_loops",
    1
)

cmd.set(
    "stick_radius",
    0.20
)

cmd.set(
    "sphere_scale",
    0.32
)

cmd.set(
    "specular",
    0.20
)


## ================================================================ ##
## Canonical AlphaFold pLDDT colors
## ================================================================ ##


cmd.set_color(
    "plddt_vhigh",
    [
        0.000000,
        0.325490,
        0.839216
    ]
)


cmd.set_color(
    "plddt_conf",
    [
        0.396078,
        0.796078,
        0.952941
    ]
)


cmd.set_color(
    "plddt_low",
    [
        1.000000,
        0.858824,
        0.074510
    ]
)


cmd.set_color(
    "plddt_vlow",
    [
        1.000000,
        0.490196,
        0.270588
    ]
)


## ================================================================ ##
## Heme styling
## ================================================================ ##


cmd.set_color(
    "heme_red",
    [
        0.720,
        0.100,
        0.100
    ]
)


cmd.set_color(
    "fe_gold",
    [
        1.000,
        0.670,
        0.000
    ]
)


## ================================================================ ##
## Load and style structures
## ================================================================ ##


qc = []


for panel in PANELS:

    panel_id = panel[
        "panel"
    ]


    obj = (
        "obj_"
        + panel_id
    )


    cmd.load(
        panel[
            "cif"
        ],
        obj
    )


    cmd.hide(
        "everything",
        obj
    )


    protein = (
        f"{{obj}} and polymer.protein"
    )


    protein_atoms = cmd.count_atoms(
        protein
    )


    ca_atoms = cmd.count_atoms(
        f"({{protein}}) and name CA"
    )


    if (
        protein_atoms == 0
        or ca_atoms == 0
    ):

        raise RuntimeError(
            f"{{FIGURE_SET}} {{panel_id}}: "
            "no protein atoms detected."
        )


    cmd.show(
        "cartoon",
        protein
    )


    ## ------------------------------------------------------------ ##
    ## pLDDT coloring
    ##
    ## AlphaFold3 pLDDT is stored in the B-factor field.
    ##
    ## IMPORTANT:
    ## Avoid >= because this PyMOL parser previously split it into
    ## '> =' and rejected the selection.
    ## ------------------------------------------------------------ ##


    sel_vlow = (
        f"({{protein}}) and b < 50"
    )


    sel_low = (
        f"({{protein}}) "
        "and b > 49.999 "
        "and b < 70"
    )


    sel_conf = (
        f"({{protein}}) "
        "and b > 69.999 "
        "and b < 90"
    )


    sel_vhigh = (
        f"({{protein}}) "
        "and b > 89.999"
    )


    cmd.color(
        "plddt_vlow",
        sel_vlow
    )


    cmd.color(
        "plddt_low",
        sel_low
    )


    cmd.color(
        "plddt_conf",
        sel_conf
    )


    cmd.color(
        "plddt_vhigh",
        sel_vhigh
    )


    n_vlow = cmd.count_atoms(
        sel_vlow
    )


    n_low = cmd.count_atoms(
        sel_low
    )


    n_conf = cmd.count_atoms(
        sel_conf
    )


    n_vhigh = cmd.count_atoms(
        sel_vhigh
    )


    colored_atoms = (
        n_vlow
        + n_low
        + n_conf
        + n_vhigh
    )


    if colored_atoms != protein_atoms:

        raise RuntimeError(
            f"{{FIGURE_SET}} {{panel_id}}: "
            f"pLDDT class atoms={{colored_atoms}}, "
            f"protein atoms={{protein_atoms}}."
        )


    ## ------------------------------------------------------------ ##
    ## Heme C
    ## ------------------------------------------------------------ ##


    heme = (
        f"{{obj}} and resn HEC"
    )


    heme_non_fe = (
        f"{{heme}} and not name FE"
    )


    heme_fe = (
        f"{{heme}} and name FE"
    )


    n_heme_fe = cmd.count_atoms(
        heme_fe
    )


    expected_hemes = int(
        panel[
            "expected_hemes"
        ]
    )


    if n_heme_fe != expected_hemes:

        raise RuntimeError(
            f"{{FIGURE_SET}} {{panel_id}}: "
            f"expected {{expected_hemes}} HEC Fe atoms, "
            f"found {{n_heme_fe}}."
        )


    cmd.show(
        "sticks",
        heme_non_fe
    )


    cmd.color(
        "heme_red",
        heme_non_fe
    )


    cmd.show(
        "spheres",
        heme_fe
    )


    cmd.color(
        "fe_gold",
        heme_fe
    )


    qc.append(
        [
            FIGURE_SET,
            panel_id,
            ca_atoms,
            expected_hemes,
            n_heme_fe,
            n_vlow,
            n_low,
            n_conf,
            n_vhigh,
        ]
    )


## ================================================================ ##
## Structural superposition for VISUAL ORIENTATION ONLY
##
## The resulting PyMOL RMSD values are intentionally ignored.
##
## The reported quantitative comparisons are the independent
## Stage-36 US-align results.
## ================================================================ ##


reference_obj = (
    "obj_A"
)


for panel in PANELS:

    panel_id = panel[
        "panel"
    ]


    if panel_id == "A":

        continue


    mobile_obj = (
        "obj_"
        + panel_id
    )


    cmd.super(

        f"{{mobile_obj}} and polymer.protein",

        f"{{reference_obj}} and polymer.protein",

        cycles=5,

        transform=1
    )


## ================================================================ ##
## Shared camera/view
##
## All aligned structures are used to determine one common framing.
## Individual panels are then rendered using that identical view.
## ================================================================ ##


cmd.enable(
    "all"
)


cmd.orient(
    "all"
)


cmd.zoom(
    "all",
    buffer=4.0
)


shared_view = cmd.get_view()


## ================================================================ ##
## Render panels
## ================================================================ ##


for panel in PANELS:

    panel_id = panel[
        "panel"
    ]


    obj = (
        "obj_"
        + panel_id
    )


    cmd.disable(
        "all"
    )


    cmd.enable(
        obj
    )


    cmd.set_view(
        shared_view
    )


    cmd.ray(
        1600,
        1250
    )


    cmd.png(

        panel[
            "output_png"
        ],

        dpi=300
    )


    print(
        "RENDERED",
        FIGURE_SET,
        panel_id,
        panel[
            "output_png"
        ]
    )


## ================================================================ ##
## Write QC directly to TSV
##
## This is inside a generated Python helper.
## Backslashes therefore need to survive the outer Stage-38 template.
## ================================================================ ##


with open(
    QC_FILE,
    "w"
) as handle:

    handle.write(
        "figure_set\\t"
        "panel\\t"
        "CA_atoms\\t"
        "expected_hemes\\t"
        "HEC_Fe_atoms\\t"
        "pLDDT_lt50_atoms\\t"
        "pLDDT_50_70_atoms\\t"
        "pLDDT_70_90_atoms\\t"
        "pLDDT_gt90_atoms\\n"
    )


    for row in qc:

        handle.write(
            "\\t".join(
                str(x)
                for x in row
            )
            + "\\n"
        )


print(
    "QC written:",
    QC_FILE,
    flush=True
)


cmd.quit()
'''


    helper.write_text(
        helper_code
    )


    ## ------------------------------------------------------------- ##
    ## Validate generated helper syntax before invoking PyMOL.
    ## ------------------------------------------------------------- ##

    try:

        compile(
            helper_code,
            str(
                helper
            ),
            "exec"
        )

    except SyntaxError as exc:

        raise RuntimeError(
            f"Generated PyMOL helper has invalid syntax:\n"
            f"{helper}\n\n"
            f"{exc}"
        ) from exc


    print()

    print(
        "=" * 100
    )

    print(
        f"RENDERING {figure_set}"
    )

    print(
        "=" * 100
    )


    result = subprocess.run(

        [
            str(
                PYMOL
            ),

            "-cq",

            "-r",

            str(
                helper
            ),
        ],

        stdout=subprocess.PIPE,

        stderr=subprocess.STDOUT,

        text=True,
    )


    log_file.write_text(
        result.stdout
    )


    print(
        result.stdout
    )


    if result.returncode != 0:

        raise RuntimeError(
            f"PyMOL failed for {figure_set}.\n"
            f"See:\n{log_file}"
        )


    ## ============================================================= ##
    ## Read QC written directly by the PyMOL helper
    ## ============================================================= ##


    figure_qc_file = (
        TABLES
        / (
            figure_set
            + "_rendering_qc.tsv"
        )
    )


    if (
        not figure_qc_file.is_file()
        or figure_qc_file.stat().st_size == 0
    ):

        raise RuntimeError(
            f"{figure_set}: PyMOL did not produce its QC table.\n"
            f"Check:\n{log_file}"
        )


    figure_qc = pd.read_csv(
        figure_qc_file,
        sep="\t"
    )


    if len(
        figure_qc
    ) != EXPECTED[
        figure_set
    ]:

        raise RuntimeError(
            f"{figure_set}: expected "
            f"{EXPECTED[figure_set]} QC rows; "
            f"found {len(figure_qc)}."
        )


    qc_rows.extend(
        figure_qc.to_dict(
            "records"
        )
    )


## ================================================================== ##
## Validate outputs
## ================================================================== ##


qc = pd.DataFrame(
    qc_rows
)


if len(
    qc
) != 15:

    raise RuntimeError(
        f"Expected 15 QC rows; "
        f"found {len(qc)}."
    )


for _, row in df.iterrows():

    png = (
        PANELS
        / row[
            "figure_set"
        ]
        / (
            f"{row['panel']}.png"
        )
    )


    if (
        not png.is_file()
        or png.stat().st_size == 0
    ):

        raise RuntimeError(
            f"Missing rendered panel:\n{png}"
        )


qc.to_csv(
    QC_OUT,
    sep="\t",
    index=False
)


## ================================================================== ##
## Summary
## ================================================================== ##


print()

print(
    "=" * 100
)

print(
    "STAGE 38 COMPLETE"
)

print(
    "=" * 100
)

print()

print(
    "Rendered panels:"
)

print(
    PANELS
)

print()

print(
    "Rendering QC:"
)

print(
    QC_OUT
)

print()

print(
    "All final panels were rendered from the original holo CIFs."
)

print(
    "Heme C ligands were retained and validated against expected counts."
)

print(
    "Protein cartoons use the standardized AlphaFold pLDDT palette."
)

print(
    "US-align protein-only PDBs were NOT used for rendering."
)
