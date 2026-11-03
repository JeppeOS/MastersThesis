#!/usr/bin/env python3

from pathlib import Path
from itertools import combinations
import os
import re
import shutil
import subprocess
import pandas as pd


## ================================================================== ##
## PURPOSE
##
## Standardized quantitative structure comparison for the four
## homologous AlphaFold3 figure sets:
##
##   - Cluster_00035
##   - Cyc2
##   - MCA0421 / Cluster_00063
##   - Cluster_00121
##
## Workflow:
##
##   selected mature holo CIF
##           ↓
##   PyMOL: retain polymer.protein only
##           ↓
##   protein-only PDB
##           ↓
##   US-align -a T
##           ↓
##   average-length-normalized TM-score + RMSD
##
## We calculate:
##
##   1. every pair within each figure set
##   2. panel A vs all other panels
##
## The latter is what will be displayed in the standardized figures.
##
## Original source CIF files are NEVER modified.
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
    / "01_manifest"
    / "standard_structure_manifest.tsv"
)


ROOT = (
    BASE
    / "structure_standardization"
    / "02_structural_comparison"
)


PDBDIR = (
    ROOT
    / "protein_only_pdb"
)


RAWDIR = (
    ROOT
    / "usalign_raw"
)


TABLEDIR = (
    ROOT
    / "tables"
)


PMLDIR = (
    ROOT
    / "pymol_scripts"
)


for directory in (
    PDBDIR,
    RAWDIR,
    TABLEDIR,
    PMLDIR,
):

    directory.mkdir(
        parents=True,
        exist_ok=True
    )


ALL_PAIRS_OUT = (
    TABLEDIR
    / "USalign_all_pairwise.tsv"
)


REFERENCE_OUT = (
    TABLEDIR
    / "USalign_reference_vs_A.tsv"
)


TM_MATRIX_OUT = (
    TABLEDIR
    / "USalign_TMscore_matrix.tsv"
)


RMSD_MATRIX_OUT = (
    TABLEDIR
    / "USalign_RMSD_matrix.tsv"
)


FINAL_MANIFEST_OUT = (
    TABLEDIR
    / "standard_structure_manifest_with_USalign.tsv"
)


## ================================================================== ##
## Read manifest
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
## Locate PyMOL
## ================================================================== ##


PYMOL = shutil.which(
    "pymol"
)


if PYMOL is None:

    raise RuntimeError(
        "Could not find 'pymol' in PATH.\n"
        "Activate the genome_viz environment and rerun."
    )


## ================================================================== ##
## Locate US-align
##
## Try common executable names and installation locations.
## ================================================================== ##


def locate_usalign():

    executable_names = [
        "USalign",
        "US-align",
        "usalign",
    ]


    for name in executable_names:

        hit = shutil.which(
            name
        )

        if hit:

            return Path(
                hit
            )


    candidates = [

        Path.home()
        / "bin"
        / "USalign",

        Path.home()
        / "software"
        / "USalign",

        Path.home()
        / "software"
        / "US-align"
        / "USalign",

        Path.home()
        / "USalign",
    ]


    for path in candidates:

        if (
            path.is_file()
            and os.access(
                path,
                os.X_OK
            )
        ):

            return path


    return None


USALIGN = locate_usalign()


## ================================================================== ##
## Filename-safe identifier
## ================================================================== ##


def safe_name(
    value
):

    return re.sub(
        r"[^A-Za-z0-9_.-]+",
        "_",
        str(
            value
        )
    ).strip(
        "_"
    )


## ================================================================== ##
## Convert holo CIF → protein-only PDB
##
## polymer.protein excludes:
##
##   - HEC ligands
##   - ions
##   - water
##   - other non-protein atoms
##
## These protein-only files are used ONLY for quantitative alignment.
## The holo CIF remains the structure used for final rendering.
## ================================================================== ##


def make_protein_only_pdb(
    source_cif,
    output_pdb,
    helper_file
):

    source_cif = Path(
        source_cif
    )


    output_pdb = Path(
        output_pdb
    )


    helper_file = Path(
        helper_file
    )


    if not source_cif.is_file():

        raise RuntimeError(
            f"Missing source CIF:\n"
            f"{source_cif}"
        )


    ## -------------------------------------------------------------- ##
    ## Build a small Python helper executed inside PyMOL.
    ##
    ## We deliberately use several possible selections because
    ## AlphaFold3 mmCIF handling differs somewhat between PyMOL builds.
    ## -------------------------------------------------------------- ##

    helper_code = f"""
from pymol import cmd


cmd.reinitialize()


source_cif = {str(source_cif)!r}

output_pdb = {str(output_pdb)!r}


cmd.load(
    source_cif,
    "afmodel"
)


print(
    "Loaded atoms:",
    cmd.count_atoms(
        "afmodel"
    )
)


candidate_selections = [

    (
        "polymer.protein",
        "afmodel and polymer.protein"
    ),

    (
        "polymer",
        "afmodel and polymer"
    ),

    (
        "non-HETATM fallback",
        (
            "afmodel "
            "and not hetatm "
            "and not solvent "
            "and not resn HEC+HEM+HEA+HEB"
        )
    ),
]


selected_name = None

selected_expression = None


for name, expression in candidate_selections:

    n_atoms = cmd.count_atoms(
        expression
    )


    n_ca = cmd.count_atoms(
        f"({{expression}}) and name CA"
    )


    print(
        f"Selection {{name}}: "
        f"{{n_atoms}} atoms, "
        f"{{n_ca}} CA atoms"
    )


    if (
        n_atoms > 0
        and n_ca > 0
    ):

        selected_name = name

        selected_expression = expression

        break


if selected_expression is None:

    raise RuntimeError(
        "No usable protein selection was found "
        "in the AlphaFold structure."
    )


print(
    "Using protein selection:",
    selected_name
)


## -------------------------------------------------------------- ##
## Create independent protein-only object.
## -------------------------------------------------------------- ##

cmd.create(
    "protein_model",
    selected_expression
)


## -------------------------------------------------------------- ##
## Defensive removal of non-protein material.
##
## Again: this only affects the temporary alignment copy,
## NOT the original holo CIF.
## -------------------------------------------------------------- ##

cmd.remove(
    "protein_model and solvent"
)


cmd.remove(
    "protein_model and hetatm"
)


cmd.remove(
    "protein_model and resn HEC+HEM+HEA+HEB"
)


final_atoms = cmd.count_atoms(
    "protein_model"
)


final_ca = cmd.count_atoms(
    "protein_model and name CA"
)


print(
    "Final protein-only object:",
    final_atoms,
    "atoms;",
    final_ca,
    "CA atoms"
)


if (
    final_atoms == 0
    or final_ca == 0
):

    raise RuntimeError(
        "Protein-only object became empty."
    )


cmd.save(
    output_pdb,
    "protein_model",
    state=1
)


print(
    "Saved:",
    output_pdb
)


cmd.quit()
"""


    helper_file.write_text(
        helper_code
    )


    result = subprocess.run(

        [
            str(
                PYMOL
            ),

            "-cq",

            "-r",

            str(
                helper_file
            ),
        ],

        stdout=subprocess.PIPE,

        stderr=subprocess.STDOUT,

        text=True,
    )


    ## -------------------------------------------------------------- ##
    ## Retain extraction log for QC.
    ## -------------------------------------------------------------- ##

    log_file = helper_file.with_suffix(
        ".log"
    )


    log_file.write_text(
        result.stdout
    )


    if result.returncode != 0:

        raise RuntimeError(
            "PyMOL failed while creating "
            "protein-only PDB:\n"
            f"{source_cif}\n\n"
            f"{result.stdout}"
        )


    if (
        not output_pdb.is_file()
        or output_pdb.stat().st_size == 0
    ):

        raise RuntimeError(
            "Protein-only PDB was not created:\n"
            f"{output_pdb}\n\n"
            "PyMOL output:\n"
            f"{result.stdout}"
        )



## ================================================================== ##
## Count residues in protein-only PDB
##
## Count unique CA-containing residues.
## ================================================================== ##


def pdb_residue_count(
    pdb
):

    residues = set()


    with open(
        pdb
    ) as handle:

        for line in handle:

            if not line.startswith(
                "ATOM"
            ):

                continue


            atom = line[
                12:16
            ].strip()


            if atom != "CA":

                continue


            chain = line[
                21:22
            ]


            resno = line[
                22:26
            ].strip()


            icode = line[
                26:27
            ]


            residues.add(
                (
                    chain,
                    resno,
                    icode
                )
            )


    return len(
        residues
    )


## ================================================================== ##
## Build protein-only PDB collection
## ================================================================== ##


protein_paths = []


print()

print(
    "=" * 100
)

print(
    "CREATING PROTEIN-ONLY STRUCTURES"
)

print(
    "=" * 100
)


for _, row in df.iterrows():

    figure = row[
        "figure_set"
    ]


    panel = row[
        "panel"
    ]


    prefix = (
        f"{figure}_{panel}"
    )


    pdb = (
        PDBDIR
        / f"{prefix}.pdb"
    )


    pml = (
        PMLDIR
        / f"{prefix}.py"
    )


    make_protein_only_pdb(

        row[
            "structure_cif"
        ],

        pdb,

        pml,
    )


    n_res = pdb_residue_count(
        pdb
    )


    if n_res == 0:

        raise RuntimeError(
            f"No protein residues found in:\n"
            f"{pdb}"
        )


    protein_paths.append(
        {
            "figure_set":
                figure,

            "panel":
                panel,

            "protein_pdb":
                str(
                    pdb
                ),

            "protein_residues":
                n_res,
        }
    )


    print(
        f"{figure:15s} "
        f"{panel} "
        f"{n_res:4d} residues "
        f"→ {pdb.name}"
    )


protein_df = pd.DataFrame(
    protein_paths
)


df = df.merge(

    protein_df,

    on=[
        "figure_set",
        "panel"
    ],

    how="left",

    validate="one_to_one"
)


## ================================================================== ##
## If US-align is missing:
##
## protein-only PDB preparation is still useful, so stop only here.
## ================================================================== ##


if USALIGN is None:

    prep_out = (
        TABLEDIR
        / "standard_structure_manifest_with_protein_PDB.tsv"
    )


    df.to_csv(
        prep_out,
        sep="\t",
        index=False
    )


    print()

    print(
        "=" * 100
    )

    print(
        "US-align NOT FOUND"
    )

    print(
        "=" * 100
    )

    print()

    print(
        "Protein-only PDB preparation completed successfully."
    )

    print()

    print(
        "Prepared manifest:"
    )

    print(
        prep_out
    )

    print()

    print(
        "No US-align comparisons were performed."
    )

    print()

    print(
        "Checked executable names:"
    )

    print(
        "  USalign"
    )

    print(
        "  US-align"
    )

    print(
        "  usalign"
    )

    print()

    print(
        "Also checked:"
    )

    print(
        "  ~/bin/USalign"
    )

    print(
        "  ~/software/USalign"
    )

    print(
        "  ~/software/US-align/USalign"
    )

    raise SystemExit(
        2
    )


print()

print(
    f"US-align executable: {USALIGN}"
)


## ================================================================== ##
## Parse US-align output
##
## Run with:
##
##     -a T
##
## so TM-score is normalized by the average length of both structures.
##
## We explicitly extract the average-length-normalized TM-score.
## ================================================================== ##


def parse_usalign_output(
    text
):

    ## -------------------------------------------------------------- ##
    ## Aligned length / RMSD / sequence identity
    ##
    ## Example:
    ##
    ## Aligned length= 281, RMSD= 2.08,
    ## Seq_ID=n_identical/n_aligned= 0.947
    ## -------------------------------------------------------------- ##

    alignment_match = re.search(

        r"Aligned length=\s*"
        r"([0-9]+)"
        r",\s*RMSD=\s*"
        r"([0-9.]+)"
        r",\s*Seq_ID=.*?=\s*"
        r"([0-9.]+)",

        text
    )


    if alignment_match is None:

        raise RuntimeError(
            "Could not parse aligned length / RMSD / "
            "sequence identity from US-align output:\n\n"
            + text
        )


    aligned_length = int(
        alignment_match.group(
            1
        )
    )


    rmsd = float(
        alignment_match.group(
            2
        )
    )


    seq_id = float(
        alignment_match.group(
            3
        )
    )


    ## -------------------------------------------------------------- ##
    ## Symmetric TM-score
    ##
    ## US-align 20241108 reports:
    ##
    ## TM-score= 0.xxxxx
    ## (if normalized by average length of two structures:
    ##  L=..., d0=...)
    ##
    ## This is the score we want because it is symmetric with respect
    ## to the two compared proteins.
    ## -------------------------------------------------------------- ##

    average_match = re.search(

        r"TM-score=\s*"
        r"([0-9.]+)"
        r"\s*"
        r"\("
        r"if normalized by average length "
        r"of two structures:",

        text,

        flags=re.IGNORECASE
    )


    if average_match is not None:

        tm_average = float(
            average_match.group(
                1
            )
        )

    else:

        ## ---------------------------------------------------------- ##
        ## More permissive fallback for other US-align versions.
        ## ---------------------------------------------------------- ##

        average_match = re.search(

            r"TM-score=\s*"
            r"([0-9.]+)"
            r".{0,100}?"
            r"average length",

            text,

            flags=(
                re.IGNORECASE
                |
                re.DOTALL
            )
        )


        if average_match is None:

            raise RuntimeError(
                "Could not identify the average-length-normalized "
                "TM-score from US-align output:\n\n"
                + text
            )


        tm_average = float(
            average_match.group(
                1
            )
        )


    ## -------------------------------------------------------------- ##
    ## Also retain both conventional directional TM-scores for QC.
    ##
    ## These are normalized separately by Structure_1 and Structure_2.
    ## -------------------------------------------------------------- ##

    tm_directional = re.findall(

        r"TM-score=\s*"
        r"([0-9.]+)"
        r"\s*"
        r"\("
        r"normalized by length of Structure_([12])",

        text,

        flags=re.IGNORECASE
    )


    tm_structure_1 = None
    tm_structure_2 = None


    for score, structure_no in tm_directional:

        if structure_no == "1":

            tm_structure_1 = float(
                score
            )

        elif structure_no == "2":

            tm_structure_2 = float(
                score
            )


    return {

        "aligned_length":
            aligned_length,

        "rmsd_A":
            rmsd,

        "seq_identity_aligned":
            seq_id,

        "tm_score":
            tm_average,

        "tm_score_norm_structure_1":
            tm_structure_1,

        "tm_score_norm_structure_2":
            tm_structure_2,
    }



## ================================================================== ##
## Run one US-align comparison
## ================================================================== ##


def run_usalign(
    row1,
    row2
):

    figure = row1[
        "figure_set"
    ]


    panel1 = row1[
        "panel"
    ]


    panel2 = row2[
        "panel"
    ]


    raw_file = (
        RAWDIR
        / (
            f"{figure}_"
            f"{panel1}_vs_{panel2}.txt"
        )
    )


    cmd = [

        str(
            USALIGN
        ),

        row1[
            "protein_pdb"
        ],

        row2[
            "protein_pdb"
        ],

        "-a",
        "T",
    ]


    result = subprocess.run(

        cmd,

        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )


    raw_file.write_text(
        result.stdout
    )


    if result.returncode != 0:

        raise RuntimeError(
            "US-align failed:\n"
            +
            " ".join(
                cmd
            )
            +
            "\n\n"
            +
            result.stdout
        )


    metrics = parse_usalign_output(
        result.stdout
    )


    return {

        "figure_set":
            figure,

        "panel_1":
            panel1,

        "panel_2":
            panel2,

        "protein_1":
            row1[
                "protein_id"
            ],

        "protein_2":
            row2[
                "protein_id"
            ],

        "length_1":
            int(
                row1[
                    "protein_residues"
                ]
            ),

        "length_2":
            int(
                row2[
                    "protein_residues"
                ]
            ),

        **metrics,

        "raw_output":
            str(
                raw_file
            ),
    }


## ================================================================== ##
## All-vs-all comparison within each figure set
## ================================================================== ##


print()

print(
    "=" * 100
)

print(
    "RUNNING US-align ALL-VS-ALL"
)

print(
    "=" * 100
)


comparison_rows = []


for figure, group in df.groupby(
    "figure_set",
    sort=False
):

    group = group.sort_values(
        "panel"
    )


    group_records = group.to_dict(
        "records"
    )


    print()

    print(
        figure
    )

    print(
        "-" * len(
            figure
        )
    )


    for row1, row2 in combinations(
        group_records,
        2
    ):

        result = run_usalign(
            row1,
            row2
        )


        comparison_rows.append(
            result
        )


        print(

            f"{result['panel_1']} vs "
            f"{result['panel_2']}  "

            f"TM={result['tm_score']:.4f}  "

            f"RMSD={result['rmsd_A']:.2f} Å  "

            f"aligned={result['aligned_length']}"
        )


comparisons = pd.DataFrame(
    comparison_rows
)


comparisons.to_csv(
    ALL_PAIRS_OUT,
    sep="\t",
    index=False
)


## ================================================================== ##
## Reference-panel comparisons
##
## Panel A is reference for every figure set.
## ================================================================== ##


reference_rows = []


for figure, group in df.groupby(
    "figure_set",
    sort=False
):

    reference = group[
        group[
            "reference_panel"
        ]
    ]


    if len(
        reference
    ) != 1:

        raise RuntimeError(
            f"{figure}: expected exactly "
            f"one reference panel."
        )


    ref_panel = str(
        reference.iloc[
            0
        ][
            "panel"
        ]
    )


    if ref_panel != "A":

        raise RuntimeError(
            f"{figure}: expected A as reference; "
            f"found {ref_panel}."
        )


    for _, row in group.iterrows():

        panel = str(
            row[
                "panel"
            ]
        )


        if panel == "A":

            reference_rows.append({

                "figure_set":
                    figure,

                "reference_panel":
                    "A",

                "comparison_panel":
                    "A",

                "tm_score_vs_A":
                    1.0,

                "rmsd_vs_A_A":
                    0.0,

                "aligned_length_vs_A":
                    int(
                        row[
                            "protein_residues"
                        ]
                    ),

                "seq_identity_vs_A":
                    1.0,
            })


            continue


        hit = comparisons[
            (
                comparisons[
                    "figure_set"
                ]
                == figure
            )
            &
            (
                (
                    (
                        comparisons[
                            "panel_1"
                        ]
                        == "A"
                    )
                    &
                    (
                        comparisons[
                            "panel_2"
                        ]
                        == panel
                    )
                )
                |
                (
                    (
                        comparisons[
                            "panel_1"
                        ]
                        == panel
                    )
                    &
                    (
                        comparisons[
                            "panel_2"
                        ]
                        == "A"
                    )
                )
            )
        ]


        if len(
            hit
        ) != 1:

            raise RuntimeError(
                f"{figure}: could not uniquely locate "
                f"A-vs-{panel} comparison."
            )


        metric = hit.iloc[
            0
        ]


        reference_rows.append({

            "figure_set":
                figure,

            "reference_panel":
                "A",

            "comparison_panel":
                panel,

            "tm_score_vs_A":
                metric[
                    "tm_score"
                ],

            "rmsd_vs_A_A":
                metric[
                    "rmsd_A"
                ],

            "aligned_length_vs_A":
                metric[
                    "aligned_length"
                ],

            "seq_identity_vs_A":
                metric[
                    "seq_identity_aligned"
                ],
        })


reference_df = pd.DataFrame(
    reference_rows
)


reference_df.to_csv(
    REFERENCE_OUT,
    sep="\t",
    index=False
)


## ================================================================== ##
## Join metrics back into standardized manifest
## ================================================================== ##


final_manifest = df.merge(

    reference_df,

    left_on=[
        "figure_set",
        "panel"
    ],

    right_on=[
        "figure_set",
        "comparison_panel"
    ],

    how="left",

    validate="one_to_one"
)


final_manifest = final_manifest.drop(
    columns=[
        "comparison_panel",
        "reference_panel_y"
    ]
).rename(
    columns={
        "reference_panel_x":
            "reference_panel"
    }
)


final_manifest.to_csv(
    FINAL_MANIFEST_OUT,
    sep="\t",
    index=False
)


## ================================================================== ##
## TM-score matrix
##
## One global matrix with labels figure:panel.
## ================================================================== ##


labels = [

    f"{row.figure_set}:{row.panel}"

    for row in df.itertuples(
        index=False
    )
]


tm_matrix = pd.DataFrame(

    float(
        "nan"
    ),

    index=labels,

    columns=labels
)


rmsd_matrix = pd.DataFrame(

    float(
        "nan"
    ),

    index=labels,

    columns=labels
)


for label in labels:

    tm_matrix.loc[
        label,
        label
    ] = 1.0


    rmsd_matrix.loc[
        label,
        label
    ] = 0.0


for row in comparisons.itertuples(
    index=False
):

    x = (
        f"{row.figure_set}:"
        f"{row.panel_1}"
    )


    y = (
        f"{row.figure_set}:"
        f"{row.panel_2}"
    )


    tm_matrix.loc[
        x,
        y
    ] = row.tm_score


    tm_matrix.loc[
        y,
        x
    ] = row.tm_score


    rmsd_matrix.loc[
        x,
        y
    ] = row.rmsd_A


    rmsd_matrix.loc[
        y,
        x
    ] = row.rmsd_A


tm_matrix.to_csv(
    TM_MATRIX_OUT,
    sep="\t"
)


rmsd_matrix.to_csv(
    RMSD_MATRIX_OUT,
    sep="\t"
)


## ================================================================== ##
## Final report
## ================================================================== ##


print()

print(
    "=" * 100
)

print(
    "REFERENCE COMPARISONS FOR FINAL FIGURES"
)

print(
    "=" * 100
)

print()


for figure, group in reference_df.groupby(
    "figure_set",
    sort=False
):

    print(
        figure
    )

    print(
        "-" * len(
            figure
        )
    )


    for _, row in group.iterrows():

        panel = row[
            "comparison_panel"
        ]


        if panel == "A":

            print(
                "A  reference"
            )

        else:

            print(

                f"{panel} vs A  "

                f"TM-score "
                f"{row['tm_score_vs_A']:.3f}  "

                f"RMSD "
                f"{row['rmsd_vs_A_A']:.2f} Å  "

                f"aligned "
                f"{int(row['aligned_length_vs_A'])}"
            )


    print()


print(
    "=" * 100
)

print(
    "OUTPUTS"
)

print(
    "=" * 100
)

print()

print(
    f"Protein-only PDBs:\n"
    f"{PDBDIR}"
)

print()

print(
    f"All pairwise comparisons:\n"
    f"{ALL_PAIRS_OUT}"
)

print()

print(
    f"Panel-A reference comparisons:\n"
    f"{REFERENCE_OUT}"
)

print()

print(
    f"TM-score matrix:\n"
    f"{TM_MATRIX_OUT}"
)

print()

print(
    f"RMSD matrix:\n"
    f"{RMSD_MATRIX_OUT}"
)

print()

print(
    f"Updated standardized manifest:\n"
    f"{FINAL_MANIFEST_OUT}"
)
