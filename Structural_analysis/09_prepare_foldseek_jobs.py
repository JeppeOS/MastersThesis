#!/usr/bin/env python3

import csv
import re
import shutil
import subprocess
from pathlib import Path


## ================================================================== ##
## PURPOSE
##
## Prepare selected AlphaFold3 best models for Foldseek.
##
## Input:
##
##   02_results/
##       selected_model_manifest.tsv
##       selected_models/<job>/best_model.cif
##
## Output:
##
##   03_foldseek/
##       protein_only_pdb/
##       foldseek_query_manifest.tsv
##       create_query_db.sh
##
## Each holo AlphaFold3 CIF is converted to a protein-only PDB:
##
##   - protein polymer retained
##   - HEC ligands removed
##   - Fe atoms removed
##   - water / ions / other heteroatoms removed
##
## The original holo CIFs remain untouched.
##
## These protein-only models can then be supplied directly to:
##
##   foldseek easy-search
##
## or converted to a reusable Foldseek database using createdb.
## ================================================================== ##


## ================================================================== ##
## Paths
## ================================================================== ##

ROOT = (
    Path.home()
    / "methanotrophs"
    / "methanotroph_project"
    / "jeppe"
    / "genome_analysis_workflow"
    / "18_structure_prediction"
    / "09_mature_heme_structure_figure"
)


RESULT_DIR = (
    ROOT
    / "02_results"
)


INPUT_MANIFEST = (
    RESULT_DIR
    / "selected_model_manifest.tsv"
)


OUTDIR = (
    ROOT
    / "03_foldseek"
)


PDB_DIR = (
    OUTDIR
    / "protein_only_pdb"
)


OUTPUT_MANIFEST = (
    OUTDIR
    / "foldseek_query_manifest.tsv"
)


CREATE_DB_SCRIPT = (
    OUTDIR
    / "create_query_db.sh"
)


PYMOL_HELPER = (
    OUTDIR
    / "extract_protein_only.py"
)


FOLDSEEK_DB = (
    OUTDIR
    / "query_db"
)


PDB_DIR.mkdir(
    parents=True,
    exist_ok=True
)


## ================================================================== ##
## Find PyMOL
## ================================================================== ##


PYMOL = shutil.which(
    "pymol"
)


if PYMOL is None:

    raise RuntimeError(
        "Could not find PyMOL in PATH.\n"
        "Activate the genome_viz environment first."
    )


print(
    "PyMOL:",
    PYMOL
)


## ================================================================== ##
## Helpers
## ================================================================== ##


def safe_id(text):
    """
    Make a stable Foldseek query identifier.
    """

    text = str(
        text
    ).strip()

    text = re.sub(
        r"[^A-Za-z0-9_.-]+",
        "_",
        text
    )

    text = re.sub(
        r"_+",
        "_",
        text
    )

    return text.strip(
        "_"
    )


def count_pdb_atoms(path):

    atom_count = 0
    hetatm_count = 0
    ca_count = 0
    residues = set()


    with open(
        path
    ) as fh:

        for line in fh:

            record = line[
                0:6
            ].strip()


            if record == "ATOM":

                atom_count += 1


                atom_name = line[
                    12:16
                ].strip()


                chain = line[
                    21:22
                ].strip()


                resid = line[
                    22:27
                ].strip()


                residues.add(
                    (
                        chain,
                        resid
                    )
                )


                if atom_name == "CA":

                    ca_count += 1


            elif record == "HETATM":

                hetatm_count += 1


    return {
        "atom_count": atom_count,
        "hetatm_count": hetatm_count,
        "ca_count": ca_count,
        "residue_count": len(
            residues
        ),
    }


## ================================================================== ##
## Read selected-model manifest
## ================================================================== ##


if not INPUT_MANIFEST.exists():

    raise FileNotFoundError(
        "Missing selected-model manifest:\n"
        f"{INPUT_MANIFEST}"
    )


with open(
    INPUT_MANIFEST,
    newline=""
) as fh:

    rows = list(
        csv.DictReader(
            fh,
            delimiter="\t"
        )
    )


if not rows:

    raise RuntimeError(
        "Selected-model manifest contains no rows."
    )


required = {
    "family",
    "protein_id",
    "job_name",
    "best_model",
}


missing = (
    required
    -
    set(
        rows[
            0
        ].keys()
    )
)


if missing:

    raise RuntimeError(
        "Selected-model manifest is missing columns:\n"
        + "\n".join(
            sorted(
                missing
            )
        )
    )


print(
    "Selected models:",
    len(
        rows
    )
)


## ================================================================== ##
## Define unique Foldseek query names
## ================================================================== ##


prepared = []


seen_ids = set()


for row in rows:

    family = safe_id(
        row[
            "family"
        ]
    )


    protein_id = safe_id(
        row[
            "protein_id"
        ]
    )


    query_id = (
        f"{family}__{protein_id}"
    )


    if query_id in seen_ids:

        raise RuntimeError(
            f"Duplicate Foldseek query ID: {query_id}"
        )


    seen_ids.add(
        query_id
    )


    source_cif = Path(
        row[
            "best_model"
        ]
    )


    if not source_cif.exists():

        raise FileNotFoundError(
            "Selected model missing:\n"
            f"{source_cif}"
        )


    output_pdb = (
        PDB_DIR
        / f"{query_id}.pdb"
    )


    prepared.append(
        {
            **row,
            "foldseek_query_id": query_id,
            "source_holo_cif": str(
                source_cif
            ),
            "protein_only_pdb": str(
                output_pdb
            ),
        }
    )


## ================================================================== ##
## Build one PyMOL helper script
##
## We perform all conversions in one PyMOL process.
## ================================================================== ##


helper_lines = [

    "from pymol import cmd",
    "",
    "cmd.reinitialize()",
    "",
]


for i, row in enumerate(
    prepared,
    start=1
):

    source = Path(
        row[
            "source_holo_cif"
        ]
    )


    dest = Path(
        row[
            "protein_only_pdb"
        ]
    )


    obj = (
        f"model_{i}"
    )


    selection = (
        f"protein_{i}"
    )


    helper_lines.extend(
        [
            (
                f'cmd.load(r"{source}", '
                f'"{obj}")'
            ),

            (
                f'cmd.select("{selection}", '
                f'"{obj} and polymer.protein")'
            ),

            (
                f'n = cmd.count_atoms("{selection}")'
            ),

            (
                f'assert n > 0, '
                f'"No protein atoms found in {source}"'
            ),

            (
                f'cmd.save(r"{dest}", '
                f'"{selection}")'
            ),

            (
                f'cmd.delete("{selection}")'
            ),

            (
                f'cmd.delete("{obj}")'
            ),

            "",
        ]
    )


helper_lines.extend(
    [
        'print("Protein-only Foldseek models written successfully.")',
        "cmd.quit()",
        "",
    ]
)


PYMOL_HELPER.write_text(
    "\n".join(
        helper_lines
    )
)


## ================================================================== ##
## Syntax-check helper
## ================================================================== ##


compile(
    PYMOL_HELPER.read_text(),
    str(
        PYMOL_HELPER
    ),
    "exec"
)


## ================================================================== ##
## Run PyMOL
## ================================================================== ##


print()
print(
    "Generating protein-only PDBs..."
)


command = [
    PYMOL,
    "-cq",
    str(
        PYMOL_HELPER
    ),
]


result = subprocess.run(
    command,
    text=True,
    capture_output=True
)


if result.stdout:

    print(
        result.stdout
    )


if result.stderr:

    print(
        result.stderr
    )


if result.returncode != 0:

    raise RuntimeError(
        "PyMOL protein extraction failed."
    )


## ================================================================== ##
## Validate generated PDB files
## ================================================================== ##


for row in prepared:

    pdb = Path(
        row[
            "protein_only_pdb"
        ]
    )


    if not pdb.exists():

        raise RuntimeError(
            "Expected Foldseek PDB was not generated:\n"
            f"{pdb}"
        )


    stats = count_pdb_atoms(
        pdb
    )


    if stats[
        "atom_count"
    ] == 0:

        raise RuntimeError(
            f"No protein ATOM records in:\n{pdb}"
        )


    if stats[
        "ca_count"
    ] == 0:

        raise RuntimeError(
            f"No C-alpha atoms found in:\n{pdb}"
        )


    if stats[
        "hetatm_count"
    ] != 0:

        raise RuntimeError(
            "Protein-only Foldseek input unexpectedly "
            f"contains HETATM records:\n{pdb}"
        )


    row[
        "protein_atom_count"
    ] = stats[
        "atom_count"
    ]


    row[
        "protein_residue_count"
    ] = stats[
        "residue_count"
    ]


    row[
        "ca_atom_count"
    ] = stats[
        "ca_count"
    ]


    row[
        "heteroatom_count"
    ] = stats[
        "hetatm_count"
    ]


## ================================================================== ##
## Write Foldseek query manifest
## ================================================================== ##


fields = [

    "family",

    "label",

    "annotation",

    "protein_id",

    "source_type",

    "job_name",

    "selected_model_index",

    "ptm",

    "expected_hec",

    "hec_records_detected",

    "foldseek_query_id",

    "source_holo_cif",

    "protein_only_pdb",

    "protein_residue_count",

    "protein_atom_count",

    "ca_atom_count",

    "heteroatom_count",
]


## Keep only columns that exist or are generated.

with open(
    OUTPUT_MANIFEST,
    "w",
    newline=""
) as fh:

    writer = csv.DictWriter(
        fh,
        fieldnames=fields,
        delimiter="\t",
        extrasaction="ignore"
    )


    writer.writeheader()


    for row in prepared:

        writer.writerow(
            row
        )


## ================================================================== ##
## Write Foldseek createdb helper
##
## This is deliberately generated but NOT executed.
## ================================================================== ##


CREATE_DB_SCRIPT.write_text(
    f'''#!/usr/bin/env bash
set -euo pipefail

ROOT="{ROOT}"
QUERY_DIR="{PDB_DIR}"
QUERY_DB="{FOLDSEEK_DB}"

echo "Building Foldseek query database"
echo "Input:  $QUERY_DIR"
echo "Output: $QUERY_DB"

foldseek createdb \\
    "$QUERY_DIR" \\
    "$QUERY_DB"

echo
echo "Foldseek query database created:"
echo "  $QUERY_DB"
'''
)


CREATE_DB_SCRIPT.chmod(
    0o755
)


## ================================================================== ##
## Final QC
## ================================================================== ##


print()
print(
    "=" * 130
)

print(
    "FOLDSEEK QUERY PREPARATION"
)

print(
    "=" * 130
)

print()


header = (
    f"{'Family':24s} "
    f"{'Protein':39s} "
    f"{'Residues':>8s} "
    f"{'Cα':>6s} "
    f"{'HETATM':>7s}"
)


print(
    header
)

print(
    "-" * len(
        header
    )
)


for row in sorted(
    prepared,
    key=lambda x: (
        x[
            "family"
        ],
        x[
            "protein_id"
        ],
    )
):

    print(
        f"{row['family']:24s} "
        f"{row['protein_id']:39s} "
        f"{row['protein_residue_count']:8d} "
        f"{row['ca_atom_count']:6d} "
        f"{row['heteroatom_count']:7d}"
    )


print()
print(
    "Prepared Foldseek structures:",
    len(
        prepared
    )
)


print()
print(
    "Protein-only PDB directory:"
)

print(
    " ",
    PDB_DIR
)


print()
print(
    "Foldseek query manifest:"
)

print(
    " ",
    OUTPUT_MANIFEST
)


print()
print(
    "Optional createdb helper:"
)

print(
    " ",
    CREATE_DB_SCRIPT
)


print()
print(
    "Original holo AlphaFold models were not modified."
)

print(
    "All Foldseek query PDBs contain protein atoms only."
)
PY
