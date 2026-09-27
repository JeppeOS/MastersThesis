#!/usr/bin/env python3

from __future__ import annotations

import json
import math
import re
import sys
import zipfile
from functools import lru_cache
from pathlib import Path

import pandas as pd

try:
    from Bio.PDB import MMCIFParser, PDBIO, Select
    from Bio.SeqUtils import seq1
except ImportError:
    sys.exit(
        "\nERROR: Biopython is required.\n"
        "Run this from a conda environment containing Biopython.\n"
    )


## ================================================================== ##
## STAGE 71
##
## QC HEME-AWARE ALPHAFOLD MODELS FOR METHYLOBACTER MULTIHEME
## CANDIDATES
##
## Expected jobs:
##
##     Cluster_00042 candidate       -> 5 HEC
##     no-match candidate            -> 7 HEC
##     Cluster_00288/00218 candidate -> 8 HEC
##
## For every AlphaFold model:
##
##     - identify all canonical CXXCH motifs
##     - identify all HEC molecules
##     - assign motifs to HEC molecules one-to-one
##     - measure:
##
##           Cys1-S -> nearest HEC atom
##           Cys2-S -> nearest HEC atom
##           His-NE2 -> HEC Fe
##
## A site passes when all three distances are <= 3.0 Å.
##
## Model selection:
##
##     highest AlphaFold ranking score
##
## The selected holo CIF remains unchanged.
## A protein-only PDB copy is written for Foldseek.
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
    / "Methylobacter_multiheme"
    / "structure_analysis"
)


AF_ROOT = (
    ROOT
    / "alphafold"
)


RAW_DIR = (
    AF_ROOT
    / "raw_data"
)


EXTRACTED_DIR = (
    AF_ROOT
    / "extracted"
)


QC_DIR = (
    AF_ROOT
    / "qc"
)


FOLDSEEK_DIR = (
    ROOT
    / "foldseek"
    / "queries"
)


OUT_ALL = (
    QC_DIR
    / "Methylobacter_multiheme_AF_heme_geometry_all_models.tsv"
)


OUT_SELECTED = (
    QC_DIR
    / "Methylobacter_multiheme_AF_selected_models.tsv"
)


OUT_SITES = (
    QC_DIR
    / "Methylobacter_multiheme_AF_selected_model_heme_sites.tsv"
)


EXPECTED_JOBS = 3
EXPECTED_MODELS_PER_JOB = 5
EXPECTED_TOTAL_HEMES = 20

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


def safe_name(text):

    return re.sub(
        r"[^A-Za-z0-9_.-]+",
        "_",
        str(text),
    )


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


def read_summary_score(cif_path, model_number):

    parent = cif_path.parent

    name = cif_path.name


    prefix = re.sub(
        rf"_model_{model_number}\.cif$",
        "",
        name,
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


def motif_information(
    sequence,
    polymer,
):

    matches = list(
        re.finditer(
            r"C..CH",
            sequence,
        )
    )


    motifs = []


    for index, match in enumerate(
        matches,
        start=1,
    ):

        i = match.start()


        if (
            i + 4
            >=
            len(
                polymer
            )
        ):

            fail(
                "CXXCH motif exceeds polymer residue list."
            )


        cys1 = polymer[
            i
        ]

        cys2 = polymer[
            i + 3
        ]

        his = polymer[
            i + 4
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
                f"Sequence/residue mismatch at motif {index}."
            )


        motifs.append(
            {
                "motif_number":
                    index,

                "sequence_start_1based":
                    i + 1,

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


def optimal_assignment(cost_matrix):

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

            if (
                used_mask
                &
                (
                    1
                    <<
                    heme_index
                )
            ):

                continue


            remainder_cost, remainder_assignment = solve(
                motif_index + 1,
                used_mask
                |
                (
                    1
                    <<
                    heme_index
                ),
            )


            total = (
                cost_matrix[
                    motif_index
                ][
                    heme_index
                ]
                +
                remainder_cost
            )


            if total < best_cost:

                best_cost = total

                best_assignment = (
                    heme_index,
                ) + remainder_assignment


        return (
            best_cost,
            best_assignment,
        )


    return solve(
        0,
        0,
    )


class ProteinOnlySelect(
    Select
):

    def __init__(
        self,
        chain_id,
    ):

        self.chain_id = chain_id


    def accept_chain(
        self,
        chain,
    ):

        return int(
            chain.id
            ==
            self.chain_id
        )


    def accept_residue(
        self,
        residue,
    ):

        return int(
            residue.id[
                0
            ]
            ==
            " "
        )


## ================================================================== ##
## Directories
## ================================================================== ##

for directory in [
    EXTRACTED_DIR,
    QC_DIR,
    FOLDSEEK_DIR,
]:

    directory.mkdir(
        parents=True,
        exist_ok=True,
    )


## ================================================================== ##
## Extract archives
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
    ) as zf:

        zf.extractall(
            destination
        )


    print(
        f"  extracted: {archive.name}"
    )


## ================================================================== ##
## Discover models
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


    job = match.group(
        "job"
    )


    heme_match = re.search(
        r"_hemec(\d+)",
        job,
        flags=re.IGNORECASE,
    )


    if heme_match is None:

        fail(
            f"Could not infer expected HEC count from job:\n{job}"
        )


    expected_hemes = int(
        heme_match.group(
            1
        )
    )


    discovered.append(
        {
            "job":
                job,

            "model_number":
                int(
                    match.group(
                        "model"
                    )
                ),

            "expected_hemes":
                expected_hemes,

            "cif":
                cif,
        }
    )


jobs = sorted(
    set(
        item[
            "job"
        ]

        for item
        in discovered
    )
)


if len(
    jobs
) != EXPECTED_JOBS:

    fail(
        f"Expected {EXPECTED_JOBS} AlphaFold jobs; "
        f"found {len(jobs)}:\n"
        +
        "\n".join(
            jobs
        )
    )


expected_counts = sorted(
    set(
        item[
            "expected_hemes"
        ]

        for item
        in discovered
    )
)


if expected_counts != [
    5,
    7,
    8,
]:

    fail(
        f"Expected heme-count set [5, 7, 8]; "
        f"found {expected_counts}."
    )


print()
print(
    f"Discovered {len(jobs)} Methylobacter jobs."
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

    job_models = sorted(
        [
            item

            for item
            in discovered

            if item[
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
        job_models
    ) != EXPECTED_MODELS_PER_JOB:

        fail(
            f"{job}: expected "
            f"{EXPECTED_MODELS_PER_JOB} models; "
            f"found {len(job_models)}."
        )


    expected_hemes = job_models[
        0
    ][
        "expected_hemes"
    ]


    print()
    print(
        f"{job}  [expected HEC={expected_hemes}]"
    )


    for item in job_models:

        model_number = item[
            "model_number"
        ]

        cif = item[
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


        motifs = motif_information(
            sequence,
            polymer,
        )


        hemes = get_hemes(
            structure
        )


        if len(
            motifs
        ) != expected_hemes:

            fail(
                f"{job} model {model_number}: "
                f"found {len(motifs)} CXXCH motifs; "
                f"expected {expected_hemes}."
            )


        if len(
            hemes
        ) != expected_hemes:

            fail(
                f"{job} model {model_number}: "
                f"found {len(hemes)} HEC residues; "
                f"expected {expected_hemes}."
            )


        geometry_matrix = []
        full_geometry = []


        for motif in motifs:

            cost_row = []
            geometry_row = []


            for heme in hemes:

                geometry = pair_geometry(
                    motif,
                    heme,
                )


                if geometry is None:

                    fail(
                        f"{job} model {model_number}: "
                        "required SG/NE2/FE atoms missing."
                    )


                cost_row.append(
                    geometry[
                        "assignment_score"
                    ]
                )


                geometry_row.append(
                    geometry
                )


            geometry_matrix.append(
                cost_row
            )

            full_geometry.append(
                geometry_row
            )


        total_cost, assignment = optimal_assignment(
            geometry_matrix
        )


        n_pass = 0
        current_sites = []


        for motif_index, heme_index in enumerate(
            assignment
        ):

            motif = motifs[
                motif_index
            ]

            heme = hemes[
                heme_index
            ]

            geometry = full_geometry[
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


            current_sites.append(
                {
                    "job":
                        job,

                    "expected_hemes":
                        expected_hemes,

                    "model_number":
                        model_number,

                    "motif_number":
                        motif[
                            "motif_number"
                        ],

                    "motif_start_1based":
                        motif[
                            "sequence_start_1based"
                        ],

                    "heme_assignment":
                        heme_index + 1,

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
            )


        ranking_score, summary_file = read_summary_score(
            cif,
            model_number,
        )


        model_rows.append(
            {
                "job":
                    job,

                "expected_hemes":
                    expected_hemes,

                "model_number":
                    model_number,

                "protein_chain":
                    chain.id,

                "sequence_length":
                    len(
                        sequence
                    ),

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
                        x[
                            "cys1_heme_distance"
                        ]

                        for x
                        in current_sites
                    ),

                "max_cys2_heme_distance":
                    max(
                        x[
                            "cys2_heme_distance"
                        ]

                        for x
                        in current_sites
                    ),

                "max_his_fe_distance":
                    max(
                        x[
                            "his_fe_distance"
                        ]

                        for x
                        in current_sites
                    ),

                "cif":
                    str(
                        cif
                    ),

                "summary_file":
                    summary_file,
            }
        )


        site_rows.extend(
            current_sites
        )


        print(
            f"  model_{model_number}: "
            f"{n_pass}/{expected_hemes} sites <=3 Å"
            +
            (
                f", ranking={ranking_score:.4f}"
                if not math.isnan(
                    ranking_score
                )
                else
                ""
            )
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


## ================================================================== ##
## Select best model per protein
## ================================================================== ##

selected_rows = []
selected_site_frames = []


for job, group in models.groupby(
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
                "all_sites_pass",
                "sites_passing_3A",
                "model_number",
            ],
            ascending=[
                False,
                False,
                False,
                True,
            ],
            kind="stable",
        )


    else:

        group = group.sort_values(
            [
                "model_number",
            ],
            kind="stable",
        )


    best = group.iloc[
        0
    ].copy()


    expected_hemes = int(
        best[
            "expected_hemes"
        ]
    )


    if int(
        best[
            "all_sites_pass"
        ]
    ) != 1:

        fail(
            f"Best model for {job} does not have "
            f"{expected_hemes}/{expected_hemes} "
            "geometrically matched heme sites."
        )


    cif = Path(
        best[
            "cif"
        ]
    )


    structure = parser.get_structure(
        f"{job}_selected",
        str(
            cif
        ),
    )


    chain, polymer = get_protein_chain(
        structure
    )


    pdb_path = (
        FOLDSEEK_DIR
        /
        (
            safe_name(
                job
            )
            +
            ".pdb"
        )
    )


    io = PDBIO()

    io.set_structure(
        structure
    )


    io.save(
        str(
            pdb_path
        ),
        ProteinOnlySelect(
            chain.id
        ),
    )


    best[
        "selected_for_foldseek"
    ] = 1

    best[
        "protein_only_pdb"
    ] = str(
        pdb_path
    )


    selected_rows.append(
        best
    )


    selected_site_frames.append(
        sites[
            (
                sites[
                    "job"
                ]
                ==
                job
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
    )


selected = pd.DataFrame(
    selected_rows
)


selected_sites = pd.concat(
    selected_site_frames,
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
    selected
) != EXPECTED_JOBS:

    fail(
        "Did not select exactly three structures."
    )


if int(
    selected[
        "expected_hemes"
    ].sum()
) != EXPECTED_TOTAL_HEMES:

    fail(
        "Selected structures do not sum to 20 expected hemes."
    )


if int(
    selected[
        "sites_passing_3A"
    ].sum()
) != EXPECTED_TOTAL_HEMES:

    fail(
        "Not all 20 selected heme sites passed geometry QC."
    )


pdb_files = sorted(
    FOLDSEEK_DIR.glob(
        "*.pdb"
    )
)


if len(
    pdb_files
) != EXPECTED_JOBS:

    fail(
        f"Expected three Foldseek PDB files; "
        f"found {len(pdb_files)}."
    )


## ================================================================== ##
## Report
## ================================================================== ##

print()
print("=" * 80)

print(
    "STAGE 71 COMPLETE"
)

print("=" * 80)

print()

print(
    "Selected models:"
)


print(
    selected[
        [
            "job",
            "expected_hemes",
            "model_number",
            "ranking_score",
            "sites_passing_3A",
            "max_cys1_heme_distance",
            "max_cys2_heme_distance",
            "max_his_fe_distance",
        ]
    ]
    .to_string(
        index=False
    )
)


print()

print(
    f"Selected structures:         {len(selected)}"
)

print(
    "Expected heme sites:        "
    f"{int(selected['expected_hemes'].sum())}"
)

print(
    "Selected sites passing QC: "
    f"{int(selected['sites_passing_3A'].sum())}/"
    f"{EXPECTED_TOTAL_HEMES}"
)

print()

print(
    "Foldseek protein-only queries:"
)


for path in pdb_files:

    print(
        f"  {path}"
    )


print()

print(
    "QC table:"
)

print(
    f"  {OUT_ALL}"
)

print()

print(
    "Selected models:"
)

print(
    f"  {OUT_SELECTED}"
)

print()

print(
    "Selected heme-site geometry:"
)

print(
    f"  {OUT_SITES}"
)
