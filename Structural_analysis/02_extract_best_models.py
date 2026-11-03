#!/usr/bin/env python3

import csv
import json
import re
import shutil
import zipfile
from pathlib import Path


## ================================================================ ##
## Paths
## ================================================================ ##

ROOT = (
    Path.home()
    / "methanotrophs"
    / "methanotroph_project"
    / "jeppe"
    / "genome_analysis_workflow"
    / "18_structure_prediction"
    / "09_mature_heme_structure_figure"
)

ZIP_DIR = ROOT / "02_results"
EXTRACT_DIR = ZIP_DIR / "extracted"
SELECTED_DIR = ZIP_DIR / "selected_models"

INPUT_MANIFEST = (
    ROOT
    / "01_inputs"
    / "mature_heme_manifest.tsv"
)

OUTPUT_MANIFEST = (
    ZIP_DIR
    / "selected_model_manifest.tsv"
)

EXTRACT_DIR.mkdir(
    parents=True,
    exist_ok=True
)

SELECTED_DIR.mkdir(
    parents=True,
    exist_ok=True
)


## ================================================================ ##
## Helpers
## ================================================================ ##

def canonical(s):

    s = s.lower()

    for prefix in (
        "fold_",
        "m10_",
        "m20_",
        "m35_",
    ):
        if s.startswith(prefix):
            pass

    return re.sub(
        r"[^a-z0-9]+",
        "",
        s
    )


def read_json(path):

    try:
        return json.loads(
            path.read_text()
        )

    except Exception:
        return {}


def model_index_from_name(path):

    name = path.name.lower()

    patterns = [
        r"model[_-]?(\d+)",
        r"summary_confidences[_-]?(\d+)",
        r"full_data[_-]?(\d+)",
    ]

    for pattern in patterns:

        m = re.search(
            pattern,
            name
        )

        if m:
            return int(
                m.group(1)
            )

    return None


def score_summary(d):

    """
    AlphaFold Server outputs vary slightly.

    Prefer ranking_score if present.
    Otherwise fall back to ptm / iptm / mean pLDDT.
    """

    candidates = [
        "ranking_score",
        "rankingScore",
        "ptm",
        "pTM",
        "iptm",
        "ipTM",
        "fraction_disordered",
    ]

    for key in candidates:

        value = d.get(key)

        if isinstance(
            value,
            (int, float)
        ):
            return float(value), key

    ## Some files contain chain-level PTM arrays.
    for key in (
        "chain_ptm",
        "chain_iptm",
    ):

        value = d.get(key)

        if (
            isinstance(value, list)
            and value
        ):

            vals = [
                float(x)
                for x in value
                if isinstance(
                    x,
                    (int, float)
                )
            ]

            if vals:
                return (
                    sum(vals) / len(vals),
                    key
                )

    return None, None


def count_hec_in_cif(path):

    text = path.read_text(
        errors="replace"
    )

    ## Count unique HEC residue occurrences conservatively.
    ##
    ## We inspect atom-site lines containing HEC and extract
    ## chain + residue number where possible.

    residues = set()

    for line in text.splitlines():

        if "HEC" not in line:
            continue

        ## CIF atom lines usually begin ATOM/HETATM.
        if not (
            line.startswith("ATOM")
            or line.startswith("HETATM")
        ):
            continue

        parts = line.split()

        ## Store whole relevant token pattern if parsing is uncertain.
        hec_positions = [
            i
            for i, x in enumerate(parts)
            if x == "HEC"
        ]

        if not hec_positions:
            continue

        i = hec_positions[0]

        neighborhood = tuple(
            parts[
                max(0, i - 2):
                min(len(parts), i + 6)
            ]
        )

        residues.add(
            neighborhood
        )

    return len(residues)


## ================================================================ ##
## Read intended jobs
## ================================================================ ##

with open(
    INPUT_MANIFEST
) as fh:

    rows = list(
        csv.DictReader(
            fh,
            delimiter="\t"
        )
    )


expected = {
    canonical(
        row["job_name"]
    ): row
    for row in rows
}


## ================================================================ ##
## Extract ZIPs
## ================================================================ ##

zip_files = sorted(
    ZIP_DIR.glob("*.zip")
)

print(
    f"ZIP archives found: {len(zip_files)}"
)

if len(zip_files) != 10:

    raise RuntimeError(
        f"Expected 10 ZIP files, found {len(zip_files)}"
    )


results = []


for zpath in zip_files:

    zstem = zpath.stem

    outdir = (
        EXTRACT_DIR
        / zstem
    )

    if outdir.exists():
        shutil.rmtree(
            outdir
        )

    outdir.mkdir(
        parents=True
    )

    with zipfile.ZipFile(
        zpath
    ) as z:

        z.extractall(
            outdir
        )


    ## ------------------------------------------------------------ ##
    ## Match ZIP to intended job
    ## ------------------------------------------------------------ ##

    zip_key = canonical(
        zstem
    )

    matches = []

    for key, row in expected.items():

        ## AlphaFold Server adds "fold_" prefix and lowercases names.
        if (
            key in zip_key
            or zip_key in key
        ):
            matches.append(
                row
            )


    if len(matches) != 1:

        ## More forgiving matching:
        ## compare after removing "fold".
        zip_key2 = zip_key.removeprefix(
            "fold"
        )

        matches = [
            row
            for key, row in expected.items()
            if (
                key == zip_key2
                or key in zip_key2
                or zip_key2 in key
            )
        ]


    if len(matches) != 1:

        raise RuntimeError(
            f"\nCould not uniquely map ZIP:\n"
            f"  {zpath.name}\n"
            f"matches={len(matches)}"
        )


    row = matches[0]

    job_name = row[
        "job_name"
    ]


    ## ------------------------------------------------------------ ##
    ## Find output files
    ## ------------------------------------------------------------ ##

    cif_files = sorted(
        outdir.rglob("*.cif")
    )

    summary_files = sorted(
        outdir.rglob(
            "*summary_confidences*.json"
        )
    )

    full_data_files = sorted(
        outdir.rglob(
            "*full_data*.json"
        )
    )


    if not cif_files:

        raise RuntimeError(
            f"No CIF files in {zpath}"
        )


    ## ------------------------------------------------------------ ##
    ## Build score table
    ## ------------------------------------------------------------ ##

    scored = []

    for sf in summary_files:

        d = read_json(
            sf
        )

        score, score_name = score_summary(
            d
        )

        idx = model_index_from_name(
            sf
        )

        scored.append(
            {
                "summary": sf,
                "index": idx,
                "score": score,
                "score_name": score_name,
                "json": d,
            }
        )


    ## ------------------------------------------------------------ ##
    ## Select best model
    ##
    ## ranking_score preferred.
    ## If scores unavailable, use model 0 if present,
    ## otherwise first CIF.
    ## ------------------------------------------------------------ ##

    scored_valid = [
        x
        for x in scored
        if x["score"] is not None
    ]


    if scored_valid:

        best_score_record = max(
            scored_valid,
            key=lambda x: x["score"]
        )

        best_index = best_score_record[
            "index"
        ]

        best_score = best_score_record[
            "score"
        ]

        best_score_name = best_score_record[
            "score_name"
        ]

    else:

        best_index = 0
        best_score = ""
        best_score_name = "fallback_model0"


    ## ------------------------------------------------------------ ##
    ## Match CIF to selected index
    ## ------------------------------------------------------------ ##

    candidate_cifs = []

    if best_index is not None:

        for cif in cif_files:

            idx = model_index_from_name(
                cif
            )

            if idx == best_index:
                candidate_cifs.append(
                    cif
                )


    if len(candidate_cifs) == 1:

        best_cif = candidate_cifs[0]

    elif len(candidate_cifs) > 1:

        ## Prefer files explicitly named model.
        model_named = [
            p
            for p in candidate_cifs
            if "model" in p.name.lower()
        ]

        best_cif = (
            model_named[0]
            if model_named
            else candidate_cifs[0]
        )

    else:

        ## Some AF server downloads contain a single model.cif
        ## without an index.
        if len(cif_files) == 1:

            best_cif = cif_files[0]

        else:

            model0 = [
                p
                for p in cif_files
                if model_index_from_name(p) == 0
            ]

            best_cif = (
                model0[0]
                if model0
                else cif_files[0]
            )


    ## ------------------------------------------------------------ ##
    ## Find corresponding JSON files
    ## ------------------------------------------------------------ ##

    best_summary = None
    best_full_data = None

    if best_index is not None:

        for p in summary_files:

            if model_index_from_name(p) == best_index:
                best_summary = p
                break

        for p in full_data_files:

            if model_index_from_name(p) == best_index:
                best_full_data = p
                break


    if (
        best_summary is None
        and len(summary_files) == 1
    ):
        best_summary = summary_files[0]


    if (
        best_full_data is None
        and len(full_data_files) == 1
    ):
        best_full_data = full_data_files[0]


    ## ------------------------------------------------------------ ##
    ## Copy selected files
    ## ------------------------------------------------------------ ##

    target_dir = (
        SELECTED_DIR
        / job_name
    )

    target_dir.mkdir(
        parents=True,
        exist_ok=True
    )


    dest_cif = (
        target_dir
        / "best_model.cif"
    )

    shutil.copy2(
        best_cif,
        dest_cif
    )


    if best_summary:

        shutil.copy2(
            best_summary,
            target_dir
            / "best_summary_confidences.json"
        )


    if best_full_data:

        shutil.copy2(
            best_full_data,
            target_dir
            / "best_full_data.json"
        )


    hec_count_cif = count_hec_in_cif(
        dest_cif
    )


    results.append(
        {
            "family": row["family"],
            "pair_no": row["pair_no"],
            "label": row["old_label"],
            "protein_id": row["protein_id"],
            "job_name": job_name,
            "expected_hec": row["cxxch_count"],
            "selected_model_index": (
                best_index
                if best_index is not None
                else ""
            ),
            "selection_score": best_score,
            "selection_metric": best_score_name,
            "n_cif_files": len(cif_files),
            "n_summary_files": len(summary_files),
            "hec_records_detected": hec_count_cif,
            "best_model": str(dest_cif),
        }
    )


## ================================================================ ##
## Write selected-model manifest
## ================================================================ ##

fields = [
    "family",
    "pair_no",
    "label",
    "protein_id",
    "job_name",
    "expected_hec",
    "selected_model_index",
    "selection_score",
    "selection_metric",
    "n_cif_files",
    "n_summary_files",
    "hec_records_detected",
    "best_model",
]


with open(
    OUTPUT_MANIFEST,
    "w",
    newline=""
) as fh:

    writer = csv.DictWriter(
        fh,
        fieldnames=fields,
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(
        results
    )


## ================================================================ ##
## Console QC
## ================================================================ ##

print()
print("=" * 120)
print("SELECTED MODELS")
print("=" * 120)

for r in sorted(
    results,
    key=lambda x: (
        x["family"],
        int(x["pair_no"])
    )
):

    print(
        f"{r['family']:10s} "
        f"{r['label']:18s} "
        f"HEC expected={str(r['expected_hec']):>2s} "
        f"CIF-detected={str(r['hec_records_detected']):>2s} "
        f"model={str(r['selected_model_index']):>2s} "
        f"{r['selection_metric']}={r['selection_score']}"
    )


print()
print(
    "Manifest:",
    OUTPUT_MANIFEST
)

print(
    "Selected models:",
    SELECTED_DIR
)

