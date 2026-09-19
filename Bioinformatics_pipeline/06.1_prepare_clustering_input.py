#!/usr/bin/env python3

from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path


PROJECT = Path.home() / "methanotrophs/methanotroph_project/jeppe/genome_analysis_workflow"

MASTER = PROJECT / "05_integrated_annotations" / "MASTER_PROTEIN_ANNOTATIONS.tsv"
OUT_DIR = PROJECT / "06_clustering"

FASTA_OUT = OUT_DIR / "exported_heme_candidates.faa"
METADATA_OUT = OUT_DIR / "clustering_input_metadata.tsv"
SUMMARY_OUT = OUT_DIR / "clustering_input_summary.tsv"


def as_int(value) -> int:
    value = str(value).strip()
    if value in {"1", "TRUE", "True", "true"}:
        return 1
    if value in {"0", "FALSE", "False", "false", ""}:
        return 0
    return int(float(value))


if not MASTER.is_file():
    raise SystemExit(f"ERROR: master annotation table not found: {MASTER}")

OUT_DIR.mkdir(parents=True, exist_ok=True)

with MASTER.open("r", encoding="utf-8", newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    rows = list(reader)
    fields = reader.fieldnames or []

required = {
    "genome",
    "protein_id",
    "findmehemes_positive",
    "number_of_hemes",
    "export_evidence",
    "fegenie_positive",
    "fegenie_HMMs",
    "signalp_prediction",
    "deeptmhmm_class",
    "deeptmhmm_sequence",
    "localization_class",
    "strong_fegenie_unannotated_soluble_periplasmic_like_candidate",
}

missing = required.difference(fields)
if missing:
    raise SystemExit(
        f"ERROR: MASTER_PROTEIN_ANNOTATIONS.tsv is missing columns: {sorted(missing)}"
    )

selected = []
seen_ids = set()

for row in rows:
    if as_int(row["findmehemes_positive"]) != 1:
        continue

    if as_int(row["export_evidence"]) != 1:
        continue

    protein_id = row["protein_id"].strip()
    genome = row["genome"].strip()

    if not protein_id or not genome:
        raise SystemExit("ERROR: selected row has an empty genome or protein_id.")

    if protein_id in seen_ids:
        raise SystemExit(
            "ERROR: protein_id is not globally unique in clustering input: "
            f"{protein_id}"
        )
    seen_ids.add(protein_id)

    sequence = row["deeptmhmm_sequence"].strip()

    # Prodigal translations in this workflow retain a terminal stop symbol.
    # Remove only the terminal stop for a conventional protein FASTA.
    if sequence.endswith("*"):
        sequence = sequence[:-1]

    if "*" in sequence:
        raise SystemExit(
            f"ERROR: internal stop character found in selected protein {protein_id}"
        )

    if not sequence:
        raise SystemExit(f"ERROR: empty sequence for selected protein {protein_id}")

    row = dict(row)
    row["clustering_sequence_length"] = str(len(sequence))
    row["clustering_sequence"] = sequence
    selected.append(row)

if not selected:
    raise SystemExit("ERROR: no exported heme-positive proteins were selected.")

# ---------------------------------------------------------------------------
# FASTA
# ---------------------------------------------------------------------------

with FASTA_OUT.open("w", encoding="utf-8") as handle:
    for row in selected:
        handle.write(f">{row['protein_id']}\n")
        sequence = row["clustering_sequence"]
        for start in range(0, len(sequence), 80):
            handle.write(sequence[start : start + 80] + "\n")

# ---------------------------------------------------------------------------
# Metadata sidecar
# ---------------------------------------------------------------------------

metadata_fields = [
    "genome",
    "protein_id",
    "length",
    "clustering_sequence_length",
    "candidate_source",
    "fegenie_positive",
    "fegenie_HMMs",
    "fegenie_categories",
    "fegenie_max_bitscore",
    "fegenie_max_cutoff",
    "findmehemes_positive",
    "number_of_hemes",
    "signalp_prediction",
    "signalp_cs_position",
    "deeptmhmm_class",
    "deeptmhmm_n_tm_helices",
    "deeptmhmm_sp_length",
    "deeptmhmm_side_after_sp",
    "deeptmhmm_c_terminal_side",
    "export_evidence",
    "localization_class",
    "fegenie_unannotated_heme_candidate",
    "fegenie_unannotated_exported_heme_candidate",
    "strong_fegenie_unannotated_soluble_periplasmic_like_candidate",
    "source_proteome",
]

metadata_fields = [f for f in metadata_fields if f in selected[0]]

with METADATA_OUT.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=metadata_fields,
        delimiter="\t",
        lineterminator="\n",
        extrasaction="ignore",
    )
    writer.writeheader()
    writer.writerows(selected)

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

n_total = len(selected)
n_genomes = len({r["genome"] for r in selected})
n_fegenie_pos = sum(as_int(r["fegenie_positive"]) for r in selected)
n_fegenie_neg = n_total - n_fegenie_pos
n_strong_unknown = sum(
    as_int(r["strong_fegenie_unannotated_soluble_periplasmic_like_candidate"])
    for r in selected
)

localization_counts = Counter(r["localization_class"] for r in selected)
signalp_counts = Counter(r["signalp_prediction"] for r in selected)
deep_counts = Counter(r["deeptmhmm_class"] for r in selected)

summary_rows = [
    ("selected_exported_heme_proteins", n_total),
    ("genomes_represented", n_genomes),
    ("fegenie_positive_selected", n_fegenie_pos),
    ("fegenie_negative_selected", n_fegenie_neg),
    ("strong_unknown_soluble_periplasmic_like_selected", n_strong_unknown),
]

for key, count in sorted(localization_counts.items()):
    summary_rows.append((f"localization::{key}", count))

for key, count in sorted(signalp_counts.items()):
    summary_rows.append((f"signalp::{key}", count))

for key, count in sorted(deep_counts.items()):
    summary_rows.append((f"deeptmhmm::{key}", count))

with SUMMARY_OUT.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(["metric", "value"])
    writer.writerows(summary_rows)

print("Clustering input preparation complete.")
print(f"Selected exported heme proteins:       {n_total}")
print(f"Genomes represented:                   {n_genomes}")
print(f"FeGenie-positive selected:             {n_fegenie_pos}")
print(f"FeGenie-negative selected:             {n_fegenie_neg}")
print(f"Strong unknown periplasmic-like:       {n_strong_unknown}")
print()
print(f"FASTA:    {FASTA_OUT}")
print(f"Metadata: {METADATA_OUT}")
print(f"Summary:  {SUMMARY_OUT}")
