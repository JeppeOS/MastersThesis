#!/usr/bin/env python3

from __future__ import annotations

import csv
import re
from collections import Counter
from pathlib import Path


PROJECT = Path.home() / "methanotrophs/methanotroph_project/jeppe/genome_analysis_workflow"

DEEPTMHMM_DIR = PROJECT / "04_deeptmhmm"
MANIFEST = DEEPTMHMM_DIR / "deeptmhmm_candidate_manifest.tsv"
RESULTS_DIR = DEEPTMHMM_DIR / "results_by_genome"

OUT_DIR = PROJECT / "05_integrated_annotations"

MASTER_OUT = OUT_DIR / "MASTER_PROTEIN_ANNOTATIONS.tsv"
UNKNOWN_EXPORTED_OUT = OUT_DIR / "unknown_exported_cytochrome_candidates.tsv"
UNKNOWN_SOLUBLE_OUT = OUT_DIR / "unknown_soluble_periplasmic_like_cytochromes.tsv"
UNKNOWN_LIPO_OUT = OUT_DIR / "unknown_lipoprotein_cytochromes.tsv"
UNKNOWN_MEMBRANE_OUT = OUT_DIR / "unknown_membrane_associated_cytochromes.tsv"
SUMMARY_OUT = OUT_DIR / "integration_summary.tsv"
CLASS_COUNTS_OUT = OUT_DIR / "localization_class_counts.tsv"


SIGNALP_EXPORT_CLASSES = {"SP", "LIPO", "TAT", "TATLIPO", "PILIN"}
SIGNALP_SOLUBLE_EXPORT_CLASSES = {"SP", "TAT"}
SIGNALP_LIPO_CLASSES = {"LIPO", "TATLIPO"}


def count_runs(text: str, char: str) -> int:
    return sum(1 for _ in re.finditer(f"{re.escape(char)}+", text))


def first_non_s(topology: str) -> str:
    stripped = topology.lstrip("S")
    return stripped[0] if stripped else ""


def last_side(topology: str) -> str:
    for char in reversed(topology):
        if char in {"I", "O"}:
            return char
    return ""


def parse_deeptmhmm_3line(path: Path):
    lines = []

    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.rstrip("\r\n")
            if line.strip():
                lines.append(line)

    if len(lines) % 3 != 0:
        raise ValueError(
            f"{path} contains {len(lines)} non-empty lines; expected a multiple of 3."
        )

    for i in range(0, len(lines), 3):
        header = lines[i]
        sequence = lines[i + 1]
        topology = lines[i + 2]

        if not header.startswith(">"):
            raise ValueError(f"Unexpected DeepTMHMM header in {path}: {header!r}")

        header_text = header[1:].strip()

        if "|" not in header_text:
            raise ValueError(f"DeepTMHMM header lacks class separator '|': {header!r}")

        protein_part, class_part = header_text.rsplit("|", 1)
        protein_id = protein_part.strip().split()[0]
        deeptmhmm_class = class_part.strip()

        sequence_clean = sequence.strip()
        topology_clean = topology.strip()

        # For these DeepTMHMM outputs, the topology string has one state for
        # every character in the emitted sequence line, including the terminal
        # '*' retained from the Prodigal translation.
        if len(sequence_clean) != len(topology_clean):
            raise ValueError(
                f"Length mismatch for {protein_id} in {path}: "
                f"sequence={len(sequence_clean)}, topology={len(topology_clean)}"
            )

        yield {
            "protein_id": protein_id,
            "deeptmhmm_class": deeptmhmm_class,
            "deeptmhmm_sequence": sequence_clean,
            "deeptmhmm_topology": topology_clean,
            "deeptmhmm_n_tm_helices": count_runs(topology_clean, "M"),
            "deeptmhmm_sp_length": len(topology_clean) - len(topology_clean.lstrip("S")),
            "deeptmhmm_side_after_sp": first_non_s(topology_clean),
            "deeptmhmm_c_terminal_side": last_side(topology_clean),
        }


def localization_class(signalp_prediction: str, deeptmhmm_class: str, n_tm: int) -> str:
    sp = signalp_prediction.strip().upper()
    dt = deeptmhmm_class.strip().upper()

    signalp_export = sp in SIGNALP_EXPORT_CLASSES
    deep_sp = dt == "SP"

    if sp in SIGNALP_LIPO_CLASSES:
        return "exported_lipoprotein_candidate"

    if sp == "PILIN":
        return "pilin_like_exported_candidate"

    if (sp in SIGNALP_SOLUBLE_EXPORT_CLASSES or deep_sp) and n_tm == 0:
        return "soluble_exported_periplasmic_like_candidate"

    if signalp_export or deep_sp:
        if n_tm == 1:
            return "exported_single_pass_membrane_candidate"
        if n_tm >= 2:
            return "exported_multipass_membrane_candidate"
        return "exported_candidate_uncertain_topology"

    if n_tm == 1:
        return "single_pass_membrane_no_export_signal"

    if n_tm >= 2:
        return "multipass_membrane_no_export_signal"

    if dt == "GLOB":
        return "globular_no_export_signal"

    return "uncertain"


def bool01(value) -> int:
    return 1 if str(value).strip() in {"1", "TRUE", "True", "true"} else 0


if not MANIFEST.is_file():
    raise SystemExit(f"ERROR: candidate manifest not found: {MANIFEST}")

if not RESULTS_DIR.is_dir():
    raise SystemExit(f"ERROR: DeepTMHMM results directory not found: {RESULTS_DIR}")

OUT_DIR.mkdir(parents=True, exist_ok=True)

with MANIFEST.open("r", encoding="utf-8", newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    manifest_rows = list(reader)
    manifest_fields = reader.fieldnames or []

required_manifest = {
    "genome",
    "protein_id",
    "length",
    "candidate_source",
    "fegenie_positive",
    "fegenie_HMMs",
    "findmehemes_positive",
    "number_of_hemes",
    "signalp_prediction",
}

missing = required_manifest.difference(manifest_fields)
if missing:
    raise SystemExit(
        f"ERROR: candidate manifest missing required columns: {sorted(missing)}"
    )

manifest_by_key = {}
for row in manifest_rows:
    key = (row["genome"].strip(), row["protein_id"].strip())
    if key in manifest_by_key:
        raise SystemExit(f"ERROR: duplicate manifest key: {key}")
    manifest_by_key[key] = row

deep_by_key = {}
result_files = sorted(RESULTS_DIR.glob("*/predicted_topologies.3line"))

if len(result_files) != 631:
    print(
        f"WARNING: found {len(result_files)} predicted_topologies.3line files; "
        "expected 631."
    )

for result_file in result_files:
    genome = result_file.parent.name
    for deep_row in parse_deeptmhmm_3line(result_file):
        key = (genome, deep_row["protein_id"])
        if key in deep_by_key:
            raise SystemExit(f"ERROR: duplicate DeepTMHMM result key: {key}")
        deep_by_key[key] = deep_row

manifest_keys = set(manifest_by_key)
deep_keys = set(deep_by_key)

missing_deep = sorted(manifest_keys - deep_keys)
unexpected_deep = sorted(deep_keys - manifest_keys)

if missing_deep:
    preview = "\n".join(f"  {g}\t{p}" for g, p in missing_deep[:20])
    raise SystemExit(
        "ERROR: candidate proteins are missing DeepTMHMM results.\n"
        f"Missing count: {len(missing_deep)}\n{preview}"
    )

if unexpected_deep:
    preview = "\n".join(f"  {g}\t{p}" for g, p in unexpected_deep[:20])
    raise SystemExit(
        "ERROR: DeepTMHMM results contain proteins not present in the candidate manifest.\n"
        f"Unexpected count: {len(unexpected_deep)}\n{preview}"
    )

final_rows = []

for row in manifest_rows:
    genome = row["genome"].strip()
    protein_id = row["protein_id"].strip()
    deep = deep_by_key[(genome, protein_id)]

    signalp_prediction = row.get("signalp_prediction", "").strip().upper()
    deeptmhmm_class = deep["deeptmhmm_class"].strip().upper()
    n_tm = int(deep["deeptmhmm_n_tm_helices"])

    signalp_export = int(signalp_prediction in SIGNALP_EXPORT_CLASSES)
    deeptmhmm_sp = int(deeptmhmm_class == "SP")
    export_evidence = int(signalp_export or deeptmhmm_sp)

    loc_class = localization_class(
        signalp_prediction=signalp_prediction,
        deeptmhmm_class=deeptmhmm_class,
        n_tm=n_tm,
    )

    fegenie_positive = bool01(row.get("fegenie_positive", "0"))
    heme_positive = bool01(row.get("findmehemes_positive", "0"))

    unknown_fegenie_cytochrome = int(
        fegenie_positive == 0 and heme_positive == 1
    )

    unknown_exported_cytochrome = int(
        unknown_fegenie_cytochrome == 1 and export_evidence == 1
    )

    strong_unknown_soluble = int(
        unknown_fegenie_cytochrome == 1
        and signalp_prediction in SIGNALP_SOLUBLE_EXPORT_CLASSES
        and deeptmhmm_class == "SP"
        and n_tm == 0
    )

    merged = dict(row)
    merged.update(deep)
    merged.update(
        {
            "signalp_export_positive": signalp_export,
            "deeptmhmm_sp_positive": deeptmhmm_sp,
            "export_evidence": export_evidence,
            "localization_class": loc_class,
            "fegenie_unannotated_heme_candidate": unknown_fegenie_cytochrome,
            "fegenie_unannotated_exported_heme_candidate": unknown_exported_cytochrome,
            "strong_fegenie_unannotated_soluble_periplasmic_like_candidate": strong_unknown_soluble,
        }
    )
    final_rows.append(merged)

deep_fields = [
    "deeptmhmm_class",
    "deeptmhmm_sequence",
    "deeptmhmm_topology",
    "deeptmhmm_n_tm_helices",
    "deeptmhmm_sp_length",
    "deeptmhmm_side_after_sp",
    "deeptmhmm_c_terminal_side",
]

derived_fields = [
    "signalp_export_positive",
    "deeptmhmm_sp_positive",
    "export_evidence",
    "localization_class",
    "fegenie_unannotated_heme_candidate",
    "fegenie_unannotated_exported_heme_candidate",
    "strong_fegenie_unannotated_soluble_periplasmic_like_candidate",
]

output_fields = manifest_fields + [
    f for f in deep_fields + derived_fields if f not in manifest_fields
]


def write_tsv(path: Path, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=output_fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


write_tsv(MASTER_OUT, final_rows)

unknown_exported = [
    r for r in final_rows
    if int(r["fegenie_unannotated_exported_heme_candidate"]) == 1
]

unknown_soluble = [
    r for r in final_rows
    if r["localization_class"] == "soluble_exported_periplasmic_like_candidate"
    and int(r["fegenie_unannotated_heme_candidate"]) == 1
]

unknown_lipo = [
    r for r in final_rows
    if r["localization_class"] == "exported_lipoprotein_candidate"
    and int(r["fegenie_unannotated_heme_candidate"]) == 1
]

unknown_membrane = [
    r for r in final_rows
    if r["localization_class"]
    in {
        "exported_single_pass_membrane_candidate",
        "exported_multipass_membrane_candidate",
    }
    and int(r["fegenie_unannotated_heme_candidate"]) == 1
]

write_tsv(UNKNOWN_EXPORTED_OUT, unknown_exported)
write_tsv(UNKNOWN_SOLUBLE_OUT, unknown_soluble)
write_tsv(UNKNOWN_LIPO_OUT, unknown_lipo)
write_tsv(UNKNOWN_MEMBRANE_OUT, unknown_membrane)

localization_counts = Counter(r["localization_class"] for r in final_rows)

with CLASS_COUNTS_OUT.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(["localization_class", "protein_count"])
    for class_name, count in sorted(localization_counts.items()):
        writer.writerow([class_name, count])

summary_rows = [
    ("total_integrated_candidates", len(final_rows)),
    ("fegenie_positive", sum(bool01(r["fegenie_positive"]) for r in final_rows)),
    ("findmehemes_positive", sum(bool01(r["findmehemes_positive"]) for r in final_rows)),
    ("signalp_export_positive", sum(int(r["signalp_export_positive"]) for r in final_rows)),
    ("deeptmhmm_sp_positive", sum(int(r["deeptmhmm_sp_positive"]) for r in final_rows)),
    (
        "fegenie_unannotated_heme_candidates",
        sum(int(r["fegenie_unannotated_heme_candidate"]) for r in final_rows),
    ),
    ("fegenie_unannotated_exported_heme_candidates", len(unknown_exported)),
    ("fegenie_unannotated_soluble_periplasmic_like_candidates", len(unknown_soluble)),
    (
        "strong_fegenie_unannotated_soluble_periplasmic_like_candidates",
        sum(
            int(r["strong_fegenie_unannotated_soluble_periplasmic_like_candidate"])
            for r in final_rows
        ),
    ),
    ("fegenie_unannotated_lipoprotein_candidates", len(unknown_lipo)),
    ("fegenie_unannotated_exported_membrane_candidates", len(unknown_membrane)),
]

with SUMMARY_OUT.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(["metric", "value"])
    writer.writerows(summary_rows)

print("Integrated annotation build complete.")
print(f"DeepTMHMM result files:                              {len(result_files)}")
print(f"Integrated candidate proteins:                       {len(final_rows)}")
print(f"FeGenie-unannotated heme candidates:                 {summary_rows[5][1]}")
print(f"FeGenie-unannotated exported heme candidates:        {len(unknown_exported)}")
print(f"FeGenie-unannotated soluble/periplasmic-like:        {len(unknown_soluble)}")
print(
    "Strong FeGenie-unannotated soluble/periplasmic-like: "
    f"{summary_rows[8][1]}"
)
print(f"FeGenie-unannotated lipoprotein candidates:          {len(unknown_lipo)}")
print(f"FeGenie-unannotated exported membrane candidates:    {len(unknown_membrane)}")
print()
print(f"Master table: {MASTER_OUT}")
print(f"Summary:      {SUMMARY_OUT}")
print(f"Class counts: {CLASS_COUNTS_OUT}")
