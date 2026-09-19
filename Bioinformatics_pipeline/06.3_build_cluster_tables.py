#!/usr/bin/env python3

from __future__ import annotations

import csv
from collections import Counter, defaultdict
from pathlib import Path


PROJECT = Path.home() / "methanotrophs/methanotroph_project/jeppe/genome_analysis_workflow"
WORK_DIR = PROJECT / "06_clustering"

METADATA = WORK_DIR / "clustering_input_metadata.tsv"
MMSEQS_TSV = WORK_DIR / "mmseqs40_cov80_cluster.tsv"

MEMBERSHIP_OUT = WORK_DIR / "cluster_membership.tsv"
SUMMARY_OUT = WORK_DIR / "cluster_summary.tsv"
GENOME_CLUSTER_OUT = WORK_DIR / "genome_cluster_pairs.tsv"
REP_MAP_OUT = WORK_DIR / "cluster_representatives.tsv"


def as_int(value) -> int:
    value = str(value).strip()
    if value in {"1", "TRUE", "True", "true"}:
        return 1
    if value in {"0", "FALSE", "False", "false", ""}:
        return 0
    return int(float(value))


for path in (METADATA, MMSEQS_TSV):
    if not path.is_file():
        raise SystemExit(f"ERROR: required file not found: {path}")

with METADATA.open("r", encoding="utf-8", newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    metadata_rows = list(reader)
    metadata_fields = reader.fieldnames or []

required_meta = {
    "genome",
    "protein_id",
    "fegenie_positive",
    "fegenie_HMMs",
    "number_of_hemes",
    "signalp_prediction",
    "deeptmhmm_class",
    "localization_class",
    "strong_fegenie_unannotated_soluble_periplasmic_like_candidate",
}

missing = required_meta.difference(metadata_fields)
if missing:
    raise SystemExit(f"ERROR: metadata missing columns: {sorted(missing)}")

meta = {}
for row in metadata_rows:
    pid = row["protein_id"].strip()
    if pid in meta:
        raise SystemExit(f"ERROR: duplicate protein_id in metadata: {pid}")
    meta[pid] = row

rep_to_members = defaultdict(list)
member_to_rep = {}

with MMSEQS_TSV.open("r", encoding="utf-8") as handle:
    for line_no, raw in enumerate(handle, start=1):
        line = raw.rstrip("\r\n")
        if not line:
            continue

        fields = line.split("\t")
        if len(fields) < 2:
            raise SystemExit(
                f"ERROR: malformed MMseqs TSV line {line_no}: {line!r}"
            )

        rep, member = fields[0], fields[1]

        if member in member_to_rep:
            raise SystemExit(
                f"ERROR: protein occurs in more than one cluster: {member}"
            )

        member_to_rep[member] = rep
        rep_to_members[rep].append(member)

metadata_ids = set(meta)
member_ids = set(member_to_rep)

missing_members = sorted(metadata_ids - member_ids)
unexpected_members = sorted(member_ids - metadata_ids)

if missing_members:
    raise SystemExit(
        f"ERROR: {len(missing_members)} metadata proteins missing from MMseqs clusters. "
        f"First: {missing_members[:10]}"
    )

if unexpected_members:
    raise SystemExit(
        f"ERROR: {len(unexpected_members)} MMseqs proteins missing from metadata. "
        f"First: {unexpected_members[:10]}"
    )

# Stable human-readable cluster IDs. Larger clusters receive lower numbers;
# ties are resolved by representative protein ID.
ordered_reps = sorted(
    rep_to_members,
    key=lambda rep: (-len(rep_to_members[rep]), rep),
)

rep_to_cluster = {
    rep: f"Cluster_{idx:05d}"
    for idx, rep in enumerate(ordered_reps, start=1)
}

membership_fields = [
    "cluster",
    "representative_protein",
] + metadata_fields

membership_rows = []

for rep in ordered_reps:
    cluster = rep_to_cluster[rep]

    for member in sorted(rep_to_members[rep]):
        row = {
            "cluster": cluster,
            "representative_protein": rep,
            **meta[member],
        }
        membership_rows.append(row)

with MEMBERSHIP_OUT.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=membership_fields,
        delimiter="\t",
        lineterminator="\n",
        extrasaction="ignore",
    )
    writer.writeheader()
    writer.writerows(membership_rows)

with REP_MAP_OUT.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(["cluster", "representative_protein"])
    for rep in ordered_reps:
        writer.writerow([rep_to_cluster[rep], rep])

summary_fields = [
    "cluster",
    "representative_protein",
    "n_proteins",
    "n_genomes",
    "mean_heme_count",
    "min_heme_count",
    "max_heme_count",
    "n_fegenie_positive",
    "n_fegenie_negative",
    "pct_fegenie_negative",
    "dominant_fegenie_HMM",
    "dominant_signalp_prediction",
    "dominant_deeptmhmm_class",
    "dominant_localization_class",
    "n_strong_unknown_soluble_periplasmic_like",
    "pct_strong_unknown_soluble_periplasmic_like",
]

summary_rows = []

for rep in ordered_reps:
    members = rep_to_members[rep]
    rows = [meta[m] for m in members]

    hemes = [int(r["number_of_hemes"]) for r in rows]
    n_fg_pos = sum(as_int(r["fegenie_positive"]) for r in rows)
    n_fg_neg = len(rows) - n_fg_pos
    n_strong = sum(
        as_int(r["strong_fegenie_unannotated_soluble_periplasmic_like_candidate"])
        for r in rows
    )

    hmm_counter = Counter()
    for r in rows:
        hmms = [x for x in r["fegenie_HMMs"].split(";") if x]
        hmm_counter.update(hmms)

    signalp_counter = Counter(r["signalp_prediction"] for r in rows)
    deep_counter = Counter(r["deeptmhmm_class"] for r in rows)
    loc_counter = Counter(r["localization_class"] for r in rows)

    dominant_hmm = (
        hmm_counter.most_common(1)[0][0]
        if hmm_counter
        else "NA"
    )

    summary_rows.append(
        {
            "cluster": rep_to_cluster[rep],
            "representative_protein": rep,
            "n_proteins": len(rows),
            "n_genomes": len({r["genome"] for r in rows}),
            "mean_heme_count": f"{sum(hemes) / len(hemes):.3f}",
            "min_heme_count": min(hemes),
            "max_heme_count": max(hemes),
            "n_fegenie_positive": n_fg_pos,
            "n_fegenie_negative": n_fg_neg,
            "pct_fegenie_negative": f"{100 * n_fg_neg / len(rows):.3f}",
            "dominant_fegenie_HMM": dominant_hmm,
            "dominant_signalp_prediction": signalp_counter.most_common(1)[0][0],
            "dominant_deeptmhmm_class": deep_counter.most_common(1)[0][0],
            "dominant_localization_class": loc_counter.most_common(1)[0][0],
            "n_strong_unknown_soluble_periplasmic_like": n_strong,
            "pct_strong_unknown_soluble_periplasmic_like": (
                f"{100 * n_strong / len(rows):.3f}"
            ),
        }
    )

with SUMMARY_OUT.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=summary_fields,
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(summary_rows)

# Presence/absence input for the later co-occurrence network.
genome_cluster_pairs = sorted(
    {
        (row["genome"], row["cluster"])
        for row in membership_rows
    }
)

with GENOME_CLUSTER_OUT.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(["genome", "cluster"])
    writer.writerows(genome_cluster_pairs)

print("Cluster table construction complete.")
print(f"Proteins:                {len(membership_rows)}")
print(f"Clusters:                {len(ordered_reps)}")
print(f"Genome-cluster pairs:    {len(genome_cluster_pairs)}")
print(f"Membership table:        {MEMBERSHIP_OUT}")
print(f"Cluster summary:         {SUMMARY_OUT}")
print(f"Genome-cluster pairs:    {GENOME_CLUSTER_OUT}")
print(f"Representative mapping:  {REP_MAP_OUT}")
