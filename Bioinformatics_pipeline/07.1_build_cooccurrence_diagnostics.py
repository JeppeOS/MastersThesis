#!/usr/bin/env python3

from __future__ import annotations

import csv
from collections import defaultdict, Counter, deque
from itertools import combinations
from pathlib import Path


PROJECT = Path.home() / "methanotrophs/methanotroph_project/jeppe/genome_analysis_workflow"

CLUSTER_DIR = PROJECT / "06_clustering"
WORK_DIR = PROJECT / "07_cooccurrence_network"

GENOME_CLUSTER_TSV = CLUSTER_DIR / "genome_cluster_pairs.tsv"
CLUSTER_SUMMARY_TSV = CLUSTER_DIR / "cluster_summary.tsv"

MIN_PREVALENCE = 5
MIN_SHARED_GENOMES = 5

RETAINED_CLUSTERS_OUT = WORK_DIR / "retained_clusters_prevalence_ge5.tsv"
PAIRWISE_OUT = WORK_DIR / "pairwise_jaccard_all_retained.tsv"
ELIGIBLE_EDGES_OUT = WORK_DIR / "pairwise_jaccard_shared_ge5.tsv"
JACCARD_BINS_OUT = WORK_DIR / "jaccard_distribution.tsv"
THRESHOLD_SENSITIVITY_OUT = WORK_DIR / "jaccard_threshold_sensitivity.tsv"
PREVALENCE_OUT = WORK_DIR / "cluster_prevalence_distribution.tsv"
SUMMARY_OUT = WORK_DIR / "cooccurrence_diagnostics_summary.tsv"


def read_tsv(path: Path):
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader), (reader.fieldnames or [])


def connected_components(nodes, edges):
    adjacency = {node: set() for node in nodes}

    for a, b in edges:
        adjacency[a].add(b)
        adjacency[b].add(a)

    visited = set()
    component_sizes = []

    for node in nodes:
        if node in visited:
            continue

        queue = deque([node])
        visited.add(node)
        size = 0

        while queue:
            current = queue.popleft()
            size += 1

            for neighbor in adjacency[current]:
                if neighbor not in visited:
                    visited.add(neighbor)
                    queue.append(neighbor)

        component_sizes.append(size)

    component_sizes.sort(reverse=True)
    return component_sizes


for required in (GENOME_CLUSTER_TSV, CLUSTER_SUMMARY_TSV):
    if not required.is_file():
        raise SystemExit(f"ERROR: required input not found: {required}")

WORK_DIR.mkdir(parents=True, exist_ok=True)

pairs_rows, pair_fields = read_tsv(GENOME_CLUSTER_TSV)
summary_rows, summary_fields = read_tsv(CLUSTER_SUMMARY_TSV)

if not {"genome", "cluster"}.issubset(pair_fields):
    raise SystemExit(
        f"ERROR: {GENOME_CLUSTER_TSV} must contain genome and cluster columns."
    )

required_summary = {
    "cluster",
    "representative_protein",
    "n_proteins",
    "n_genomes",
    "mean_heme_count",
    "pct_fegenie_negative",
    "dominant_fegenie_HMM",
    "dominant_signalp_prediction",
    "dominant_deeptmhmm_class",
    "dominant_localization_class",
    "pct_strong_unknown_soluble_periplasmic_like",
}
missing_summary = required_summary.difference(summary_fields)
if missing_summary:
    raise SystemExit(
        "ERROR: cluster_summary.tsv missing columns: "
        f"{sorted(missing_summary)}"
    )

cluster_to_genomes = defaultdict(set)
seen_pairs = set()

for row in pairs_rows:
    genome = row["genome"].strip()
    cluster = row["cluster"].strip()

    if not genome or not cluster:
        raise SystemExit("ERROR: empty genome or cluster in genome_cluster_pairs.tsv")

    key = (genome, cluster)
    if key in seen_pairs:
        raise SystemExit(f"ERROR: duplicate genome-cluster pair: {key}")

    seen_pairs.add(key)
    cluster_to_genomes[cluster].add(genome)

summary_by_cluster = {}

for row in summary_rows:
    cluster = row["cluster"].strip()

    if cluster in summary_by_cluster:
        raise SystemExit(f"ERROR: duplicate cluster in cluster_summary.tsv: {cluster}")

    summary_by_cluster[cluster] = row

pair_clusters = set(cluster_to_genomes)
summary_clusters = set(summary_by_cluster)

missing_summary_clusters = sorted(pair_clusters - summary_clusters)
missing_pair_clusters = sorted(summary_clusters - pair_clusters)

if missing_summary_clusters:
    raise SystemExit(
        "ERROR: clusters in genome_cluster_pairs.tsv missing from cluster_summary.tsv. "
        f"First: {missing_summary_clusters[:10]}"
    )

if missing_pair_clusters:
    raise SystemExit(
        "ERROR: clusters in cluster_summary.tsv missing from genome_cluster_pairs.tsv. "
        f"First: {missing_pair_clusters[:10]}"
    )

for cluster, genomes in cluster_to_genomes.items():
    observed = len(genomes)
    expected = int(summary_by_cluster[cluster]["n_genomes"])

    if observed != expected:
        raise SystemExit(
            f"ERROR: n_genomes mismatch for {cluster}: "
            f"pairs={observed}, cluster_summary={expected}"
        )

prevalence_counter = Counter(len(v) for v in cluster_to_genomes.values())

with PREVALENCE_OUT.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(["n_genomes", "n_clusters"])
    for prevalence, count in sorted(prevalence_counter.items()):
        writer.writerow([prevalence, count])

retained = sorted(
    cluster
    for cluster, genomes in cluster_to_genomes.items()
    if len(genomes) >= MIN_PREVALENCE
)

with RETAINED_CLUSTERS_OUT.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=summary_fields,
        delimiter="\t",
        lineterminator="\n",
        extrasaction="ignore",
    )
    writer.writeheader()

    for cluster in retained:
        writer.writerow(summary_by_cluster[cluster])

pair_fields_out = [
    "source",
    "target",
    "source_n_genomes",
    "target_n_genomes",
    "shared_genomes",
    "union_genomes",
    "jaccard",
    "min_prevalence_containment",
]

all_pair_rows = []
eligible_rows = []

for source, target in combinations(retained, 2):
    a = cluster_to_genomes[source]
    b = cluster_to_genomes[target]

    shared = len(a & b)
    union = len(a | b)

    jaccard = shared / union if union else 0.0
    containment = shared / min(len(a), len(b)) if min(len(a), len(b)) else 0.0

    row = {
        "source": source,
        "target": target,
        "source_n_genomes": len(a),
        "target_n_genomes": len(b),
        "shared_genomes": shared,
        "union_genomes": union,
        "jaccard": f"{jaccard:.6f}",
        "min_prevalence_containment": f"{containment:.6f}",
    }

    all_pair_rows.append(row)

    if shared >= MIN_SHARED_GENOMES:
        eligible_rows.append(row)

with PAIRWISE_OUT.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=pair_fields_out,
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(all_pair_rows)

with ELIGIBLE_EDGES_OUT.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=pair_fields_out,
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(eligible_rows)

bin_edges = [i / 20 for i in range(21)]
bin_counts = [0 for _ in range(20)]

eligible_values = [float(r["jaccard"]) for r in eligible_rows]

for value in eligible_values:
    if value >= 1.0:
        idx = 19
    else:
        idx = int(value / 0.05)
        idx = max(0, min(idx, 19))
    bin_counts[idx] += 1

with JACCARD_BINS_OUT.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(["jaccard_lower", "jaccard_upper", "n_pairs"])

    for idx, count in enumerate(bin_counts):
        lower = bin_edges[idx]
        upper = bin_edges[idx + 1]
        writer.writerow([f"{lower:.2f}", f"{upper:.2f}", count])

thresholds = [
    0.10,
    0.15,
    0.20,
    0.25,
    0.30,
    0.35,
    0.40,
    0.45,
    0.50,
    0.55,
    0.60,
    0.65,
    0.70,
    0.75,
    0.80,
]

sensitivity_rows = []

for threshold in thresholds:
    threshold_edges = [
        (row["source"], row["target"])
        for row in eligible_rows
        if float(row["jaccard"]) >= threshold
    ]

    nodes_with_edges = set()
    for a, b in threshold_edges:
        nodes_with_edges.add(a)
        nodes_with_edges.add(b)

    components = connected_components(nodes_with_edges, threshold_edges)

    n_nodes = len(nodes_with_edges)
    n_edges = len(threshold_edges)
    possible_edges = n_nodes * (n_nodes - 1) / 2 if n_nodes >= 2 else 0
    density = n_edges / possible_edges if possible_edges else 0.0

    sensitivity_rows.append(
        {
            "jaccard_threshold": f"{threshold:.2f}",
            "min_shared_genomes": MIN_SHARED_GENOMES,
            "nodes_with_edges": n_nodes,
            "edges": n_edges,
            "connected_components": len(components),
            "largest_component_nodes": components[0] if components else 0,
            "network_density": f"{density:.6f}",
        }
    )

with THRESHOLD_SENSITIVITY_OUT.open("w", encoding="utf-8", newline="") as handle:
    fields = [
        "jaccard_threshold",
        "min_shared_genomes",
        "nodes_with_edges",
        "edges",
        "connected_components",
        "largest_component_nodes",
        "network_density",
    ]
    writer = csv.DictWriter(
        handle,
        fieldnames=fields,
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(sensitivity_rows)

n_all_clusters = len(cluster_to_genomes)
n_retained = len(retained)
n_removed = n_all_clusters - n_retained
n_genomes = len({genome for genome, _ in seen_pairs})
n_all_pairs = len(all_pair_rows)
n_eligible = len(eligible_rows)

summary_metrics = [
    ("input_clusters", n_all_clusters),
    ("input_genomes", n_genomes),
    ("input_genome_cluster_pairs", len(seen_pairs)),
    ("minimum_cluster_prevalence", MIN_PREVALENCE),
    ("retained_clusters", n_retained),
    ("removed_low_prevalence_clusters", n_removed),
    ("all_retained_cluster_pairs", n_all_pairs),
    ("minimum_shared_genomes_for_edge_eligibility", MIN_SHARED_GENOMES),
    ("pairs_with_shared_genomes_ge5", n_eligible),
]

with SUMMARY_OUT.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(["metric", "value"])
    writer.writerows(summary_metrics)

print("Co-occurrence diagnostics complete.")
print(f"Input clusters:                     {n_all_clusters}")
print(f"Input genomes:                      {n_genomes}")
print(f"Genome-cluster pairs:               {len(seen_pairs)}")
print(f"Clusters retained (>= {MIN_PREVALENCE} genomes):      {n_retained}")
print(f"Clusters removed (< {MIN_PREVALENCE} genomes):       {n_removed}")
print(f"All retained cluster pairs:         {n_all_pairs}")
print(f"Pairs sharing >= {MIN_SHARED_GENOMES} genomes:         {n_eligible}")
print()
print(f"Retained clusters:      {RETAINED_CLUSTERS_OUT}")
print(f"All pairwise Jaccard:   {PAIRWISE_OUT}")
print(f"Eligible pair table:    {ELIGIBLE_EDGES_OUT}")
print(f"Jaccard distribution:   {JACCARD_BINS_OUT}")
print(f"Threshold sensitivity:  {THRESHOLD_SENSITIVITY_OUT}")
print(f"Summary:                {SUMMARY_OUT}")
