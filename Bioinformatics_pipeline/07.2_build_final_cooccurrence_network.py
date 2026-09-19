#!/usr/bin/env python3

from __future__ import annotations

import csv
from collections import defaultdict, deque
from pathlib import Path


PROJECT = Path.home() / "methanotrophs/methanotroph_project/jeppe/genome_analysis_workflow"
WORK_DIR = PROJECT / "07_cooccurrence_network"

RETAINED_CLUSTERS = WORK_DIR / "retained_clusters_prevalence_ge5.tsv"
ELIGIBLE_PAIRS = WORK_DIR / "pairwise_jaccard_shared_ge5.tsv"

JACCARD_THRESHOLD = 0.25
MIN_SHARED_GENOMES = 5

TAG = "jaccard025_shared5"

EDGES_OUT = WORK_DIR / f"cytoscape_edges_{TAG}.tsv"
NODES_OUT = WORK_DIR / f"cytoscape_nodes_{TAG}.tsv"
ALL_RETAINED_OUT = WORK_DIR / f"all_retained_nodes_{TAG}.tsv"
ISOLATES_OUT = WORK_DIR / f"isolated_retained_clusters_{TAG}.tsv"
SUMMARY_OUT = WORK_DIR / f"network_summary_{TAG}.tsv"


def read_tsv(path: Path):
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader), (reader.fieldnames or [])


def component_map(nodes, adjacency):
    visited = set()
    components = []

    for node in sorted(nodes):
        if node in visited:
            continue

        queue = deque([node])
        visited.add(node)
        members = []

        while queue:
            current = queue.popleft()
            members.append(current)

            for neighbor in sorted(adjacency[current]):
                if neighbor not in visited:
                    visited.add(neighbor)
                    queue.append(neighbor)

        components.append(sorted(members))

    components.sort(key=lambda members: (-len(members), members[0]))

    node_to_component = {}
    component_sizes = {}

    for idx, members in enumerate(components, start=1):
        component_id = f"Component_{idx:03d}"
        component_sizes[component_id] = len(members)

        for node in members:
            node_to_component[node] = component_id

    return node_to_component, component_sizes


for required in (RETAINED_CLUSTERS, ELIGIBLE_PAIRS):
    if not required.is_file():
        raise SystemExit(f"ERROR: required input not found: {required}")

node_rows, node_fields = read_tsv(RETAINED_CLUSTERS)
pair_rows, pair_fields = read_tsv(ELIGIBLE_PAIRS)

required_node_fields = {
    "cluster",
    "representative_protein",
    "n_proteins",
    "n_genomes",
    "mean_heme_count",
    "n_fegenie_positive",
    "n_fegenie_negative",
    "pct_fegenie_negative",
    "dominant_fegenie_HMM",
    "dominant_signalp_prediction",
    "dominant_deeptmhmm_class",
    "dominant_localization_class",
    "n_strong_unknown_soluble_periplasmic_like",
    "pct_strong_unknown_soluble_periplasmic_like",
}

missing_nodes = required_node_fields.difference(node_fields)
if missing_nodes:
    raise SystemExit(
        f"ERROR: retained node table missing columns: {sorted(missing_nodes)}"
    )

required_pair_fields = {
    "source",
    "target",
    "source_n_genomes",
    "target_n_genomes",
    "shared_genomes",
    "union_genomes",
    "jaccard",
    "min_prevalence_containment",
}

missing_pairs = required_pair_fields.difference(pair_fields)
if missing_pairs:
    raise SystemExit(
        f"ERROR: pair table missing columns: {sorted(missing_pairs)}"
    )

node_by_cluster = {}

for row in node_rows:
    cluster = row["cluster"].strip()

    if cluster in node_by_cluster:
        raise SystemExit(f"ERROR: duplicate retained cluster: {cluster}")

    node_by_cluster[cluster] = row

final_edges = []

for row in pair_rows:
    shared = int(row["shared_genomes"])
    jaccard = float(row["jaccard"])

    if shared < MIN_SHARED_GENOMES:
        continue

    if jaccard < JACCARD_THRESHOLD:
        continue

    source = row["source"].strip()
    target = row["target"].strip()

    if source not in node_by_cluster or target not in node_by_cluster:
        raise SystemExit(
            f"ERROR: edge references non-retained node: {source}, {target}"
        )

    edge = dict(row)
    edge["interaction"] = "cooccurs"
    edge["weight"] = f"{jaccard:.6f}"
    final_edges.append(edge)

final_edges.sort(
    key=lambda row: (
        -float(row["jaccard"]),
        -int(row["shared_genomes"]),
        row["source"],
        row["target"],
    )
)

network_nodes = set()
adjacency = defaultdict(set)
weighted_degree = defaultdict(float)

for row in final_edges:
    source = row["source"]
    target = row["target"]
    weight = float(row["jaccard"])

    network_nodes.add(source)
    network_nodes.add(target)

    adjacency[source].add(target)
    adjacency[target].add(source)

    weighted_degree[source] += weight
    weighted_degree[target] += weight

isolates = sorted(set(node_by_cluster) - network_nodes)

node_to_component, component_sizes = component_map(network_nodes, adjacency)

largest_component = max(component_sizes.values()) if component_sizes else 0

derived_node_fields = [
    "degree",
    "weighted_degree",
    "component",
    "component_size",
    "is_connected_at_threshold",
]

output_node_fields = node_fields + [
    field for field in derived_node_fields if field not in node_fields
]

network_node_rows = []

for cluster in sorted(network_nodes):
    row = dict(node_by_cluster[cluster])
    component = node_to_component[cluster]

    row.update(
        {
            "degree": len(adjacency[cluster]),
            "weighted_degree": f"{weighted_degree[cluster]:.6f}",
            "component": component,
            "component_size": component_sizes[component],
            "is_connected_at_threshold": 1,
        }
    )

    network_node_rows.append(row)

all_retained_rows = []

for cluster in sorted(node_by_cluster):
    row = dict(node_by_cluster[cluster])

    if cluster in network_nodes:
        component = node_to_component[cluster]
        row.update(
            {
                "degree": len(adjacency[cluster]),
                "weighted_degree": f"{weighted_degree[cluster]:.6f}",
                "component": component,
                "component_size": component_sizes[component],
                "is_connected_at_threshold": 1,
            }
        )
    else:
        row.update(
            {
                "degree": 0,
                "weighted_degree": "0.000000",
                "component": "",
                "component_size": 0,
                "is_connected_at_threshold": 0,
            }
        )

    all_retained_rows.append(row)

edge_output_fields = [
    "source",
    "interaction",
    "target",
    "weight",
    "jaccard",
    "shared_genomes",
    "union_genomes",
    "source_n_genomes",
    "target_n_genomes",
    "min_prevalence_containment",
]

with EDGES_OUT.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=edge_output_fields,
        delimiter="\t",
        lineterminator="\n",
        extrasaction="ignore",
    )
    writer.writeheader()
    writer.writerows(final_edges)

with NODES_OUT.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=output_node_fields,
        delimiter="\t",
        lineterminator="\n",
        extrasaction="ignore",
    )
    writer.writeheader()
    writer.writerows(network_node_rows)

with ALL_RETAINED_OUT.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=output_node_fields,
        delimiter="\t",
        lineterminator="\n",
        extrasaction="ignore",
    )
    writer.writeheader()
    writer.writerows(all_retained_rows)

with ISOLATES_OUT.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=node_fields,
        delimiter="\t",
        lineterminator="\n",
        extrasaction="ignore",
    )
    writer.writeheader()
    for cluster in isolates:
        writer.writerow(node_by_cluster[cluster])

possible_edges = (
    len(network_nodes) * (len(network_nodes) - 1) / 2
    if len(network_nodes) >= 2
    else 0
)

network_density = (
    len(final_edges) / possible_edges
    if possible_edges
    else 0.0
)

summary_rows = [
    ("jaccard_threshold", JACCARD_THRESHOLD),
    ("minimum_shared_genomes", MIN_SHARED_GENOMES),
    ("retained_clusters_before_edge_filter", len(node_by_cluster)),
    ("network_nodes", len(network_nodes)),
    ("isolated_retained_clusters", len(isolates)),
    ("network_edges", len(final_edges)),
    ("connected_components", len(component_sizes)),
    ("largest_component_nodes", largest_component),
    ("network_density", f"{network_density:.6f}"),
]

with SUMMARY_OUT.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(["metric", "value"])
    writer.writerows(summary_rows)

print("Final co-occurrence network build complete.")
print(f"Jaccard threshold:             {JACCARD_THRESHOLD}")
print(f"Minimum shared genomes:       {MIN_SHARED_GENOMES}")
print(f"Retained clusters available:  {len(node_by_cluster)}")
print(f"Network nodes:                {len(network_nodes)}")
print(f"Network edges:                {len(final_edges)}")
print(f"Connected components:         {len(component_sizes)}")
print(f"Largest component:            {largest_component}")
print(f"Isolated retained clusters:   {len(isolates)}")
print()
print(f"Cytoscape edges:      {EDGES_OUT}")
print(f"Cytoscape nodes:      {NODES_OUT}")
print(f"All retained nodes:   {ALL_RETAINED_OUT}")
print(f"Isolates:             {ISOLATES_OUT}")
print(f"Summary:              {SUMMARY_OUT}")
