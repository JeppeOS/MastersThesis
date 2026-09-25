#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import numpy as np


## ================================================================== ##
## Configuration
## ================================================================== ##

WF = (
    Path.home()
    / "methanotrophs"
    / "methanotroph_project"
    / "jeppe"
    / "genome_analysis_workflow"
)

TARGET_CLUSTER = "Cluster_00063"
TARGET_MODULE = "Module_06"

FOCAL_OTHER_MODULES = {
    "Module_10",
    "Module_20",
    "Module_35",
}

N_METHYLOCOCCUS = 5
N_METHYLOBACTER = 5
N_OTHER = 2


## ================================================================== ##
## Paths
## ================================================================== ##

FRAMEWORK = (
    WF
    / "16_visualization"
    / "16A_module06_gene_map_framework"
)

OUTDIR = (
    WF
    / "16_visualization"
    / "16A3_module06_representative_selection"
)

OUTDIR.mkdir(
    parents=True,
    exist_ok=True
)

GENE_MAP = (
    FRAMEWORK
    / "Module_06_Cluster_00063_gene_map_data.tsv"
)

CLUSTER_PROTEINS = (
    FRAMEWORK
    / "cluster00063_proteins.tsv"
)

PROTEIN_MODULES = (
    WF
    / "09_module_occurrence"
    / "protein_module_membership.tsv"
)

OBSERVABILITY = (
    WF
    / "12_neighborhood_extraction"
    / "focal_context_observability.tsv"
)

TAXONOMY = (
    WF
    / "16_visualization"
    / "globdb_r226_taxonomy.tsv"
)


## ================================================================== ##
## Helpers
## ================================================================== ##

def as_num(x):
    return pd.to_numeric(x, errors="coerce")


def rank_group(genus):

    genus = "" if pd.isna(genus) else str(genus)

    if genus == "Methylococcus":
        return "Methylococcus"

    if genus.startswith("Methylobacter"):
        return "Methylobacter"

    return "Other"


def rank_candidates(df):

    return df.sort_values(
        [
            "n_other_focal_modules",
            "focal_module_completeness_sum",
            "n_cluster00063_copies",
            "best_full_20kb",
            "best_full_10kb",
            "best_local_heme",
            "best_local_module06",
            "genome",
        ],
        ascending=[
            False,
            False,
            False,
            False,
            False,
            False,
            False,
            True,
        ],
    )


## ================================================================== ##
## Read data
## ================================================================== ##

genes = pd.read_csv(
    GENE_MAP,
    sep="\t",
    low_memory=False
)

cluster_proteins = pd.read_csv(
    CLUSTER_PROTEINS,
    sep="\t",
    low_memory=False
)

pm = pd.read_csv(
    PROTEIN_MODULES,
    sep="\t",
    low_memory=False
)

obs = pd.read_csv(
    OBSERVABILITY,
    sep="\t",
    low_memory=False
)


## ================================================================== ##
## Taxonomy
## ================================================================== ##

tax = pd.read_csv(
    TAXONOMY,
    sep="\t",
    header=None,
    dtype=str
)

tax = tax.iloc[:, :2].copy()
tax.columns = [
    "genome",
    "taxonomy",
]


def get_rank(s, prefix):

    if pd.isna(s):
        return np.nan

    for item in str(s).split(";"):

        if item.startswith(prefix):
            return item[len(prefix):]

    return np.nan


tax["family"] = tax["taxonomy"].map(
    lambda z: get_rank(z, "f__")
)

tax["genus"] = tax["taxonomy"].map(
    lambda z: get_rank(z, "g__")
)

tax["species"] = tax["taxonomy"].map(
    lambda z: get_rank(z, "s__")
)


## ================================================================== ##
## Genome-level Cluster_00063 copy number
## ================================================================== ##

copies = (
    cluster_proteins
    .groupby("genome")
    .agg(
        n_cluster00063_copies=(
            "protein_id",
            "nunique"
        )
    )
    .reset_index()
)


## ================================================================== ##
## Presence in Modules 10 / 20 / 35
## ================================================================== ##

other_presence = (
    pm[
        pm["module"].isin(
            FOCAL_OTHER_MODULES
        )
    ]
    [["genome", "module"]]
    .drop_duplicates()
)

presence_wide = (
    other_presence
    .assign(present=1)
    .pivot(
        index="genome",
        columns="module",
        values="present"
    )
    .fillna(0)
    .reset_index()
)

for module in sorted(FOCAL_OTHER_MODULES):

    if module not in presence_wide.columns:
        presence_wide[module] = 0

presence_wide["n_other_focal_modules"] = (
    presence_wide[
        sorted(FOCAL_OTHER_MODULES)
    ]
    .sum(axis=1)
)


## ================================================================== ##
## Optional completeness values for Modules 10 / 20 / 35
## ================================================================== ##

completeness_file = (
    WF
    / "09_module_occurrence"
    / "genome_module_completeness_matrix.tsv"
)

completeness = pd.DataFrame(
    {"genome": copies["genome"]}
)

if completeness_file.exists():

    comp = pd.read_csv(
        completeness_file,
        sep="\t",
        low_memory=False
    )

    if "genome" in comp.columns:

        comp_cols = [
            m for m in sorted(FOCAL_OTHER_MODULES)
            if m in comp.columns
        ]

        if comp_cols:

            for c in comp_cols:
                comp[c] = as_num(comp[c])

            comp[
                "focal_module_completeness_sum"
            ] = (
                comp[comp_cols]
                .fillna(0)
                .sum(axis=1)
            )

            completeness = comp[
                [
                    "genome",
                    "focal_module_completeness_sum",
                ]
            ].copy()


if (
    "focal_module_completeness_sum"
    not in completeness.columns
):
    completeness[
        "focal_module_completeness_sum"
    ] = 0


## ================================================================== ##
## Candidate focal-locus quality
## ================================================================== ##

obs63 = obs[
    obs["focal_cluster"]
    == TARGET_CLUSTER
].copy()

for c in [
    "full_10kb_context_observable",
    "full_20kb_context_observable",
    "n_same_module_other_proteins_in_envelope",
]:

    obs63[c] = (
        as_num(obs63[c])
        .fillna(0)
    )


## Count highlighted local heme proteins per focal locus.

heme_counts = (
    genes[
        genes["display_class"]
        == "findmehemes_other"
    ]
    .groupby(
        [
            "genome",
            "focal_protein_id",
        ]
    )
    .agg(
        n_local_heme=(
            "neighbor_protein_id",
            "nunique"
        )
    )
    .reset_index()
)


module_counts = (
    genes[
        genes["display_class"]
        == "module06_other"
    ]
    .groupby(
        [
            "genome",
            "focal_protein_id",
        ]
    )
    .agg(
        n_local_module06=(
            "neighbor_protein_id",
            "nunique"
        )
    )
    .reset_index()
)


loci = (
    obs63
    .merge(
        heme_counts,
        on=[
            "genome",
            "focal_protein_id",
        ],
        how="left"
    )
    .merge(
        module_counts,
        on=[
            "genome",
            "focal_protein_id",
        ],
        how="left"
    )
)

loci[
    [
        "n_local_heme",
        "n_local_module06",
    ]
] = (
    loci[
        [
            "n_local_heme",
            "n_local_module06",
        ]
    ]
    .fillna(0)
)


## ================================================================== ##
## Pick best Cluster_00063 locus within each genome
##
## Priority:
## 1. complete ±20-kb context
## 2. complete ±10-kb context
## 3. more other heme proteins
## 4. more other Module_06 proteins
## ================================================================== ##

loci = loci.sort_values(
    [
        "genome",
        "full_20kb_context_observable",
        "full_10kb_context_observable",
        "n_local_heme",
        "n_local_module06",
        "n_genes_within_20kb_window",
        "focal_protein_id",
    ],
    ascending=[
        True,
        False,
        False,
        False,
        False,
        False,
        True,
    ],
)

best_locus = (
    loci
    .drop_duplicates(
        subset=["genome"],
        keep="first"
    )
    .copy()
)


## ================================================================== ##
## Genome-level candidate table
## ================================================================== ##

candidate = (
    copies
    .merge(
        tax[
            [
                "genome",
                "family",
                "genus",
                "species",
            ]
        ],
        on="genome",
        how="left"
    )
    .merge(
        presence_wide,
        on="genome",
        how="left"
    )
    .merge(
        completeness,
        on="genome",
        how="left"
    )
    .merge(
        best_locus[
            [
                "genome",
                "focal_protein_id",
                "contig",
                "full_10kb_context_observable",
                "full_20kb_context_observable",
                "n_local_heme",
                "n_local_module06",
            ]
        ],
        on="genome",
        how="left"
    )
)


for module in sorted(FOCAL_OTHER_MODULES):

    candidate[module] = (
        as_num(candidate[module])
        .fillna(0)
        .astype(int)
    )


candidate[
    "n_other_focal_modules"
] = (
    candidate[
        sorted(FOCAL_OTHER_MODULES)
    ]
    .sum(axis=1)
)

candidate[
    "focal_module_completeness_sum"
] = (
    as_num(
        candidate[
            "focal_module_completeness_sum"
        ]
    )
    .fillna(0)
)

candidate[
    "best_full_10kb"
] = (
    as_num(
        candidate[
            "full_10kb_context_observable"
        ]
    )
    .fillna(0)
)

candidate[
    "best_full_20kb"
] = (
    as_num(
        candidate[
            "full_20kb_context_observable"
        ]
    )
    .fillna(0)
)

candidate[
    "best_local_heme"
] = (
    as_num(
        candidate["n_local_heme"]
    )
    .fillna(0)
)

candidate[
    "best_local_module06"
] = (
    as_num(
        candidate["n_local_module06"]
    )
    .fillna(0)
)


candidate[
    "selection_group"
] = (
    candidate["genus"]
    .map(rank_group)
)


## ================================================================== ##
## Rank within requested biological groups
## ================================================================== ##

mc = rank_candidates(
    candidate[
        candidate["selection_group"]
        == "Methylococcus"
    ]
).copy()

mb = rank_candidates(
    candidate[
        candidate["selection_group"]
        == "Methylobacter"
    ]
).copy()

other = rank_candidates(
    candidate[
        candidate["selection_group"]
        == "Other"
    ]
).copy()


selected_mc = mc.head(
    N_METHYLOCOCCUS
)

selected_mb = mb.head(
    N_METHYLOBACTER
)

selected_other = other.head(
    N_OTHER
)


selected = pd.concat(
    [
        selected_mc,
        selected_mb,
        selected_other,
    ],
    ignore_index=True
)


if len(selected) != 12:

    raise RuntimeError(
        f"Expected 12 selected genomes, "
        f"found {len(selected)}."
    )


## ================================================================== ##
## Figure order
##
## Methylococcus first, then Methylobacter, then other genera.
## Within groups retain ranking.
## ================================================================== ##

selected[
    "selection_order"
] = range(
    1,
    len(selected) + 1
)


selected[
    "selected_for_figure"
] = 1


candidate = candidate.merge(
    selected[
        [
            "genome",
            "selection_order",
            "selected_for_figure",
        ]
    ],
    on="genome",
    how="left"
)

candidate[
    "selected_for_figure"
] = (
    candidate[
        "selected_for_figure"
    ]
    .fillna(0)
    .astype(int)
)


## ================================================================== ##
## Selected gene-map rows
## ================================================================== ##

selected_key = selected[
    [
        "genome",
        "focal_protein_id",
        "selection_order",
        "family",
        "genus",
        "species",
    ]
].copy()


selected_genes = genes.merge(
    selected_key,
    on=[
        "genome",
        "focal_protein_id",
    ],
    how="inner",
    suffixes=(
        "",
        "_selected"
    )
)


## ================================================================== ##
## Orient every locus so Cluster_00063 points to the right
##
## Coordinates are relative to focal midpoint.
## ================================================================== ##

for c in [
    "neighbor_start",
    "neighbor_end",
    "focal_start",
    "focal_end",
]:

    selected_genes[c] = as_num(
        selected_genes[c]
    )


focal_mid = (
    selected_genes["focal_start"]
    +
    selected_genes["focal_end"]
) / 2


forward_focal = (
    selected_genes["focal_strand"]
    == "+"
)


selected_genes[
    "plot_start_bp"
] = np.where(
    forward_focal,
    selected_genes["neighbor_start"]
        - focal_mid,
    focal_mid
        - selected_genes["neighbor_end"],
)


selected_genes[
    "plot_end_bp"
] = np.where(
    forward_focal,
    selected_genes["neighbor_end"]
        - focal_mid,
    focal_mid
        - selected_genes["neighbor_start"],
)


selected_genes[
    "plot_strand"
] = np.where(
    forward_focal,
    selected_genes["neighbor_strand"],
    np.where(
        selected_genes["neighbor_strand"]
        == "+",
        "-",
        "+"
    )
)


## Keep the same plotted scale as previous figures:
## ±20 kb around the focal protein.

selected_genes = selected_genes[
    (
        selected_genes["plot_end_bp"]
        >= -20000
    )
    &
    (
        selected_genes["plot_start_bp"]
        <= 20000
    )
].copy()


## ================================================================== ##
## Taxonomy display label
## ================================================================== ##

def taxonomy_display(row):

    genus = row.get("genus", "")
    species = row.get("species", "")

    genus = (
        ""
        if pd.isna(genus)
        else str(genus)
    )

    species = (
        ""
        if pd.isna(species)
        else str(species)
    )

    if species:
        tax_label = species
    elif genus:
        tax_label = genus
    else:
        tax_label = row["genome"]

    return (
        f"{tax_label} | "
        f"{row['genome']}"
    )


selected_genes[
    "taxonomy_display"
] = selected_genes.apply(
    taxonomy_display,
    axis=1
)

selected_genes[
    "taxonomy_family"
] = selected_genes[
    "family"
]

selected_genes[
    "taxonomy_genus"
] = selected_genes[
    "genus"
]

selected_genes[
    "taxonomy_species"
] = selected_genes[
    "species"
]

selected_genes[
    "protein_id"
] = selected_genes[
    "neighbor_protein_id"
]

selected_genes[
    "neighbor_mmseq_cluster"
] = selected_genes[
    "neighbor_cluster"
]


## Unique row/region ID
selected_genes[
    "plot_region_id"
] = (
    selected_genes["genome"]
    + "|"
    + selected_genes["focal_protein_id"]
)


## ================================================================== ##
## Save
## ================================================================== ##

candidate_out = (
    OUTDIR
    / "Module_06_representative_candidate_table.tsv"
)

selected_out = (
    OUTDIR
    / "Module_06_selected_regions.tsv"
)

gene_out = (
    OUTDIR
    / "Module_06_selected_gene_map_data.tsv"
)


candidate.sort_values(
    [
        "selection_group",
        "n_other_focal_modules",
        "focal_module_completeness_sum",
    ],
    ascending=[
        True,
        False,
        False,
    ],
).to_csv(
    candidate_out,
    sep="\t",
    index=False
)


selected.to_csv(
    selected_out,
    sep="\t",
    index=False
)


selected_genes.to_csv(
    gene_out,
    sep="\t",
    index=False
)


## ================================================================== ##
## Report
## ================================================================== ##

print()
print("=" * 120)
print("MODULE_06 — 12 GENOMES SELECTED FOR GENE MAP")
print("=" * 120)

show_cols = [
    "selection_order",
    "genome",
    "genus",
    "species",
    "n_cluster00063_copies",
    "Module_10",
    "Module_20",
    "Module_35",
    "n_other_focal_modules",
    "focal_module_completeness_sum",
    "focal_protein_id",
    "best_full_20kb",
    "best_local_heme",
    "best_local_module06",
]

print(
    selected[
        show_cols
    ].to_string(
        index=False
    )
)


print()
print("=" * 120)
print("SELECTION COMPOSITION")
print("=" * 120)

print(
    selected[
        "selection_group"
    ]
    .value_counts()
    .to_string()
)


print()
print("=" * 120)
print("SELECTED GENE-MAP DATA")
print("=" * 120)

print(
    "Regions:",
    selected_genes[
        "plot_region_id"
    ].nunique()
)

print(
    "Genes:",
    len(selected_genes)
)

print()
print(
    selected_genes[
        "display_class"
    ]
    .value_counts()
    .to_string()
)


print()
print("Written:")
print(candidate_out)
print(selected_out)
print(gene_out)

