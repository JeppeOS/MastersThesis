#!/usr/bin/env python3

from pathlib import Path
import pandas as pd


## ================================================================ ##
## Paths
## ================================================================ ##

PROJECT = (
    Path.home()
    / "methanotrophs"
    / "methanotroph_project"
    / "jeppe"
)

WF = PROJECT / "genome_analysis_workflow"
VIZ = WF / "16_visualization"

OUTDIR = VIZ / "16A_module06_gene_map_framework"
OUTDIR.mkdir(parents=True, exist_ok=True)

TARGET_MODULE = "Module_06"
TARGET_CLUSTER = "Cluster_00063"

FILES = {
    "cluster_membership":
        WF / "06_clustering" / "cluster_membership.tsv",

    "protein_module_membership":
        WF / "09_module_occurrence" / "protein_module_membership.tsv",

    "taxonomy":
        WF / "16_visualization" / "globdb_r226_taxonomy.tsv",

    "module_protein_coordinates":
        WF / "10_genomic_coordinates" / "module_protein_coordinates.tsv",

    "focal_module_occurrences":
        WF / "12_neighborhood_extraction" / "focal_module_occurrences.tsv",

    "observed_neighborhood_genes":
        WF / "12_neighborhood_extraction" / "observed_neighborhood_genes.tsv",

    "focal_context_observability":
        WF / "12_neighborhood_extraction" / "focal_context_observability.tsv",

    "module_member_pair_observations":
        WF / "12_neighborhood_extraction" / "module_member_pair_observations.tsv",
}


## ================================================================ ##
## Helpers
## ================================================================ ##

def read_tsv(path):
    return pd.read_csv(path, sep="\t", dtype=str)

def detect_col(df, candidates):
    for c in candidates:
        if c in df.columns:
            return c
    return None

def add_taxonomy(df, genome_col="genome"):
    tax = read_tsv(FILES["taxonomy"])

    if tax.shape[1] < 2:
        raise RuntimeError("taxonomy file does not have >=2 columns")

    tax = tax.iloc[:, :2].copy()
    tax.columns = ["genome", "taxonomy"]

    def get_rank(taxonomy, prefix):
        if pd.isna(taxonomy):
            return None
        for part in taxonomy.split(";"):
            if part.startswith(prefix):
                return part[len(prefix):]
        return None

    tax["family"] = tax["taxonomy"].apply(lambda z: get_rank(z, "f__"))
    tax["genus"] = tax["taxonomy"].apply(lambda z: get_rank(z, "g__"))
    tax["species"] = tax["taxonomy"].apply(lambda z: get_rank(z, "s__"))

    if genome_col not in df.columns:
        return df

    return df.merge(
        tax[["genome", "family", "genus", "species", "taxonomy"]],
        left_on=genome_col,
        right_on="genome",
        how="left"
    )

def save(df, name):
    path = OUTDIR / name
    df.to_csv(path, sep="\t", index=False)
    print(path)
    return path

def filter_by_values(df, values_by_column):
    mask = pd.Series(False, index=df.index)

    for col, values in values_by_column.items():
        if col in df.columns and values:
            mask |= df[col].isin(values)

    return df[mask].copy()

def occurrence_columns(df):
    return [c for c in df.columns if "occurrence" in c.lower()]

def possible_protein_columns(df):
    return [
        c for c in df.columns
        if (
            "protein" in c.lower()
            or c.lower() in {"member", "representative_protein"}
        )
    ]


## ================================================================ ##
## Load primary tables
## ================================================================ ##

cm = read_tsv(FILES["cluster_membership"])
pm = read_tsv(FILES["protein_module_membership"])

print("=" * 100)
print("PRIMARY INPUTS")
print("=" * 100)
print("cluster_membership:", FILES["cluster_membership"])
print("protein_module_membership:", FILES["protein_module_membership"])
print()

print("cluster_membership columns:")
print(cm.columns.tolist())
print()

print("protein_module_membership columns:")
print(pm.columns.tolist())
print()


## ================================================================ ##
## Standardize key columns
## ================================================================ ##

cm_cluster_col = detect_col(cm, ["cluster"])
cm_genome_col = detect_col(cm, ["genome"])
cm_protein_col = detect_col(cm, ["protein_id"])
cm_length_col = detect_col(cm, ["length"])
cm_heme_col = detect_col(cm, ["number_of_hemes"])

pm_module_col = detect_col(pm, ["module"])
pm_cluster_col = detect_col(pm, ["cluster"])
pm_genome_col = detect_col(pm, ["genome"])
pm_protein_col = detect_col(pm, ["protein_id"])


## ================================================================ ##
## Cluster_00063 proteins
## ================================================================ ##

cluster63 = cm[
    cm[cm_cluster_col] == TARGET_CLUSTER
].copy()

cluster63 = cluster63.drop_duplicates(subset=[cm_protein_col])

cluster63 = add_taxonomy(cluster63, genome_col=cm_genome_col)

print("=" * 100)
print(f"{TARGET_CLUSTER} BASIC")
print("=" * 100)
print("Proteins:", cluster63[cm_protein_col].nunique())
print("Genomes:", cluster63[cm_genome_col].nunique())
print()


## ================================================================ ##
## Module_06 proteins
## ================================================================ ##

module06 = pm[
    pm[pm_module_col] == TARGET_MODULE
].copy()

module06 = module06.drop_duplicates(subset=[pm_protein_col])

module06 = add_taxonomy(module06, genome_col=pm_genome_col)

print("=" * 100)
print(f"{TARGET_MODULE} BASIC")
print("=" * 100)
print("Proteins:", module06[pm_protein_col].nunique())
print("Genomes:", module06[pm_genome_col].nunique())
print()


## ================================================================ ##
## Save direct protein tables
## ================================================================ ##

save(cluster63, "cluster00063_proteins.tsv")
save(module06, "module06_proteins.tsv")


## ================================================================ ##
## Cluster_00063 summary
## ================================================================ ##

cluster63_num = cluster63.copy()
if cm_length_col in cluster63_num.columns:
    cluster63_num[cm_length_col] = pd.to_numeric(cluster63_num[cm_length_col], errors="coerce")
if cm_heme_col in cluster63_num.columns:
    cluster63_num[cm_heme_col] = pd.to_numeric(cluster63_num[cm_heme_col], errors="coerce")

summary = pd.DataFrame({
    "cluster": [TARGET_CLUSTER],
    "module": [TARGET_MODULE],
    "n_proteins": [cluster63[cm_protein_col].nunique()],
    "n_genomes": [cluster63[cm_genome_col].nunique()],
    "min_length": [cluster63_num[cm_length_col].min() if cm_length_col else None],
    "max_length": [cluster63_num[cm_length_col].max() if cm_length_col else None],
    "mean_length": [cluster63_num[cm_length_col].mean() if cm_length_col else None],
    "median_length": [cluster63_num[cm_length_col].median() if cm_length_col else None],
    "min_hemes": [cluster63_num[cm_heme_col].min() if cm_heme_col else None],
    "max_hemes": [cluster63_num[cm_heme_col].max() if cm_heme_col else None],
    "median_hemes": [cluster63_num[cm_heme_col].median() if cm_heme_col else None],
})

save(summary, "cluster00063_summary.tsv")

print("=" * 100)
print("CLUSTER_00063 SUMMARY")
print("=" * 100)
print(summary.to_string(index=False))
print()


## ================================================================ ##
## Module_06 cluster composition
## ================================================================ ##

modcomp = (
    module06
    .groupby(pm_cluster_col, dropna=False)
    .agg(
        n_proteins=(pm_protein_col, "nunique"),
        n_genomes=(pm_genome_col, "nunique"),
    )
    .reset_index()
    .sort_values(["n_proteins", "n_genomes"], ascending=[False, False])
)

save(modcomp, "module06_cluster_composition.tsv")

print("=" * 100)
print("MODULE_06 CLUSTER COMPOSITION")
print("=" * 100)
print(modcomp.to_string(index=False))
print()


## ================================================================ ##
## Cluster_00063 genome table
## ================================================================ ##

genome_summary = (
    cluster63
    .groupby(
        [cm_genome_col, "family", "genus", "species"],
        dropna=False
    )
    .agg(
        n_cluster00063_proteins=(cm_protein_col, "nunique"),
        protein_ids=(cm_protein_col, lambda z: ";".join(sorted(set(z.astype(str))))),
        lengths=(cm_length_col, lambda z: ";".join(str(x) for x in z if pd.notna(x))) if cm_length_col else (cm_protein_col, lambda z: ""),
        heme_counts=(cm_heme_col, lambda z: ";".join(str(x) for x in z if pd.notna(x))) if cm_heme_col else (cm_protein_col, lambda z: ""),
    )
    .reset_index()
    .sort_values(["genus", cm_genome_col], na_position="last")
)

save(genome_summary, "cluster00063_genome_summary.tsv")

print("=" * 100)
print("GENOMES CARRYING CLUSTER_00063")
print("=" * 100)
print(genome_summary.to_string(index=False))
print()


## ================================================================ ##
## Representative selection candidates
## ================================================================ ##

# Add total number of Module_06 proteins per genome
module06_copy = (
    module06
    .groupby(pm_genome_col, dropna=False)
    .agg(
        n_module06_proteins=(pm_protein_col, "nunique")
    )
    .reset_index()
)

rep_candidates = genome_summary.merge(
    module06_copy,
    left_on=cm_genome_col,
    right_on=pm_genome_col,
    how="left"
)

if pm_genome_col != cm_genome_col and pm_genome_col in rep_candidates.columns:
    rep_candidates = rep_candidates.drop(columns=[pm_genome_col])

rep_candidates = rep_candidates.sort_values(
    ["genus", "n_cluster00063_proteins", "n_module06_proteins", cm_genome_col],
    ascending=[True, False, False, True],
    na_position="last"
)

save(rep_candidates, "module06_representative_selection_candidates.tsv")


## ================================================================ ##
## Stage 10 / 12 extraction
## ================================================================ ##

cluster63_proteins = set(cluster63[cm_protein_col].dropna())
module06_proteins = set(module06[pm_protein_col].dropna())
module06_genomes = set(module06[pm_genome_col].dropna())
cluster63_genomes = set(cluster63[cm_genome_col].dropna())

print("=" * 100)
print("STAGE 10 / 12 FILTERING")
print("=" * 100)

## --------------------------------------------------------------- ##
## 10_genomic_coordinates/module_protein_coordinates.tsv
## --------------------------------------------------------------- ##

coords = read_tsv(FILES["module_protein_coordinates"])
print("module_protein_coordinates columns:")
print(coords.columns.tolist())
print()

coord_values = {}
for col in coords.columns:
    low = col.lower()

    if low == "module":
        coord_values[col] = {TARGET_MODULE}
    elif low == "cluster":
        coord_values[col] = {TARGET_CLUSTER}
    elif low == "genome":
        coord_values[col] = module06_genomes | cluster63_genomes
    elif "protein" in low or low in {"member", "representative_protein"}:
        coord_values[col] = module06_proteins | cluster63_proteins

coords_subset = filter_by_values(coords, coord_values)
coords_subset = add_taxonomy(coords_subset, genome_col=detect_col(coords_subset, ["genome"]))
save(coords_subset, "module06_coordinates.tsv")


## --------------------------------------------------------------- ##
## 12_neighborhood_extraction/focal_module_occurrences.tsv
## --------------------------------------------------------------- ##

focal_occ = read_tsv(FILES["focal_module_occurrences"])
print("focal_module_occurrences columns:")
print(focal_occ.columns.tolist())
print()

focal_values = {}
for col in focal_occ.columns:
    low = col.lower()

    if low == "module":
        focal_values[col] = {TARGET_MODULE}
    elif low == "cluster":
        focal_values[col] = {TARGET_CLUSTER}
    elif low == "genome":
        focal_values[col] = module06_genomes | cluster63_genomes
    elif "protein" in low or low in {"member", "representative_protein"}:
        focal_values[col] = module06_proteins | cluster63_proteins

focal_occ_subset = filter_by_values(focal_occ, focal_values)
focal_occ_subset = add_taxonomy(focal_occ_subset, genome_col=detect_col(focal_occ_subset, ["genome"]))
save(focal_occ_subset, "module06_focal_module_occurrences.tsv")

occ_cols = occurrence_columns(focal_occ_subset)
occ_ids = set()
for c in occ_cols:
    occ_ids |= set(focal_occ_subset[c].dropna())

print("Occurrence columns in focal_module_occurrences subset:", occ_cols)
print("Unique occurrence IDs captured:", len(occ_ids))
print()


## --------------------------------------------------------------- ##
## 12_neighborhood_extraction/observed_neighborhood_genes.tsv
## --------------------------------------------------------------- ##

obs = read_tsv(FILES["observed_neighborhood_genes"])
print("observed_neighborhood_genes columns:")
print(obs.columns.tolist())
print()

obs_values = {}
for col in obs.columns:
    low = col.lower()

    if low == "module":
        obs_values[col] = {TARGET_MODULE}
    elif low == "cluster":
        obs_values[col] = {TARGET_CLUSTER}
    elif low == "genome":
        obs_values[col] = module06_genomes | cluster63_genomes
    elif "occurrence" in low:
        obs_values[col] = occ_ids
    elif "protein" in low or low in {"member", "representative_protein"}:
        obs_values[col] = module06_proteins | cluster63_proteins

obs_subset = filter_by_values(obs, obs_values)
obs_subset = add_taxonomy(obs_subset, genome_col=detect_col(obs_subset, ["genome"]))
save(obs_subset, "module06_observed_neighborhood_genes.tsv")


## --------------------------------------------------------------- ##
## Cluster_00063-only neighborhood subset
## --------------------------------------------------------------- ##

cluster63_occ_ids = set()

for c in occ_cols:
    if c in focal_occ_subset.columns:
        sub = focal_occ_subset.copy()

        protein_hits = pd.Series(False, index=sub.index)
        for pcol in possible_protein_columns(sub):
            protein_hits |= sub[pcol].isin(cluster63_proteins)

        cluster_hits = pd.Series(False, index=sub.index)
        if "cluster" in sub.columns:
            cluster_hits |= sub["cluster"].eq(TARGET_CLUSTER)

        keep = protein_hits | cluster_hits

        if keep.any():
            cluster63_occ_ids |= set(sub.loc[keep, c].dropna())

cluster63_obs_values = {}
for col in obs.columns:
    low = col.lower()

    if low == "cluster":
        cluster63_obs_values[col] = {TARGET_CLUSTER}
    elif low == "genome":
        cluster63_obs_values[col] = cluster63_genomes
    elif "occurrence" in low:
        cluster63_obs_values[col] = cluster63_occ_ids
    elif "protein" in low or low in {"member", "representative_protein"}:
        cluster63_obs_values[col] = cluster63_proteins

cluster63_obs = filter_by_values(obs, cluster63_obs_values)
cluster63_obs = add_taxonomy(cluster63_obs, genome_col=detect_col(cluster63_obs, ["genome"]))
save(cluster63_obs, "cluster00063_observed_neighborhood_genes.tsv")


## --------------------------------------------------------------- ##
## 12_neighborhood_extraction/focal_context_observability.tsv
## --------------------------------------------------------------- ##

ctx = read_tsv(FILES["focal_context_observability"])
print("focal_context_observability columns:")
print(ctx.columns.tolist())
print()

ctx_values = {}
for col in ctx.columns:
    low = col.lower()

    if low == "module":
        ctx_values[col] = {TARGET_MODULE}
    elif low == "cluster":
        ctx_values[col] = {TARGET_CLUSTER}
    elif low == "genome":
        ctx_values[col] = module06_genomes | cluster63_genomes
    elif "occurrence" in low:
        ctx_values[col] = occ_ids
    elif "protein" in low or low in {"member", "representative_protein"}:
        ctx_values[col] = module06_proteins | cluster63_proteins

ctx_subset = filter_by_values(ctx, ctx_values)
ctx_subset = add_taxonomy(ctx_subset, genome_col=detect_col(ctx_subset, ["genome"]))
save(ctx_subset, "module06_focal_context_observability.tsv")


## --------------------------------------------------------------- ##
## 12_neighborhood_extraction/module_member_pair_observations.tsv
## --------------------------------------------------------------- ##

pair = read_tsv(FILES["module_member_pair_observations"])
print("module_member_pair_observations columns:")
print(pair.columns.tolist())
print()

pair_values = {}
for col in pair.columns:
    low = col.lower()

    if low == "module":
        pair_values[col] = {TARGET_MODULE}
    elif low == "cluster":
        pair_values[col] = {TARGET_CLUSTER}
    elif low == "genome":
        pair_values[col] = module06_genomes | cluster63_genomes
    elif "occurrence" in low:
        pair_values[col] = occ_ids
    elif "protein" in low or low in {"member", "representative_protein"}:
        pair_values[col] = module06_proteins | cluster63_proteins

pair_subset = filter_by_values(pair, pair_values)
pair_subset = add_taxonomy(pair_subset, genome_col=detect_col(pair_subset, ["genome"]))
save(pair_subset, "module06_module_member_pair_observations.tsv")


## ================================================================ ##
## Final report
## ================================================================ ##

print()
print("=" * 100)
print("OUTPUT FILES WRITTEN")
print("=" * 100)

for p in sorted(OUTDIR.glob("*.tsv")):
    print(p)

print()
print("Done.")
