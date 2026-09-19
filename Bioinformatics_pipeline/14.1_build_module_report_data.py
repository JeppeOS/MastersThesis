#!/usr/bin/env python3

from pathlib import Path
from collections import Counter
import sys

import numpy as np
import pandas as pd


## ================================================================== ##
## STAGE 14A - BUILD INTEGRATED MODULE REPORT DATA
##
## PURPOSE
## -------
##
## Assemble all authoritative upstream evidence into one reproducible
## report-data package for each MCL-defined module.
##
##
## IMPORTANT HIERARCHY
## -------------------
##
## Module
##   -> contains MMseqs2 sequence families
##
## MMseqs2 family
##   -> remains the primary protein-family identity
##
## Neighborhood evidence
##   -> Stage 13A-D
##
## Functional / protein evidence
##   -> FeGenie
##   -> FindMeHemes
##   -> SignalP
##   -> DeepTMHMM
##   -> heme counts
##   -> GlobDB COG/product
##
##
## This stage DOES NOT:
##
##   * assume a module is an operon
##   * force one linear module architecture
##   * impose a conserved/not-conserved percentage cutoff
##   * replace MMseqs2 family identity with annotation labels
##
##
## Instead it builds the evidence package needed for the later
## module-level architecture/report builder.
##
##
## PERCENTAGE RULE
## ---------------
##
## Existing numerator / denominator / support fields are preserved.
##
## A value such as:
##
##     100% (1/1)
##
## remains visibly weaker support than:
##
##     100% (44/44)
##
##
## FEGENIE / PROTEIN EVIDENCE
## --------------------------
##
## Protein evidence is recovered from the Stage-12 observed-neighborhood
## catalogue, where protein-level annotations from the integrated
## annotation pipeline were already propagated.
##
## This lets Stage 14 integrate FeGenie etc. WITHOUT changing the
## Stage-13 neighborhood definitions.
## ================================================================== ##


## ================================================================== ##
## Paths
## ================================================================== ##

HERE = Path(__file__).resolve().parent
WORKFLOW = HERE.parent.parent


## Core sequence-family / module definitions ##

MODULE_MEMBERSHIP = (
    WORKFLOW
    / "07_cooccurrence_network"
    / "module_membership.tsv"
)


CLUSTER_SUMMARY = (
    WORKFLOW
    / "06_clustering"
    / "cluster_summary.tsv"
)


## Module occurrence / completeness ##

PROTEIN_MODULE_MEMBERSHIP = (
    WORKFLOW
    / "09_module_occurrence"
    / "protein_module_membership.tsv"
)


GENOME_MODULE_SUMMARY_ALL = (
    WORKFLOW
    / "09_module_occurrence"
    / "genome_module_summary_all.tsv"
)


MODULE_COMPLETENESS_SUMMARY = (
    WORKFLOW
    / "09_module_occurrence"
    / "module_completeness_summary.tsv"
)


## Stage 12 protein evidence / genomic catalogue ##

OBSERVED_NEIGHBORHOODS = (
    WORKFLOW
    / "12_neighborhood_extraction"
    / "observed_neighborhood_genes.tsv"
)


## Stage 13A ##

PAIR_COLOCALIZATION = (
    WORKFLOW
    / "13_neighborhood_conservation"
    / "13A_module_colocalization"
    / "module_cluster_pair_colocalization.tsv"
)


MODULE_COLOCALIZATION = (
    WORKFLOW
    / "13_neighborhood_conservation"
    / "13A_module_colocalization"
    / "module_colocalization_summary.tsv"
)


## Stage 13C2 ##

LOCAL_INTERVENING_GENES = (
    WORKFLOW
    / "13_neighborhood_conservation"
    / "13C_annotation_context"
    / "local_module_pair_intervening_genes.tsv"
)


LOCAL_INTERVENING_PATTERNS = (
    WORKFLOW
    / "13_neighborhood_conservation"
    / "13C_annotation_context"
    / "local_module_pair_intervening_pattern_summary.tsv"
)


LOCAL_PAIR_SUMMARY = (
    WORKFLOW
    / "13_neighborhood_conservation"
    / "13C_annotation_context"
    / "local_module_pair_intervening_pair_summary.tsv"
)


## Stage 13D ##

FOCAL_POSITION_CONSENSUS = (
    WORKFLOW
    / "13_neighborhood_conservation"
    / "13D_consensus_neighborhoods"
    / "focal_consensus_position_summary.tsv"
)


MMSEQ_CONSENSUS = (
    WORKFLOW
    / "13_neighborhood_conservation"
    / "13D_consensus_neighborhoods"
    / "focal_mmseqs_consensus_associations.tsv"
)


COG_CONSENSUS = (
    WORKFLOW
    / "13_neighborhood_conservation"
    / "13D_consensus_neighborhoods"
    / "focal_cog_consensus_associations.tsv"
)


## ================================================================== ##
## Output structure
## ================================================================== ##

GLOBAL_DIR = (
    HERE
    / "global"
)


MODULE_DIR = (
    HERE
    / "modules"
)


GLOBAL_DIR.mkdir(
    parents=True,
    exist_ok=True,
)


MODULE_DIR.mkdir(
    parents=True,
    exist_ok=True,
)


OUT_OVERVIEW = (
    GLOBAL_DIR
    / "module_report_overview.tsv"
)


OUT_CLUSTER_EVIDENCE = (
    GLOBAL_DIR
    / "module_member_cluster_evidence.tsv"
)


OUT_PROTEIN_EVIDENCE = (
    GLOBAL_DIR
    / "module_protein_evidence.tsv"
)


OUT_INTERVENING_EVIDENCE = (
    GLOBAL_DIR
    / "local_intervening_gene_integrated_evidence.tsv"
)


OUT_INDEX = (
    GLOBAL_DIR
    / "module_report_index.tsv"
)


OUT_MANIFEST = (
    GLOBAL_DIR
    / "module_report_manifest.tsv"
)


OUT_QC = (
    GLOBAL_DIR
    / "module_report_data_qc.tsv"
)


## ================================================================== ##
## Expected authoritative upstream values
## ================================================================== ##

EXPECTED_MODULES = 35
EXPECTED_MODULE_CLUSTERS = 156
EXPECTED_MODULE_PROTEINS = 10_537

EXPECTED_POSITION_ROWS = 6_396
EXPECTED_MMSEQ_CONSENSUS_ROWS = 3_050
EXPECTED_COG_CONSENSUS_ROWS = 107_728

EXPECTED_LOCAL_INTERVENING_GENES = 6_950

EXPECTED_PAIR_OBSERVATIONS_STAGE12 = 34_037
EXPECTED_LOCAL_WITHIN_20KB = 2_757


## ================================================================== ##
## Utilities
## ================================================================== ##

def fail(message):

    print(
        f"\nERROR: {message}",
        file=sys.stderr,
    )

    sys.exit(1)


def require_file(path):

    if not path.exists():

        fail(
            f"Required file does not exist:\n"
            f"{path}"
        )


def read_tsv(path):

    require_file(
        path
    )

    return pd.read_csv(
        path,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )


def require_columns(
    df,
    columns,
    source,
):

    missing = (
        set(columns)
        -
        set(df.columns)
    )

    if missing:

        fail(
            f"{source} missing required column(s): "
            +
            ", ".join(
                sorted(missing)
            )
        )


def as_numeric(
    series,
):

    return pd.to_numeric(
        series.replace(
            "",
            np.nan,
        ),
        errors="coerce",
    )


def pct(
    numerator,
    denominator,
):

    if denominator == 0:

        return np.nan

    return (
        100.0
        *
        numerator
        /
        denominator
    )


def distribution_string(
    values,
):

    cleaned = [
        str(value).strip()

        for value in values

        if str(value).strip() != ""
    ]


    if len(
        cleaned
    ) == 0:

        return ""


    counts = Counter(
        cleaned
    )


    ordered = sorted(
        counts.items(),
        key=lambda x: (
            -x[1],
            x[0],
        ),
    )


    return "; ".join(
        f"{value}:{count}"

        for value, count
        in ordered
    )


def prefix_nonkey_columns(
    df,
    key,
    prefix,
):

    rename = {
        column:
            (
                column

                if column == key

                else
                prefix + column
            )

        for column
        in df.columns
    }


    return df.rename(
        columns=rename
    )


def ensure_one_row_per(
    df,
    key,
    source,
):

    if df[
        key
    ].duplicated().any():

        fail(
            f"{source} does not contain "
            f"exactly one row per {key}."
        )


def subset_module(
    df,
    module,
):

    if "module" in df.columns:

        return df[
            df[
                "module"
            ]
            ==
            module
        ].copy()


    if "focal_module" in df.columns:

        return df[
            df[
                "focal_module"
            ]
            ==
            module
        ].copy()


    fail(
        "Cannot determine module column "
        "for per-module table."
    )


def write_table(
    df,
    path,
):

    df.to_csv(
        path,
        sep="\t",
        index=False,
        na_rep="",
        float_format="%.6f",
    )


## ================================================================== ##
## Start
## ================================================================== ##

print("=" * 80)
print("STAGE 14A - BUILD INTEGRATED MODULE REPORT DATA")
print("=" * 80)


## ================================================================== ##
## 1. Read authoritative module membership
## ================================================================== ##

print()
print("Reading authoritative MCL module membership...")


membership = read_tsv(
    MODULE_MEMBERSHIP
)


require_columns(
    membership,
    [
        "cluster",
        "module",
    ],
    MODULE_MEMBERSHIP.name,
)


if membership[
    "cluster"
].duplicated().any():

    fail(
        "A MMseqs2 cluster occurs in more than "
        "one MCL module."
    )


if len(
    membership
) != EXPECTED_MODULE_CLUSTERS:

    fail(
        f"Expected "
        f"{EXPECTED_MODULE_CLUSTERS} "
        f"module clusters but found "
        f"{len(membership)}."
    )


if membership[
    "module"
].nunique() != EXPECTED_MODULES:

    fail(
        f"Expected "
        f"{EXPECTED_MODULES} "
        f"modules but found "
        f"{membership['module'].nunique()}."
    )


modules = sorted(
    membership[
        "module"
    ].unique()
)


print(
    f"  Modules:          "
    f"{len(modules):,}"
)

print(
    f"  Member clusters:  "
    f"{len(membership):,}"
)


## ================================================================== ##
## 2. Read cluster summary
## ================================================================== ##

print()
print("Reading MMseqs2 cluster summary...")


cluster_summary = read_tsv(
    CLUSTER_SUMMARY
)


require_columns(
    cluster_summary,
    [
        "cluster",
    ],
    CLUSTER_SUMMARY.name,
)


if cluster_summary[
    "cluster"
].duplicated().any():

    fail(
        "cluster_summary.tsv contains "
        "duplicate cluster IDs."
    )


member_clusters = membership.merge(
    cluster_summary,
    on="cluster",
    how="left",
    validate="one_to_one",
    indicator=True,
)


if (
    member_clusters[
        "_merge"
    ]
    !=
    "both"
).any():

    fail(
        "At least one MCL cluster is missing "
        "from cluster_summary.tsv."
    )


member_clusters = member_clusters.drop(
    columns="_merge"
)


## ================================================================== ##
## 3. Read module protein membership
## ================================================================== ##

print()
print("Reading module protein membership...")


protein_membership = read_tsv(
    PROTEIN_MODULE_MEMBERSHIP
)


require_columns(
    protein_membership,
    [
        "genome",
        "protein_id",
        "cluster",
        "module",
    ],
    PROTEIN_MODULE_MEMBERSHIP.name,
)


if len(
    protein_membership
) != EXPECTED_MODULE_PROTEINS:

    fail(
        f"Expected "
        f"{EXPECTED_MODULE_PROTEINS:,} "
        f"module proteins but found "
        f"{len(protein_membership):,}."
    )


if protein_membership[
    [
        "genome",
        "protein_id",
    ]
].duplicated().any():

    fail(
        "Duplicate genome + protein_id "
        "module-protein keys."
    )


## Validate protein cluster -> module mapping. ##

protein_membership_check = protein_membership.merge(
    membership.rename(
        columns={
            "module":
                "expected_module",
        }
    ),
    on="cluster",
    how="left",
    validate="many_to_one",
)


if (
    protein_membership_check[
        "module"
    ]
    !=
    protein_membership_check[
        "expected_module"
    ]
).any():

    fail(
        "Protein module assignments disagree "
        "with authoritative cluster -> module membership."
    )


print(
    f"  Module proteins:  "
    f"{len(protein_membership):,}"
)

print(
    f"  Genomes:          "
    f"{protein_membership['genome'].nunique():,}"
)


## ================================================================== ##
## 4. Recover invariant protein evidence from Stage 12
##
## Stage-12 neighborhood rows repeat the same neighbor protein in
## multiple focal neighborhoods.
##
## We retain ONLY invariant protein properties here.
## ================================================================== ##

print()
print(
    "Recovering integrated protein evidence "
    "from Stage-12 neighborhood catalogue..."
)


neigh = read_tsv(
    OBSERVED_NEIGHBORHOODS
)


require_columns(
    neigh,
    [
        "genome",
        "neighbor_protein_id",

        "neighbor_start",
        "neighbor_end",
        "neighbor_strand",
        "neighbor_gene_rank",

        "neighbor_annotation_match_type",
        "neighbor_annotation_accepted",

        "neighbor_globdb_cog",
        "neighbor_globdb_gene",
        "neighbor_globdb_product",

        "neighbor_cluster",
        "neighbor_module",

        "neighbor_candidate_source",

        "neighbor_fegenie_positive",
        "neighbor_fegenie_HMMs",
        "neighbor_fegenie_categories",

        "neighbor_findmehemes_positive",
        "neighbor_number_of_hemes",

        "neighbor_signalp_prediction",

        "neighbor_deeptmhmm_class",
        "neighbor_deeptmhmm_n_tm_helices",

        "neighbor_export_evidence",
        "neighbor_localization_class",

        "is_focal",
    ],
    OBSERVED_NEIGHBORHOODS.name,
)


evidence_source_columns = [
    "genome",
    "neighbor_protein_id",

    "neighbor_start",
    "neighbor_end",
    "neighbor_strand",
    "neighbor_gene_rank",

    "neighbor_annotation_match_type",
    "neighbor_annotation_accepted",

    "neighbor_globdb_cog",
    "neighbor_globdb_gene",
    "neighbor_globdb_product",

    "neighbor_cluster",
    "neighbor_module",

    "neighbor_candidate_source",

    "neighbor_fegenie_positive",
    "neighbor_fegenie_HMMs",
    "neighbor_fegenie_categories",

    "neighbor_findmehemes_positive",
    "neighbor_number_of_hemes",

    "neighbor_signalp_prediction",

    "neighbor_deeptmhmm_class",
    "neighbor_deeptmhmm_n_tm_helices",

    "neighbor_export_evidence",
    "neighbor_localization_class",
]


evidence = neigh[
    evidence_source_columns
].copy()


rename_evidence = {
    "neighbor_protein_id":
        "protein_id",

    "neighbor_start":
        "start",

    "neighbor_end":
        "end",

    "neighbor_strand":
        "strand",

    "neighbor_gene_rank":
        "gene_rank",

    "neighbor_annotation_match_type":
        "annotation_match_type",

    "neighbor_annotation_accepted":
        "annotation_accepted",

    "neighbor_globdb_cog":
        "globdb_cog",

    "neighbor_globdb_gene":
        "globdb_gene",

    "neighbor_globdb_product":
        "globdb_product",

    "neighbor_cluster":
        "cluster",

    "neighbor_module":
        "module",

    "neighbor_candidate_source":
        "candidate_source",

    "neighbor_fegenie_positive":
        "fegenie_positive",

    "neighbor_fegenie_HMMs":
        "fegenie_HMMs",

    "neighbor_fegenie_categories":
        "fegenie_categories",

    "neighbor_findmehemes_positive":
        "findmehemes_positive",

    "neighbor_number_of_hemes":
        "number_of_hemes",

    "neighbor_signalp_prediction":
        "signalp_prediction",

    "neighbor_deeptmhmm_class":
        "deeptmhmm_class",

    "neighbor_deeptmhmm_n_tm_helices":
        "deeptmhmm_n_tm_helices",

    "neighbor_export_evidence":
        "export_evidence",

    "neighbor_localization_class":
        "localization_class",
}


evidence = evidence.rename(
    columns=rename_evidence
)


## Remove exact duplicate evidence rows produced because the same
## protein is observed around multiple focal genes. ##

evidence_distinct = evidence.drop_duplicates()


## After removing identical rows, each genome + protein must have only
## one invariant evidence record. ##

if evidence_distinct[
    [
        "genome",
        "protein_id",
    ]
].duplicated().any():

    duplicated_keys = (
        evidence_distinct.loc[
            evidence_distinct[
                [
                    "genome",
                    "protein_id",
                ]
            ].duplicated(
                keep=False
            ),
            [
                "genome",
                "protein_id",
            ],
        ]
        .drop_duplicates()
    )


    fail(
        f"{len(duplicated_keys):,} proteins have "
        f"conflicting supposedly invariant evidence "
        f"across Stage-12 neighborhood rows."
    )


evidence = evidence_distinct.copy()


print(
    f"  Unique proteins represented in Stage 12: "
    f"{len(evidence):,}"
)


## ================================================================== ##
## 5. Integrated evidence for all 10,537 module proteins
## ================================================================== ##

print()
print("Building module-protein evidence table...")


## Prefix evidence fields to avoid ambiguity with authoritative
## module membership fields. ##

evidence_for_join = evidence.rename(
    columns={
        column:
            (
                column

                if column in [
                    "genome",
                    "protein_id",
                ]

                else
                "evidence_" + column
            )

        for column
        in evidence.columns
    }
)


module_proteins = protein_membership.merge(
    evidence_for_join,
    on=[
        "genome",
        "protein_id",
    ],
    how="left",
    validate="one_to_one",
    indicator=True,
)


if (
    module_proteins[
        "_merge"
    ]
    !=
    "both"
).any():

    fail(
        "At least one module protein is absent "
        "from the Stage-12 protein evidence catalogue."
    )


module_proteins = module_proteins.drop(
    columns="_merge"
)


## Validate evidence-level cluster/module where present. ##

if "evidence_cluster" in module_proteins.columns:

    bad = module_proteins[
        (
            module_proteins[
                "evidence_cluster"
            ]
            !=
            ""
        )
        &
        (
            module_proteins[
                "evidence_cluster"
            ]
            !=
            module_proteins[
                "cluster"
            ]
        )
    ]


    if len(
        bad
    ) > 0:

        fail(
            "Stage-12 evidence cluster disagrees "
            "with module protein membership."
        )


if "evidence_module" in module_proteins.columns:

    bad = module_proteins[
        (
            module_proteins[
                "evidence_module"
            ]
            !=
            ""
        )
        &
        (
            module_proteins[
                "evidence_module"
            ]
            !=
            module_proteins[
                "module"
            ]
        )
    ]


    if len(
        bad
    ) > 0:

        fail(
            "Stage-12 evidence module disagrees "
            "with module protein membership."
        )


## ================================================================== ##
## 6. Aggregate protein evidence per MMseqs2 family
##
## This complements cluster_summary.tsv with report-oriented fields.
## ================================================================== ##

print()
print("Summarizing protein evidence by MMseqs2 family...")


cluster_evidence_rows = []


for (
    module,
    cluster,
), group in module_proteins.groupby(
    [
        "module",
        "cluster",
    ],
    sort=True,
):

    row = {
        "module":
            module,

        "cluster":
            cluster,

        "report_n_proteins":
            len(
                group
            ),

        "report_n_genomes":
            group[
                "genome"
            ]
            .nunique(),
    }


    ## -------------------------------------------------------------- ##
    ## Boolean/count evidence
    ## -------------------------------------------------------------- ##

    boolean_fields = [
        (
            "evidence_fegenie_positive",
            "report_n_fegenie_positive",
            "report_pct_fegenie_positive",
        ),

        (
            "evidence_findmehemes_positive",
            "report_n_findmehemes_positive",
            "report_pct_findmehemes_positive",
        ),

        (
            "evidence_export_evidence",
            "report_n_export_evidence",
            "report_pct_export_evidence",
        ),

        (
            "evidence_annotation_accepted",
            "report_n_globdb_annotation_accepted",
            "report_pct_globdb_annotation_accepted",
        ),
    ]


    for (
        source,
        count_name,
        pct_name,
    ) in boolean_fields:

        if source not in group.columns:

            continue


        numeric = as_numeric(
            group[
                source
            ]
        )


        informative = int(
            numeric.notna().sum()
        )


        positive = int(
            (
                numeric
                ==
                1
            )
            .sum()
        )


        row[
            count_name
        ] = positive


        row[
            count_name
            +
            "_informative"
        ] = informative


        row[
            pct_name
        ] = pct(
            positive,
            informative,
        )


    ## -------------------------------------------------------------- ##
    ## Heme count
    ## -------------------------------------------------------------- ##

    if "evidence_number_of_hemes" in group.columns:

        hemes = as_numeric(
            group[
                "evidence_number_of_hemes"
            ]
        )


        if hemes.notna().any():

            row[
                "report_mean_heme_count"
            ] = hemes.mean()


            row[
                "report_min_heme_count"
            ] = hemes.min()


            row[
                "report_max_heme_count"
            ] = hemes.max()


    ## -------------------------------------------------------------- ##
    ## Detailed categorical compositions
    ##
    ## These retain the original evidence rather than selecting one
    ## annotation as "truth".
    ## -------------------------------------------------------------- ##

    distribution_fields = [
        (
            "evidence_fegenie_HMMs",
            "report_fegenie_HMM_distribution",
        ),

        (
            "evidence_fegenie_categories",
            "report_fegenie_category_distribution",
        ),

        (
            "evidence_signalp_prediction",
            "report_signalp_distribution",
        ),

        (
            "evidence_deeptmhmm_class",
            "report_deeptmhmm_class_distribution",
        ),

        (
            "evidence_localization_class",
            "report_localization_distribution",
        ),

        (
            "evidence_candidate_source",
            "report_candidate_source_distribution",
        ),

        (
            "evidence_globdb_cog",
            "report_globdb_cog_distribution",
        ),

        (
            "evidence_globdb_product",
            "report_globdb_product_distribution",
        ),
    ]


    for source, target in distribution_fields:

        if source in group.columns:

            row[
                target
            ] = distribution_string(
                group[
                    source
                ]
            )


    cluster_evidence_rows.append(
        row
    )


cluster_evidence = pd.DataFrame(
    cluster_evidence_rows
)


if len(
    cluster_evidence
) != EXPECTED_MODULE_CLUSTERS:

    fail(
        f"Expected "
        f"{EXPECTED_MODULE_CLUSTERS} "
        f"cluster evidence rows but produced "
        f"{len(cluster_evidence)}."
    )


member_clusters = member_clusters.merge(
    cluster_evidence,
    on=[
        "module",
        "cluster",
    ],
    how="left",
    validate="one_to_one",
)


## ================================================================== ##
## 7. Read occurrence / completeness evidence
## ================================================================== ##

print()
print("Reading module occurrence/completeness evidence...")


genome_module = read_tsv(
    GENOME_MODULE_SUMMARY_ALL
)


require_columns(
    genome_module,
    [
        "genome",
        "module",
    ],
    GENOME_MODULE_SUMMARY_ALL.name,
)


module_completeness = read_tsv(
    MODULE_COMPLETENESS_SUMMARY
)


require_columns(
    module_completeness,
    [
        "module",
    ],
    MODULE_COMPLETENESS_SUMMARY.name,
)


ensure_one_row_per(
    module_completeness,
    "module",
    MODULE_COMPLETENESS_SUMMARY.name,
)


## ================================================================== ##
## 8. Read Stage 13A physical association evidence
## ================================================================== ##

print()
print("Reading module colocalization evidence...")


pair_coloc = read_tsv(
    PAIR_COLOCALIZATION
)


require_columns(
    pair_coloc,
    [
        "module",
        "cluster_a",
        "cluster_b",
    ],
    PAIR_COLOCALIZATION.name,
)


module_coloc = read_tsv(
    MODULE_COLOCALIZATION
)


require_columns(
    module_coloc,
    [
        "module",
    ],
    MODULE_COLOCALIZATION.name,
)


ensure_one_row_per(
    module_coloc,
    "module",
    MODULE_COLOCALIZATION.name,
)


## ================================================================== ##
## 9. Read Stage 13D consensus evidence
## ================================================================== ##

print()
print("Reading focal consensus-neighborhood evidence...")


positions = read_tsv(
    FOCAL_POSITION_CONSENSUS
)


require_columns(
    positions,
    [
        "focal_cluster",
        "focal_module",
        "oriented_gene_offset",
    ],
    FOCAL_POSITION_CONSENSUS.name,
)


if len(
    positions
) != EXPECTED_POSITION_ROWS:

    fail(
        f"Expected "
        f"{EXPECTED_POSITION_ROWS:,} "
        f"consensus-position rows but found "
        f"{len(positions):,}."
    )


mmseq_consensus = read_tsv(
    MMSEQ_CONSENSUS
)


require_columns(
    mmseq_consensus,
    [
        "focal_cluster",
        "focal_module",
        "neighbor_cluster",
        "direction",
    ],
    MMSEQ_CONSENSUS.name,
)


if len(
    mmseq_consensus
) != EXPECTED_MMSEQ_CONSENSUS_ROWS:

    fail(
        f"Expected "
        f"{EXPECTED_MMSEQ_CONSENSUS_ROWS:,} "
        f"MMseq consensus rows but found "
        f"{len(mmseq_consensus):,}."
    )


cog_consensus = read_tsv(
    COG_CONSENSUS
)


require_columns(
    cog_consensus,
    [
        "focal_cluster",
        "focal_module",
        "neighbor_cog",
        "direction",
    ],
    COG_CONSENSUS.name,
)


if len(
    cog_consensus
) != EXPECTED_COG_CONSENSUS_ROWS:

    fail(
        f"Expected "
        f"{EXPECTED_COG_CONSENSUS_ROWS:,} "
        f"COG consensus rows but found "
        f"{len(cog_consensus):,}."
    )


## ================================================================== ##
## 10. Read authoritative local intervening-gene context
## ================================================================== ##

print()
print("Reading local intervening-gene evidence...")


local_genes = read_tsv(
    LOCAL_INTERVENING_GENES
)


require_columns(
    local_genes,
    [
        "genome",
        "module",
        "cluster_a",
        "cluster_b",
        "intervening_protein_id",
    ],
    LOCAL_INTERVENING_GENES.name,
)


if len(
    local_genes
) != EXPECTED_LOCAL_INTERVENING_GENES:

    fail(
        f"Expected "
        f"{EXPECTED_LOCAL_INTERVENING_GENES:,} "
        f"local intervening genes but found "
        f"{len(local_genes):,}."
    )


local_patterns = read_tsv(
    LOCAL_INTERVENING_PATTERNS
)


require_columns(
    local_patterns,
    [
        "module",
        "cluster_a",
        "cluster_b",
    ],
    LOCAL_INTERVENING_PATTERNS.name,
)


local_pair_summary = read_tsv(
    LOCAL_PAIR_SUMMARY
)


require_columns(
    local_pair_summary,
    [
        "module",
        "cluster_a",
        "cluster_b",
    ],
    LOCAL_PAIR_SUMMARY.name,
)


## ================================================================== ##
## 11. Integrate protein-level evidence into intervening genes
##
## THIS is where FeGenie etc. can finally contribute to X/Y genes
## in module reports without changing Stage-13 neighborhood calling.
## ================================================================== ##

print()
print(
    "Integrating FeGenie / heme / localization evidence "
    "for local intervening genes..."
)


intervening_lookup = evidence_for_join.rename(
    columns={
        "protein_id":
            "intervening_protein_id",
    }
)


intervening_evidence = local_genes.merge(
    intervening_lookup,
    on=[
        "genome",
        "intervening_protein_id",
    ],
    how="left",
    validate="many_to_one",
    indicator=True,
)


n_unmatched_intervening = int(
    (
        intervening_evidence[
            "_merge"
        ]
        !=
        "both"
    )
    .sum()
)


if n_unmatched_intervening > 0:

    fail(
        f"{n_unmatched_intervening:,} "
        f"local intervening genes could not be "
        f"recovered from the Stage-12 protein evidence catalogue."
    )


intervening_evidence = intervening_evidence.drop(
    columns="_merge"
)


## ================================================================== ##
## 12. Cross-check focal cluster -> module identities
## ================================================================== ##

expected_module_lookup = dict(
    zip(
        membership[
            "cluster"
        ],
        membership[
            "module"
        ],
    )
)


for (
    table,
    cluster_column,
    module_column,
    label,
) in [
    (
        positions,
        "focal_cluster",
        "focal_module",
        "position consensus",
    ),

    (
        mmseq_consensus,
        "focal_cluster",
        "focal_module",
        "MMseq consensus",
    ),

    (
        cog_consensus,
        "focal_cluster",
        "focal_module",
        "COG consensus",
    ),
]:

    expected = table[
        cluster_column
    ].map(
        expected_module_lookup
    )


    bad = (
        expected
        !=
        table[
            module_column
        ]
    )


    if bad.any():

        fail(
            f"{int(bad.sum()):,} "
            f"{label} rows disagree with "
            f"authoritative cluster -> module membership."
        )


## ================================================================== ##
## 13. Build one-row module overview
## ================================================================== ##

print()
print("Building integrated module overview...")


overview_rows = []


for module in modules:

    module_clusters = membership[
        membership[
            "module"
        ]
        ==
        module
    ][
        "cluster"
    ]


    module_protein_rows = module_proteins[
        module_proteins[
            "module"
        ]
        ==
        module
    ]


    module_pair_rows = pair_coloc[
        pair_coloc[
            "module"
        ]
        ==
        module
    ]


    module_position_rows = positions[
        positions[
            "focal_module"
        ]
        ==
        module
    ]


    module_mmseq_rows = mmseq_consensus[
        mmseq_consensus[
            "focal_module"
        ]
        ==
        module
    ]


    module_cog_rows = cog_consensus[
        cog_consensus[
            "focal_module"
        ]
        ==
        module
    ]


    module_local_genes = intervening_evidence[
        intervening_evidence[
            "module"
        ]
        ==
        module
    ]


    cluster_ev = cluster_evidence[
        cluster_evidence[
            "module"
        ]
        ==
        module
    ]


    if "report_n_fegenie_positive" in cluster_ev.columns:

        n_fegenie_clusters = int(
            (
                cluster_ev[
                    "report_n_fegenie_positive"
                ]
                >
                0
            )
            .sum()
        )

    else:

        n_fegenie_clusters = 0


    overview_rows.append(
        {
            "module":
                module,

            "n_member_clusters":
                len(
                    module_clusters
                ),

            "n_module_proteins":
                len(
                    module_protein_rows
                ),

            "n_genomes_with_any_module_protein":
                module_protein_rows[
                    "genome"
                ]
                .nunique(),

            "n_member_clusters_with_fegenie_positive_protein":
                n_fegenie_clusters,

            "n_pair_colocalization_rows":
                len(
                    module_pair_rows
                ),

            "n_focal_consensus_position_rows":
                len(
                    module_position_rows
                ),

            "n_mmseq_directional_consensus_rows":
                len(
                    module_mmseq_rows
                ),

            "n_cog_directional_consensus_rows":
                len(
                    module_cog_rows
                ),

            "n_local_intervening_gene_rows":
                len(
                    module_local_genes
                ),
        }
    )


overview = pd.DataFrame(
    overview_rows
)


## Merge authoritative Stage-09 completeness summary. ##

overview = overview.merge(
    prefix_nonkey_columns(
        module_completeness,
        "module",
        "occurrence_",
    ),
    on="module",
    how="left",
    validate="one_to_one",
)


## Merge Stage-13A module colocalization summary. ##

overview = overview.merge(
    prefix_nonkey_columns(
        module_coloc,
        "module",
        "colocalization_",
    ),
    on="module",
    how="left",
    validate="one_to_one",
)


## ================================================================== ##
## 14. Write global integrated tables
## ================================================================== ##

print()
print("Writing global report-data tables...")


write_table(
    overview,
    OUT_OVERVIEW,
)


write_table(
    member_clusters,
    OUT_CLUSTER_EVIDENCE,
)


write_table(
    module_proteins,
    OUT_PROTEIN_EVIDENCE,
)


write_table(
    intervening_evidence,
    OUT_INTERVENING_EVIDENCE,
)


## ================================================================== ##
## 15. Build per-module report packages
## ================================================================== ##

print()
print("Writing per-module report packages...")


manifest_rows = []
index_rows = []


for module in modules:

    outdir = (
        MODULE_DIR
        / module
    )


    outdir.mkdir(
        parents=True,
        exist_ok=True,
    )


    tables = {
        "00_module_overview.tsv":
            overview[
                overview[
                    "module"
                ]
                ==
                module
            ].copy(),

        "01_member_clusters.tsv":
            member_clusters[
                member_clusters[
                    "module"
                ]
                ==
                module
            ].copy(),

        "02_genome_occurrence_completeness.tsv":
            genome_module[
                genome_module[
                    "module"
                ]
                ==
                module
            ].copy(),

        "03_module_proteins.tsv":
            module_proteins[
                module_proteins[
                    "module"
                ]
                ==
                module
            ].copy(),

        "04_pair_colocalization.tsv":
            pair_coloc[
                pair_coloc[
                    "module"
                ]
                ==
                module
            ].copy(),

        "05_focal_consensus_positions.tsv":
            positions[
                positions[
                    "focal_module"
                ]
                ==
                module
            ].copy(),

        "06_mmseq_consensus_associations.tsv":
            mmseq_consensus[
                mmseq_consensus[
                    "focal_module"
                ]
                ==
                module
            ].copy(),

        "07_cog_consensus_associations.tsv":
            cog_consensus[
                cog_consensus[
                    "focal_module"
                ]
                ==
                module
            ].copy(),

        "08_local_pair_summary.tsv":
            local_pair_summary[
                local_pair_summary[
                    "module"
                ]
                ==
                module
            ].copy(),

        "09_local_intervening_patterns.tsv":
            local_patterns[
                local_patterns[
                    "module"
                ]
                ==
                module
            ].copy(),

        "10_local_intervening_genes.tsv":
            local_genes[
                local_genes[
                    "module"
                ]
                ==
                module
            ].copy(),

        "11_local_intervening_gene_integrated_evidence.tsv":
            intervening_evidence[
                intervening_evidence[
                    "module"
                ]
                ==
                module
            ].copy(),
    }


    for filename, table in tables.items():

        path = (
            outdir
            / filename
        )


        write_table(
            table,
            path,
        )


        manifest_rows.append(
            {
                "module":
                    module,

                "table":
                    filename,

                "relative_path":
                    str(
                        path.relative_to(
                            HERE
                        )
                    ),

                "n_rows":
                    len(
                        table
                    ),

                "n_columns":
                    len(
                        table.columns
                    ),
            }
        )


    ## -------------------------------------------------------------- ##
    ## Small human-readable README for every module
    ## -------------------------------------------------------------- ##

    ov = tables[
        "00_module_overview.tsv"
    ].iloc[0]


    readme = f"""\
{module} - Stage 14A report-data package
{'=' * 72}

Primary identity
----------------
MMseqs2 sequence-family clusters remain the primary protein-family unit.

Module interpretation
---------------------
The MCL module represents genome-level co-occurrence structure. It is NOT
assumed to be an operon, biochemical complex, or single linear gene cassette.

Evidence included
-----------------
01_member_clusters.tsv
    MMseqs2 family statistics plus integrated FeGenie/heme/localization
    summaries for the module member families.

02_genome_occurrence_completeness.tsv
    Genome-level module occurrence and completeness.

03_module_proteins.tsv
    Protein-level evidence for all proteins belonging to module families.

04_pair_colocalization.tsv
    Stage-13A pairwise physical association evidence.

05_focal_consensus_positions.tsv
    Evidence-preserving exact oriented neighborhoods around every focal
    module family.

06_mmseq_consensus_associations.tsv
    Directional MMseqs2-family neighborhood relationships, including exact
    positions, cumulative saturation radius, physical windows, and support.

07_cog_consensus_associations.tsv
    COG-based functional context around focal sequence families.

08_local_pair_summary.tsv
    Local <=20-kb pair-level architecture summary.

09_local_intervening_patterns.tsv
    Repeated patterns between locally associated module families.

10_local_intervening_genes.tsv
    Actual Prodigal genes occupying local intervening positions.

11_local_intervening_gene_integrated_evidence.tsv
    Same intervening genes enriched with FeGenie, FindMeHemes, SignalP,
    DeepTMHMM, heme, GlobDB, and localization evidence where available.

Evidence-strength principle
---------------------------
Percentages must be interpreted with their numerator and denominator.
For example, 100% (1/1) is not equivalent in evidential strength to
100% (44/44).

Module size
-----------
Member MMseqs2 clusters: {ov['n_member_clusters']}
Module proteins:         {ov['n_module_proteins']}
Genomes represented:     {ov['n_genomes_with_any_module_protein']}

This package is an evidence source for the later integrated module report.
"""


    readme_path = (
        outdir
        / "README.txt"
    )


    readme_path.write_text(
        readme
    )


    manifest_rows.append(
        {
            "module":
                module,

            "table":
                "README.txt",

            "relative_path":
                str(
                    readme_path.relative_to(
                        HERE
                    )
                ),

            "n_rows":
                "",

            "n_columns":
                "",
        }
    )


    index_rows.append(
        {
            "module":
                module,

            "report_directory":
                str(
                    outdir.relative_to(
                        HERE
                    )
                ),

            "n_member_clusters":
                int(
                    ov[
                        "n_member_clusters"
                    ]
                ),

            "n_module_proteins":
                int(
                    ov[
                        "n_module_proteins"
                    ]
                ),

            "n_genomes_with_any_module_protein":
                int(
                    ov[
                        "n_genomes_with_any_module_protein"
                    ]
                ),

            "n_pair_colocalization_rows":
                int(
                    ov[
                        "n_pair_colocalization_rows"
                    ]
                ),

            "n_focal_consensus_position_rows":
                int(
                    ov[
                        "n_focal_consensus_position_rows"
                    ]
                ),

            "n_mmseq_directional_consensus_rows":
                int(
                    ov[
                        "n_mmseq_directional_consensus_rows"
                    ]
                ),

            "n_cog_directional_consensus_rows":
                int(
                    ov[
                        "n_cog_directional_consensus_rows"
                    ]
                ),

            "n_local_intervening_gene_rows":
                int(
                    ov[
                        "n_local_intervening_gene_rows"
                    ]
                ),
        }
    )


manifest = pd.DataFrame(
    manifest_rows
)


index = pd.DataFrame(
    index_rows
)


write_table(
    manifest,
    OUT_MANIFEST,
)


write_table(
    index,
    OUT_INDEX,
)


## ================================================================== ##
## 16. Cross-module split QC
## ================================================================== ##

print()
print("Running report-package QC...")


if index[
    "n_member_clusters"
].sum() != EXPECTED_MODULE_CLUSTERS:

    fail(
        "Per-module cluster counts do not sum "
        "to 156."
    )


if index[
    "n_module_proteins"
].sum() != EXPECTED_MODULE_PROTEINS:

    fail(
        "Per-module protein counts do not sum "
        "to 10,537."
    )


if index[
    "n_focal_consensus_position_rows"
].sum() != EXPECTED_POSITION_ROWS:

    fail(
        "Per-module position rows do not sum "
        "to 6,396."
    )


if index[
    "n_mmseq_directional_consensus_rows"
].sum() != EXPECTED_MMSEQ_CONSENSUS_ROWS:

    fail(
        "Per-module MMseq consensus rows "
        "do not sum to 3,050."
    )


if index[
    "n_cog_directional_consensus_rows"
].sum() != EXPECTED_COG_CONSENSUS_ROWS:

    fail(
        "Per-module COG consensus rows "
        "do not sum to 107,728."
    )


if index[
    "n_local_intervening_gene_rows"
].sum() != EXPECTED_LOCAL_INTERVENING_GENES:

    fail(
        "Per-module local intervening genes "
        "do not sum to 6,950."
    )


## Every focal cluster must have exactly 41 position rows. ##

position_cluster_counts = (
    positions
    .groupby(
        "focal_cluster"
    )
    .size()
)


if (
    position_cluster_counts
    !=
    41
).any():

    fail(
        "At least one focal MMseqs2 cluster "
        "does not have 41 position-summary rows."
    )


## ================================================================== ##
## 17. Extra integrated-evidence QC
## ================================================================== ##

## Module proteins with FeGenie information. ##

if "evidence_fegenie_positive" in module_proteins.columns:

    fegenie_numeric = as_numeric(
        module_proteins[
            "evidence_fegenie_positive"
        ]
    )


    n_fegenie_positive = int(
        (
            fegenie_numeric
            ==
            1
        )
        .sum()
    )

else:

    n_fegenie_positive = 0


## Intervening genes with any FeGenie hit. ##

if "evidence_fegenie_positive" in intervening_evidence.columns:

    intervening_fegenie = as_numeric(
        intervening_evidence[
            "evidence_fegenie_positive"
        ]
    )


    n_intervening_fegenie_positive = int(
        (
            intervening_fegenie
            ==
            1
        )
        .sum()
    )

else:

    n_intervening_fegenie_positive = 0


## Intervening genes with FindMeHemes evidence. ##

if "evidence_findmehemes_positive" in intervening_evidence.columns:

    intervening_fmh = as_numeric(
        intervening_evidence[
            "evidence_findmehemes_positive"
        ]
    )


    n_intervening_findmehemes_positive = int(
        (
            intervening_fmh
            ==
            1
        )
        .sum()
    )

else:

    n_intervening_findmehemes_positive = 0


## ================================================================== ##
## 18. QC table
## ================================================================== ##

qc = pd.DataFrame(
    [
        [
            "mcl_modules",
            len(
                modules
            ),
        ],

        [
            "module_clusters",
            len(
                membership
            ),
        ],

        [
            "module_proteins",
            len(
                module_proteins
            ),
        ],

        [
            "module_proteins_with_fegenie_positive",
            n_fegenie_positive,
        ],

        [
            "consensus_position_rows",
            len(
                positions
            ),
        ],

        [
            "mmseq_consensus_rows",
            len(
                mmseq_consensus
            ),
        ],

        [
            "cog_consensus_rows",
            len(
                cog_consensus
            ),
        ],

        [
            "local_intervening_gene_rows",
            len(
                intervening_evidence
            ),
        ],

        [
            "local_intervening_genes_with_fegenie_positive",
            n_intervening_fegenie_positive,
        ],

        [
            "local_intervening_genes_with_findmehemes_positive",
            n_intervening_findmehemes_positive,
        ],

        [
            "module_report_directories",
            len(
                index
            ),
        ],

        [
            "protein_evidence_mapping_qc_pass",
            1,
        ],

        [
            "intervening_evidence_mapping_qc_pass",
            1,
        ],

        [
            "cluster_module_mapping_qc_pass",
            1,
        ],

        [
            "consensus_module_mapping_qc_pass",
            1,
        ],

        [
            "per_module_split_totals_qc_pass",
            1,
        ],
    ],
    columns=[
        "metric",
        "value",
    ],
)


write_table(
    qc,
    OUT_QC,
)


## ================================================================== ##
## 19. Terminal summary
## ================================================================== ##

print()
print("Integrated module-report data summary")


print(
    f"  Modules:                              "
    f"{len(modules):,}"
)

print(
    f"  Member MMseqs2 clusters:              "
    f"{len(membership):,}"
)

print(
    f"  Module proteins:                      "
    f"{len(module_proteins):,}"
)

print(
    f"  Module proteins FeGenie-positive:     "
    f"{n_fegenie_positive:,}"
)

print(
    f"  Focal consensus-position rows:        "
    f"{len(positions):,}"
)

print(
    f"  MMseq directional consensus rows:     "
    f"{len(mmseq_consensus):,}"
)

print(
    f"  COG directional consensus rows:       "
    f"{len(cog_consensus):,}"
)

print(
    f"  Local intervening genes:              "
    f"{len(intervening_evidence):,}"
)

print(
    f"  Intervening genes FeGenie-positive:   "
    f"{n_intervening_fegenie_positive:,}"
)

print(
    f"  Intervening genes FindMeHemes-positive:"
    f" {n_intervening_findmehemes_positive:,}"
)

print(
    f"  Per-module report packages:           "
    f"{len(index):,}"
)


## ================================================================== ##
## 20. Spotlight Module 20 / Module 35 integrated evidence
## ================================================================== ##

for spotlight_module in [
    "Module_20",
    "Module_35",
]:

    print()
    print(
        f"{spotlight_module} report-data spotlight"
    )


    ov = overview[
        overview[
            "module"
        ]
        ==
        spotlight_module
    ].iloc[0]


    print(
        f"  Member clusters:              "
        f"{ov['n_member_clusters']}"
    )

    print(
        f"  Module proteins:              "
        f"{ov['n_module_proteins']}"
    )

    print(
        f"  Genomes represented:          "
        f"{ov['n_genomes_with_any_module_protein']}"
    )

    print(
        f"  Pair colocalization rows:     "
        f"{ov['n_pair_colocalization_rows']}"
    )

    print(
        f"  Local intervening genes:      "
        f"{ov['n_local_intervening_gene_rows']}"
    )


    spot_clusters = member_clusters[
        member_clusters[
            "module"
        ]
        ==
        spotlight_module
    ]


    display_columns = [
        column

        for column in [
            "cluster",
            "representative_protein",
            "n_proteins",
            "n_genomes",
            "mean_heme_count",
            "dominant_fegenie_HMM",
            "report_fegenie_HMM_distribution",
            "report_signalp_distribution",
            "report_deeptmhmm_class_distribution",
            "report_localization_distribution",
        ]

        if column
        in spot_clusters.columns
    ]


    print()
    print(
        spot_clusters[
            display_columns
        ]
        .to_string(
            index=False
        )
    )


print()
print("=" * 80)
print("SUCCESS")
print("=" * 80)


print(
    f"Global overview:           "
    f"{OUT_OVERVIEW}"
)

print(
    f"Member-cluster evidence:   "
    f"{OUT_CLUSTER_EVIDENCE}"
)

print(
    f"Module-protein evidence:   "
    f"{OUT_PROTEIN_EVIDENCE}"
)

print(
    f"Intervening-gene evidence: "
    f"{OUT_INTERVENING_EVIDENCE}"
)

print(
    f"Report index:              "
    f"{OUT_INDEX}"
)

print(
    f"Manifest:                  "
    f"{OUT_MANIFEST}"
)

print(
    f"QC:                        "
    f"{OUT_QC}"
)

print(
    f"Per-module packages:       "
    f"{MODULE_DIR}"
)
