#!/usr/bin/env python3

from __future__ import annotations

import hashlib
import math
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd


## ================================================================== ##
## STAGE 58
##
## BUILD STABLE LOCAL-FAMILY + ARCHITECTURE TABLES
##
## Purpose
## -------
##
## Convert the raw Stage-57 MMseqs2 clustering into a stable
## family-level representation of the 21 priority neighborhoods.
##
## Core principles:
##
##   1. Actual Prodigal genes remain the biological units within
##      individual neighborhoods.
##
##   2. MMseqs2 families describe sequence homology across
##      neighborhoods.
##
##   3. Family identity is defined from the exact SET OF MEMBERS,
##      not from whichever protein MMseqs happened to choose as
##      representative.
##
##   4. No recurrent-family threshold is imposed here.
##
##   5. No genes are removed because they lack annotation.
##
## Outputs include:
##
##   - stable MG_local_### membership
##   - gene-level neighborhood table with local-family membership
##   - overall local-family summary
##   - focal-gene -> local-family crosswalk
##   - region x family presence table
##   - family recurrence within each focal group
##   - ordered family architecture for each focal region
##
## The resulting gene-level table is intended to become the input
## for later gene-map visualization.
## ================================================================== ##


## ================================================================== ##
## Paths
## ================================================================== ##

PROJECT = (
    Path.home()
    / "methanotrophs"
    / "methanotroph_project"
    / "jeppe"
    / "metagenome"
)


WORK_DIR = (
    PROJECT
    / "functional_analysis"
    / "methylococcales_neighborhoods"
)


LOCAL_DIR = (
    WORK_DIR
    / "local_clustering"
)


MMSEQS_CLUSTERS = (
    LOCAL_DIR
    / "MG_local40_cov80_cluster.tsv"
)


PROTEIN_METADATA = (
    WORK_DIR
    / "priority_neighborhood_protein_metadata_annotated.tsv"
)


NEIGHBORHOODS = (
    WORK_DIR
    / "priority_observed_neighborhood_genes_annotated.tsv"
)


FOCAL_OCCURRENCES = (
    WORK_DIR
    / "priority_focal_occurrences.tsv"
)


OUT_MEMBERSHIP = (
    LOCAL_DIR
    / "MG_local_family_membership.tsv"
)


OUT_FAMILY_SUMMARY = (
    LOCAL_DIR
    / "MG_local_family_summary.tsv"
)


OUT_FOCAL_CROSSWALK = (
    LOCAL_DIR
    / "priority_focal_local_family_crosswalk.tsv"
)


OUT_REGION_PRESENCE = (
    LOCAL_DIR
    / "MG_local_family_region_presence.tsv"
)


OUT_GROUP_SUMMARY = (
    LOCAL_DIR
    / "MG_local_family_by_focal_group.tsv"
)


OUT_ARCHITECTURES = (
    LOCAL_DIR
    / "priority_region_family_architectures.tsv"
)


OUT_GENE_TABLE = (
    WORK_DIR
    / "priority_observed_neighborhood_genes_local_families.tsv"
)


OUT_QC = (
    LOCAL_DIR
    / "local_family_architecture_qc.tsv"
)


## ================================================================== ##
## Locked upstream expectations
## ================================================================== ##

EXPECTED_PROTEINS = 845

EXPECTED_NEIGHBORHOOD_ROWS = 845

EXPECTED_FOCALS = 21


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

    if not path.is_file():

        fail(
            f"Required file does not exist:\n{path}"
        )


def clean(value):

    if value is None:

        return ""

    if pd.isna(value):

        return ""

    return str(value).strip()


def as_int(value):

    value = clean(value)

    if value == "":

        return 0

    try:

        return int(
            float(
                value
            )
        )

    except Exception:

        return 0


def as_float(value):

    value = clean(value)

    if value == "":

        return np.nan

    try:

        return float(
            value
        )

    except Exception:

        return np.nan


def dominant_nonempty(values):

    values = [
        clean(value)

        for value
        in values

        if clean(value)
    ]


    if not values:

        return {
            "value":
                "",

            "support":
                0,

            "total":
                0,

            "fraction":
                np.nan,

            "distribution":
                "",
        }


    counts = Counter(
        values
    )


    ordered = sorted(
        counts.items(),
        key=lambda item: (
            -item[1],
            item[0],
        ),
    )


    value, support = (
        ordered[0]
    )


    total = len(
        values
    )


    return {
        "value":
            value,

        "support":
            support,

        "total":
            total,

        "fraction":
            support / total,

        "distribution":
            ";".join(
                f"{key}:{count}"

                for key, count
                in ordered
            ),
    }


def semicolon_tokens(values):

    tokens = []


    for value in values:

        value = clean(value)

        if not value:

            continue


        for token in value.split(";"):

            token = token.strip()

            if token:

                tokens.append(
                    token
                )


    return tokens


def safe_sd(values):

    values = [
        float(value)

        for value
        in values

        if not pd.isna(
            value
        )
    ]


    if len(values) <= 1:

        return 0.0


    return float(
        np.std(
            values,
            ddof=1,
        )
    )


def stable_member_signature(
    group,
):

    members = sorted(
        (
            clean(row.genome)
            +
            "|"
            +
            clean(row.protein_id)
        )

        for row
        in group.itertuples(
            index=False
        )
    )


    return "||".join(
        members
    )


def stable_member_hash(
    signature,
):

    return hashlib.sha256(
        signature.encode(
            "utf-8"
        )
    ).hexdigest()[
        :16
    ]


def make_focal_group(row):

    cyc2 = as_int(
        row.get(
            "priority_cyc2",
            0,
        )
    )


    multiheme = as_int(
        row.get(
            "priority_ge5_hemes",
            0,
        )
    )


    top_cluster = clean(
        row.get(
            "globdb_top_scoring_cluster",
            ""
        )
    )


    focal_id = clean(
        row.get(
            "focal_id",
            ""
        )
    )


    if cyc2 == 1:

        if top_cluster:

            return (
                "Cyc2__"
                +
                top_cluster
            )

        return (
            "Cyc2__NoMatch__"
            +
            focal_id
        )


    if multiheme == 1:

        if top_cluster:

            return (
                "Multiheme__"
                +
                top_cluster
            )

        ## -------------------------------------------- ##
        ## Do not combine unrelated no-match proteins
        ## into a pseudo-family.
        ## -------------------------------------------- ##

        return (
            "Multiheme__NoMatch__"
            +
            focal_id
        )


    return (
        "Other__"
        +
        focal_id
    )


## ================================================================== ##
## Start
## ================================================================== ##

print("=" * 80)

print(
    "STAGE 58 - BUILD LOCAL FAMILY ARCHITECTURE TABLES"
)

print("=" * 80)


for path in (
    MMSEQS_CLUSTERS,
    PROTEIN_METADATA,
    NEIGHBORHOODS,
    FOCAL_OCCURRENCES,
):

    require_file(
        path
    )


LOCAL_DIR.mkdir(
    parents=True,
    exist_ok=True,
)


## ================================================================== ##
## 1. Read raw MMseqs cluster membership
##
## easy-cluster cluster TSV:
##
##     representative    member
##
## ================================================================== ##

print()
print("Reading Stage-57 MMseqs membership...")


raw_membership = pd.read_csv(
    MMSEQS_CLUSTERS,
    sep="\t",
    names=[
        "mmseqs_representative",
        "neighborhood_protein_id",
    ],
    header=None,
    dtype=str,
    keep_default_na=False,
)


if len(
    raw_membership
) != EXPECTED_PROTEINS:

    fail(
        f"Expected {EXPECTED_PROTEINS} MMseqs assignments "
        f"but found {len(raw_membership)}."
    )


if raw_membership[
    "neighborhood_protein_id"
].duplicated().any():

    duplicates = (
        raw_membership.loc[
            raw_membership[
                "neighborhood_protein_id"
            ].duplicated(
                keep=False
            ),
            "neighborhood_protein_id",
        ]
        .tolist()
    )


    fail(
        "One or more proteins have multiple MMseqs "
        "cluster assignments.\n"
        +
        "\n".join(
            duplicates[
                :20
            ]
        )
    )


n_raw_families = raw_membership[
    "mmseqs_representative"
].nunique()


print(
    f"  Protein assignments: "
    f"{len(raw_membership):,}"
)

print(
    f"  Raw MMseqs families: "
    f"{n_raw_families:,}"
)


## ================================================================== ##
## 2. Read neighborhood protein metadata
## ================================================================== ##

print()
print("Reading neighborhood protein metadata...")


metadata = pd.read_csv(
    PROTEIN_METADATA,
    sep="\t",
    dtype=str,
    keep_default_na=False,
)


required_metadata = {
    "neighborhood_protein_id",
    "genome",
    "protein_id",
}


missing = (
    required_metadata
    -
    set(
        metadata.columns
    )
)


if missing:

    fail(
        "Protein metadata is missing:\n"
        +
        "\n".join(
            sorted(
                missing
            )
        )
    )


if len(
    metadata
) != EXPECTED_PROTEINS:

    fail(
        f"Expected {EXPECTED_PROTEINS} metadata rows "
        f"but found {len(metadata)}."
    )


if metadata[
    "neighborhood_protein_id"
].duplicated().any():

    fail(
        "Duplicate neighborhood_protein_id values "
        "in protein metadata."
    )


if (
    set(
        raw_membership[
            "neighborhood_protein_id"
        ]
    )
    !=
    set(
        metadata[
            "neighborhood_protein_id"
        ]
    )
):

    fail(
        "MMseqs membership and protein metadata "
        "contain different protein ID sets."
    )


membership = raw_membership.merge(
    metadata,
    on="neighborhood_protein_id",
    how="left",
    validate="one_to_one",
)


## ================================================================== ##
## 3. Read actual-gene neighborhood observations
## ================================================================== ##

print()
print("Reading annotated actual-gene neighborhoods...")


genes = pd.read_csv(
    NEIGHBORHOODS,
    sep="\t",
    dtype=str,
    keep_default_na=False,
)


required_gene_columns = {
    "focal_id",
    "genome",

    "focal_protein_id",

    "neighbor_protein_id",

    "oriented_gene_offset",
    "oriented_midpoint_offset_bp",

    "neighbor_strand",
    "oriented_neighbor_strand",

    "is_focal",

    "neighbor_annotation_transfer_status",
    "neighbor_transferred_cog",
    "neighbor_transferred_product",
}


missing = (
    required_gene_columns
    -
    set(
        genes.columns
    )
)


if missing:

    fail(
        "Annotated neighborhood table is missing:\n"
        +
        "\n".join(
            sorted(
                missing
            )
        )
    )


if len(
    genes
) != EXPECTED_NEIGHBORHOOD_ROWS:

    fail(
        f"Expected {EXPECTED_NEIGHBORHOOD_ROWS} "
        f"neighborhood rows but found "
        f"{len(genes)}."
    )


## Stage 52 established that all 845 observations represent
## different actual genes.

if genes[
    [
        "genome",
        "neighbor_protein_id",
    ]
].duplicated().any():

    fail(
        "Expected each actual neighborhood gene to "
        "occur once, but duplicate genome + "
        "neighbor_protein_id keys were found."
    )


genes[
    "oriented_gene_offset"
] = pd.to_numeric(
    genes[
        "oriented_gene_offset"
    ],
    errors="raise",
).astype(int)


genes[
    "oriented_midpoint_offset_bp"
] = pd.to_numeric(
    genes[
        "oriented_midpoint_offset_bp"
    ],
    errors="raise",
)


genes[
    "is_focal"
] = pd.to_numeric(
    genes[
        "is_focal"
    ],
    errors="raise",
).astype(int)


## ================================================================== ##
## 4. Map raw MMseq families onto the actual genes
## ================================================================== ##

membership_key = membership[
    [
        "genome",
        "protein_id",
        "neighborhood_protein_id",
        "mmseqs_representative",
    ]
].rename(
    columns={
        "protein_id":
            "neighbor_protein_id",
    }
)


genes = genes.merge(
    membership_key,
    on=[
        "genome",
        "neighbor_protein_id",
    ],
    how="left",
    validate="one_to_one",
)


if genes[
    "mmseqs_representative"
].eq("").any():

    fail(
        "At least one actual neighborhood gene lacks "
        "an MMseqs family assignment."
    )


if genes[
    "mmseqs_representative"
].isna().any():

    fail(
        "At least one actual neighborhood gene lacks "
        "an MMseqs family assignment."
    )


## ================================================================== ##
## 5. Read and classify the 21 focal loci
## ================================================================== ##

print()
print("Building focal-group definitions...")


focals = pd.read_csv(
    FOCAL_OCCURRENCES,
    sep="\t",
    dtype=str,
    keep_default_na=False,
)


if len(
    focals
) != EXPECTED_FOCALS:

    fail(
        f"Expected {EXPECTED_FOCALS} focal rows "
        f"but found {len(focals)}."
    )


required_focal_columns = {
    "focal_id",
    "genome",
    "protein_id",

    "priority_cyc2",
    "priority_ge5_hemes",

    "number_of_hemes",

    "globdb_top_scoring_cluster",
}


missing = (
    required_focal_columns
    -
    set(
        focals.columns
    )
)


if missing:

    fail(
        "Focal occurrence table is missing:\n"
        +
        "\n".join(
            sorted(
                missing
            )
        )
    )


if focals[
    "focal_id"
].duplicated().any():

    fail(
        "Duplicate focal_id values."
    )


focals[
    "focal_group"
] = focals.apply(
    make_focal_group,
    axis=1,
)


group_sizes = (
    focals
    .groupby(
        "focal_group",
        as_index=False,
    )
    .agg(
        n_focal_regions_in_group=(
            "focal_id",
            "nunique",
        ),

        n_genomes_in_group=(
            "genome",
            "nunique",
        ),
    )
)


print()
print("Focal groups:")


for row in group_sizes.itertuples(
    index=False
):

    print(
        f"  {row.focal_group:<40} "
        f"{row.n_focal_regions_in_group:>2} regions"
    )


focal_group_lookup = focals[
    [
        "focal_id",
        "focal_group",
    ]
]


genes = genes.merge(
    focal_group_lookup,
    on="focal_id",
    how="left",
    validate="many_to_one",
)


if genes[
    "focal_group"
].eq("").any():

    fail(
        "One or more neighborhood rows lack "
        "a focal-group assignment."
    )


## ================================================================== ##
## 6. Compute biological metrics for RAW MMseq families
##
## These metrics are used only to give deterministic, readable
## MG_local_### ordering.
##
## Family identity itself remains the exact member set.
## ================================================================== ##

raw_family_rows = []


for raw_family, group in genes.groupby(
    "mmseqs_representative",
    sort=False,
):

    member_group = membership[
        membership[
            "mmseqs_representative"
        ]
        ==
        raw_family
    ]


    signature = stable_member_signature(
        member_group
    )


    raw_family_rows.append(
        {
            "mmseqs_representative":
                raw_family,

            "n_proteins":
                len(
                    member_group
                ),

            "n_focal_regions":
                group[
                    "focal_id"
                ].nunique(),

            "n_genomes":
                group[
                    "genome"
                ].nunique(),

            "n_focal_proteins":
                int(
                    group[
                        "is_focal"
                    ].sum()
                ),

            "member_signature":
                signature,

            "member_set_hash":
                stable_member_hash(
                    signature
                ),
        }
    )


raw_family_summary = pd.DataFrame(
    raw_family_rows
)


if len(
    raw_family_summary
) != n_raw_families:

    fail(
        "Raw family-summary row count does not "
        "equal number of MMseqs families."
    )


## ================================================================== ##
## 7. Assign stable MG_local_### IDs
##
## Ordering:
##
##   1. most focal regions
##   2. most genomes
##   3. most proteins
##   4. exact member signature
##
## Thus a different MMseqs representative will not change family
## identity or ordering if the actual cluster membership is unchanged.
## ================================================================== ##

raw_family_summary = (
    raw_family_summary
    .sort_values(
        [
            "n_focal_regions",
            "n_genomes",
            "n_proteins",
            "member_signature",
        ],
        ascending=[
            False,
            False,
            False,
            True,
        ],
        kind="stable",
    )
    .reset_index(
        drop=True
    )
)


raw_family_summary[
    "local_family_id"
] = [
    f"MG_local_{index:03d}"

    for index
    in range(
        1,
        len(
            raw_family_summary
        )
        + 1,
    )
]


family_id_map = dict(
    zip(
        raw_family_summary[
            "mmseqs_representative"
        ],
        raw_family_summary[
            "local_family_id"
        ],
    )
)


genes[
    "local_family_id"
] = genes[
    "mmseqs_representative"
].map(
    family_id_map
)


membership[
    "local_family_id"
] = membership[
    "mmseqs_representative"
].map(
    family_id_map
)


if genes[
    "local_family_id"
].isna().any():

    fail(
        "Stable local-family mapping failed "
        "for neighborhood genes."
    )


if membership[
    "local_family_id"
].isna().any():

    fail(
        "Stable local-family mapping failed "
        "for membership table."
    )


## ================================================================== ##
## 8. Focal protein -> stable local-family crosswalk
##
## This is one of the most important QC/interpretation tables.
## ================================================================== ##

focal_gene_rows = genes[
    genes[
        "is_focal"
    ]
    ==
    1
][
    [
        "focal_id",
        "genome",
        "focal_protein_id",
        "focal_group",
        "local_family_id",
        "mmseqs_representative",
        "neighborhood_protein_id",
    ]
].copy()


if len(
    focal_gene_rows
) != EXPECTED_FOCALS:

    fail(
        f"Expected {EXPECTED_FOCALS} focal self rows "
        f"but found {len(focal_gene_rows)}."
    )


focal_meta = focals.copy()


focal_meta = focal_meta.rename(
    columns={
        "protein_id":
            "focal_protein_id",
    }
)


## Avoid duplicate focal_group / genome fields from the self rows.

focal_meta_keep = [
    column

    for column
    in focal_meta.columns

    if column
    not in {
        "genome",
        "focal_group",
    }
]


focal_crosswalk = focal_gene_rows.merge(
    focal_meta[
        focal_meta_keep
    ],
    on=[
        "focal_id",
        "focal_protein_id",
    ],
    how="left",
    validate="one_to_one",
)


focal_crosswalk = (
    focal_crosswalk
    .sort_values(
        [
            "focal_group",
            "genome",
            "focal_protein_id",
        ],
        kind="stable",
    )
    .reset_index(
        drop=True
    )
)


focal_crosswalk.to_csv(
    OUT_FOCAL_CROSSWALK,
    sep="\t",
    index=False,
    na_rep="",
)


## ================================================================== ##
## 9. Write stable family membership
## ================================================================== ##

membership = (
    membership
    .sort_values(
        [
            "local_family_id",
            "genome",
            "protein_id",
        ],
        kind="stable",
    )
    .reset_index(
        drop=True
    )
)


membership.to_csv(
    OUT_MEMBERSHIP,
    sep="\t",
    index=False,
    na_rep="",
)


## ================================================================== ##
## 10. Build overall local-family summary
## ================================================================== ##

print()
print("Summarizing local families...")


family_rows = []


for local_family_id, group in genes.groupby(
    "local_family_id",
    sort=False,
):

    raw_family = group[
        "mmseqs_representative"
    ].iloc[0]


    raw_info = raw_family_summary[
        raw_family_summary[
            "local_family_id"
        ]
        ==
        local_family_id
    ].iloc[0]


    positions = pd.to_numeric(
        group[
            "oriented_midpoint_offset_bp"
        ],
        errors="coerce",
    )


    strand_summary = dominant_nonempty(
        group[
            "oriented_neighbor_strand"
        ]
    )


    cog_summary = dominant_nonempty(
        group[
            "neighbor_transferred_cog"
        ]
    )


    product_summary = dominant_nonempty(
        group[
            "neighbor_transferred_product"
        ]
    )


    transfer_summary = dominant_nonempty(
        group[
            "neighbor_annotation_transfer_status"
        ]
    )


    fe_tokens = semicolon_tokens(
        group.get(
            "neighbor_fegenie_HMMs",
            pd.Series(
                dtype=str
            ),
        )
    )


    fe_summary = dominant_nonempty(
        fe_tokens
    )


    focal_rows = group[
        group[
            "is_focal"
        ]
        ==
        1
    ]


    focal_groups = sorted(
        set(
            focal_rows[
                "focal_group"
            ]
        )
    )


    region_groups = sorted(
        set(
            group[
                "focal_group"
            ]
        )
    )


    ## -------------------------------------------------------------- ##
    ## Functional-evidence counts
    ## -------------------------------------------------------------- ##

    fe_positive = (
        pd.to_numeric(
            group.get(
                "neighbor_fegenie_positive",
                0,
            ),
            errors="coerce",
        )
        .fillna(0)
        .astype(int)
    )


    hemes = (
        pd.to_numeric(
            group.get(
                "neighbor_number_of_hemes",
                0,
            ),
            errors="coerce",
        )
        .fillna(0)
    )


    findmehemes = (
        pd.to_numeric(
            group.get(
                "neighbor_findmehemes_positive",
                0,
            ),
            errors="coerce",
        )
        .fillna(0)
        .astype(int)
    )


    integrated = (
        pd.to_numeric(
            group.get(
                "neighbor_is_integrated_candidate",
                0,
            ),
            errors="coerce",
        )
        .fillna(0)
        .astype(int)
    )


    ## -------------------------------------------------------------- ##
    ## Conservative family-level display annotation
    ## -------------------------------------------------------------- ##

    if product_summary[
        "value"
    ]:

        annotation_label = (
            product_summary[
                "value"
            ]
        )

        annotation_source = (
            "transferred_product"
        )


    elif cog_summary[
        "value"
    ]:

        annotation_label = (
            cog_summary[
                "value"
            ]
        )

        annotation_source = (
            "transferred_COG"
        )


    elif fe_summary[
        "value"
    ]:

        annotation_label = (
            fe_summary[
                "value"
            ]
        )

        annotation_source = (
            "FeGenie"
        )


    else:

        annotation_label = ""

        annotation_source = ""


    family_rows.append(
        {
            "local_family_id":
                local_family_id,

            "mmseqs_representative":
                raw_family,

            "member_set_hash":
                raw_info[
                    "member_set_hash"
                ],

            "n_proteins":
                len(
                    group
                ),

            "n_focal_regions":
                group[
                    "focal_id"
                ].nunique(),

            "n_genomes":
                group[
                    "genome"
                ].nunique(),

            "n_focal_proteins":
                int(
                    group[
                        "is_focal"
                    ].sum()
                ),

            "contains_priority_focal":
                int(
                    group[
                        "is_focal"
                    ].sum()
                    >
                    0
                ),

            "focal_groups_represented":
                ";".join(
                    focal_groups
                ),

            "n_focal_groups_represented":
                len(
                    focal_groups
                ),

            "neighborhood_groups_present":
                ";".join(
                    region_groups
                ),

            "n_neighborhood_groups_present":
                len(
                    region_groups
                ),

            "median_oriented_midpoint_bp":
                float(
                    positions.median()
                ),

            "min_oriented_midpoint_bp":
                float(
                    positions.min()
                ),

            "max_oriented_midpoint_bp":
                float(
                    positions.max()
                ),

            "position_span_bp":
                float(
                    positions.max()
                    -
                    positions.min()
                ),

            "position_sd_bp":
                safe_sd(
                    positions.dropna()
                ),

            "dominant_oriented_strand":
                strand_summary[
                    "value"
                ],

            "dominant_oriented_strand_support":
                strand_summary[
                    "support"
                ],

            "oriented_strand_observations":
                strand_summary[
                    "total"
                ],

            "oriented_strand_conservation_fraction":
                (
                    strand_summary[
                        "fraction"
                    ]
                ),

            "n_integrated_candidates":
                int(
                    integrated.sum()
                ),

            "n_fegenie_positive":
                int(
                    fe_positive.sum()
                ),

            "dominant_fegenie_HMM":
                fe_summary[
                    "value"
                ],

            "dominant_fegenie_HMM_support":
                fe_summary[
                    "support"
                ],

            "fegenie_HMM_distribution":
                fe_summary[
                    "distribution"
                ],

            "n_findmehemes_positive":
                int(
                    findmehemes.sum()
                ),

            "mean_heme_count":
                float(
                    hemes.mean()
                ),

            "max_heme_count":
                int(
                    hemes.max()
                ),

            "n_with_transferred_COG":
                int(
                    (
                        group[
                            "neighbor_transferred_cog"
                        ]
                        .str.strip()
                        !=
                        ""
                    ).sum()
                ),

            "dominant_transferred_COG":
                cog_summary[
                    "value"
                ],

            "dominant_COG_support":
                cog_summary[
                    "support"
                ],

            "COG_annotated_members":
                cog_summary[
                    "total"
                ],

            "dominant_COG_fraction_among_annotated":
                cog_summary[
                    "fraction"
                ],

            "COG_distribution":
                cog_summary[
                    "distribution"
                ],

            "n_with_transferred_product":
                int(
                    (
                        group[
                            "neighbor_transferred_product"
                        ]
                        .str.strip()
                        !=
                        ""
                    ).sum()
                ),

            "dominant_transferred_product":
                product_summary[
                    "value"
                ],

            "dominant_product_support":
                product_summary[
                    "support"
                ],

            "product_annotated_members":
                product_summary[
                    "total"
                ],

            "dominant_product_fraction_among_annotated":
                product_summary[
                    "fraction"
                ],

            "annotation_transfer_status_distribution":
                transfer_summary[
                    "distribution"
                ],

            "family_annotation_label":
                annotation_label,

            "family_annotation_source":
                annotation_source,

            "focal_regions_present":
                ";".join(
                    sorted(
                        set(
                            group[
                                "focal_id"
                            ]
                        )
                    )
                ),

            "member_protein_ids":
                ";".join(
                    sorted(
                        set(
                            group[
                                "neighbor_protein_id"
                            ]
                        )
                    )
                ),
        }
    )


family_summary = pd.DataFrame(
    family_rows
)


## Sort explicitly by local-family number.

family_summary[
    "_family_number"
] = family_summary[
    "local_family_id"
].str.extract(
    r"(\d+)$"
)[
    0
].astype(int)


family_summary = (
    family_summary
    .sort_values(
        "_family_number",
        kind="stable",
    )
    .drop(
        columns=[
            "_family_number",
        ]
    )
    .reset_index(
        drop=True
    )
)


family_summary.to_csv(
    OUT_FAMILY_SUMMARY,
    sep="\t",
    index=False,
    na_rep="",
    float_format="%.6f",
)


## ================================================================== ##
## 11. Region x local-family presence table
##
## One row = one focal region x one local family.
##
## Paralogue copy number is retained.
## ================================================================== ##

print()
print("Building focal-region family presence table...")


presence_rows = []


for (
    focal_id,
    local_family_id,
), group in genes.groupby(
    [
        "focal_id",
        "local_family_id",
    ],
    sort=False,
):

    focal_group = group[
        "focal_group"
    ].iloc[0]


    positions = pd.to_numeric(
        group[
            "oriented_midpoint_offset_bp"
        ],
        errors="coerce",
    )


    offsets = pd.to_numeric(
        group[
            "oriented_gene_offset"
        ],
        errors="coerce",
    )


    strand_summary = dominant_nonempty(
        group[
            "oriented_neighbor_strand"
        ]
    )


    presence_rows.append(
        {
            "focal_id":
                focal_id,

            "focal_group":
                focal_group,

            "genome":
                group[
                    "genome"
                ].iloc[0],

            "focal_protein_id":
                group[
                    "focal_protein_id"
                ].iloc[0],

            "local_family_id":
                local_family_id,

            "n_gene_copies":
                len(
                    group
                ),

            "contains_focal_gene":
                int(
                    (
                        group[
                            "is_focal"
                        ]
                        ==
                        1
                    ).any()
                ),

            "median_oriented_gene_offset":
                float(
                    offsets.median()
                ),

            "min_oriented_gene_offset":
                int(
                    offsets.min()
                ),

            "max_oriented_gene_offset":
                int(
                    offsets.max()
                ),

            "median_oriented_midpoint_bp":
                float(
                    positions.median()
                ),

            "min_oriented_midpoint_bp":
                float(
                    positions.min()
                ),

            "max_oriented_midpoint_bp":
                float(
                    positions.max()
                ),

            "dominant_oriented_strand":
                strand_summary[
                    "value"
                ],

            "strand_support_fraction":
                strand_summary[
                    "fraction"
                ],

            "member_proteins":
                ";".join(
                    sorted(
                        group[
                            "neighbor_protein_id"
                        ]
                    )
                ),
        }
    )


region_presence = pd.DataFrame(
    presence_rows
)


region_presence = region_presence.merge(
    family_summary[
        [
            "local_family_id",
            "family_annotation_label",
            "family_annotation_source",

            "dominant_transferred_COG",
            "dominant_transferred_product",

            "n_fegenie_positive",
            "n_findmehemes_positive",
            "max_heme_count",
        ]
    ],
    on="local_family_id",
    how="left",
    validate="many_to_one",
)


region_presence = (
    region_presence
    .sort_values(
        [
            "focal_group",
            "focal_id",
            "median_oriented_midpoint_bp",
            "local_family_id",
        ],
        kind="stable",
    )
    .reset_index(
        drop=True
    )
)


region_presence.to_csv(
    OUT_REGION_PRESENCE,
    sep="\t",
    index=False,
    na_rep="",
    float_format="%.6f",
)


## ================================================================== ##
## 12. Family recurrence within each focal group
##
## This is the main comparison table.
##
## Example:
##
##     Multiheme__Cluster_00121
##
## denominator = 4 focal regions
##
## A family found in all four gets:
##
##     n_regions_with_family = 4
##     region_prevalence     = 1.0
##
## No censoring correction is applied here. Stage-51 observability
## remains available separately and should be considered when a
## family is absent near a contig boundary.
## ================================================================== ##

print()
print("Summarizing family recurrence by focal group...")


group_rows = []


for focal_group, group_info in group_sizes.set_index(
    "focal_group"
).iterrows():

    n_total_regions = int(
        group_info[
            "n_focal_regions_in_group"
        ]
    )


    group_presence = region_presence[
        region_presence[
            "focal_group"
        ]
        ==
        focal_group
    ]


    for local_family_id, fam in group_presence.groupby(
        "local_family_id",
        sort=False,
    ):

        n_regions = fam[
            "focal_id"
        ].nunique()


        n_gene_copies = int(
            fam[
                "n_gene_copies"
            ].sum()
        )


        positions = pd.to_numeric(
            fam[
                "median_oriented_midpoint_bp"
            ],
            errors="coerce",
        )


        family_meta = family_summary[
            family_summary[
                "local_family_id"
            ]
            ==
            local_family_id
        ].iloc[0]


        group_rows.append(
            {
                "focal_group":
                    focal_group,

                "n_focal_regions_in_group":
                    n_total_regions,

                "local_family_id":
                    local_family_id,

                "n_regions_with_family":
                    n_regions,

                "region_prevalence":
                    (
                        n_regions
                        /
                        n_total_regions
                    ),

                "n_genomes_with_family":
                    fam[
                        "genome"
                    ].nunique(),

                "n_gene_copies":
                    n_gene_copies,

                "mean_copies_per_positive_region":
                    (
                        n_gene_copies
                        /
                        n_regions
                    ),

                "median_region_position_bp":
                    float(
                        positions.median()
                    ),

                "min_region_position_bp":
                    float(
                        positions.min()
                    ),

                "max_region_position_bp":
                    float(
                        positions.max()
                    ),

                "region_position_span_bp":
                    float(
                        positions.max()
                        -
                        positions.min()
                    ),

                "n_regions_where_family_is_focal":
                    int(
                        fam[
                            "contains_focal_gene"
                        ].sum()
                    ),

                "family_annotation_label":
                    family_meta[
                        "family_annotation_label"
                    ],

                "family_annotation_source":
                    family_meta[
                        "family_annotation_source"
                    ],

                "dominant_transferred_COG":
                    family_meta[
                        "dominant_transferred_COG"
                    ],

                "dominant_transferred_product":
                    family_meta[
                        "dominant_transferred_product"
                    ],

                "n_fegenie_positive_members":
                    family_meta[
                        "n_fegenie_positive"
                    ],

                "n_findmehemes_positive_members":
                    family_meta[
                        "n_findmehemes_positive"
                    ],

                "max_heme_count":
                    family_meta[
                        "max_heme_count"
                    ],

                "focal_regions_with_family":
                    ";".join(
                        sorted(
                            fam[
                                "focal_id"
                            ].unique()
                        )
                    ),
            }
        )


group_summary = pd.DataFrame(
    group_rows
)


group_summary = (
    group_summary
    .sort_values(
        [
            "focal_group",
            "region_prevalence",
            "n_regions_with_family",
            "local_family_id",
        ],
        ascending=[
            True,
            False,
            False,
            True,
        ],
        kind="stable",
    )
    .reset_index(
        drop=True
    )
)


group_summary.to_csv(
    OUT_GROUP_SUMMARY,
    sep="\t",
    index=False,
    na_rep="",
    float_format="%.6f",
)


## ================================================================== ##
## 13. Build ordered family architecture for every focal region
##
## The signature retains paralogues.
##
## Example:
##
##     MG_local_022 >
##     MG_local_004 >
##     MG_local_004 >
##     MG_local_011
##
## Two copies of the same family therefore remain two genes.
## ================================================================== ##

print()
print("Building ordered region architecture signatures...")


annotation_lookup = dict(
    zip(
        family_summary[
            "local_family_id"
        ],
        family_summary[
            "family_annotation_label"
        ],
    )
)


architecture_rows = []


for focal_id, group in genes.groupby(
    "focal_id",
    sort=False,
):

    group = (
        group
        .sort_values(
            [
                "oriented_midpoint_offset_bp",
                "oriented_gene_offset",
                "neighbor_protein_id",
            ],
            kind="stable",
        )
        .reset_index(
            drop=True
        )
    )


    families = group[
        "local_family_id"
    ].tolist()


    proteins = group[
        "neighbor_protein_id"
    ].tolist()


    offsets = [
        str(
            int(value)
        )

        for value
        in group[
            "oriented_gene_offset"
        ]
    ]


    annotations = [
        annotation_lookup.get(
            family,
            ""
        )

        for family
        in families
    ]


    architecture_rows.append(
        {
            "focal_id":
                focal_id,

            "focal_group":
                group[
                    "focal_group"
                ].iloc[0],

            "genome":
                group[
                    "genome"
                ].iloc[0],

            "focal_protein_id":
                group[
                    "focal_protein_id"
                ].iloc[0],

            "n_genes":
                len(
                    group
                ),

            "n_local_families":
                group[
                    "local_family_id"
                ].nunique(),

            "n_paralogue_excess":
                (
                    len(
                        group
                    )
                    -
                    group[
                        "local_family_id"
                    ].nunique()
                ),

            "ordered_local_family_signature":
                ">".join(
                    families
                ),

            "ordered_protein_ids":
                ">".join(
                    proteins
                ),

            "ordered_gene_offsets":
                ">".join(
                    offsets
                ),

            "ordered_family_annotations":
                ">".join(
                    annotations
                ),
        }
    )


architectures = pd.DataFrame(
    architecture_rows
)


architectures = (
    architectures
    .sort_values(
        [
            "focal_group",
            "genome",
            "focal_protein_id",
        ],
        kind="stable",
    )
    .reset_index(
        drop=True
    )
)


architectures.to_csv(
    OUT_ARCHITECTURES,
    sep="\t",
    index=False,
    na_rep="",
)


## ================================================================== ##
## 14. Add family-level annotation to every actual gene
##
## This is the file we will use later for visualization.
## ================================================================== ##

family_annotation_columns = [
    "local_family_id",

    "member_set_hash",

    "n_proteins",
    "n_focal_regions",
    "n_genomes",

    "contains_priority_focal",
    "focal_groups_represented",

    "family_annotation_label",
    "family_annotation_source",

    "dominant_transferred_COG",
    "dominant_COG_fraction_among_annotated",

    "dominant_transferred_product",
    "dominant_product_fraction_among_annotated",

    "n_fegenie_positive",
    "dominant_fegenie_HMM",

    "n_findmehemes_positive",
    "mean_heme_count",
    "max_heme_count",
]


gene_output = genes.merge(
    family_summary[
        family_annotation_columns
    ],
    on="local_family_id",
    how="left",
    validate="many_to_one",
    suffixes=(
        "",
        "_family",
    ),
)


if len(
    gene_output
) != EXPECTED_NEIGHBORHOOD_ROWS:

    fail(
        "Family-summary merge changed gene row count."
    )


if gene_output[
    "family_annotation_label"
].isna().any():

    fail(
        "At least one actual gene failed family-summary mapping."
    )


gene_output = (
    gene_output
    .sort_values(
        [
            "focal_group",
            "focal_id",
            "oriented_midpoint_offset_bp",
            "oriented_gene_offset",
            "neighbor_protein_id",
        ],
        kind="stable",
    )
    .reset_index(
        drop=True
    )
)


gene_output.to_csv(
    OUT_GENE_TABLE,
    sep="\t",
    index=False,
    na_rep="",
    float_format="%.6f",
)


## ================================================================== ##
## 15. QC and useful headline summaries
## ================================================================== ##

n_stable_families = family_summary[
    "local_family_id"
].nunique()


n_singletons = int(
    (
        family_summary[
            "n_proteins"
        ]
        ==
        1
    ).sum()
)


n_recurrent_regions = int(
    (
        family_summary[
            "n_focal_regions"
        ]
        >=
        2
    ).sum()
)


n_unannotated_families = int(
    (
        family_summary[
            "family_annotation_label"
        ]
        .str.strip()
        ==
        ""
    ).sum()
)


n_unannotated_recurrent = int(
    (
        (
            family_summary[
                "family_annotation_label"
            ]
            .str.strip()
            ==
            ""
        )
        &
        (
            family_summary[
                "n_focal_regions"
            ]
            >=
            2
        )
    ).sum()
)


## -------------------------------------------------------------- ##
## Do all focal proteins belonging to the same old top-family
## group map consistently to the same local family?
##
## This is descriptive QC, not a forced requirement.
## -------------------------------------------------------------- ##

focal_group_family_counts = (
    focal_crosswalk
    .groupby(
        "focal_group",
        as_index=False,
    )
    .agg(
        n_focal_proteins=(
            "focal_id",
            "nunique",
        ),

        n_local_families_among_focals=(
            "local_family_id",
            "nunique",
        ),

        focal_local_families=(
            "local_family_id",
            lambda x:
                ";".join(
                    sorted(
                        set(
                            x
                        )
                    )
                ),
        ),
    )
)


print()
print("Focal-group sequence-family check")


for row in focal_group_family_counts.itertuples(
    index=False
):

    print(
        f"  {row.focal_group:<40} "
        f"{row.n_focal_proteins:>2} focals -> "
        f"{row.n_local_families_among_focals:>2} local families "
        f"({row.focal_local_families})"
    )


qc_rows = [
    (
        "input_proteins",
        len(
            membership
        ),
    ),

    (
        "raw_MMseqs_families",
        n_raw_families,
    ),

    (
        "stable_local_families",
        n_stable_families,
    ),

    (
        "family_singletons",
        n_singletons,
    ),

    (
        "families_present_in_ge2_focal_regions",
        n_recurrent_regions,
    ),

    (
        "families_without_display_annotation",
        n_unannotated_families,
    ),

    (
        "recurrent_families_without_display_annotation",
        n_unannotated_recurrent,
    ),

    (
        "focal_regions",
        focals[
            "focal_id"
        ].nunique(),
    ),

    (
        "focal_groups",
        focals[
            "focal_group"
        ].nunique(),
    ),

    (
        "gene_level_output_rows",
        len(
            gene_output
        ),
    ),
]


qc = pd.DataFrame(
    qc_rows,
    columns=[
        "metric",
        "value",
    ],
)


qc.to_csv(
    OUT_QC,
    sep="\t",
    index=False,
)


## ================================================================== ##
## Final report
## ================================================================== ##

print()
print("=" * 80)
print("STAGE 58 COMPLETE")
print("=" * 80)

print(
    f"Actual neighborhood genes:             "
    f"{len(gene_output):,}"
)

print(
    f"Stable local families:                 "
    f"{n_stable_families:,}"
)

print(
    f"Singleton local families:              "
    f"{n_singletons:,}"
)

print(
    f"Families in >=2 focal regions:         "
    f"{n_recurrent_regions:,}"
)

print(
    f"Unannotated families:                  "
    f"{n_unannotated_families:,}"
)

print(
    f"Recurrent unannotated families:        "
    f"{n_unannotated_recurrent:,}"
)

print(
    f"Focal groups:                          "
    f"{focals['focal_group'].nunique():,}"
)

print()
print(
    f"Membership:\n  "
    f"{OUT_MEMBERSHIP}"
)

print(
    f"Family summary:\n  "
    f"{OUT_FAMILY_SUMMARY}"
)

print(
    f"Focal-family crosswalk:\n  "
    f"{OUT_FOCAL_CROSSWALK}"
)

print(
    f"Region presence:\n  "
    f"{OUT_REGION_PRESENCE}"
)

print(
    f"Family x focal-group summary:\n  "
    f"{OUT_GROUP_SUMMARY}"
)

print(
    f"Region architectures:\n  "
    f"{OUT_ARCHITECTURES}"
)

print(
    f"Gene-level table for later visualization:\n  "
    f"{OUT_GENE_TABLE}"
)

print(
    f"QC:\n  "
    f"{OUT_QC}"
)
