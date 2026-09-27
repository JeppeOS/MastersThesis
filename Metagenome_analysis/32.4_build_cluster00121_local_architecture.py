#!/usr/bin/env python3

from __future__ import annotations

import hashlib
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd


## ================================================================== ##
## STAGE 65
##
## BUILD COMBINED LOCAL-FAMILY ARCHITECTURE FOR CLUSTER_00121
##
## Inputs
## ------
##
## 23 focal regions:
##
##     19 GlobDB
##      4 metagenome
##
## 686 actual neighborhood proteins
##
## MMseqs clustering:
##
##     40% minimum identity
##     80% coverage
##     cov-mode 0
##
##
## Goals
## -----
##
## 1. Convert raw MMseqs representatives into stable:
##
##        C121_local_###
##
##    family identifiers.
##
## 2. Map these families back to every actual gene.
##
## 3. Quantify family recurrence across:
##
##        GlobDB loci
##        metagenome loci
##        all 23 loci
##
## 4. Summarize COG/product evidence.
##
## 5. Build ordered architecture signatures.
##
## 6. Identify recurrent unannotated families that may merit
##    targeted structural follow-up.
##
##
## IMPORTANT
## ---------
##
## The focal Cluster_00121 family itself is NOT forced into one family.
## We report whatever the sequence clustering actually produced.
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


WORK = (
    PROJECT
    / "comparative_analysis"
    / "Cluster_00121"
)


LOCAL_DIR = (
    WORK
    / "local_clustering"
)


NEIGHBORHOODS = (
    WORK
    / "Cluster_00121_combined_neighborhood_genes.tsv"
)


PROTEIN_METADATA = (
    WORK
    / "Cluster_00121_combined_neighborhood_protein_metadata.tsv"
)


FOCAL_REGIONS = (
    WORK
    / "Cluster_00121_focal_regions.tsv"
)


RAW_CLUSTERS = (
    LOCAL_DIR
    / "C121_local40_cov80_cluster.tsv"
)


## ================================================================== ##
## Outputs
## ================================================================== ##

OUT_MEMBERSHIP = (
    LOCAL_DIR
    / "C121_local_family_membership.tsv"
)


OUT_SUMMARY = (
    LOCAL_DIR
    / "C121_local_family_summary.tsv"
)


OUT_REGION_PRESENCE = (
    LOCAL_DIR
    / "C121_local_family_region_presence.tsv"
)


OUT_DATASET_COMPARISON = (
    LOCAL_DIR
    / "C121_local_family_dataset_comparison.tsv"
)


OUT_ARCHITECTURES = (
    LOCAL_DIR
    / "C121_region_architectures.tsv"
)


OUT_GENES = (
    WORK
    / "Cluster_00121_combined_neighborhood_genes_local_families.tsv"
)


OUT_UNKNOWN = (
    LOCAL_DIR
    / "C121_recurrent_unknown_families.tsv"
)


OUT_FOCAL_CROSSWALK = (
    LOCAL_DIR
    / "C121_focal_family_crosswalk.tsv"
)


OUT_QC = (
    LOCAL_DIR
    / "C121_local_architecture_qc.tsv"
)


## ================================================================== ##
## Constants
## ================================================================== ##

EXPECTED_PROTEINS = 686

EXPECTED_RAW_FAMILIES = 559

EXPECTED_REGIONS = 23

EXPECTED_GLOBDB_REGIONS = 19

EXPECTED_METAGENOME_REGIONS = 4


## ================================================================== ##
## Helpers
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
            f"Required file missing:\n{path}"
        )


def clean(value):

    if value is None:

        return ""

    if pd.isna(value):

        return ""

    return str(
        value
    ).strip()


def ordered_unique(values):

    seen = set()

    result = []


    for value in values:

        value = clean(
            value
        )


        if not value:

            continue


        if value in seen:

            continue


        seen.add(
            value
        )

        result.append(
            value
        )


    return result


def join_unique(
    values,
    separator=";",
):

    return separator.join(
        ordered_unique(
            values
        )
    )


def dominant_nonempty(values):

    values = [
        clean(
            value
        )

        for value
        in values

        if clean(
            value
        )
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


    value, support = ordered[
        0
    ]


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


def member_signature(values):

    members = sorted(
        set(
            clean(
                value
            )

            for value
            in values

            if clean(
                value
            )
        )
    )


    signature = "|".join(
        members
    )


    digest = hashlib.sha256(
        signature.encode(
            "utf-8"
        )
    ).hexdigest()[
        :12
    ]


    return signature, digest


def informative_product(value):

    value = clean(
        value
    )


    if not value:

        return False


    lower = value.lower()


    generic_exact = {
        "hypothetical protein",
        "uncharacterized protein",
        "unknown protein",
        "predicted protein",
        "conserved hypothetical protein",
        "putative protein",
    }


    if lower in generic_exact:

        return False


    generic_fragments = [
        "hypothetical protein",
        "uncharacterized protein",
        "protein of unknown function",
        "unknown function",
    ]


    if any(
        fragment in lower

        for fragment
        in generic_fragments
    ):

        return False


    return True


## ================================================================== ##
## Validate inputs
## ================================================================== ##

print("=" * 80)

print(
    "STAGE 65 - CLUSTER_00121 LOCAL ARCHITECTURE"
)

print("=" * 80)


for path in [
    NEIGHBORHOODS,
    PROTEIN_METADATA,
    FOCAL_REGIONS,
    RAW_CLUSTERS,
]:

    require_file(
        path
    )


## ================================================================== ##
## 1. Read actual-gene neighborhood table
## ================================================================== ##

print()
print(
    "Reading combined neighborhood genes..."
)


genes = pd.read_csv(
    NEIGHBORHOODS,
    sep="\t",
    dtype=str,
    keep_default_na=False,
)


required_gene_columns = {
    "source_dataset",
    "focal_id",
    "genome",

    "focal_protein_id",
    "neighbor_protein_id",

    "combined_protein_id",

    "oriented_gene_offset",
    "oriented_midpoint_offset_bp",
    "oriented_neighbor_strand",

    "is_focal",

    "comparison_cog",
    "comparison_product",
    "comparison_annotation_source",

    "fegenie_positive",
    "fegenie_HMMs",

    "findmehemes_positive",
    "number_of_hemes",
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
        "Combined neighborhood table is missing:\n"
        +
        "\n".join(
            sorted(
                missing
            )
        )
    )


if genes[
    "focal_id"
].nunique() != EXPECTED_REGIONS:

    fail(
        f"Expected {EXPECTED_REGIONS} focal regions; "
        f"found {genes['focal_id'].nunique()}."
    )


dataset_region_counts = (
    genes[
        [
            "source_dataset",
            "focal_id",
        ]
    ]
    .drop_duplicates()
    [
        "source_dataset"
    ]
    .value_counts()
)


if int(
    dataset_region_counts.get(
        "GlobDB",
        0,
    )
) != EXPECTED_GLOBDB_REGIONS:

    fail(
        "Unexpected number of GlobDB regions."
    )


if int(
    dataset_region_counts.get(
        "Metagenome",
        0,
    )
) != EXPECTED_METAGENOME_REGIONS:

    fail(
        "Unexpected number of metagenome regions."
    )


## ================================================================== ##
## 2. Read protein metadata
## ================================================================== ##

print()
print(
    "Reading unique protein metadata..."
)


protein_metadata = pd.read_csv(
    PROTEIN_METADATA,
    sep="\t",
    dtype=str,
    keep_default_na=False,
)


if len(
    protein_metadata
) != EXPECTED_PROTEINS:

    fail(
        f"Expected {EXPECTED_PROTEINS} unique proteins; "
        f"found {len(protein_metadata)}."
    )


if protein_metadata[
    "combined_protein_id"
].duplicated().any():

    fail(
        "Duplicate combined_protein_id values."
    )


## ================================================================== ##
## 3. Read raw MMseqs clustering
##
## easy-cluster TSV:
##
##     representative    member
##
## ================================================================== ##

print()
print(
    "Reading Stage-64 MMseqs clustering..."
)


raw = pd.read_csv(
    RAW_CLUSTERS,
    sep="\t",
    header=None,
    names=[
        "raw_representative",
        "combined_protein_id",
    ],
    dtype=str,
    keep_default_na=False,
)


if len(
    raw
) != EXPECTED_PROTEINS:

    fail(
        f"Expected {EXPECTED_PROTEINS} MMseqs assignments; "
        f"found {len(raw)}."
    )


if raw[
    "combined_protein_id"
].duplicated().any():

    fail(
        "A protein has multiple raw MMseqs assignments."
    )


n_raw_families = raw[
    "raw_representative"
].nunique()


if n_raw_families != EXPECTED_RAW_FAMILIES:

    fail(
        f"Expected {EXPECTED_RAW_FAMILIES} raw families; "
        f"found {n_raw_families}."
    )


protein_ids = set(
    protein_metadata[
        "combined_protein_id"
    ]
)


cluster_ids = set(
    raw[
        "combined_protein_id"
    ]
)


if protein_ids != cluster_ids:

    fail(
        "Protein IDs in metadata and MMseqs clustering "
        "do not match exactly."
    )


print(
    f"  Protein assignments: "
    f"{len(raw)}"
)

print(
    f"  Raw families:        "
    f"{n_raw_families}"
)


## ================================================================== ##
## 4. Build protein-level evidence table
##
## A physical protein should occur once in this particular dataset,
## but we aggregate defensively rather than assuming that forever.
## ================================================================== ##

protein_evidence_rows = []


for combined_id, group in genes.groupby(
    "combined_protein_id",
    sort=False,
):

    first = group.iloc[
        0
    ]


    protein_evidence_rows.append(
        {
            "combined_protein_id":
                combined_id,

            "source_dataset":
                clean(
                    first[
                        "source_dataset"
                    ]
                ),

            "genome":
                clean(
                    first[
                        "genome"
                    ]
                ),

            "protein_id":
                clean(
                    first[
                        "neighbor_protein_id"
                    ]
                ),

            "comparison_cog":
                join_unique(
                    group[
                        "comparison_cog"
                    ]
                ),

            "comparison_product":
                join_unique(
                    group[
                        "comparison_product"
                    ],
                    separator=" || ",
                ),

            "comparison_annotation_source":
                join_unique(
                    group[
                        "comparison_annotation_source"
                    ]
                ),

            "fegenie_positive":
                int(
                    pd.to_numeric(
                        group[
                            "fegenie_positive"
                        ],
                        errors="coerce",
                    )
                    .fillna(0)
                    .max()
                    ==
                    1
                ),

            "fegenie_HMMs":
                join_unique(
                    group[
                        "fegenie_HMMs"
                    ]
                ),

            "findmehemes_positive":
                int(
                    pd.to_numeric(
                        group[
                            "findmehemes_positive"
                        ],
                        errors="coerce",
                    )
                    .fillna(0)
                    .max()
                    ==
                    1
                ),

            "number_of_hemes":
                pd.to_numeric(
                    group[
                        "number_of_hemes"
                    ],
                    errors="coerce",
                ).max(),

            "n_focal_regions_observed":
                group[
                    "focal_id"
                ].nunique(),
        }
    )


protein_evidence = pd.DataFrame(
    protein_evidence_rows
)


if len(
    protein_evidence
) != EXPECTED_PROTEINS:

    fail(
        f"Expected {EXPECTED_PROTEINS} protein-evidence rows; "
        f"found {len(protein_evidence)}."
    )


## ================================================================== ##
## 5. Attach raw families to proteins
## ================================================================== ##

membership = (
    raw
    .merge(
        protein_metadata,
        on="combined_protein_id",
        how="left",
        validate="one_to_one",
    )
    .merge(
        protein_evidence,
        on=[
            "combined_protein_id",
            "source_dataset",
            "genome",
            "protein_id",
        ],
        how="left",
        validate="one_to_one",
    )
)


## ================================================================== ##
## 6. Calculate raw-family statistics needed for deterministic IDs
## ================================================================== ##

family_sort_rows = []


for raw_family, group in membership.groupby(
    "raw_representative",
    sort=False,
):

    members = group[
        "combined_protein_id"
    ].tolist()


    region_subset = (
        genes[
            genes[
                "combined_protein_id"
            ].isin(
                members
            )
        ]
        [
            [
                "focal_id",
                "source_dataset",
            ]
        ]
        .drop_duplicates()
    )


    signature, digest = member_signature(
        members
    )


    family_sort_rows.append(
        {
            "raw_representative":
                raw_family,

            "n_proteins":
                len(
                    group
                ),

            "n_genomes":
                group[
                    "genome"
                ].nunique(),

            "n_focal_regions":
                region_subset[
                    "focal_id"
                ].nunique(),

            "n_GlobDB_regions":
                region_subset[
                    region_subset[
                        "source_dataset"
                    ]
                    ==
                    "GlobDB"
                ]
                [
                    "focal_id"
                ].nunique(),

            "n_metagenome_regions":
                region_subset[
                    region_subset[
                        "source_dataset"
                    ]
                    ==
                    "Metagenome"
                ]
                [
                    "focal_id"
                ].nunique(),

            "member_signature":
                signature,

            "member_signature_sha256_12":
                digest,
        }
    )


family_sort = pd.DataFrame(
    family_sort_rows
)


## Stable IDs:
##
## highest focal-region prevalence first,
## then genome prevalence,
## then family size,
## then exact member signature.
family_sort = (
    family_sort
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


family_sort[
    "local_family_id"
] = [
    f"C121_local_{i:03d}"

    for i in range(
        1,
        len(
            family_sort
        )
        + 1,
    )
]


raw_to_stable = dict(
    zip(
        family_sort[
            "raw_representative"
        ],
        family_sort[
            "local_family_id"
        ],
    )
)


membership[
    "local_family_id"
] = membership[
    "raw_representative"
].map(
    raw_to_stable
)


if membership[
    "local_family_id"
].isna().any():

    fail(
        "Failed to assign one or more stable family IDs."
    )


## ================================================================== ##
## 7. Family summaries
## ================================================================== ##

print()
print(
    "Summarizing stable local families..."
)


summary_rows = []


for family_id, group in membership.groupby(
    "local_family_id",
    sort=False,
):

    raw_family = clean(
        group.iloc[
            0
        ][
            "raw_representative"
        ]
    )


    sort_meta = family_sort[
        family_sort[
            "local_family_id"
        ]
        ==
        family_id
    ].iloc[
        0
    ]


    globdb_group = group[
        group[
            "source_dataset"
        ]
        ==
        "GlobDB"
    ]


    mg_group = group[
        group[
            "source_dataset"
        ]
        ==
        "Metagenome"
    ]


    cog = dominant_nonempty(
        group[
            "comparison_cog"
        ]
    )


    product = dominant_nonempty(
        group[
            "comparison_product"
        ]
    )


    informative_products = [
        clean(
            value
        )

        for value
        in group[
            "comparison_product"
        ]

        if informative_product(
            value
        )
    ]


    informative_product_summary = dominant_nonempty(
        informative_products
    )


    hmm = dominant_nonempty(
        group[
            "fegenie_HMMs"
        ]
    )


    n_regions = int(
        sort_meta[
            "n_focal_regions"
        ]
    )


    n_globdb_regions = int(
        sort_meta[
            "n_GlobDB_regions"
        ]
    )


    n_mg_regions = int(
        sort_meta[
            "n_metagenome_regions"
        ]
    )


    if (
        n_globdb_regions > 0
        and
        n_mg_regions > 0
    ):

        dataset_distribution = (
            "shared_GlobDB_metagenome"
        )


    elif n_globdb_regions > 0:

        dataset_distribution = (
            "GlobDB_only"
        )


    else:

        dataset_distribution = (
            "metagenome_only"
        )


    has_cog = bool(
        cog[
            "value"
        ]
    )


    has_informative_product = bool(
        informative_product_summary[
            "value"
        ]
    )


    functionally_unannotated = int(
        not has_cog
        and
        not has_informative_product
    )


    recurrent = int(
        n_regions
        >=
        2
    )


    recurrent_unknown = int(
        recurrent
        ==
        1
        and
        functionally_unannotated
        ==
        1
    )


    if recurrent_unknown:

        if dataset_distribution == "shared_GlobDB_metagenome":

            structural_followup_class = (
                "shared_recurrent_unknown"
            )


        elif dataset_distribution == "GlobDB_only":

            structural_followup_class = (
                "GlobDB_recurrent_unknown"
            )


        else:

            structural_followup_class = (
                "metagenome_recurrent_unknown"
            )


    else:

        structural_followup_class = ""


    hemes = pd.to_numeric(
        group[
            "number_of_hemes"
        ],
        errors="coerce",
    )


    summary_rows.append(
        {
            "local_family_id":
                family_id,

            "raw_representative":
                raw_family,

            "member_signature_sha256_12":
                sort_meta[
                    "member_signature_sha256_12"
                ],

            "n_proteins":
                len(
                    group
                ),

            "n_genomes":
                group[
                    "genome"
                ].nunique(),

            "n_focal_regions":
                n_regions,

            "prevalence_all_23":
                n_regions
                /
                EXPECTED_REGIONS,

            "n_GlobDB_proteins":
                len(
                    globdb_group
                ),

            "n_metagenome_proteins":
                len(
                    mg_group
                ),

            "n_GlobDB_regions":
                n_globdb_regions,

            "prevalence_GlobDB_19":
                n_globdb_regions
                /
                EXPECTED_GLOBDB_REGIONS,

            "n_metagenome_regions":
                n_mg_regions,

            "prevalence_metagenome_4":
                n_mg_regions
                /
                EXPECTED_METAGENOME_REGIONS,

            "dataset_distribution":
                dataset_distribution,

            "dominant_cog":
                cog[
                    "value"
                ],

            "dominant_cog_support":
                cog[
                    "support"
                ],

            "dominant_cog_fraction_annotated":
                cog[
                    "fraction"
                ],

            "cog_distribution":
                cog[
                    "distribution"
                ],

            "dominant_product":
                product[
                    "value"
                ],

            "dominant_informative_product":
                informative_product_summary[
                    "value"
                ],

            "all_annotation_sources":
                join_unique(
                    group[
                        "comparison_annotation_source"
                    ]
                ),

            "n_fegenie_positive":
                int(
                    pd.to_numeric(
                        group[
                            "fegenie_positive"
                        ],
                        errors="coerce",
                    )
                    .fillna(0)
                    .sum()
                ),

            "dominant_fegenie_HMM":
                hmm[
                    "value"
                ],

            "n_findmehemes_positive":
                int(
                    pd.to_numeric(
                        group[
                            "findmehemes_positive"
                        ],
                        errors="coerce",
                    )
                    .fillna(0)
                    .sum()
                ),

            "mean_heme_count":
                (
                    float(
                        hemes.mean()
                    )
                    if hemes.notna().any()
                    else np.nan
                ),

            "max_heme_count":
                (
                    float(
                        hemes.max()
                    )
                    if hemes.notna().any()
                    else np.nan
                ),

            "recurrent_family":
                recurrent,

            "functionally_unannotated":
                functionally_unannotated,

            "recurrent_unknown_family":
                recurrent_unknown,

            "structural_followup_class":
                structural_followup_class,
        }
    )


family_summary = pd.DataFrame(
    summary_rows
)


## Preserve deterministic stable-ID ordering.
family_summary[
    "_family_number"
] = (
    family_summary[
        "local_family_id"
    ]
    .str.extract(
        r"(\d+)$"
    )[
        0
    ]
    .astype(int)
)


family_summary = (
    family_summary
    .sort_values(
        "_family_number"
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
    OUT_SUMMARY,
    sep="\t",
    index=False,
    na_rep="",
    float_format="%.6f",
)


## ================================================================== ##
## 8. Write protein membership
## ================================================================== ##

membership = (
    membership
    .sort_values(
        [
            "local_family_id",
            "source_dataset",
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
## 9. Map stable families back to actual-gene neighborhoods
## ================================================================== ##

genes_local = genes.merge(
    membership[
        [
            "combined_protein_id",
            "local_family_id",
            "raw_representative",
        ]
    ],
    on="combined_protein_id",
    how="left",
    validate="many_to_one",
)


if genes_local[
    "local_family_id"
].isna().any():

    fail(
        "One or more neighborhood genes lack a local family."
    )


genes_local.to_csv(
    OUT_GENES,
    sep="\t",
    index=False,
    na_rep="",
)


## ================================================================== ##
## 10. Family x region presence
## ================================================================== ##

presence_rows = []


for (
    focal_id,
    local_family_id,
), group in genes_local.groupby(
    [
        "focal_id",
        "local_family_id",
    ],
    sort=False,
):

    first = group.iloc[
        0
    ]


    offsets = pd.to_numeric(
        group[
            "oriented_gene_offset"
        ],
        errors="coerce",
    )


    midpoint_offsets = pd.to_numeric(
        group[
            "oriented_midpoint_offset_bp"
        ],
        errors="coerce",
    )


    presence_rows.append(
        {
            "focal_id":
                focal_id,

            "source_dataset":
                first[
                    "source_dataset"
                ],

            "genome":
                first[
                    "genome"
                ],

            "focal_protein_id":
                first[
                    "focal_protein_id"
                ],

            "local_family_id":
                local_family_id,

            "n_family_genes":
                len(
                    group
                ),

            "min_oriented_gene_offset":
                offsets.min(),

            "max_oriented_gene_offset":
                offsets.max(),

            "mean_oriented_gene_offset":
                offsets.mean(),

            "min_oriented_midpoint_offset_bp":
                midpoint_offsets.min(),

            "max_oriented_midpoint_offset_bp":
                midpoint_offsets.max(),

            "contains_focal_gene":
                int(
                    pd.to_numeric(
                        group[
                            "is_focal"
                        ],
                        errors="coerce",
                    )
                    .fillna(0)
                    .max()
                    ==
                    1
                ),
        }
    )


region_presence = pd.DataFrame(
    presence_rows
)


region_presence.to_csv(
    OUT_REGION_PRESENCE,
    sep="\t",
    index=False,
    na_rep="",
    float_format="%.6f",
)


## ================================================================== ##
## 11. Dataset-comparison table
##
## This is the useful table for asking:
##
##     Is this family conserved in GlobDB?
##     Is it also recovered in our MAGs?
## ================================================================== ##

dataset_comparison = family_summary[
    [
        "local_family_id",

        "n_proteins",
        "n_genomes",
        "n_focal_regions",

        "n_GlobDB_regions",
        "prevalence_GlobDB_19",

        "n_metagenome_regions",
        "prevalence_metagenome_4",

        "dataset_distribution",

        "dominant_cog",
        "dominant_informative_product",

        "n_fegenie_positive",
        "n_findmehemes_positive",

        "mean_heme_count",
        "max_heme_count",

        "recurrent_family",
        "functionally_unannotated",
        "recurrent_unknown_family",

        "structural_followup_class",
    ]
].copy()


dataset_comparison.to_csv(
    OUT_DATASET_COMPARISON,
    sep="\t",
    index=False,
    na_rep="",
    float_format="%.6f",
)


## ================================================================== ##
## 12. Ordered architecture signatures
##
## Family identity is retained WITH PARALOGUES.
##
## Example:
##
## C121_local_004(+) ->
## C121_local_001(+) ->
## C121_local_004(+)
##
## We do not collapse repeated genes.
## ================================================================== ##

architecture_rows = []


for focal_id, group in genes_local.groupby(
    "focal_id",
    sort=False,
):

    first = group.iloc[
        0
    ]


    ordered = group.copy()


    ordered[
        "_offset"
    ] = pd.to_numeric(
        ordered[
            "oriented_midpoint_offset_bp"
        ],
        errors="coerce",
    )


    ordered = ordered.sort_values(
        [
            "_offset",
            "local_family_id",
            "neighbor_protein_id",
        ],
        kind="stable",
    )


    family_architecture = " -> ".join(
        (
            clean(
                row.local_family_id
            )
            +
            "("
            +
            clean(
                row.oriented_neighbor_strand
            )
            +
            ")"
        )

        for row
        in ordered.itertuples(
            index=False
        )
    )


    cog_architecture = " -> ".join(
        (
            clean(
                row.comparison_cog
            )
            if clean(
                row.comparison_cog
            )
            else
            "NA"
        )

        for row
        in ordered.itertuples(
            index=False
        )
    )


    architecture_rows.append(
        {
            "focal_id":
                focal_id,

            "source_dataset":
                first[
                    "source_dataset"
                ],

            "genome":
                first[
                    "genome"
                ],

            "focal_protein_id":
                first[
                    "focal_protein_id"
                ],

            "n_genes":
                len(
                    ordered
                ),

            "n_local_families":
                ordered[
                    "local_family_id"
                ].nunique(),

            "family_architecture":
                family_architecture,

            "cog_architecture":
                cog_architecture,
        }
    )


architectures = pd.DataFrame(
    architecture_rows
)


architectures.to_csv(
    OUT_ARCHITECTURES,
    sep="\t",
    index=False,
)


## ================================================================== ##
## 13. Focal-family crosswalk
##
## Critical sanity check:
##
## Do the 19 old + 4 new Cluster_00121 focal proteins all remain
## in one common local sequence family?
##
## We REPORT rather than force the answer.
## ================================================================== ##

focal_crosswalk = (
    genes_local[
        pd.to_numeric(
            genes_local[
                "is_focal"
            ],
            errors="coerce",
        )
        .fillna(0)
        .astype(int)
        ==
        1
    ]
    [
        [
            "focal_id",
            "source_dataset",
            "genome",
            "focal_protein_id",
            "combined_protein_id",
            "local_family_id",
            "number_of_hemes",
            "comparison_cog",
            "comparison_product",
        ]
    ]
    .drop_duplicates()
    .sort_values(
        [
            "local_family_id",
            "source_dataset",
            "genome",
        ],
        kind="stable",
    )
    .reset_index(
        drop=True
    )
)


if len(
    focal_crosswalk
) != EXPECTED_REGIONS:

    fail(
        f"Expected {EXPECTED_REGIONS} focal crosswalk rows; "
        f"found {len(focal_crosswalk)}."
    )


focal_crosswalk.to_csv(
    OUT_FOCAL_CROSSWALK,
    sep="\t",
    index=False,
)


## ================================================================== ##
## 14. Recurrent unknowns for later structural consideration
## ================================================================== ##

unknown = (
    family_summary[
        family_summary[
            "recurrent_unknown_family"
        ]
        ==
        1
    ]
    .copy()
)


## Put shared GlobDB+MAG unknowns first because they are immediately
## relevant to the comparative question.
distribution_order = {
    "shared_GlobDB_metagenome":
        0,

    "GlobDB_only":
        1,

    "metagenome_only":
        2,
}


unknown[
    "_distribution_order"
] = unknown[
    "dataset_distribution"
].map(
    distribution_order
).fillna(
    99
)


unknown = (
    unknown
    .sort_values(
        [
            "_distribution_order",
            "n_focal_regions",
            "n_genomes",
            "n_proteins",
            "local_family_id",
        ],
        ascending=[
            True,
            False,
            False,
            False,
            True,
        ],
        kind="stable",
    )
    .drop(
        columns=[
            "_distribution_order",
        ]
    )
    .reset_index(
        drop=True
    )
)


unknown.to_csv(
    OUT_UNKNOWN,
    sep="\t",
    index=False,
    na_rep="",
    float_format="%.6f",
)


## ================================================================== ##
## 15. Headline QC
## ================================================================== ##

n_singletons = int(
    (
        family_summary[
            "n_proteins"
        ]
        ==
        1
    ).sum()
)


n_recurrent = int(
    (
        family_summary[
            "n_focal_regions"
        ]
        >=
        2
    ).sum()
)


n_shared = int(
    (
        family_summary[
            "dataset_distribution"
        ]
        ==
        "shared_GlobDB_metagenome"
    ).sum()
)


n_shared_recurrent = int(
    (
        (
            family_summary[
                "dataset_distribution"
            ]
            ==
            "shared_GlobDB_metagenome"
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


n_recurrent_unknown = int(
    family_summary[
        "recurrent_unknown_family"
    ].sum()
)


n_shared_recurrent_unknown = int(
    (
        family_summary[
            "structural_followup_class"
        ]
        ==
        "shared_recurrent_unknown"
    ).sum()
)


n_focal_families = focal_crosswalk[
    "local_family_id"
].nunique()


qc = pd.DataFrame(
    [
        (
            "neighborhood_proteins",
            EXPECTED_PROTEINS,
        ),

        (
            "stable_local_families",
            len(
                family_summary
            ),
        ),

        (
            "singleton_families",
            n_singletons,
        ),

        (
            "families_in_ge2_regions",
            n_recurrent,
        ),

        (
            "families_shared_GlobDB_metagenome",
            n_shared,
        ),

        (
            "shared_families_in_ge2_regions",
            n_shared_recurrent,
        ),

        (
            "recurrent_unknown_families",
            n_recurrent_unknown,
        ),

        (
            "shared_recurrent_unknown_families",
            n_shared_recurrent_unknown,
        ),

        (
            "distinct_local_families_among_23_focal_proteins",
            n_focal_families,
        ),
    ],
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
## Terminal report
## ================================================================== ##

print()
print(
    "Focal Cluster_00121 sequence-family check"
)


focal_family_counts = (
    focal_crosswalk
    .groupby(
        "local_family_id"
    )
    .agg(
        n_focals=(
            "focal_id",
            "nunique",
        ),

        n_GlobDB=(
            "source_dataset",
            lambda x:
                int(
                    (
                        x
                        ==
                        "GlobDB"
                    ).sum()
                ),
        ),

        n_metagenome=(
            "source_dataset",
            lambda x:
                int(
                    (
                        x
                        ==
                        "Metagenome"
                    ).sum()
                ),
        ),
    )
    .reset_index()
    .sort_values(
        [
            "n_focals",
            "local_family_id",
        ],
        ascending=[
            False,
            True,
        ],
    )
)


for row in focal_family_counts.itertuples(
    index=False
):

    print(
        f"  {row.local_family_id:<18} "
        f"{row.n_focals:>2} focals "
        f"(GlobDB={row.n_GlobDB}, "
        f"MAG={row.n_metagenome})"
    )


print()
print("=" * 80)

print(
    "STAGE 65 COMPLETE"
)

print("=" * 80)

print()
print(
    f"Neighborhood proteins:             "
    f"{EXPECTED_PROTEINS}"
)

print(
    f"Stable local families:             "
    f"{len(family_summary)}"
)

print(
    f"Singleton families:                "
    f"{n_singletons}"
)

print(
    f"Families in >=2 focal regions:     "
    f"{n_recurrent}"
)

print()
print(
    f"Shared GlobDB + MAG families:      "
    f"{n_shared}"
)

print(
    f"Shared recurrent families:         "
    f"{n_shared_recurrent}"
)

print()
print(
    f"Recurrent unknown families:        "
    f"{n_recurrent_unknown}"
)

print(
    f"Shared recurrent unknown families: "
    f"{n_shared_recurrent_unknown}"
)

print()
print(
    f"Local families among 23 focals:    "
    f"{n_focal_families}"
)

print()
print("Family summary:")

print(
    f"  {OUT_SUMMARY}"
)

print()
print("GlobDB vs MAG comparison:")

print(
    f"  {OUT_DATASET_COMPARISON}"
)

print()
print("Recurrent unknowns:")

print(
    f"  {OUT_UNKNOWN}"
)

print()
print("Gene-level map input:")

print(
    f"  {OUT_GENES}"
)

print()
print("Focal-family crosswalk:")

print(
    f"  {OUT_FOCAL_CROSSWALK}"
)

print()
print("Region architectures:")

print(
    f"  {OUT_ARCHITECTURES}"
)

print()
print("QC:")

print(
    f"  {OUT_QC}"
)
