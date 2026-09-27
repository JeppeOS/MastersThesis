#!/usr/bin/env python3

from __future__ import annotations

import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd


## ================================================================== ##
## STAGE 60
##
## INTEGRATE DIRECT COG20 ANNOTATION
##
## Purpose
## -------
##
## Integrate direct anvi'o COG20 annotation of the 845 priority-
## neighborhood proteins with:
##
##   - conservative GlobDB homology-transfer annotations;
##   - stable MG_local_### sequence families;
##   - actual-gene neighborhood observations;
##   - focal-group family summaries.
##
##
## Important:
##
## anvi-run-ncbi-cogs writes multiple SOURCE rows per protein:
##
##   COG20_FUNCTION
##   COG20_CATEGORY
##   COG20_PATHWAY
##
## Only COG20_FUNCTION rows represent direct COG assignments.
##
## CATEGORY and PATHWAY information are retained separately.
##
##
## No existing annotation is overwritten.
##
## Protein-level comparison classes:
##
##   direct_and_transfer_agree
##   direct_and_transfer_disagree
##   direct_only
##   transfer_only
##   neither
##
##
## Family-level direct COG summaries are descriptive.
## No family is forcibly reannotated here.
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


COG_DIR = (
    WORK_DIR
    / "direct_cog20"
)


LOCAL_DIR = (
    WORK_DIR
    / "local_clustering"
)


DIRECT_COG = (
    COG_DIR
    / "priority_neighborhood_proteins_COG20.tsv"
)


PROTEIN_METADATA = (
    WORK_DIR
    / "priority_neighborhood_protein_metadata_annotated.tsv"
)


LOCAL_MEMBERSHIP = (
    LOCAL_DIR
    / "MG_local_family_membership.tsv"
)


LOCAL_FAMILY_SUMMARY = (
    LOCAL_DIR
    / "MG_local_family_summary.tsv"
)


LOCAL_GROUP_SUMMARY = (
    LOCAL_DIR
    / "MG_local_family_by_focal_group.tsv"
)


GENE_TABLE = (
    WORK_DIR
    / "priority_observed_neighborhood_genes_local_families.tsv"
)


## ================================================================== ##
## Outputs
## ================================================================== ##

OUT_PROTEIN_COG = (
    COG_DIR
    / "direct_COG20_protein_annotations.tsv"
)


OUT_COMPARISON = (
    COG_DIR
    / "COG20_vs_GlobDB_transfer.tsv"
)


OUT_METADATA = (
    WORK_DIR
    / "priority_neighborhood_protein_metadata_COG20.tsv"
)


OUT_MEMBERSHIP = (
    LOCAL_DIR
    / "MG_local_family_membership_COG20.tsv"
)


OUT_FAMILY_SUMMARY = (
    LOCAL_DIR
    / "MG_local_family_summary_COG20.tsv"
)


OUT_GROUP_SUMMARY = (
    LOCAL_DIR
    / "MG_local_family_by_focal_group_COG20.tsv"
)


OUT_GENE_TABLE = (
    WORK_DIR
    / "priority_observed_neighborhood_genes_local_families_COG20.tsv"
)


OUT_QC = (
    COG_DIR
    / "COG20_integration_qc.tsv"
)


EXPECTED_PROTEINS = 845

EXPECTED_FAMILIES = 718

EXPECTED_GENE_ROWS = 845


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

    return str(
        value
    ).strip()


def ordered_unique(values):

    result = []

    seen = set()


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


    dominant, support = ordered[
        0
    ]


    total = len(
        values
    )


    return {
        "value":
            dominant,

        "support":
            support,

        "total":
            total,

        "fraction":
            support / total,

        "distribution":
            ";".join(
                f"{value}:{count}"

                for value, count
                in ordered
            ),
    }


## ================================================================== ##
## Start
## ================================================================== ##

print("=" * 80)

print(
    "STAGE 60 - INTEGRATE DIRECT COG20 ANNOTATIONS"
)

print("=" * 80)


for path in (
    DIRECT_COG,
    PROTEIN_METADATA,
    LOCAL_MEMBERSHIP,
    LOCAL_FAMILY_SUMMARY,
    LOCAL_GROUP_SUMMARY,
    GENE_TABLE,
):

    require_file(
        path
    )


## ================================================================== ##
## 1. Read direct anvi'o output
## ================================================================== ##

print()
print("Reading direct COG20 output...")


direct = pd.read_csv(
    DIRECT_COG,
    sep="\t",
    dtype=str,
    keep_default_na=False,
)


required_direct = {
    "accession",
    "e_value",
    "function",
    "gene_callers_id",
    "source",
}


missing = (
    required_direct
    -
    set(
        direct.columns
    )
)


if missing:

    fail(
        "Direct COG20 output is missing:\n"
        +
        "\n".join(
            sorted(
                missing
            )
        )
    )


direct[
    "e_value_numeric"
] = pd.to_numeric(
    direct[
        "e_value"
    ],
    errors="coerce",
)


print(
    f"  Annotation rows: "
    f"{len(direct):,}"
)

print()
print("  Rows by source:")


for source, count in (
    direct[
        "source"
    ]
    .value_counts()
    .sort_index()
    .items()
):

    print(
        f"    {source:<20} "
        f"{count:>6,}"
    )


## ================================================================== ##
## 2. Read authoritative 845-protein metadata
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
    "transfer_status",
    "transferred_cog",
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
        f"Expected {EXPECTED_PROTEINS} proteins but "
        f"found {len(metadata)}."
    )


if metadata[
    "neighborhood_protein_id"
].duplicated().any():

    fail(
        "Duplicate neighborhood_protein_id values."
    )


query_ids = set(
    metadata[
        "neighborhood_protein_id"
    ]
)


unexpected_ids = (
    set(
        direct[
            "gene_callers_id"
        ]
    )
    -
    query_ids
)


if unexpected_ids:

    fail(
        f"{len(unexpected_ids)} anvi'o IDs are not "
        f"present in the 845-protein metadata."
    )


## ================================================================== ##
## 3. Separate FUNCTION / CATEGORY / PATHWAY rows
## ================================================================== ##

function_rows = direct[
    direct[
        "source"
    ]
    ==
    "COG20_FUNCTION"
].copy()


category_rows = direct[
    direct[
        "source"
    ]
    ==
    "COG20_CATEGORY"
].copy()


pathway_rows = direct[
    direct[
        "source"
    ]
    ==
    "COG20_PATHWAY"
].copy()


print()
print(
    f"  COG20_FUNCTION rows: "
    f"{len(function_rows):,}"
)

print(
    f"  COG20_CATEGORY rows: "
    f"{len(category_rows):,}"
)

print(
    f"  COG20_PATHWAY rows:  "
    f"{len(pathway_rows):,}"
)


## ================================================================== ##
## 4. Aggregate direct COG function calls per protein
##
## We retain all function calls.
##
## The primary direct COG is the function row with the lowest
## e-value, with accession/function as deterministic tie-breakers.
## ================================================================== ##

protein_rows = []


for protein_id in metadata[
    "neighborhood_protein_id"
]:

    funcs = function_rows[
        function_rows[
            "gene_callers_id"
        ]
        ==
        protein_id
    ].copy()


    cats = category_rows[
        category_rows[
            "gene_callers_id"
        ]
        ==
        protein_id
    ].copy()


    paths = pathway_rows[
        pathway_rows[
            "gene_callers_id"
        ]
        ==
        protein_id
    ].copy()


    ## -------------------------------------------------------------- ##
    ## Direct function calls
    ## -------------------------------------------------------------- ##

    if len(
        funcs
    ) > 0:

        funcs = (
            funcs
            .sort_values(
                [
                    "e_value_numeric",
                    "accession",
                    "function",
                ],
                ascending=[
                    True,
                    True,
                    True,
                ],
                kind="stable",
                na_position="last",
            )
            .reset_index(
                drop=True
            )
        )


        primary = funcs.iloc[
            0
        ]


        primary_cog = clean(
            primary[
                "accession"
            ]
        )


        primary_function = clean(
            primary[
                "function"
            ]
        )


        primary_evalue = clean(
            primary[
                "e_value"
            ]
        )


        direct_cogs = join_unique(
            funcs[
                "accession"
            ]
        )


        direct_functions = join_unique(
            funcs[
                "function"
            ],
            separator=" || ",
        )


        n_function_calls = (
            funcs[
                [
                    "accession",
                    "function",
                ]
            ]
            .drop_duplicates()
            .shape[0]
        )


    else:

        primary_cog = ""

        primary_function = ""

        primary_evalue = ""

        direct_cogs = ""

        direct_functions = ""

        n_function_calls = 0


    ## -------------------------------------------------------------- ##
    ## Category rows
    ##
    ## Example accession:
    ##
    ##     H!!!I
    ##
    ## Keep exactly as anvi'o reports it.
    ## -------------------------------------------------------------- ##

    category_accessions = join_unique(
        cats[
            "accession"
        ]
    )


    category_functions = join_unique(
        cats[
            "function"
        ],
        separator=" || ",
    )


    ## -------------------------------------------------------------- ##
    ## Pathway rows
    ## -------------------------------------------------------------- ##

    pathway_accessions = join_unique(
        paths[
            "accession"
        ]
    )


    pathway_functions = join_unique(
        paths[
            "function"
        ],
        separator=" || ",
    )


    protein_rows.append(
        {
            "neighborhood_protein_id":
                protein_id,

            "direct_COG20_positive":
                int(
                    n_function_calls > 0
                ),

            "direct_COG20_primary_cog":
                primary_cog,

            "direct_COG20_primary_function":
                primary_function,

            "direct_COG20_primary_evalue":
                primary_evalue,

            "direct_COG20_n_function_calls":
                n_function_calls,

            "direct_COG20_all_cogs":
                direct_cogs,

            "direct_COG20_all_functions":
                direct_functions,

            "direct_COG20_category_codes":
                category_accessions,

            "direct_COG20_category_functions":
                category_functions,

            "direct_COG20_pathway_accessions":
                pathway_accessions,

            "direct_COG20_pathways":
                pathway_functions,
        }
    )


protein_cog = pd.DataFrame(
    protein_rows
)


if len(
    protein_cog
) != EXPECTED_PROTEINS:

    fail(
        "Protein-level direct COG table has "
        "incorrect row count."
    )


protein_cog.to_csv(
    OUT_PROTEIN_COG,
    sep="\t",
    index=False,
    na_rep="",
)


## ================================================================== ##
## 5. Compare direct COG20 with transferred GlobDB COG
## ================================================================== ##

print()
print("Comparing direct and transferred COG annotations...")


comparison = metadata[
    [
        "neighborhood_protein_id",
        "genome",
        "protein_id",

        "transfer_status",
        "transferred_cog",
        "transferred_product",
    ]
].merge(
    protein_cog,
    on="neighborhood_protein_id",
    how="left",
    validate="one_to_one",
)


def compare_annotation(row):

    transferred = clean(
        row[
            "transferred_cog"
        ]
    )


    direct_all = {
        value.strip()

        for value
        in clean(
            row[
                "direct_COG20_all_cogs"
            ]
        ).split(";")

        if value.strip()
    }


    direct_present = bool(
        direct_all
    )


    transfer_present = bool(
        transferred
    )


    if (
        direct_present
        and
        transfer_present
    ):

        if transferred in direct_all:

            return (
                "direct_and_transfer_agree"
            )

        return (
            "direct_and_transfer_disagree"
        )


    if direct_present:

        return (
            "direct_only"
        )


    if transfer_present:

        return (
            "transfer_only"
        )


    return (
        "neither"
    )


comparison[
    "COG_annotation_comparison"
] = comparison.apply(
    compare_annotation,
    axis=1,
)


comparison[
    "primary_direct_equals_transfer"
] = (
    (
        comparison[
            "direct_COG20_primary_cog"
        ]
        ==
        comparison[
            "transferred_cog"
        ]
    )
    &
    (
        comparison[
            "transferred_cog"
        ]
        .str.strip()
        !=
        ""
    )
).astype(int)


comparison.to_csv(
    OUT_COMPARISON,
    sep="\t",
    index=False,
    na_rep="",
)


comparison_counts = (
    comparison[
        "COG_annotation_comparison"
    ]
    .value_counts()
)


## ================================================================== ##
## 6. Merge direct COG onto protein metadata
## ================================================================== ##

direct_for_merge = protein_cog.copy()


annotated_metadata = metadata.merge(
    direct_for_merge,
    on="neighborhood_protein_id",
    how="left",
    validate="one_to_one",
)


annotated_metadata = annotated_metadata.merge(
    comparison[
        [
            "neighborhood_protein_id",
            "COG_annotation_comparison",
            "primary_direct_equals_transfer",
        ]
    ],
    on="neighborhood_protein_id",
    how="left",
    validate="one_to_one",
)


annotated_metadata.to_csv(
    OUT_METADATA,
    sep="\t",
    index=False,
    na_rep="",
)


## ================================================================== ##
## 7. Merge direct COG onto stable local-family membership
## ================================================================== ##

membership = pd.read_csv(
    LOCAL_MEMBERSHIP,
    sep="\t",
    dtype=str,
    keep_default_na=False,
)


if len(
    membership
) != EXPECTED_PROTEINS:

    fail(
        f"Expected {EXPECTED_PROTEINS} local-family "
        f"membership rows but found {len(membership)}."
    )


required_membership = {
    "neighborhood_protein_id",
    "local_family_id",
}


missing = (
    required_membership
    -
    set(
        membership.columns
    )
)


if missing:

    fail(
        "Local-family membership is missing:\n"
        +
        "\n".join(
            sorted(
                missing
            )
        )
    )


membership = membership.merge(
    protein_cog,
    on="neighborhood_protein_id",
    how="left",
    validate="one_to_one",
)


membership = membership.merge(
    comparison[
        [
            "neighborhood_protein_id",
            "COG_annotation_comparison",
            "primary_direct_equals_transfer",
        ]
    ],
    on="neighborhood_protein_id",
    how="left",
    validate="one_to_one",
)


membership.to_csv(
    OUT_MEMBERSHIP,
    sep="\t",
    index=False,
    na_rep="",
)


## ================================================================== ##
## 8. Build local-family direct COG20 summaries
##
## Count direct COG occurrence by actual protein.
##
## All COG20_FUNCTION accessions are considered, not only the
## primary hit.
## ================================================================== ##

print()
print("Summarizing direct COG20 evidence by MG_local family...")


family_summary = pd.read_csv(
    LOCAL_FAMILY_SUMMARY,
    sep="\t",
    dtype=str,
    keep_default_na=False,
)


if len(
    family_summary
) != EXPECTED_FAMILIES:

    fail(
        f"Expected {EXPECTED_FAMILIES} local families "
        f"but found {len(family_summary)}."
    )


## -------------------------------------------------------------- ##
## Build protein x COG membership.
## -------------------------------------------------------------- ##

protein_to_family = membership[
    [
        "neighborhood_protein_id",
        "local_family_id",
    ]
]


direct_cog_members = function_rows[
    [
        "gene_callers_id",
        "accession",
        "function",
        "e_value_numeric",
    ]
].rename(
    columns={
        "gene_callers_id":
            "neighborhood_protein_id",

        "accession":
            "direct_cog",

        "function":
            "direct_cog_function",
    }
)


direct_cog_members = (
    direct_cog_members
    .drop_duplicates(
        subset=[
            "neighborhood_protein_id",
            "direct_cog",
        ]
    )
    .merge(
        protein_to_family,
        on="neighborhood_protein_id",
        how="left",
        validate="many_to_one",
    )
)


if direct_cog_members[
    "local_family_id"
].isna().any():

    fail(
        "A direct COG20 annotation could not be "
        "mapped to an MG_local family."
    )


## -------------------------------------------------------------- ##
## Stable COG -> function dictionary from direct annotations.
## -------------------------------------------------------------- ##

cog_function_rows = []


for cog, group in direct_cog_members.groupby(
    "direct_cog",
    sort=False,
):

    summary = dominant_nonempty(
        group[
            "direct_cog_function"
        ]
    )


    cog_function_rows.append(
        {
            "direct_cog":
                cog,

            "direct_cog_function":
                summary[
                    "value"
                ],
        }
    )


cog_function_map = dict(
    zip(
        pd.DataFrame(
            cog_function_rows
        )[
            "direct_cog"
        ],
        pd.DataFrame(
            cog_function_rows
        )[
            "direct_cog_function"
        ],
    )
) if cog_function_rows else {}


family_direct_rows = []


for family_id in family_summary[
    "local_family_id"
]:

    members = membership[
        membership[
            "local_family_id"
        ]
        ==
        family_id
    ]


    n_members = len(
        members
    )


    annotated_members = members[
        members[
            "direct_COG20_positive"
        ].astype(str)
        ==
        "1"
    ]


    n_direct = len(
        annotated_members
    )


    cog_hits = direct_cog_members[
        direct_cog_members[
            "local_family_id"
        ]
        ==
        family_id
    ]


    if len(
        cog_hits
    ) > 0:

        counts = (
            cog_hits[
                "direct_cog"
            ]
            .value_counts()
        )


        max_support = int(
            counts.iloc[
                0
            ]
        )


        tied = sorted(
            counts[
                counts
                ==
                max_support
            ]
            .index
            .tolist()
        )


        dominant_cog = tied[
            0
        ]


        dominant_function = (
            cog_function_map.get(
                dominant_cog,
                ""
            )
        )


        support_all = (
            max_support
            /
            n_members
        )


        support_annotated = (
            max_support
            /
            n_direct
            if n_direct > 0
            else np.nan
        )


        distribution = ";".join(
            f"{cog}:{count}"

            for cog, count
            in sorted(
                counts.items(),
                key=lambda x: (
                    -x[1],
                    x[0],
                ),
            )
        )


        all_cogs = ";".join(
            sorted(
                counts.index
            )
        )


    else:

        dominant_cog = ""

        dominant_function = ""

        max_support = 0

        support_all = 0.0

        support_annotated = np.nan

        distribution = ""

        all_cogs = ""


    ## -------------------------------------------------------------- ##
    ## Descriptive family status.
    ##
    ## This does NOT overwrite the existing family annotation.
    ## -------------------------------------------------------------- ##

    if n_direct == 0:

        status = (
            "no_direct_COG20"
        )


    elif n_members == 1:

        status = (
            "singleton_direct_COG20"
        )


    elif (
        max_support
        ==
        n_members
    ):

        status = (
            "unanimous_all_members"
        )


    elif (
        support_all
        >=
        0.50
        and
        max_support
        >=
        2
    ):

        status = (
            "majority_family_support"
        )


    elif (
        support_annotated
        == 1.0
    ):

        status = (
            "unanimous_annotated_subset"
        )


    else:

        status = (
            "mixed_or_sparse"
        )


    family_direct_rows.append(
        {
            "local_family_id":
                family_id,

            "direct_COG20_annotated_members":
                n_direct,

            "direct_COG20_member_coverage_fraction":
                n_direct
                /
                n_members,

            "direct_COG20_dominant_cog":
                dominant_cog,

            "direct_COG20_dominant_function":
                dominant_function,

            "direct_COG20_dominant_support":
                max_support,

            "direct_COG20_dominant_fraction_all_members":
                support_all,

            "direct_COG20_dominant_fraction_annotated_members":
                support_annotated,

            "direct_COG20_all_cogs":
                all_cogs,

            "direct_COG20_distribution":
                distribution,

            "direct_COG20_family_status":
                status,
        }
    )


family_direct = pd.DataFrame(
    family_direct_rows
)


family_summary_cog = family_summary.merge(
    family_direct,
    on="local_family_id",
    how="left",
    validate="one_to_one",
)


## ================================================================== ##
## 9. Flag previously unannotated families rescued by direct COG20
## ================================================================== ##

family_summary_cog[
    "preexisting_family_annotation"
] = (
    family_summary_cog[
        "family_annotation_label"
    ]
    .str.strip()
    !=
    ""
).astype(int)


family_summary_cog[
    "direct_COG20_rescues_unannotated_family"
] = (
    (
        family_summary_cog[
            "preexisting_family_annotation"
        ]
        ==
        0
    )
    &
    (
        family_summary_cog[
            "direct_COG20_dominant_cog"
        ]
        .str.strip()
        !=
        ""
    )
).astype(int)


family_summary_cog[
    "recurrent_family"
] = (
    pd.to_numeric(
        family_summary_cog[
            "n_focal_regions"
        ],
        errors="coerce",
    )
    >=
    2
).astype(int)


family_summary_cog[
    "direct_COG20_rescues_recurrent_unannotated_family"
] = (
    (
        family_summary_cog[
            "direct_COG20_rescues_unannotated_family"
        ]
        ==
        1
    )
    &
    (
        family_summary_cog[
            "recurrent_family"
        ]
        ==
        1
    )
).astype(int)


family_summary_cog.to_csv(
    OUT_FAMILY_SUMMARY,
    sep="\t",
    index=False,
    na_rep="",
    float_format="%.6f",
)


## ================================================================== ##
## 10. Attach family-level COG evidence to focal-group summary
## ================================================================== ##

group_summary = pd.read_csv(
    LOCAL_GROUP_SUMMARY,
    sep="\t",
    dtype=str,
    keep_default_na=False,
)


family_direct_columns = [
    "local_family_id",

    "direct_COG20_annotated_members",
    "direct_COG20_member_coverage_fraction",

    "direct_COG20_dominant_cog",
    "direct_COG20_dominant_function",

    "direct_COG20_dominant_support",
    "direct_COG20_dominant_fraction_all_members",
    "direct_COG20_dominant_fraction_annotated_members",

    "direct_COG20_all_cogs",
    "direct_COG20_distribution",

    "direct_COG20_family_status",

    "direct_COG20_rescues_unannotated_family",
    "direct_COG20_rescues_recurrent_unannotated_family",
]


group_summary = group_summary.merge(
    family_summary_cog[
        family_direct_columns
    ],
    on="local_family_id",
    how="left",
    validate="many_to_one",
)


group_summary.to_csv(
    OUT_GROUP_SUMMARY,
    sep="\t",
    index=False,
    na_rep="",
)


## ================================================================== ##
## 11. Attach direct protein-level COG20 to actual genes
##
## This becomes the next gene-level visualization input.
## ================================================================== ##

genes = pd.read_csv(
    GENE_TABLE,
    sep="\t",
    dtype=str,
    keep_default_na=False,
)


if len(
    genes
) != EXPECTED_GENE_ROWS:

    fail(
        f"Expected {EXPECTED_GENE_ROWS} gene rows "
        f"but found {len(genes)}."
    )


if "neighborhood_protein_id" not in genes.columns:

    fail(
        "Gene-level table lacks neighborhood_protein_id."
    )


genes_cog = genes.merge(
    protein_cog,
    on="neighborhood_protein_id",
    how="left",
    validate="one_to_one",
)


genes_cog = genes_cog.merge(
    comparison[
        [
            "neighborhood_protein_id",
            "COG_annotation_comparison",
            "primary_direct_equals_transfer",
        ]
    ],
    on="neighborhood_protein_id",
    how="left",
    validate="one_to_one",
)


genes_cog.to_csv(
    OUT_GENE_TABLE,
    sep="\t",
    index=False,
    na_rep="",
)


## ================================================================== ##
## 12. QC / headline results
## ================================================================== ##

n_direct_positive = int(
    protein_cog[
        "direct_COG20_positive"
    ].sum()
)


n_no_direct = (
    EXPECTED_PROTEINS
    -
    n_direct_positive
)


n_agree = int(
    comparison_counts.get(
        "direct_and_transfer_agree",
        0,
    )
)


n_disagree = int(
    comparison_counts.get(
        "direct_and_transfer_disagree",
        0,
    )
)


n_direct_only = int(
    comparison_counts.get(
        "direct_only",
        0,
    )
)


n_transfer_only = int(
    comparison_counts.get(
        "transfer_only",
        0,
    )
)


n_neither = int(
    comparison_counts.get(
        "neither",
        0,
    )
)


n_rescued_families = int(
    family_summary_cog[
        "direct_COG20_rescues_unannotated_family"
    ].sum()
)


n_rescued_recurrent = int(
    family_summary_cog[
        "direct_COG20_rescues_recurrent_unannotated_family"
    ].sum()
)


n_recurrent_unannotated_before = int(
    (
        (
            family_summary_cog[
                "preexisting_family_annotation"
            ]
            ==
            0
        )
        &
        (
            family_summary_cog[
                "recurrent_family"
            ]
            ==
            1
        )
    ).sum()
)


n_recurrent_still_without_cog = int(
    (
        (
            family_summary_cog[
                "preexisting_family_annotation"
            ]
            ==
            0
        )
        &
        (
            family_summary_cog[
                "recurrent_family"
            ]
            ==
            1
        )
        &
        (
            family_summary_cog[
                "direct_COG20_dominant_cog"
            ]
            .str.strip()
            ==
            ""
        )
    ).sum()
)


qc = pd.DataFrame(
    [
        (
            "proteins",
            EXPECTED_PROTEINS,
        ),

        (
            "direct_COG20_positive_proteins",
            n_direct_positive,
        ),

        (
            "direct_COG20_negative_proteins",
            n_no_direct,
        ),

        (
            "direct_and_transfer_agree",
            n_agree,
        ),

        (
            "direct_and_transfer_disagree",
            n_disagree,
        ),

        (
            "direct_only",
            n_direct_only,
        ),

        (
            "transfer_only",
            n_transfer_only,
        ),

        (
            "neither",
            n_neither,
        ),

        (
            "local_families",
            len(
                family_summary_cog
            ),
        ),

        (
            "families_rescued_by_direct_COG20",
            n_rescued_families,
        ),

        (
            "recurrent_unannotated_families_before_COG20",
            n_recurrent_unannotated_before,
        ),

        (
            "recurrent_unannotated_families_rescued_by_COG20",
            n_rescued_recurrent,
        ),

        (
            "recurrent_unannotated_families_still_without_COG20",
            n_recurrent_still_without_cog,
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
## Final report
## ================================================================== ##

print()
print("=" * 80)
print("STAGE 60 COMPLETE")
print("=" * 80)

print()
print("Protein-level direct COG20")

print(
    f"  Direct COG20 positive:           "
    f"{n_direct_positive:>3}/{EXPECTED_PROTEINS}"
)

print(
    f"  No direct COG20 call:            "
    f"{n_no_direct:>3}/{EXPECTED_PROTEINS}"
)


print()
print("Direct COG20 vs GlobDB transfer")

print(
    f"  Agree:                           "
    f"{n_agree:>3}"
)

print(
    f"  Disagree:                        "
    f"{n_disagree:>3}"
)

print(
    f"  Direct only:                     "
    f"{n_direct_only:>3}"
)

print(
    f"  Transfer only:                   "
    f"{n_transfer_only:>3}"
)

print(
    f"  Neither:                         "
    f"{n_neither:>3}"
)


print()
print("Local-family annotation rescue")

print(
    f"  Families rescued by direct COG: "
    f"{n_rescued_families:>3}"
)

print(
    f"  Recurrent unannotated before:    "
    f"{n_recurrent_unannotated_before:>3}"
)

print(
    f"  Recurrent rescued by COG20:      "
    f"{n_rescued_recurrent:>3}"
)

print(
    f"  Recurrent still without COG20:   "
    f"{n_recurrent_still_without_cog:>3}"
)


print()
print("Outputs:")

print(
    f"  Protein COG20 annotations:\n"
    f"    {OUT_PROTEIN_COG}"
)

print(
    f"  Direct-vs-transfer comparison:\n"
    f"    {OUT_COMPARISON}"
)

print(
    f"  Updated protein metadata:\n"
    f"    {OUT_METADATA}"
)

print(
    f"  Updated family membership:\n"
    f"    {OUT_MEMBERSHIP}"
)

print(
    f"  Updated family summary:\n"
    f"    {OUT_FAMILY_SUMMARY}"
)

print(
    f"  Updated focal-group summary:\n"
    f"    {OUT_GROUP_SUMMARY}"
)

print(
    f"  Updated actual-gene table:\n"
    f"    {OUT_GENE_TABLE}"
)

print(
    f"  QC:\n"
    f"    {OUT_QC}"
)
