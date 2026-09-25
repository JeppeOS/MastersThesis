#!/usr/bin/env python3

from __future__ import annotations

import sys
from pathlib import Path
from collections import Counter

import pandas as pd


## ================================================================== ##
## STAGE 56
##
## CONSERVATIVE GLOBDB ANNOTATION TRANSFER
##
## Query:
##     845 proteins from the 21 priority metagenome neighborhoods
##
## Reference:
##     2,002,656 proteins from the original 631-genome GlobDB
##     Methylococcales workflow
##
## Existing GlobDB annotations:
##
##     globdb_cog
##     globdb_gene
##     globdb_product
##
## Transfer policy
## ---------------
##
## HIGH CONFIDENCE
##
##     pident >= 70%
##     qcov   >= 0.80
##     tcov   >= 0.80
##
##     -> transfer consensus COG where supported
##     -> transfer gene/product from best qualifying annotated hit
##
##
## FAMILY LEVEL
##
##     pident >= 50%
##     qcov   >= 0.80
##     tcov   >= 0.80
##
##     -> transfer consensus COG only
##     -> do NOT transfer gene/product
##
##
## NO TRANSFER
##
##     no suitable accepted annotated reference homolog
##
##
## Important:
## ---------
##
## Raw best-hit information is retained even when no annotation is
## transferred.
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


TRANSFER_DIR = (
    WORK_DIR
    / "globdb_annotation_transfer"
)


QUERY_METADATA = (
    WORK_DIR
    / "priority_neighborhood_protein_metadata.tsv"
)


HITS_FILE = (
    TRANSFER_DIR
    / "neighborhood_vs_GlobDB631.tsv"
)


REFERENCE_METADATA = (
    TRANSFER_DIR
    / "GlobDB631_annotation_metadata.tsv"
)


NEIGHBORHOODS = (
    WORK_DIR
    / "priority_observed_neighborhood_genes.tsv"
)


OUT_TRANSFER = (
    TRANSFER_DIR
    / "neighborhood_GlobDB_annotation_transfer.tsv"
)


OUT_CONSIDERED_HITS = (
    TRANSFER_DIR
    / "neighborhood_annotation_candidate_hits.tsv"
)


OUT_METADATA_ANNOTATED = (
    WORK_DIR
    / "priority_neighborhood_protein_metadata_annotated.tsv"
)


OUT_NEIGHBORHOODS_ANNOTATED = (
    WORK_DIR
    / "priority_observed_neighborhood_genes_annotated.tsv"
)


OUT_QC = (
    TRANSFER_DIR
    / "annotation_transfer_qc.tsv"
)


## ================================================================== ##
## Locked query expectation
## ================================================================== ##

EXPECTED_QUERIES = 845


## ================================================================== ##
## Thresholds
## ================================================================== ##

HIGH_IDENTITY = 70.0

FAMILY_IDENTITY = 50.0

MIN_QCOV = 0.80

MIN_TCOV = 0.80


## Minimum fraction of qualifying COG-bearing homologs supporting
## the dominant COG.
##
## High-confidence transfer can tolerate modest diversity.
##
## Family-level transfer is deliberately stricter.

HIGH_COG_CONSENSUS = 0.60

FAMILY_COG_CONSENSUS = 0.70


## ================================================================== ##
## Utilities
## ================================================================== ##

def fail(message):

    print(
        f"\nERROR: {message}",
        file=sys.stderr,
    )

    sys.exit(1)


def clean(value):

    if pd.isna(value):

        return ""

    return str(value).strip()


def accepted_flag(value):

    value = clean(
        value
    ).lower()


    return value in {
        "1",
        "true",
        "t",
        "yes",
        "y",
    }


def has_annotation(row):

    return bool(
        clean(
            row.get(
                "globdb_cog",
                ""
            )
        )
        or
        clean(
            row.get(
                "globdb_gene",
                ""
            )
        )
        or
        clean(
            row.get(
                "globdb_product",
                ""
            )
        )
    )


def dominant_value(
    values,
):

    values = [
        clean(value)

        for value
        in values

        if clean(
            value
        )
    ]


    if not values:

        return (
            "",
            0,
            0,
            0.0,
            "",
        )


    counts = Counter(
        values
    )


    ordered = sorted(
        counts.items(),
        key=lambda item: (
            -item[1],
            item[0],
        )
    )


    dominant, support = (
        ordered[0]
    )


    total = len(
        values
    )


    fraction = (
        support
        /
        total
    )


    distribution = ";".join(
        f"{value}:{count}"

        for value, count
        in ordered
    )


    return (
        dominant,
        support,
        total,
        fraction,
        distribution,
    )


## ================================================================== ##
## Start
## ================================================================== ##

print("=" * 80)
print("STAGE 56 - CONSERVATIVE GLOBDB ANNOTATION TRANSFER")
print("=" * 80)


for path in (
    QUERY_METADATA,
    HITS_FILE,
    REFERENCE_METADATA,
    NEIGHBORHOODS,
):

    if not path.is_file():

        fail(
            f"Missing required file:\n{path}"
        )


## ================================================================== ##
## 1. Read 845-query metadata
## ================================================================== ##

print()
print("Reading metagenome neighborhood-protein metadata...")


queries = pd.read_csv(
    QUERY_METADATA,
    sep="\t",
    dtype=str,
    keep_default_na=False,
)


required_queries = {
    "neighborhood_protein_id",
    "genome",
    "protein_id",
}


missing = (
    required_queries
    -
    set(
        queries.columns
    )
)


if missing:

    fail(
        "Query metadata is missing:\n"
        +
        "\n".join(
            sorted(
                missing
            )
        )
    )


if len(
    queries
) != EXPECTED_QUERIES:

    fail(
        f"Expected {EXPECTED_QUERIES} neighborhood "
        f"proteins but found {len(queries)}."
    )


if queries[
    "neighborhood_protein_id"
].duplicated().any():

    fail(
        "Duplicate neighborhood_protein_id values."
    )


query_ids = set(
    queries[
        "neighborhood_protein_id"
    ]
)


print(
    f"  Queries: {len(queries):,}"
)


## ================================================================== ##
## 2. Read MMseqs hits
## ================================================================== ##

print()
print("Reading MMseqs GlobDB search results...")


hit_columns = [
    "query",
    "target",
    "pident",
    "alnlen",
    "qcov",
    "tcov",
    "evalue",
    "bits",
]


hits = pd.read_csv(
    HITS_FILE,
    sep="\t",
    names=hit_columns,
    header=None,
    dtype={
        "query":
            str,

        "target":
            str,
    },
)


for column in [
    "pident",
    "alnlen",
    "qcov",
    "tcov",
    "evalue",
    "bits",
]:

    hits[
        column
    ] = pd.to_numeric(
        hits[
            column
        ],
        errors="raise",
    )


unexpected_queries = (
    set(
        hits[
            "query"
        ]
    )
    -
    query_ids
)


if unexpected_queries:

    fail(
        f"{len(unexpected_queries)} MMseqs query IDs "
        f"are absent from query metadata."
    )


target_ids = set(
    hits[
        "target"
    ]
)


print(
    f"  Search rows:          "
    f"{len(hits):,}"
)

print(
    f"  Queries with >=1 hit: "
    f"{hits['query'].nunique():,}"
)

print(
    f"  Unique target hits:   "
    f"{len(target_ids):,}"
)


## ================================================================== ##
## 3. Load ONLY reference metadata needed for targets actually hit
##
## Avoid loading all 2 million annotation strings into memory.
## ================================================================== ##

print()
print("Recovering GlobDB metadata for hit proteins...")


reference_header = pd.read_csv(
    REFERENCE_METADATA,
    sep="\t",
    nrows=0,
)


available_reference_columns = list(
    reference_header.columns
)


reference_columns = [
    "protein_id",
    "globdb_cog",
    "globdb_gene",
    "globdb_product",
]


for optional in [
    "annotation_match_type",
    "annotation_accepted",
    "reciprocal_overlap",
    "has_globdb_annotation",
]:

    if optional in available_reference_columns:

        reference_columns.append(
            optional
        )


reference_parts = []


for chunk in pd.read_csv(
    REFERENCE_METADATA,
    sep="\t",
    dtype=str,
    keep_default_na=False,
    usecols=reference_columns,
    chunksize=200_000,
):

    keep = chunk[
        "protein_id"
    ].isin(
        target_ids
    )


    if keep.any():

        reference_parts.append(
            chunk.loc[
                keep
            ].copy()
        )


if not reference_parts:

    fail(
        "No hit proteins could be recovered from "
        "reference metadata."
    )


reference = pd.concat(
    reference_parts,
    ignore_index=True,
)


if reference[
    "protein_id"
].duplicated().any():

    fail(
        "Duplicate protein_id values in recovered "
        "GlobDB reference metadata."
    )


missing_target_metadata = (
    target_ids
    -
    set(
        reference[
            "protein_id"
        ]
    )
)


if missing_target_metadata:

    fail(
        f"{len(missing_target_metadata)} MMseqs target IDs "
        f"lack GlobDB metadata."
    )


reference = reference.rename(
    columns={
        "protein_id":
            "target"
    }
)


print(
    f"  Hit proteins with metadata: "
    f"{len(reference):,}"
)


## ================================================================== ##
## 4. Attach reference annotations to every hit
## ================================================================== ##

hits = hits.merge(
    reference,
    on="target",
    how="left",
    validate="many_to_one",
)


hits[
    "reference_has_annotation"
] = hits.apply(
    has_annotation,
    axis=1,
).astype(int)


## -------------------------------------------------------------- ##
## Respect old Stage-11 annotation acceptance.
##
## If annotation_accepted exists and contains meaningful accepted
## entries, use it.
##
## Otherwise fall back to presence/absence of GlobDB annotation.
## -------------------------------------------------------------- ##

if "annotation_accepted" in hits.columns:

    hits[
        "reference_annotation_accepted"
    ] = hits[
        "annotation_accepted"
    ].apply(
        lambda value:
            int(
                accepted_flag(
                    value
                )
            )
    )


    n_accepted = int(
        hits[
            "reference_annotation_accepted"
        ].sum()
    )


    if n_accepted == 0:

        print(
            "  WARNING: annotation_accepted contains no "
            "accepted hit rows; falling back to annotation presence."
        )


        hits[
            "reference_annotation_accepted"
        ] = hits[
            "reference_has_annotation"
        ]

else:

    hits[
        "reference_annotation_accepted"
    ] = hits[
        "reference_has_annotation"
    ]


hits[
    "usable_annotation_hit"
] = (
    (
        hits[
            "reference_has_annotation"
        ]
        ==
        1
    )
    &
    (
        hits[
            "reference_annotation_accepted"
        ]
        ==
        1
    )
).astype(int)


## ================================================================== ##
## 5. Mark similarity tiers
## ================================================================== ##

hits[
    "passes_bidirectional_80pct_coverage"
] = (
    (
        hits[
            "qcov"
        ]
        >=
        MIN_QCOV
    )
    &
    (
        hits[
            "tcov"
        ]
        >=
        MIN_TCOV
    )
).astype(int)


hits[
    "passes_high_confidence_threshold"
] = (
    (
        hits[
            "usable_annotation_hit"
        ]
        ==
        1
    )
    &
    (
        hits[
            "pident"
        ]
        >=
        HIGH_IDENTITY
    )
    &
    (
        hits[
            "passes_bidirectional_80pct_coverage"
        ]
        ==
        1
    )
).astype(int)


hits[
    "passes_family_threshold"
] = (
    (
        hits[
            "usable_annotation_hit"
        ]
        ==
        1
    )
    &
    (
        hits[
            "pident"
        ]
        >=
        FAMILY_IDENTITY
    )
    &
    (
        hits[
            "passes_bidirectional_80pct_coverage"
        ]
        ==
        1
    )
).astype(int)


## ================================================================== ##
## 6. Save the potentially transfer-relevant hit evidence
## ================================================================== ##

considered_hits = hits[
    (
        hits[
            "passes_family_threshold"
        ]
        ==
        1
    )
].copy()


considered_hits = (
    considered_hits
    .sort_values(
        [
            "query",
            "bits",
            "pident",
            "qcov",
            "tcov",
            "target",
        ],
        ascending=[
            True,
            False,
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


considered_hits.to_csv(
    OUT_CONSIDERED_HITS,
    sep="\t",
    index=False,
    na_rep="",
)


## ================================================================== ##
## 7. Build one conservative transfer row per query
## ================================================================== ##

print()
print("Building conservative annotation transfers...")


transfer_rows = []


hit_groups = {
    query:
        group.copy()

    for query, group
    in hits.groupby(
        "query",
        sort=False,
    )
}


for query_row in queries.itertuples(
    index=False
):

    query = query_row.neighborhood_protein_id


    group = hit_groups.get(
        query,
        pd.DataFrame(
            columns=hits.columns
        ),
    )


    ## -------------------------------------------------------------- ##
    ## Best raw hit regardless of whether annotation exists
    ## -------------------------------------------------------------- ##

    if len(
        group
    ) > 0:

        best_overall = (
            group
            .sort_values(
                [
                    "bits",
                    "pident",
                    "qcov",
                    "tcov",
                    "evalue",
                    "target",
                ],
                ascending=[
                    False,
                    False,
                    False,
                    False,
                    True,
                    True,
                ],
                kind="stable",
            )
            .iloc[0]
        )

    else:

        best_overall = None


    ## -------------------------------------------------------------- ##
    ## Best accepted annotated hit
    ## -------------------------------------------------------------- ##

    annotated = group[
        group[
            "usable_annotation_hit"
        ]
        ==
        1
    ].copy()


    if len(
        annotated
    ) > 0:

        best_annotated = (
            annotated
            .sort_values(
                [
                    "bits",
                    "pident",
                    "qcov",
                    "tcov",
                    "evalue",
                    "target",
                ],
                ascending=[
                    False,
                    False,
                    False,
                    False,
                    True,
                    True,
                ],
                kind="stable",
            )
            .iloc[0]
        )

    else:

        best_annotated = None


    ## -------------------------------------------------------------- ##
    ## Candidate tier
    ## -------------------------------------------------------------- ##

    high = group[
        group[
            "passes_high_confidence_threshold"
        ]
        ==
        1
    ].copy()


    family = group[
        group[
            "passes_family_threshold"
        ]
        ==
        1
    ].copy()


    if len(
        high
    ) > 0:

        status = (
            "high_confidence"
        )

        pool = high

        cog_consensus_required = (
            HIGH_COG_CONSENSUS
        )


    elif len(
        family
    ) > 0:

        status = (
            "family_level"
        )

        pool = family

        cog_consensus_required = (
            FAMILY_COG_CONSENSUS
        )


    else:

        status = (
            "no_transfer"
        )

        pool = group.iloc[
            0:0
        ].copy()

        cog_consensus_required = None


    ## -------------------------------------------------------------- ##
    ## Best hit inside selected transfer tier
    ## -------------------------------------------------------------- ##

    if len(
        pool
    ) > 0:

        transfer_best = (
            pool
            .sort_values(
                [
                    "bits",
                    "pident",
                    "qcov",
                    "tcov",
                    "evalue",
                    "target",
                ],
                ascending=[
                    False,
                    False,
                    False,
                    False,
                    True,
                    True,
                ],
                kind="stable",
            )
            .iloc[0]
        )

    else:

        transfer_best = None


    ## -------------------------------------------------------------- ##
    ## COG consensus
    ## -------------------------------------------------------------- ##

    (
        dominant_cog,
        cog_support,
        n_cog_hits,
        cog_fraction,
        cog_distribution,
    ) = dominant_value(
        pool[
            "globdb_cog"
        ].tolist()
        if (
            len(pool) > 0
            and
            "globdb_cog"
            in pool.columns
        )
        else
        []
    )


    transferred_cog = ""

    cog_basis = ""


    if (
        dominant_cog
        and
        cog_consensus_required is not None
        and
        cog_fraction
        >=
        cog_consensus_required
    ):

        transferred_cog = (
            dominant_cog
        )

        cog_basis = (
            "consensus"
        )


    ## Extremely close top hit can rescue COG when the qualifying
    ## homolog set contains heterogeneous annotation strings.

    elif (
        status == "high_confidence"
        and
        transfer_best is not None
        and
        clean(
            transfer_best.get(
                "globdb_cog",
                ""
            )
        )
        and
        float(
            transfer_best[
                "pident"
            ]
        )
        >= 90.0
        and
        float(
            transfer_best[
                "qcov"
            ]
        )
        >= 0.90
        and
        float(
            transfer_best[
                "tcov"
            ]
        )
        >= 0.90
    ):

        transferred_cog = clean(
            transfer_best[
                "globdb_cog"
            ]
        )

        cog_basis = (
            "very_high_identity_top_hit"
        )


    ## -------------------------------------------------------------- ##
    ## Gene / product
    ##
    ## Only high-confidence tier gets a specific product transfer.
    ##
    ## Prefer the best hit carrying the consensus COG when possible.
    ## -------------------------------------------------------------- ##

    transferred_gene = ""

    transferred_product = ""

    gene_product_source = ""


    if status == "high_confidence":

        source_pool = pool


        if transferred_cog:

            same_cog = pool[
                pool[
                    "globdb_cog"
                ]
                ==
                transferred_cog
            ]


            if len(
                same_cog
            ) > 0:

                source_pool = (
                    same_cog
                )


        if len(
            source_pool
        ) > 0:

            source_hit = (
                source_pool
                .sort_values(
                    [
                        "bits",
                        "pident",
                        "qcov",
                        "tcov",
                        "evalue",
                        "target",
                    ],
                    ascending=[
                        False,
                        False,
                        False,
                        False,
                        True,
                        True,
                    ],
                    kind="stable",
                )
                .iloc[0]
            )


            transferred_gene = clean(
                source_hit.get(
                    "globdb_gene",
                    ""
                )
            )


            transferred_product = clean(
                source_hit.get(
                    "globdb_product",
                    ""
                )
            )


            gene_product_source = clean(
                source_hit[
                    "target"
                ]
            )


    ## -------------------------------------------------------------- ##
    ## Assemble result
    ## -------------------------------------------------------------- ##

    row = {
        "neighborhood_protein_id":
            query,

        "genome":
            query_row.genome,

        "protein_id":
            query_row.protein_id,

        "transfer_status":
            status,

        "n_MMseqs_hits":
            len(
                group
            ),

        "n_usable_annotated_hits":
            len(
                annotated
            ),

        "n_high_confidence_hits":
            len(
                high
            ),

        "n_family_level_hits":
            len(
                family
            ),

        "transferred_cog":
            transferred_cog,

        "cog_transfer_basis":
            cog_basis,

        "cog_supporting_hits":
            cog_support,

        "cog_bearing_hits_in_transfer_pool":
            n_cog_hits,

        "cog_support_fraction":
            (
                f"{cog_fraction:.6f}"
                if n_cog_hits > 0
                else ""
            ),

        "cog_distribution":
            cog_distribution,

        "transferred_gene":
            transferred_gene,

        "transferred_product":
            transferred_product,

        "gene_product_source_target":
            gene_product_source,
    }


    ## -------------------------------------------------------------- ##
    ## Raw best overall hit
    ## -------------------------------------------------------------- ##

    if best_overall is not None:

        row.update(
            {
                "best_overall_target":
                    clean(
                        best_overall[
                            "target"
                        ]
                    ),

                "best_overall_pident":
                    best_overall[
                        "pident"
                    ],

                "best_overall_qcov":
                    best_overall[
                        "qcov"
                    ],

                "best_overall_tcov":
                    best_overall[
                        "tcov"
                    ],

                "best_overall_evalue":
                    best_overall[
                        "evalue"
                    ],

                "best_overall_bits":
                    best_overall[
                        "bits"
                    ],
            }
        )

    else:

        row.update(
            {
                "best_overall_target":
                    "",

                "best_overall_pident":
                    "",

                "best_overall_qcov":
                    "",

                "best_overall_tcov":
                    "",

                "best_overall_evalue":
                    "",

                "best_overall_bits":
                    "",
            }
        )


    ## -------------------------------------------------------------- ##
    ## Best accepted annotated hit
    ## -------------------------------------------------------------- ##

    if best_annotated is not None:

        row.update(
            {
                "best_annotated_target":
                    clean(
                        best_annotated[
                            "target"
                        ]
                    ),

                "best_annotated_pident":
                    best_annotated[
                        "pident"
                    ],

                "best_annotated_qcov":
                    best_annotated[
                        "qcov"
                    ],

                "best_annotated_tcov":
                    best_annotated[
                        "tcov"
                    ],

                "best_annotated_cog":
                    clean(
                        best_annotated.get(
                            "globdb_cog",
                            ""
                        )
                    ),

                "best_annotated_gene":
                    clean(
                        best_annotated.get(
                            "globdb_gene",
                            ""
                        )
                    ),

                "best_annotated_product":
                    clean(
                        best_annotated.get(
                            "globdb_product",
                            ""
                        )
                    ),
            }
        )

    else:

        row.update(
            {
                "best_annotated_target":
                    "",

                "best_annotated_pident":
                    "",

                "best_annotated_qcov":
                    "",

                "best_annotated_tcov":
                    "",

                "best_annotated_cog":
                    "",

                "best_annotated_gene":
                    "",

                "best_annotated_product":
                    "",
            }
        )


    transfer_rows.append(
        row
    )


transfer = pd.DataFrame(
    transfer_rows
)


## ================================================================== ##
## 8. Strict transfer-table QC
## ================================================================== ##

if len(
    transfer
) != EXPECTED_QUERIES:

    fail(
        f"Expected {EXPECTED_QUERIES} transfer rows "
        f"but created {len(transfer)}."
    )


if transfer[
    "neighborhood_protein_id"
].duplicated().any():

    fail(
        "Duplicate query IDs in transfer table."
    )


if set(
    transfer[
        "neighborhood_protein_id"
    ]
) != query_ids:

    fail(
        "Transfer-table query IDs differ from "
        "original 845-query set."
    )


## ================================================================== ##
## 9. Write transfer table
## ================================================================== ##

transfer = (
    transfer
    .sort_values(
        [
            "transfer_status",
            "genome",
            "protein_id",
        ],
        kind="stable",
    )
    .reset_index(
        drop=True
    )
)


transfer.to_csv(
    OUT_TRANSFER,
    sep="\t",
    index=False,
    na_rep="",
)


## ================================================================== ##
## 10. Attach transfer results to the 845-protein metadata
## ================================================================== ##

transfer_for_merge = transfer.drop(
    columns=[
        "genome",
        "protein_id",
    ]
)


annotated_metadata = queries.merge(
    transfer_for_merge,
    on="neighborhood_protein_id",
    how="left",
    validate="one_to_one",
)


if len(
    annotated_metadata
) != EXPECTED_QUERIES:

    fail(
        "Metadata annotation join changed row count."
    )


annotated_metadata.to_csv(
    OUT_METADATA_ANNOTATED,
    sep="\t",
    index=False,
    na_rep="",
)


## ================================================================== ##
## 11. Attach transferred annotation to the 845 neighborhood rows
## ================================================================== ##

neighborhoods = pd.read_csv(
    NEIGHBORHOODS,
    sep="\t",
    dtype=str,
    keep_default_na=False,
)


required_neighborhood = {
    "genome",
    "neighbor_protein_id",
}


missing = (
    required_neighborhood
    -
    set(
        neighborhoods.columns
    )
)


if missing:

    fail(
        "Observed-neighborhood table is missing:\n"
        +
        "\n".join(
            sorted(
                missing
            )
        )
    )


gene_transfer = transfer[
    [
        "genome",
        "protein_id",
        "transfer_status",

        "transferred_cog",
        "cog_transfer_basis",
        "cog_support_fraction",

        "transferred_gene",
        "transferred_product",

        "best_annotated_target",
        "best_annotated_pident",
        "best_annotated_qcov",
        "best_annotated_tcov",
        "best_annotated_cog",
        "best_annotated_gene",
        "best_annotated_product",
    ]
].rename(
    columns={
        "protein_id":
            "neighbor_protein_id",

        "transfer_status":
            "neighbor_annotation_transfer_status",

        "transferred_cog":
            "neighbor_transferred_cog",

        "cog_transfer_basis":
            "neighbor_cog_transfer_basis",

        "cog_support_fraction":
            "neighbor_cog_support_fraction",

        "transferred_gene":
            "neighbor_transferred_gene",

        "transferred_product":
            "neighbor_transferred_product",

        "best_annotated_target":
            "neighbor_best_reference_target",

        "best_annotated_pident":
            "neighbor_best_reference_pident",

        "best_annotated_qcov":
            "neighbor_best_reference_qcov",

        "best_annotated_tcov":
            "neighbor_best_reference_tcov",

        "best_annotated_cog":
            "neighbor_best_reference_cog",

        "best_annotated_gene":
            "neighbor_best_reference_gene",

        "best_annotated_product":
            "neighbor_best_reference_product",
    }
)


annotated_neighborhoods = neighborhoods.merge(
    gene_transfer,
    on=[
        "genome",
        "neighbor_protein_id",
    ],
    how="left",
    validate="many_to_one",
)


if len(
    annotated_neighborhoods
) != len(
    neighborhoods
):

    fail(
        "Neighborhood annotation join changed row count."
    )


if annotated_neighborhoods[
    "neighbor_annotation_transfer_status"
].eq("").any():

    fail(
        "One or more neighborhood genes failed "
        "annotation-transfer mapping."
    )


annotated_neighborhoods.to_csv(
    OUT_NEIGHBORHOODS_ANNOTATED,
    sep="\t",
    index=False,
    na_rep="",
)


## ================================================================== ##
## 12. QC
## ================================================================== ##

status_counts = (
    transfer[
        "transfer_status"
    ]
    .value_counts()
)


n_high = int(
    status_counts.get(
        "high_confidence",
        0,
    )
)


n_family = int(
    status_counts.get(
        "family_level",
        0,
    )
)


n_none = int(
    status_counts.get(
        "no_transfer",
        0,
    )
)


n_cog = int(
    (
        transfer[
            "transferred_cog"
        ]
        !=
        ""
    ).sum()
)


n_gene = int(
    (
        transfer[
            "transferred_gene"
        ]
        !=
        ""
    ).sum()
)


n_product = int(
    (
        transfer[
            "transferred_product"
        ]
        !=
        ""
    ).sum()
)


n_no_mmseqs_hit = int(
    (
        transfer[
            "n_MMseqs_hits"
        ]
        ==
        0
    ).sum()
)


n_hits_but_no_usable_annotation = int(
    (
        (
            transfer[
                "n_MMseqs_hits"
            ]
            >
            0
        )
        &
        (
            transfer[
                "n_usable_annotated_hits"
            ]
            ==
            0
        )
    ).sum()
)


qc = pd.DataFrame(
    [
        (
            "queries",
            len(
                transfer
            ),
        ),

        (
            "raw_MMseqs_rows",
            len(
                hits
            ),
        ),

        (
            "queries_with_no_MMseqs_hit",
            n_no_mmseqs_hit,
        ),

        (
            "queries_with_hits_but_no_usable_annotation",
            n_hits_but_no_usable_annotation,
        ),

        (
            "high_confidence_transfers",
            n_high,
        ),

        (
            "family_level_transfers",
            n_family,
        ),

        (
            "no_transfer",
            n_none,
        ),

        (
            "queries_with_transferred_COG",
            n_cog,
        ),

        (
            "queries_with_transferred_gene",
            n_gene,
        ),

        (
            "queries_with_transferred_product",
            n_product,
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
print("STAGE 56 COMPLETE")
print("=" * 80)

print(
    f"Queries:                         "
    f"{len(transfer):,}"
)

print(
    f"High-confidence transfers:       "
    f"{n_high:,}"
)

print(
    f"Family-level transfers:          "
    f"{n_family:,}"
)

print(
    f"No transfer:                     "
    f"{n_none:,}"
)

print()
print(
    f"Transferred COG:                 "
    f"{n_cog:,}"
)

print(
    f"Transferred gene name:           "
    f"{n_gene:,}"
)

print(
    f"Transferred product:             "
    f"{n_product:,}"
)

print()
print(
    f"No MMseqs hit:                   "
    f"{n_no_mmseqs_hit:,}"
)

print(
    f"Hits but no usable annotation:   "
    f"{n_hits_but_no_usable_annotation:,}"
)

print()
print(
    f"Transfer table:\n  "
    f"{OUT_TRANSFER}"
)

print(
    f"Annotated protein metadata:\n  "
    f"{OUT_METADATA_ANNOTATED}"
)

print(
    f"Annotated neighborhoods:\n  "
    f"{OUT_NEIGHBORHOODS_ANNOTATED}"
)

print(
    f"Candidate-hit evidence:\n  "
    f"{OUT_CONSIDERED_HITS}"
)

print(
    f"QC:\n  "
    f"{OUT_QC}"
)
