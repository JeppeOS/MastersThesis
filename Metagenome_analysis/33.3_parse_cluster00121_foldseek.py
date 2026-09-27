#!/usr/bin/env python3

from pathlib import Path
import io
import sys
import tarfile

import numpy as np
import pandas as pd


## ================================================================== ##
## STAGE 69
##
## PARSE CLUSTER_00121 FOLDSEEK WEB-SERVER RESULTS
##
## Inputs:
##
##     5 downloaded .gz archives
##
## Each archive is expected to contain:
##
##     alis_afdb-swissprot.m8
##     alis_afdb50.m8
##     alis_cath50.m8
##     alis_pdb100.m8
##
## Outputs:
##
##     all hits
##     top 20 hits/query/database
##     best hit/query/database
##     recurrent structural targets across the five queries
##
## No biological hit-selection rule is imposed here.
## ================================================================== ##


PROJECT = (
    Path.home()
    / "methanotrophs"
    / "methanotroph_project"
    / "jeppe"
    / "metagenome"
)


ROOT = (
    PROJECT
    / "comparative_analysis"
    / "Cluster_00121"
    / "structure_analysis"
    / "foldseek"
)


ARCHIVE_DIR = (
    ROOT
    / "raw_archives"
)


TABLE_DIR = (
    ROOT
    / "tables"
)


TABLE_DIR.mkdir(
    parents=True,
    exist_ok=True,
)


OUT_ALL = (
    TABLE_DIR
    / "Cluster_00121_foldseek_all_hits.tsv"
)


OUT_TOP20 = (
    TABLE_DIR
    / "Cluster_00121_foldseek_top20_per_query_database.tsv"
)


OUT_BEST = (
    TABLE_DIR
    / "Cluster_00121_foldseek_best_hit_per_query_database.tsv"
)


OUT_RECURRENT = (
    TABLE_DIR
    / "Cluster_00121_foldseek_recurrent_targets.tsv"
)


EXPECTED_ARCHIVES = 5


## ================================================================== ##
## Foldseek web-server M8 schema
## ================================================================== ##

COLUMNS = [
    "query",
    "target",
    "fident",
    "alnlen",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "tstart",
    "tend",
    "prob",
    "evalue",
    "bits",
    "qlen",
    "tlen",
    "qaln",
    "taln",
    "target_coordinates",
    "target_sequence",
    "taxid",
    "taxname",
]


DATABASE_FILES = {

    "afdb_swissprot":
        "alis_afdb-swissprot.m8",

    "afdb50":
        "alis_afdb50.m8",

    "cath50":
        "alis_cath50.m8",

    "pdb100":
        "alis_pdb100.m8",
}


## ================================================================== ##
## Helpers
## ================================================================== ##

def fail(message):

    print(
        f"\nERROR: {message}",
        file=sys.stderr,
    )

    sys.exit(1)


def split_target(value):

    value = str(
        value
    ).strip()


    if not value:

        return (
            "",
            "",
        )


    parts = value.split(
        None,
        1,
    )


    accession = parts[
        0
    ]


    description = (
        parts[
            1
        ]
        if len(
            parts
        )
        == 2
        else
        ""
    )


    return (
        accession,
        description,
    )


def find_archive_member(
    tar,
    required_basename,
):

    matches = [
        member

        for member
        in tar.getmembers()

        if (
            member.isfile()
            and
            Path(
                member.name
            ).name
            ==
            required_basename
        )
    ]


    if len(
        matches
    ) != 1:

        fail(
            f"Expected exactly one "
            f"{required_basename} in archive; "
            f"found {len(matches)}."
        )


    return matches[
        0
    ]


## ================================================================== ##
## Discover archives
## ================================================================== ##

archives = sorted(
    ARCHIVE_DIR.glob(
        "*.gz"
    )
)


print("=" * 80)

print(
    "STAGE 69 - PARSE CLUSTER_00121 FOLDSEEK"
)

print("=" * 80)

print()

print(
    f"Foldseek archives found: {len(archives)}"
)


if len(
    archives
) != EXPECTED_ARCHIVES:

    fail(
        f"Expected {EXPECTED_ARCHIVES} archives; "
        f"found {len(archives)}."
    )


## ================================================================== ##
## Parse all four databases from all five archives
## ================================================================== ##

all_tables = []


for archive in archives:

    query_id = archive.name.removesuffix(
        ".gz"
    )


    print()
    print(
        f"Processing: {query_id}"
    )


    with tarfile.open(
        archive,
        mode="r:gz",
    ) as tar:

        for database, filename in DATABASE_FILES.items():

            member = find_archive_member(
                tar,
                filename,
            )


            handle = tar.extractfile(
                member
            )


            if handle is None:

                fail(
                    f"Could not read {filename} "
                    f"from {archive.name}."
                )


            raw = handle.read()


            if len(
                raw
            ) == 0:

                print(
                    f"  {database:<16} 0 hits"
                )

                continue


            df = pd.read_csv(
                io.BytesIO(
                    raw
                ),
                sep="\t",
                header=None,
                names=COLUMNS,
                dtype={
                    "query":
                        str,

                    "target":
                        str,

                    "taxid":
                        str,

                    "taxname":
                        str,
                },
            )


            ## ---------------------------------------------------- ##
            ## Numeric columns
            ## ---------------------------------------------------- ##

            numeric_columns = [
                "fident",
                "alnlen",
                "mismatch",
                "gapopen",
                "qstart",
                "qend",
                "tstart",
                "tend",
                "prob",
                "evalue",
                "bits",
                "qlen",
                "tlen",
            ]


            for column in numeric_columns:

                df[
                    column
                ] = pd.to_numeric(
                    df[
                        column
                    ],
                    errors="coerce",
                )


            ## ---------------------------------------------------- ##
            ## Preserve Foldseek order as database rank.
            ##
            ## Each downloaded archive contains one query.
            ## ---------------------------------------------------- ##

            df.insert(
                0,
                "query_id",
                query_id,
            )


            df.insert(
                1,
                "database",
                database,
            )


            df.insert(
                2,
                "database_rank",
                np.arange(
                    1,
                    len(
                        df
                    )
                    + 1,
                ),
            )


            ## ---------------------------------------------------- ##
            ## Parse target accession / description
            ## ---------------------------------------------------- ##

            parsed = df[
                "target"
            ].apply(
                split_target
            )


            df[
                "target_accession"
            ] = [
                item[
                    0
                ]

                for item
                in parsed
            ]


            df[
                "target_description"
            ] = [
                item[
                    1
                ]

                for item
                in parsed
            ]


            ## ---------------------------------------------------- ##
            ## Structural coverage
            ##
            ## Coordinate span is used because Foldseek alnlen can
            ## include alignment gaps.
            ## ---------------------------------------------------- ##

            df[
                "query_aligned_span"
            ] = (
                (
                    df[
                        "qend"
                    ]
                    -
                    df[
                        "qstart"
                    ]
                )
                .abs()
                +
                1
            )


            df[
                "target_aligned_span"
            ] = (
                (
                    df[
                        "tend"
                    ]
                    -
                    df[
                        "tstart"
                    ]
                )
                .abs()
                +
                1
            )


            df[
                "query_coverage"
            ] = (
                df[
                    "query_aligned_span"
                ]
                /
                df[
                    "qlen"
                ]
            )


            df[
                "target_coverage"
            ] = (
                df[
                    "target_aligned_span"
                ]
                /
                df[
                    "tlen"
                ]
            )


            all_tables.append(
                df
            )


            print(
                f"  {database:<16} "
                f"{len(df):>6} hits"
            )


## ================================================================== ##
## Combine
## ================================================================== ##

if not all_tables:

    fail(
        "No Foldseek hits were parsed."
    )


x = pd.concat(
    all_tables,
    ignore_index=True,
)


## ================================================================== ##
## Basic QC
## ================================================================== ##

queries = sorted(
    x[
        "query_id"
    ].unique()
)


databases = sorted(
    x[
        "database"
    ].unique()
)


if len(
    queries
) != EXPECTED_ARCHIVES:

    fail(
        f"Expected results for five queries; "
        f"found {len(queries)}."
    )


missing_pairs = []


for query in queries:

    for database in DATABASE_FILES:

        n = len(
            x[
                (
                    x[
                        "query_id"
                    ]
                    ==
                    query
                )
                &
                (
                    x[
                        "database"
                    ]
                    ==
                    database
                )
            ]
        )


        if n == 0:

            missing_pairs.append(
                (
                    query,
                    database,
                )
            )


if missing_pairs:

    print()
    print(
        "WARNING: query/database combinations "
        "with zero hits:"
    )

    for query, database in missing_pairs:

        print(
            f"  {query}\t{database}"
        )


## ================================================================== ##
## Write all hits
## ================================================================== ##

x.to_csv(
    OUT_ALL,
    sep="\t",
    index=False,
    na_rep="",
    float_format="%.6g",
)


## ================================================================== ##
## Top 20 per query/database
## ================================================================== ##

top20 = (
    x[
        x[
            "database_rank"
        ]
        <=
        20
    ]
    .copy()
)


top20.to_csv(
    OUT_TOP20,
    sep="\t",
    index=False,
    na_rep="",
    float_format="%.6g",
)


## ================================================================== ##
## Best Foldseek hit per query/database
## ================================================================== ##

best = (
    x
    .sort_values(
        [
            "query_id",
            "database",
            "database_rank",
        ],
        kind="stable",
    )
    .drop_duplicates(
        [
            "query_id",
            "database",
        ],
        keep="first",
    )
    .reset_index(
        drop=True
    )
)


best.to_csv(
    OUT_BEST,
    sep="\t",
    index=False,
    na_rep="",
    float_format="%.6g",
)


## ================================================================== ##
## Recurrent targets across independent queries
##
## A target occurring in all/most of the five queries is particularly
## useful for deciding whether the structural interpretation converges.
## ================================================================== ##

recurrent = (
    x
    .groupby(
        [
            "database",
            "target_accession",
            "target_description",
        ],
        dropna=False,
    )
    .agg(
        n_queries=(
            "query_id",
            "nunique",
        ),

        median_database_rank=(
            "database_rank",
            "median",
        ),

        best_database_rank=(
            "database_rank",
            "min",
        ),

        median_probability=(
            "prob",
            "median",
        ),

        minimum_evalue=(
            "evalue",
            "min",
        ),

        median_query_coverage=(
            "query_coverage",
            "median",
        ),

        median_target_coverage=(
            "target_coverage",
            "median",
        ),

        median_fident=(
            "fident",
            "median",
        ),

        maximum_bits=(
            "bits",
            "max",
        ),
    )
    .reset_index()
)


recurrent = recurrent.sort_values(
    [
        "n_queries",
        "median_database_rank",
        "median_probability",
        "minimum_evalue",
    ],
    ascending=[
        False,
        True,
        False,
        True,
    ],
    kind="stable",
)


recurrent.to_csv(
    OUT_RECURRENT,
    sep="\t",
    index=False,
    na_rep="",
    float_format="%.6g",
)


## ================================================================== ##
## Terminal summaries
## ================================================================== ##

print()
print("=" * 80)

print(
    "STAGE 69 COMPLETE"
)

print("=" * 80)

print()

print(
    f"Queries parsed:    {len(queries)}"
)

print(
    f"Total hits:        {len(x):,}"
)

print()


## ------------------------------------------------------------------ ##
## Best hit from each database/query
## ------------------------------------------------------------------ ##

display_columns = [
    "query_id",
    "database",
    "target_accession",
    "target_description",
    "prob",
    "evalue",
    "fident",
    "query_coverage",
    "target_coverage",
]


print(
    "BEST HIT PER QUERY / DATABASE"
)

print(
    best[
        display_columns
    ]
    .to_string(
        index=False,
        max_colwidth=65,
    )
)


## ------------------------------------------------------------------ ##
## Recurrent targets in >=4/5 queries
## ------------------------------------------------------------------ ##

strong_recurrent = recurrent[
    recurrent[
        "n_queries"
    ]
    >=
    4
].copy()


print()
print(
    "TARGETS RECOVERED IN >=4/5 QUERIES"
)

if len(
    strong_recurrent
) == 0:

    print(
        "  None."
    )


else:

    print(
        strong_recurrent[
            [
                "database",
                "target_accession",
                "target_description",
                "n_queries",
                "median_database_rank",
                "median_probability",
                "median_query_coverage",
                "median_target_coverage",
                "median_fident",
            ]
        ]
        .head(
            40
        )
        .to_string(
            index=False,
            max_colwidth=65,
        )
    )


## ------------------------------------------------------------------ ##
## PDB100 top hits
## ------------------------------------------------------------------ ##

print()
print(
    "TOP PDB100 HIT FOR EACH QUERY"
)


pdb_best = best[
    best[
        "database"
    ]
    ==
    "pdb100"
]


print(
    pdb_best[
        display_columns
    ]
    .to_string(
        index=False,
        max_colwidth=80,
    )
)


## ------------------------------------------------------------------ ##
## Swiss-Prot top hits
## ------------------------------------------------------------------ ##

print()
print(
    "TOP AFDB/SWISS-PROT HIT FOR EACH QUERY"
)


sp_best = best[
    best[
        "database"
    ]
    ==
    "afdb_swissprot"
]


print(
    sp_best[
        display_columns
    ]
    .to_string(
        index=False,
        max_colwidth=80,
    )
)


print()
print(
    "Outputs:"
)

print(
    f"  {OUT_ALL}"
)

print(
    f"  {OUT_TOP20}"
)

print(
    f"  {OUT_BEST}"
)

print(
    f"  {OUT_RECURRENT}"
)
