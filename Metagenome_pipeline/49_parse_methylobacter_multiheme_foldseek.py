#!/usr/bin/env python3

from pathlib import Path
import io
import sys
import tarfile

import numpy as np
import pandas as pd


## ================================================================== ##
## STAGE 72
##
## PARSE FOLDSEEK RESULTS FOR THREE METHYLOBACTER MULTIHEME
## CANDIDATES
##
## Candidates are intentionally treated independently:
##
##     mb_c00042_b14sb1
##         5-heme
##         Cluster_00042
##
##     mb_nomatch_b15sb3
##         7-heme
##         no GlobDB cluster match
##
##     mb_c00288_b16sb1
##         8-heme
##         Cluster_00288 / secondary Cluster_00218
##
## Outputs:
##
##     - all Foldseek hits
##     - top 20 per query/database
##     - best hit per query/database
##     - cytochrome/redox-relevant hits
##     - compact candidate summary
##
## No functional identity is assigned automatically.
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
    / "Methylobacter_multiheme"
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
    / "Methylobacter_multiheme_foldseek_all_hits.tsv"
)


OUT_TOP20 = (
    TABLE_DIR
    / "Methylobacter_multiheme_foldseek_top20_per_query_database.tsv"
)


OUT_BEST = (
    TABLE_DIR
    / "Methylobacter_multiheme_foldseek_best_hit_per_query_database.tsv"
)


OUT_REDox = (
    TABLE_DIR
    / "Methylobacter_multiheme_foldseek_cytochrome_redox_hits.tsv"
)


OUT_SUMMARY = (
    TABLE_DIR
    / "Methylobacter_multiheme_foldseek_candidate_summary.tsv"
)


EXPECTED_ARCHIVES = 3


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
## Candidate metadata
## ================================================================== ##

CANDIDATES = {

    "mb_c00042_b14sb1": {
        "taxon":
            "Methylobacter_A",

        "heme_count":
            5,

        "globdb_mapping":
            "Cluster_00042",

        "caution":
            "",
    },

    "mb_nomatch_b15sb3": {
        "taxon":
            "Methylobacter_C sp002256465",

        "heme_count":
            7,

        "globdb_mapping":
            "no_match",

        "caution":
            "",
    },

    "mb_c00288_b16sb1": {
        "taxon":
            "Methylobacter_C sp002256465",

        "heme_count":
            8,

        "globdb_mapping":
            "Cluster_00288; secondary Cluster_00218",

        "caution":
            "MAG contamination 17.61%; interpret locus cautiously",
    },
}


## ================================================================== ##
## Cytochrome / EET-related description terms
##
## This is deliberately broad.
##
## A keyword match means:
##
##     "worth inspecting"
##
## not:
##
##     "this protein has this function"
## ================================================================== ##

REDOX_PATTERN = (
    r"cytochrome|"
    r"MtrA|MtrB|MtrC|MtrD|MtrE|MtrF|"
    r"MtoA|MtoB|"
    r"OmcA|OmcB|OmcS|OmcZ|"
    r"decaheme|multiheme|multi-heme|"
    r"electron transfer|"
    r"electron-transfer|"
    r"outer membrane.*heme|"
    r"heme.*outer membrane"
)


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
    "STAGE 72 - PARSE METHYLOBACTER MULTIHEME FOLDSEEK"
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


archive_ids = sorted(
    archive.name.removesuffix(
        ".gz"
    )

    for archive
    in archives
)


if set(
    archive_ids
) != set(
    CANDIDATES
):

    fail(
        "Archive names do not match expected candidates.\n\n"
        f"Observed: {archive_ids}\n"
        f"Expected: {sorted(CANDIDATES)}"
    )


## ================================================================== ##
## Parse all results
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
            ## Query / database identity
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


            ## Preserve Foldseek result ordering.
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
            ## Candidate metadata
            ## ---------------------------------------------------- ##

            metadata = CANDIDATES[
                query_id
            ]


            df[
                "taxon"
            ] = metadata[
                "taxon"
            ]


            df[
                "heme_count"
            ] = metadata[
                "heme_count"
            ]


            df[
                "globdb_mapping"
            ] = metadata[
                "globdb_mapping"
            ]


            df[
                "caution"
            ] = metadata[
                "caution"
            ]


            ## ---------------------------------------------------- ##
            ## Target accession / description
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
            ## Structural alignment coverage
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

top20 = x[
    x[
        "database_rank"
    ]
    <=
    20
].copy()


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
## Cytochrome / redox-related structural hits
## ================================================================== ##

redox = x[
    x[
        "target_description"
    ]
    .fillna("")
    .str.contains(
        REDOX_PATTERN,
        case=False,
        regex=True,
    )
].copy()


redox = redox.sort_values(
    [
        "query_id",
        "database",
        "database_rank",
    ],
    kind="stable",
)


redox.to_csv(
    OUT_REDox,
    sep="\t",
    index=False,
    na_rep="",
    float_format="%.6g",
)


## ================================================================== ##
## Compact candidate summary
##
## For each candidate:
##
##     - top AFDB50 hit
##     - top Swiss-Prot hit
##     - top PDB100 hit
##     - highest-ranked AFDB50 cytochrome/redox hit
##
## No combined score is introduced.
## ================================================================== ##

summary_rows = []


for query_id in sorted(
    CANDIDATES
):

    metadata = CANDIDATES[
        query_id
    ]


    row = {
        "query_id":
            query_id,

        "taxon":
            metadata[
                "taxon"
            ],

        "heme_count":
            metadata[
                "heme_count"
            ],

        "globdb_mapping":
            metadata[
                "globdb_mapping"
            ],

        "caution":
            metadata[
                "caution"
            ],
    }


    for database in [
        "afdb50",
        "afdb_swissprot",
        "pdb100",
    ]:

        hit = best[
            (
                best[
                    "query_id"
                ]
                ==
                query_id
            )
            &
            (
                best[
                    "database"
                ]
                ==
                database
            )
        ]


        prefix = database


        if len(
            hit
        ) == 1:

            h = hit.iloc[
                0
            ]


            row[
                f"{prefix}_target"
            ] = h[
                "target_accession"
            ]


            row[
                f"{prefix}_description"
            ] = h[
                "target_description"
            ]


            row[
                f"{prefix}_prob"
            ] = h[
                "prob"
            ]


            row[
                f"{prefix}_evalue"
            ] = h[
                "evalue"
            ]


            row[
                f"{prefix}_fident"
            ] = h[
                "fident"
            ]


            row[
                f"{prefix}_qcov"
            ] = h[
                "query_coverage"
            ]


            row[
                f"{prefix}_tcov"
            ] = h[
                "target_coverage"
            ]


    candidate_redox = redox[
        (
            redox[
                "query_id"
            ]
            ==
            query_id
        )
        &
        (
            redox[
                "database"
            ]
            ==
            "afdb50"
        )
    ].copy()


    candidate_redox = candidate_redox.sort_values(
        [
            "database_rank",
            "evalue",
        ],
        ascending=[
            True,
            True,
        ],
        kind="stable",
    )


    if len(
        candidate_redox
    ):

        h = candidate_redox.iloc[
            0
        ]


        row[
            "best_redox_rank"
        ] = h[
            "database_rank"
        ]


        row[
            "best_redox_target"
        ] = h[
            "target_accession"
        ]


        row[
            "best_redox_description"
        ] = h[
            "target_description"
        ]


        row[
            "best_redox_prob"
        ] = h[
            "prob"
        ]


        row[
            "best_redox_evalue"
        ] = h[
            "evalue"
        ]


        row[
            "best_redox_fident"
        ] = h[
            "fident"
        ]


        row[
            "best_redox_qcov"
        ] = h[
            "query_coverage"
        ]


        row[
            "best_redox_tcov"
        ] = h[
            "target_coverage"
        ]


    summary_rows.append(
        row
    )


summary = pd.DataFrame(
    summary_rows
)


summary.to_csv(
    OUT_SUMMARY,
    sep="\t",
    index=False,
    na_rep="",
    float_format="%.6g",
)


## ================================================================== ##
## Terminal report
## ================================================================== ##

display_columns = [
    "database_rank",
    "target_accession",
    "target_description",
    "prob",
    "evalue",
    "fident",
    "query_coverage",
    "target_coverage",
]


print()
print("=" * 80)

print(
    "STAGE 72 COMPLETE"
)

print("=" * 80)

print()

print(
    f"Candidates parsed: {len(CANDIDATES)}"
)

print(
    f"Total hits:        {len(x):,}"
)

print(
    f"Redox/cytochrome hits: {len(redox):,}"
)


for query_id in sorted(
    CANDIDATES
):

    metadata = CANDIDATES[
        query_id
    ]


    print()
    print("=" * 100)

    print(
        f"{query_id} | "
        f"{metadata['taxon']} | "
        f"{metadata['heme_count']}-heme | "
        f"{metadata['globdb_mapping']}"
    )

    print("=" * 100)


    ## ------------------------------------------------------------ ##
    ## Best hits by major database
    ## ------------------------------------------------------------ ##

    qbest = best[
        (
            best[
                "query_id"
            ]
            ==
            query_id
        )
        &
        (
            best[
                "database"
            ]
            .isin(
                [
                    "afdb50",
                    "afdb_swissprot",
                    "pdb100",
                ]
            )
        )
    ].copy()


    print()
    print(
        "BEST HIT PER DATABASE"
    )


    print(
        qbest[
            [
                "database",
                "target_accession",
                "target_description",
                "prob",
                "evalue",
                "fident",
                "query_coverage",
                "target_coverage",
            ]
        ]
        .to_string(
            index=False,
            max_colwidth=80,
        )
    )


    ## ------------------------------------------------------------ ##
    ## Top AFDB50 redox/cytochrome matches
    ## ------------------------------------------------------------ ##

    qr = redox[
        (
            redox[
                "query_id"
            ]
            ==
            query_id
        )
        &
        (
            redox[
                "database"
            ]
            ==
            "afdb50"
        )
    ].copy()


    print()
    print(
        "TOP AFDB50 CYTOCHROME / REDOX-RELATED HITS"
    )


    if len(
        qr
    ) == 0:

        print(
            "  None."
        )


    else:

        print(
            qr[
                display_columns
            ]
            .head(
                20
            )
            .to_string(
                index=False,
                max_colwidth=85,
            )
        )


    ## ------------------------------------------------------------ ##
    ## PDB cytochrome/redox matches
    ## ------------------------------------------------------------ ##

    qp = redox[
        (
            redox[
                "query_id"
            ]
            ==
            query_id
        )
        &
        (
            redox[
                "database"
            ]
            ==
            "pdb100"
        )
    ].copy()


    print()
    print(
        "TOP PDB100 CYTOCHROME / REDOX-RELATED HITS"
    )


    if len(
        qp
    ) == 0:

        print(
            "  None."
        )


    else:

        print(
            qp[
                display_columns
            ]
            .head(
                10
            )
            .to_string(
                index=False,
                max_colwidth=85,
            )
        )


print()
print(
    "Outputs:"
)

for path in [
    OUT_ALL,
    OUT_TOP20,
    OUT_BEST,
    OUT_REDox,
    OUT_SUMMARY,
]:

    print(
        f"  {path}"
    )
