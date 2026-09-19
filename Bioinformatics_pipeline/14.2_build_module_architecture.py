#!/usr/bin/env python3

from pathlib import Path
from itertools import combinations
from collections import Counter, defaultdict
import sys

import numpy as np
import pandas as pd


## ================================================================== ##
## STAGE 14B - SUPPORTED MODULE ARCHITECTURE
##
## PURPOSE
## -------
##
## Integrate:
##
##   * MCL module membership
##   * genome-level co-presence
##   * physical colocalization
##   * adjacency
##   * local <=5 / <=10 / <=20 kb evidence
##   * closest-pair orientation
##   * transcription-normalized gene order
##   * Stage-13D focal-family consensus directions
##
## into one architecture graph for every MCL module.
##
##
## IMPORTANT
## ---------
##
## No arbitrary percentage cutoff is used.
##
## Every possible within-module cluster pair is retained.
##
## Local architecture components are defined only by observed
## <=20-kb relationships.
##
## This does NOT imply that an entire MCL module is an operon.
## ================================================================== ##


## ================================================================== ##
## Paths
## ================================================================== ##

HERE = Path(__file__).resolve().parent
WORKFLOW = HERE.parent.parent


MODULE_MEMBERSHIP = (
    WORKFLOW
    / "07_cooccurrence_network"
    / "module_membership.tsv"
)


PAIR_OBSERVATIONS = (
    WORKFLOW
    / "12_neighborhood_extraction"
    / "module_member_pair_observations.tsv"
)


MMSEQ_CONSENSUS = (
    WORKFLOW
    / "13_neighborhood_conservation"
    / "13D_consensus_neighborhoods"
    / "focal_mmseqs_consensus_associations.tsv"
)


CLUSTER_EVIDENCE = (
    WORKFLOW
    / "14_module_reports"
    / "14A_report_data"
    / "global"
    / "module_member_cluster_evidence.tsv"
)


## ================================================================== ##
## Output structure
## ================================================================== ##

GLOBAL_DIR = HERE / "global"
MODULE_DIR = HERE / "modules"

GLOBAL_DIR.mkdir(
    parents=True,
    exist_ok=True,
)

MODULE_DIR.mkdir(
    parents=True,
    exist_ok=True,
)


OUT_NODES = (
    GLOBAL_DIR
    / "module_architecture_nodes.tsv"
)


OUT_EDGES = (
    GLOBAL_DIR
    / "module_architecture_edges.tsv"
)


OUT_SUMMARY = (
    GLOBAL_DIR
    / "module_architecture_summary.tsv"
)


OUT_QC = (
    GLOBAL_DIR
    / "module_architecture_qc.tsv"
)


## ================================================================== ##
## Locked upstream expectations
## ================================================================== ##

EXPECTED_MODULES = 35
EXPECTED_CLUSTERS = 156
EXPECTED_POSSIBLE_PAIRS = 560

EXPECTED_PAIR_OBSERVATIONS = 34_037

EXPECTED_SAME_CONTIG_OBSERVATIONS = 8_080
EXPECTED_WITHIN_5KB_OBSERVATIONS = 2_209
EXPECTED_WITHIN_10KB_OBSERVATIONS = 2_564
EXPECTED_WITHIN_20KB_OBSERVATIONS = 2_757
EXPECTED_ADJACENT_OBSERVATIONS = 947

EXPECTED_COPRESENT_PAIRS = 546
EXPECTED_SAME_CONTIG_PAIRS = 474
EXPECTED_WITHIN_10KB_PAIRS = 155
EXPECTED_ADJACENT_PAIRS = 44


## ================================================================== ##
## Utilities
## ================================================================== ##

def fail(message):

    print(
        f"\nERROR: {message}",
        file=sys.stderr,
    )

    sys.exit(1)


def read_tsv(path):

    if not path.exists():

        fail(
            f"Required file does not exist:\n"
            f"{path}"
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
            f"{source} missing column(s): "
            +
            ", ".join(
                sorted(missing)
            )
        )


def numeric(
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


def support(
    numerator,
    denominator,
):

    return (
        f"{int(numerator)}"
        f"/"
        f"{int(denominator)}"
    )


def distribution_string(
    counter,
):

    if len(
        counter
    ) == 0:

        return ""


    ordered = sorted(
        counter.items(),
        key=lambda x: (
            -x[1],
            x[0],
        ),
    )


    return "; ".join(
        f"{label}:{count}"

        for label, count
        in ordered
    )


def dominant_counter_value(
    counter,
):

    if len(
        counter
    ) == 0:

        return (
            "",
            0,
            0,
        )


    maximum = max(
        counter.values()
    )


    winners = sorted(
        label

        for label, count
        in counter.items()

        if count == maximum
    )


    if len(
        winners
    ) == 1:

        return (
            winners[0],
            maximum,
            1,
        )


    return (
        "mixed_tie",
        maximum,
        len(
            winners
        ),
    )


def canonical_pair(
    cluster_x,
    cluster_y,
):

    if cluster_x <= cluster_y:

        return (
            cluster_x,
            cluster_y,
        )

    return (
        cluster_y,
        cluster_x,
    )


def broad_orientation(
    value,
):

    value = str(
        value
    ).strip()


    if value.startswith(
        "codirectional"
    ):

        return "codirectional"


    if value == "convergent":

        return "convergent"


    if value == "divergent":

        return "divergent"


    return ""


def transcriptional_order(
    left_cluster,
    right_cluster,
    left_strand,
    right_strand,
):

    ## -------------------------------------------------------------- ##
    ## Only meaningful for codirectional pairs.
    ##
    ## genomic:
    ##
    ##   +  left -> right
    ##   -  right -> left
    ##
    ## -------------------------------------------------------------- ##

    if (
        left_strand == "+"
        and
        right_strand == "+"
    ):

        return (
            f"{left_cluster}"
            f"->{right_cluster}"
        )


    if (
        left_strand == "-"
        and
        right_strand == "-"
    ):

        return (
            f"{right_cluster}"
            f"->{left_cluster}"
        )


    return ""


def write_tsv(
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
## 1. Module membership
## ================================================================== ##

print("=" * 80)
print("STAGE 14B - SUPPORTED MODULE ARCHITECTURE")
print("=" * 80)

print()
print("Reading authoritative module membership...")


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


if len(
    membership
) != EXPECTED_CLUSTERS:

    fail(
        f"Expected "
        f"{EXPECTED_CLUSTERS} "
        f"module clusters; found "
        f"{len(membership)}."
    )


if membership[
    "module"
].nunique() != EXPECTED_MODULES:

    fail(
        "Unexpected module count."
    )


if membership[
    "cluster"
].duplicated().any():

    fail(
        "A cluster belongs to multiple modules."
    )


cluster_to_module = dict(
    zip(
        membership[
            "cluster"
        ],
        membership[
            "module"
        ],
    )
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
## 2. Generate complete possible within-module pair universe
## ================================================================== ##

print()
print("Generating all possible within-module cluster pairs...")


pair_rows = []


for module in modules:

    clusters = sorted(
        membership.loc[
            membership[
                "module"
            ]
            ==
            module,
            "cluster",
        ]
    )


    for cluster_a, cluster_b in combinations(
        clusters,
        2,
    ):

        pair_rows.append(
            {
                "module":
                    module,

                "cluster_a":
                    cluster_a,

                "cluster_b":
                    cluster_b,
            }
        )


all_pairs = pd.DataFrame(
    pair_rows
)


if len(
    all_pairs
) != EXPECTED_POSSIBLE_PAIRS:

    fail(
        f"Expected "
        f"{EXPECTED_POSSIBLE_PAIRS} "
        f"possible within-module pairs; found "
        f"{len(all_pairs)}."
    )


print(
    f"  Possible module-cluster pairs: "
    f"{len(all_pairs):,}"
)


## ================================================================== ##
## 3. Read Stage-12 pair observations
## ================================================================== ##

print()
print("Reading genome-level module-pair observations...")


obs = read_tsv(
    PAIR_OBSERVATIONS
)


require_columns(
    obs,
    [
        "genome",
        "module",
        "cluster_a",
        "cluster_b",

        "any_same_contig",
        "any_adjacent",
        "any_within_5kb",
        "any_within_10kb",
        "any_within_20kb",

        "closest_contig",
        "closest_intergenic_gap_bp",
        "closest_genes_between",

        "closest_left_cluster",
        "closest_right_cluster",
        "closest_left_strand",
        "closest_right_strand",
        "closest_orientation_class",
    ],
    PAIR_OBSERVATIONS.name,
)


if len(
    obs
) != EXPECTED_PAIR_OBSERVATIONS:

    fail(
        f"Expected "
        f"{EXPECTED_PAIR_OBSERVATIONS:,} "
        f"pair observations; found "
        f"{len(obs):,}."
    )


## Canonicalize cluster_a / cluster_b defensively. ##

canonical = obs.apply(
    lambda row:
        canonical_pair(
            row[
                "cluster_a"
            ],
            row[
                "cluster_b"
            ],
        ),
    axis=1,
)


obs[
    "cluster_a"
] = [
    value[0]

    for value
    in canonical
]


obs[
    "cluster_b"
] = [
    value[1]

    for value
    in canonical
]


## One genome x module x unordered family pair. ##

if obs[
    [
        "genome",
        "module",
        "cluster_a",
        "cluster_b",
    ]
].duplicated().any():

    fail(
        "Duplicate genome/module/cluster-pair "
        "observations detected."
    )


## Validate module identities. ##

expected_a = obs[
    "cluster_a"
].map(
    cluster_to_module
)


expected_b = obs[
    "cluster_b"
].map(
    cluster_to_module
)


if (
    expected_a
    !=
    obs[
        "module"
    ]
).any():

    fail(
        "cluster_a module identity mismatch."
    )


if (
    expected_b
    !=
    obs[
        "module"
    ]
).any():

    fail(
        "cluster_b module identity mismatch."
    )


flag_columns = [
    "any_same_contig",
    "any_adjacent",
    "any_within_5kb",
    "any_within_10kb",
    "any_within_20kb",
]


for column in flag_columns:

    obs[
        column
    ] = pd.to_numeric(
        obs[
            column
        ],
        errors="raise",
    ).astype(int)


obs[
    "closest_intergenic_gap_bp"
] = numeric(
    obs[
        "closest_intergenic_gap_bp"
    ]
)


obs[
    "closest_genes_between"
] = numeric(
    obs[
        "closest_genes_between"
    ]
)


## Physical hierarchy QC. ##

if (
    obs[
        "any_adjacent"
    ]
    >
    obs[
        "any_within_5kb"
    ]
).any():

    fail(
        "Adjacency occurs without <=5-kb flag."
    )


if (
    obs[
        "any_within_5kb"
    ]
    >
    obs[
        "any_within_10kb"
    ]
).any():

    fail(
        "<=5-kb flag occurs without <=10-kb flag."
    )


if (
    obs[
        "any_within_10kb"
    ]
    >
    obs[
        "any_within_20kb"
    ]
).any():

    fail(
        "<=10-kb flag occurs without <=20-kb flag."
    )


if (
    obs[
        "any_within_20kb"
    ]
    >
    obs[
        "any_same_contig"
    ]
).any():

    fail(
        "<=20-kb flag occurs without same-contig flag."
    )


## ================================================================== ##
## 4. Genome-level pair statistics
## ================================================================== ##

print()
print("Summarizing physical association evidence...")


PAIR_KEYS = [
    "module",
    "cluster_a",
    "cluster_b",
]


pair_stats = (
    obs
    .groupby(
        PAIR_KEYS,
        as_index=False,
    )
    .agg(
        n_genomes_both=(
            "genome",
            "nunique",
        ),

        n_genomes_same_contig=(
            "any_same_contig",
            "sum",
        ),

        n_genomes_adjacent=(
            "any_adjacent",
            "sum",
        ),

        n_genomes_within_5kb=(
            "any_within_5kb",
            "sum",
        ),

        n_genomes_within_10kb=(
            "any_within_10kb",
            "sum",
        ),

        n_genomes_within_20kb=(
            "any_within_20kb",
            "sum",
        ),
    )
)


edges = all_pairs.merge(
    pair_stats,
    on=PAIR_KEYS,
    how="left",
    validate="one_to_one",
)


count_columns = [
    "n_genomes_both",
    "n_genomes_same_contig",
    "n_genomes_adjacent",
    "n_genomes_within_5kb",
    "n_genomes_within_10kb",
    "n_genomes_within_20kb",
]


for column in count_columns:

    edges[
        column
    ] = (
        edges[
            column
        ]
        .fillna(
            0
        )
        .astype(int)
    )


## Support strings and percentages. ##

for label, numerator_column in [
    (
        "same_contig",
        "n_genomes_same_contig",
    ),
    (
        "adjacent",
        "n_genomes_adjacent",
    ),
    (
        "within_5kb",
        "n_genomes_within_5kb",
    ),
    (
        "within_10kb",
        "n_genomes_within_10kb",
    ),
    (
        "within_20kb",
        "n_genomes_within_20kb",
    ),
]:

    edges[
        f"{label}_support"
    ] = edges.apply(
        lambda row:
            support(
                row[
                    numerator_column
                ],
                row[
                    "n_genomes_both"
                ],
            )
            if row[
                "n_genomes_both"
            ] > 0
            else
            "0/0",
        axis=1,
    )


    edges[
        f"pct_{label}_of_copresent"
    ] = edges.apply(
        lambda row:
            pct(
                row[
                    numerator_column
                ],
                row[
                    "n_genomes_both"
                ],
            ),
        axis=1,
    )


## ================================================================== ##
## 5. Tightest observed physical scale
##
## This is descriptive only.
## It is NOT a significance or conservation threshold.
## ================================================================== ##

def tightest_scale(row):

    if row[
        "n_genomes_adjacent"
    ] > 0:

        return "adjacent"


    if row[
        "n_genomes_within_5kb"
    ] > 0:

        return "<=5kb"


    if row[
        "n_genomes_within_10kb"
    ] > 0:

        return "<=10kb"


    if row[
        "n_genomes_within_20kb"
    ] > 0:

        return "<=20kb"


    if row[
        "n_genomes_same_contig"
    ] > 0:

        return "same_contig_only"


    if row[
        "n_genomes_both"
    ] > 0:

        return "copresence_only"


    return "never_copresent"


edges[
    "tightest_observed_scale"
] = edges.apply(
    tightest_scale,
    axis=1,
)


edges[
    "has_local_20kb_architecture_edge"
] = (
    edges[
        "n_genomes_within_20kb"
    ]
    >
    0
).astype(int)


## ================================================================== ##
## 6. Local orientation and transcriptional order
##
## Only <=20-kb observations contribute to module architecture.
## ================================================================== ##

print()
print("Summarizing local pair orientation and gene order...")


local_obs = obs[
    obs[
        "any_within_20kb"
    ]
    ==
    1
].copy()


local_obs[
    "broad_orientation"
] = local_obs[
    "closest_orientation_class"
].map(
    broad_orientation
)


local_obs[
    "transcriptional_order"
] = local_obs.apply(
    lambda row:
        transcriptional_order(
            row[
                "closest_left_cluster"
            ],
            row[
                "closest_right_cluster"
            ],
            row[
                "closest_left_strand"
            ],
            row[
                "closest_right_strand"
            ],
        )
        if row[
            "broad_orientation"
        ]
        ==
        "codirectional"
        else
        "",
    axis=1,
)


orientation_rows = []


for (
    module,
    cluster_a,
    cluster_b,
), group in local_obs.groupby(
    PAIR_KEYS,
    sort=True,
):

    orientation_counter = Counter(
        value

        for value
        in group[
            "broad_orientation"
        ]

        if value != ""
    )


    (
        dominant_orientation,
        dominant_orientation_n,
        dominant_orientation_ties,
    ) = dominant_counter_value(
        orientation_counter
    )


    order_counter = Counter(
        value

        for value
        in group[
            "transcriptional_order"
        ]

        if value != ""
    )


    (
        dominant_order,
        dominant_order_n,
        dominant_order_ties,
    ) = dominant_counter_value(
        order_counter
    )


    n_local = len(
        group
    )


    n_codirectional = int(
        (
            group[
                "broad_orientation"
            ]
            ==
            "codirectional"
        )
        .sum()
    )


    n_convergent = int(
        (
            group[
                "broad_orientation"
            ]
            ==
            "convergent"
        )
        .sum()
    )


    n_divergent = int(
        (
            group[
                "broad_orientation"
            ]
            ==
            "divergent"
        )
        .sum()
    )


    gaps = group[
        "closest_intergenic_gap_bp"
    ].dropna()


    genes_between = group[
        "closest_genes_between"
    ].dropna()


    orientation_rows.append(
        {
            "module":
                module,

            "cluster_a":
                cluster_a,

            "cluster_b":
                cluster_b,

            "n_local_orientation_observations":
                n_local,

            "n_local_codirectional":
                n_codirectional,

            "n_local_convergent":
                n_convergent,

            "n_local_divergent":
                n_divergent,

            "local_orientation_distribution":
                distribution_string(
                    orientation_counter
                ),

            "dominant_local_orientation":
                dominant_orientation,

            "n_dominant_local_orientation":
                dominant_orientation_n,

            "n_dominant_local_orientation_ties":
                dominant_orientation_ties,

            "dominant_local_orientation_support":
                support(
                    dominant_orientation_n,
                    n_local,
                ),

            "pct_dominant_local_orientation":
                pct(
                    dominant_orientation_n,
                    n_local,
                ),

            "local_transcriptional_order_distribution":
                distribution_string(
                    order_counter
                ),

            "dominant_transcriptional_order":
                dominant_order,

            "n_dominant_transcriptional_order":
                dominant_order_n,

            "n_dominant_transcriptional_order_ties":
                dominant_order_ties,

            "dominant_transcriptional_order_support":
                (
                    support(
                        dominant_order_n,
                        n_codirectional,
                    )
                    if n_codirectional > 0
                    else
                    "0/0"
                ),

            "pct_dominant_transcriptional_order_among_codirectional":
                pct(
                    dominant_order_n,
                    n_codirectional,
                ),

            "median_local_intergenic_gap_bp":
                (
                    gaps.median()
                    if len(
                        gaps
                    ) > 0
                    else
                    np.nan
                ),

            "median_local_genes_between":
                (
                    genes_between.median()
                    if len(
                        genes_between
                    ) > 0
                    else
                    np.nan
                ),
        }
    )


orientation = pd.DataFrame(
    orientation_rows
)


edges = edges.merge(
    orientation,
    on=PAIR_KEYS,
    how="left",
    validate="one_to_one",
)


orientation_count_columns = [
    "n_local_orientation_observations",
    "n_local_codirectional",
    "n_local_convergent",
    "n_local_divergent",
    "n_dominant_local_orientation",
    "n_dominant_local_orientation_ties",
    "n_dominant_transcriptional_order",
    "n_dominant_transcriptional_order_ties",
]


for column in orientation_count_columns:

    edges[
        column
    ] = (
        edges[
            column
        ]
        .fillna(
            0
        )
        .astype(int)
    )


## ================================================================== ##
## 7. Stage-13D ordered focal consensus
## ================================================================== ##

print()
print("Reading focal-family consensus directions...")


mmseq = read_tsv(
    MMSEQ_CONSENSUS
)


require_columns(
    mmseq,
    [
        "focal_cluster",
        "focal_module",
        "neighbor_cluster",
        "direction",

        "direction_has_gene_window_hits",

        "occurrence_saturation_gene_radius",
        "occurrence_saturation_n_positive",
        "occurrence_saturation_n_informative",
        "occurrence_saturation_support",
        "occurrence_saturation_pct",

        "genome_saturation_gene_radius",
        "genome_saturation_n_positive",
        "genome_saturation_n_informative",
        "genome_saturation_support",
        "genome_saturation_pct",

        "best_exact_gene_offset",
        "best_exact_occurrence_support",
        "best_exact_pct_occurrences",
    ],
    MMSEQ_CONSENSUS.name,
)


numeric_consensus_columns = [
    "direction_has_gene_window_hits",

    "occurrence_saturation_gene_radius",
    "occurrence_saturation_n_positive",
    "occurrence_saturation_n_informative",
    "occurrence_saturation_pct",

    "genome_saturation_gene_radius",
    "genome_saturation_n_positive",
    "genome_saturation_n_informative",
    "genome_saturation_pct",

    "best_exact_gene_offset",
    "best_exact_pct_occurrences",
]


for column in numeric_consensus_columns:

    mmseq[
        column
    ] = numeric(
        mmseq[
            column
        ]
    )


## Keep only relationships between distinct member families of the
## same MCL module. ##

mmseq[
    "expected_neighbor_module"
] = mmseq[
    "neighbor_cluster"
].map(
    cluster_to_module
)


mmseq_module = mmseq[
    (
        mmseq[
            "expected_neighbor_module"
        ]
        ==
        mmseq[
            "focal_module"
        ]
    )
    &
    (
        mmseq[
            "neighbor_cluster"
        ]
        !=
        mmseq[
            "focal_cluster"
        ]
    )
].copy()


## ================================================================== ##
## 8. Collapse upstream/downstream rows to one ordered-family summary
##
## Primary direction is selected by:
##
##   1. larger number of positive genomes
##   2. larger number of positive focal occurrences
##
## If both directions contain hits, that fact remains explicit.
## ================================================================== ##

ordered_rows = []


for (
    focal_module,
    focal_cluster,
    neighbor_cluster,
), group in mmseq_module.groupby(
    [
        "focal_module",
        "focal_cluster",
        "neighbor_cluster",
    ],
    sort=True,
):

    if group[
        "direction"
    ].duplicated().any():

        fail(
            "Duplicate direction row in MMseq consensus."
        )


    by_direction = {
        row.direction:
            row

        for row
        in group.itertuples(
            index=False
        )
    }


    upstream = by_direction.get(
        "upstream"
    )


    downstream = by_direction.get(
        "downstream"
    )


    def get_value(
        row,
        name,
        default=0,
    ):

        if row is None:

            return default


        value = getattr(
            row,
            name,
        )


        if pd.isna(
            value
        ):

            return default


        return value


    up_hit = int(
        get_value(
            upstream,
            "direction_has_gene_window_hits",
            0,
        )
        ==
        1
    )


    down_hit = int(
        get_value(
            downstream,
            "direction_has_gene_window_hits",
            0,
        )
        ==
        1
    )


    up_genomes = int(
        get_value(
            upstream,
            "genome_saturation_n_positive",
            0,
        )
    )


    down_genomes = int(
        get_value(
            downstream,
            "genome_saturation_n_positive",
            0,
        )
    )


    up_occ = int(
        get_value(
            upstream,
            "occurrence_saturation_n_positive",
            0,
        )
    )


    down_occ = int(
        get_value(
            downstream,
            "occurrence_saturation_n_positive",
            0,
        )
    )


    directions_with_hits = (
        up_hit
        +
        down_hit
    )


    if directions_with_hits == 0:

        primary_direction = ""


    elif (
        up_genomes,
        up_occ,
    ) > (
        down_genomes,
        down_occ,
    ):

        primary_direction = "upstream"


    elif (
        down_genomes,
        down_occ,
    ) > (
        up_genomes,
        up_occ,
    ):

        primary_direction = "downstream"


    else:

        primary_direction = "mixed_tie"


    if primary_direction == "upstream":

        primary = upstream


    elif primary_direction == "downstream":

        primary = downstream


    else:

        primary = None


    ordered_rows.append(
        {
            "module":
                focal_module,

            "focal_cluster":
                focal_cluster,

            "neighbor_cluster":
                neighbor_cluster,

            "n_consensus_directions_with_hits":
                directions_with_hits,

            "both_consensus_directions_have_hits":
                int(
                    directions_with_hits
                    ==
                    2
                ),

            "upstream_positive_genomes":
                up_genomes,

            "downstream_positive_genomes":
                down_genomes,

            "upstream_positive_occurrences":
                up_occ,

            "downstream_positive_occurrences":
                down_occ,

            "primary_consensus_direction":
                primary_direction,

            "primary_genome_saturation_radius":
                (
                    get_value(
                        primary,
                        "genome_saturation_gene_radius",
                        np.nan,
                    )
                    if primary is not None
                    else
                    np.nan
                ),

            "primary_genome_saturation_support":
                (
                    get_value(
                        primary,
                        "genome_saturation_support",
                        "",
                    )
                    if primary is not None
                    else
                    ""
                ),

            "primary_genome_saturation_pct":
                (
                    get_value(
                        primary,
                        "genome_saturation_pct",
                        np.nan,
                    )
                    if primary is not None
                    else
                    np.nan
                ),

            "primary_occurrence_saturation_radius":
                (
                    get_value(
                        primary,
                        "occurrence_saturation_gene_radius",
                        np.nan,
                    )
                    if primary is not None
                    else
                    np.nan
                ),

            "primary_occurrence_saturation_support":
                (
                    get_value(
                        primary,
                        "occurrence_saturation_support",
                        "",
                    )
                    if primary is not None
                    else
                    ""
                ),

            "primary_occurrence_saturation_pct":
                (
                    get_value(
                        primary,
                        "occurrence_saturation_pct",
                        np.nan,
                    )
                    if primary is not None
                    else
                    np.nan
                ),

            "primary_best_exact_gene_offset":
                (
                    get_value(
                        primary,
                        "best_exact_gene_offset",
                        np.nan,
                    )
                    if primary is not None
                    else
                    np.nan
                ),

            "primary_best_exact_occurrence_support":
                (
                    get_value(
                        primary,
                        "best_exact_occurrence_support",
                        "",
                    )
                    if primary is not None
                    else
                    ""
                ),

            "primary_best_exact_pct_occurrences":
                (
                    get_value(
                        primary,
                        "best_exact_pct_occurrences",
                        np.nan,
                    )
                    if primary is not None
                    else
                    np.nan
                ),
        }
    )


ordered = pd.DataFrame(
    ordered_rows
)


## ================================================================== ##
## 9. Attach reciprocal focal perspectives to unordered module edges
## ================================================================== ##

ordered[
    "pair_a"
] = ordered.apply(
    lambda row:
        canonical_pair(
            row[
                "focal_cluster"
            ],
            row[
                "neighbor_cluster"
            ],
        )[0],
    axis=1,
)


ordered[
    "pair_b"
] = ordered.apply(
    lambda row:
        canonical_pair(
            row[
                "focal_cluster"
            ],
            row[
                "neighbor_cluster"
            ],
        )[1],
    axis=1,
)


ordered[
    "perspective"
] = np.where(
    ordered[
        "focal_cluster"
    ]
    ==
    ordered[
        "pair_a"
    ],
    "a",
    "b",
)


ordered_fields = [
    "n_consensus_directions_with_hits",
    "both_consensus_directions_have_hits",

    "upstream_positive_genomes",
    "downstream_positive_genomes",

    "upstream_positive_occurrences",
    "downstream_positive_occurrences",

    "primary_consensus_direction",

    "primary_genome_saturation_radius",
    "primary_genome_saturation_support",
    "primary_genome_saturation_pct",

    "primary_occurrence_saturation_radius",
    "primary_occurrence_saturation_support",
    "primary_occurrence_saturation_pct",

    "primary_best_exact_gene_offset",
    "primary_best_exact_occurrence_support",
    "primary_best_exact_pct_occurrences",
]


for perspective in [
    "a",
    "b",
]:

    tmp = ordered[
        ordered[
            "perspective"
        ]
        ==
        perspective
    ].copy()


    tmp = tmp.rename(
        columns={
            "pair_a":
                "cluster_a",

            "pair_b":
                "cluster_b",

            **{
                column:
                    f"{perspective}_{column}"

                for column
                in ordered_fields
            },
        }
    )


    keep = (
        [
            "module",
            "cluster_a",
            "cluster_b",
        ]
        +
        [
            f"{perspective}_{column}"

            for column
            in ordered_fields
        ]
    )


    edges = edges.merge(
        tmp[
            keep
        ],
        on=PAIR_KEYS,
        how="left",
        validate="one_to_one",
    )


## ================================================================== ##
## 10. Infer reciprocal consensus geometry
##
## For transcription-normalized focal directions:
##
## A sees B downstream + B sees A upstream
##     => codirectional A -> B
##
## A sees B upstream + B sees A downstream
##     => codirectional B -> A
##
## both see each other downstream
##     => convergent
##
## both see each other upstream
##     => divergent
##
## If either ordered relationship has meaningful hits on BOTH sides,
## geometry is conservatively labelled mixed.
## ================================================================== ##

def infer_consensus_geometry(
    row,
):

    a_direction = str(
        row.get(
            "a_primary_consensus_direction",
            "",
        )
    )


    b_direction = str(
        row.get(
            "b_primary_consensus_direction",
            "",
        )
    )


    a_both = row.get(
        "a_both_consensus_directions_have_hits",
        0,
    )


    b_both = row.get(
        "b_both_consensus_directions_have_hits",
        0,
    )


    if pd.isna(
        a_both
    ):

        a_both = 0


    if pd.isna(
        b_both
    ):

        b_both = 0


    if (
        a_direction in [
            "",
            "nan",
        ]
        or
        b_direction in [
            "",
            "nan",
        ]
    ):

        return "incomplete_consensus"


    if (
        a_direction == "mixed_tie"
        or
        b_direction == "mixed_tie"
        or
        int(
            a_both
        )
        ==
        1
        or
        int(
            b_both
        )
        ==
        1
    ):

        return "mixed_directional_context"


    if (
        a_direction == "downstream"
        and
        b_direction == "upstream"
    ):

        return (
            "codirectional:"
            f"{row['cluster_a']}"
            f"->{row['cluster_b']}"
        )


    if (
        a_direction == "upstream"
        and
        b_direction == "downstream"
    ):

        return (
            "codirectional:"
            f"{row['cluster_b']}"
            f"->{row['cluster_a']}"
        )


    if (
        a_direction == "downstream"
        and
        b_direction == "downstream"
    ):

        return "convergent"


    if (
        a_direction == "upstream"
        and
        b_direction == "upstream"
    ):

        return "divergent"


    return "unresolved"


edges[
    "consensus_geometry"
] = edges.apply(
    infer_consensus_geometry,
    axis=1,
)


def consensus_broad_geometry(
    value,
):

    value = str(
        value
    )


    if value.startswith(
        "codirectional:"
    ):

        return "codirectional"


    if value in [
        "convergent",
        "divergent",
    ]:

        return value


    return ""


edges[
    "consensus_broad_orientation"
] = edges[
    "consensus_geometry"
].map(
    consensus_broad_geometry
)


## ================================================================== ##
## 11. Compare independent Stage-12 and Stage-13D geometry evidence
## ================================================================== ##

def geometry_agreement(
    row,
):

    local = str(
        row.get(
            "dominant_local_orientation",
            "",
        )
    )


    consensus = str(
        row.get(
            "consensus_broad_orientation",
            "",
        )
    )


    if (
        local == ""
        or
        local == "mixed_tie"
        or
        consensus == ""
    ):

        return ""


    return int(
        local
        ==
        consensus
    )


edges[
    "local_vs_consensus_orientation_agreement"
] = edges.apply(
    geometry_agreement,
    axis=1,
)


def order_agreement(
    row,
):

    local_order = str(
        row.get(
            "dominant_transcriptional_order",
            "",
        )
    )


    geometry = str(
        row.get(
            "consensus_geometry",
            "",
        )
    )


    if (
        local_order == ""
        or
        local_order == "mixed_tie"
        or
        not geometry.startswith(
            "codirectional:"
        )
    ):

        return ""


    consensus_order = geometry.split(
        ":",
        1,
    )[1]


    return int(
        local_order
        ==
        consensus_order
    )


edges[
    "local_vs_consensus_order_agreement"
] = edges.apply(
    order_agreement,
    axis=1,
)


## ================================================================== ##
## 12. Build local architecture components
##
## An architecture edge exists if the two module families are observed
## within <=20 kb in at least one genome.
##
## No recurrence cutoff is imposed.
## ================================================================== ##

print()
print("Building module-level local architecture components...")


adjacency = defaultdict(
    set
)


for row in edges.itertuples(
    index=False
):

    if row.n_genomes_within_20kb > 0:

        adjacency[
            (
                row.module,
                row.cluster_a,
            )
        ].add(
            row.cluster_b
        )


        adjacency[
            (
                row.module,
                row.cluster_b,
            )
        ].add(
            row.cluster_a
        )


component_rows = []


for module in modules:

    clusters = sorted(
        membership.loc[
            membership[
                "module"
            ]
            ==
            module,
            "cluster",
        ]
    )


    visited = set()
    components = []


    for start in clusters:

        if start in visited:

            continue


        stack = [
            start
        ]


        component = []


        while stack:

            current = stack.pop()


            if current in visited:

                continue


            visited.add(
                current
            )


            component.append(
                current
            )


            neighbors = adjacency.get(
                (
                    module,
                    current,
                ),
                set(),
            )


            for neighbor in sorted(
                neighbors,
                reverse=True,
            ):

                if neighbor not in visited:

                    stack.append(
                        neighbor
                    )


        components.append(
            sorted(
                component
            )
        )


    ## Largest components first, then lexical deterministic tie-break. ##

    components = sorted(
        components,
        key=lambda component: (
            -len(
                component
            ),
            component[0],
        ),
    )


    for index, component in enumerate(
        components,
        start=1,
    ):

        component_id = (
            f"{module}_LC"
            f"{index:02d}"
        )


        for cluster in component:

            degree = len(
                adjacency.get(
                    (
                        module,
                        cluster,
                    ),
                    set(),
                )
            )


            component_rows.append(
                {
                    "module":
                        module,

                    "cluster":
                        cluster,

                    "local_component":
                        component_id,

                    "local_component_size":
                        len(
                            component
                        ),

                    "local_architecture_degree":
                        degree,

                    "has_local_20kb_edge":
                        int(
                            degree
                            >
                            0
                        ),
                }
            )


components = pd.DataFrame(
    component_rows
)


## ================================================================== ##
## 13. Build architecture node table
## ================================================================== ##

print("Building architecture nodes...")


cluster_evidence = read_tsv(
    CLUSTER_EVIDENCE
)


require_columns(
    cluster_evidence,
    [
        "module",
        "cluster",
    ],
    CLUSTER_EVIDENCE.name,
)


if cluster_evidence[
    "cluster"
].duplicated().any():

    fail(
        "Duplicate cluster in Stage-14A cluster evidence."
    )


nodes = membership.merge(
    components,
    on=[
        "module",
        "cluster",
    ],
    how="left",
    validate="one_to_one",
)


nodes = nodes.merge(
    cluster_evidence,
    on=[
        "module",
        "cluster",
    ],
    how="left",
    validate="one_to_one",
)


## ================================================================== ##
## 14. Module-level architecture summary
## ================================================================== ##

print("Building module architecture summaries...")


summary_rows = []


for module in modules:

    module_nodes = nodes[
        nodes[
            "module"
        ]
        ==
        module
    ]


    module_edges = edges[
        edges[
            "module"
        ]
        ==
        module
    ]


    component_sizes = (
        module_nodes[
            [
                "local_component",
                "local_component_size",
            ]
        ]
        .drop_duplicates()
    )


    summary_rows.append(
        {
            "module":
                module,

            "n_member_clusters":
                len(
                    module_nodes
                ),

            "n_possible_cluster_pairs":
                len(
                    module_edges
                ),

            "n_pairs_copresent":
                int(
                    (
                        module_edges[
                            "n_genomes_both"
                        ]
                        >
                        0
                    )
                    .sum()
                ),

            "n_pairs_same_contig":
                int(
                    (
                        module_edges[
                            "n_genomes_same_contig"
                        ]
                        >
                        0
                    )
                    .sum()
                ),

            "n_pairs_within_20kb":
                int(
                    (
                        module_edges[
                            "n_genomes_within_20kb"
                        ]
                        >
                        0
                    )
                    .sum()
                ),

            "n_pairs_within_10kb":
                int(
                    (
                        module_edges[
                            "n_genomes_within_10kb"
                        ]
                        >
                        0
                    )
                    .sum()
                ),

            "n_pairs_within_5kb":
                int(
                    (
                        module_edges[
                            "n_genomes_within_5kb"
                        ]
                        >
                        0
                    )
                    .sum()
                ),

            "n_pairs_with_adjacency":
                int(
                    (
                        module_edges[
                            "n_genomes_adjacent"
                        ]
                        >
                        0
                    )
                    .sum()
                ),

            "n_genome_pair_observations":
                int(
                    module_edges[
                        "n_genomes_both"
                    ]
                    .sum()
                ),

            "n_genome_pair_observations_within_20kb":
                int(
                    module_edges[
                        "n_genomes_within_20kb"
                    ]
                    .sum()
                ),

            "n_local_architecture_components":
                len(
                    component_sizes
                ),

            "largest_local_component_size":
                int(
                    component_sizes[
                        "local_component_size"
                    ]
                    .max()
                ),

            "n_clusters_with_local_edge":
                int(
                    module_nodes[
                        "has_local_20kb_edge"
                    ]
                    .sum()
                ),

            "n_clusters_without_local_edge":
                int(
                    (
                        module_nodes[
                            "has_local_20kb_edge"
                        ]
                        ==
                        0
                    )
                    .sum()
                ),

            "n_edges_with_resolved_consensus_geometry":
                int(
                    module_edges[
                        "consensus_geometry"
                    ]
                    .str.startswith(
                        "codirectional:"
                    )
                    .sum()
                    +
                    module_edges[
                        "consensus_geometry"
                    ]
                    .isin(
                        [
                            "convergent",
                            "divergent",
                        ]
                    )
                    .sum()
                ),

            "n_edges_local_consensus_orientation_agree":
                int(
                    (
                        module_edges[
                            "local_vs_consensus_orientation_agreement"
                        ]
                        ==
                        1
                    )
                    .sum()
                ),

            "n_edges_local_consensus_orientation_disagree":
                int(
                    (
                        module_edges[
                            "local_vs_consensus_orientation_agreement"
                        ]
                        ==
                        0
                    )
                    .sum()
                ),
        }
    )


summary = pd.DataFrame(
    summary_rows
)


## ================================================================== ##
## 15. Global QC
## ================================================================== ##

print()
print("Running architecture QC...")


checks = {
    "possible_pairs":
        (
            len(
                edges
            ),
            EXPECTED_POSSIBLE_PAIRS,
        ),

    "pair_observations":
        (
            int(
                edges[
                    "n_genomes_both"
                ]
                .sum()
            ),
            EXPECTED_PAIR_OBSERVATIONS,
        ),

    "same_contig_observations":
        (
            int(
                edges[
                    "n_genomes_same_contig"
                ]
                .sum()
            ),
            EXPECTED_SAME_CONTIG_OBSERVATIONS,
        ),

    "within_5kb_observations":
        (
            int(
                edges[
                    "n_genomes_within_5kb"
                ]
                .sum()
            ),
            EXPECTED_WITHIN_5KB_OBSERVATIONS,
        ),

    "within_10kb_observations":
        (
            int(
                edges[
                    "n_genomes_within_10kb"
                ]
                .sum()
            ),
            EXPECTED_WITHIN_10KB_OBSERVATIONS,
        ),

    "within_20kb_observations":
        (
            int(
                edges[
                    "n_genomes_within_20kb"
                ]
                .sum()
            ),
            EXPECTED_WITHIN_20KB_OBSERVATIONS,
        ),

    "adjacent_observations":
        (
            int(
                edges[
                    "n_genomes_adjacent"
                ]
                .sum()
            ),
            EXPECTED_ADJACENT_OBSERVATIONS,
        ),

    "copresent_pairs":
        (
            int(
                (
                    edges[
                        "n_genomes_both"
                    ]
                    >
                    0
                )
                .sum()
            ),
            EXPECTED_COPRESENT_PAIRS,
        ),

    "same_contig_pairs":
        (
            int(
                (
                    edges[
                        "n_genomes_same_contig"
                    ]
                    >
                    0
                )
                .sum()
            ),
            EXPECTED_SAME_CONTIG_PAIRS,
        ),

    "within_10kb_pairs":
        (
            int(
                (
                    edges[
                        "n_genomes_within_10kb"
                    ]
                    >
                    0
                )
                .sum()
            ),
            EXPECTED_WITHIN_10KB_PAIRS,
        ),

    "adjacent_pairs":
        (
            int(
                (
                    edges[
                        "n_genomes_adjacent"
                    ]
                    >
                    0
                )
                .sum()
            ),
            EXPECTED_ADJACENT_PAIRS,
        ),
}


for label, (
    observed,
    expected,
) in checks.items():

    if observed != expected:

        fail(
            f"QC failed for {label}: "
            f"observed {observed:,}, "
            f"expected {expected:,}."
        )


if len(
    nodes
) != EXPECTED_CLUSTERS:

    fail(
        "Architecture node count is not 156."
    )


if nodes[
    "local_component"
].eq(
    ""
).any():

    fail(
        "A node lacks a local-component assignment."
    )


## ================================================================== ##
## 16. Sort
## ================================================================== ##

nodes = (
    nodes
    .sort_values(
        [
            "module",
            "local_component",
            "cluster",
        ],
        kind="stable",
    )
    .reset_index(
        drop=True
    )
)


edges = (
    edges
    .sort_values(
        [
            "module",
            "cluster_a",
            "cluster_b",
        ],
        kind="stable",
    )
    .reset_index(
        drop=True
    )
)


summary = (
    summary
    .sort_values(
        "module",
        kind="stable",
    )
    .reset_index(
        drop=True
    )
)


## ================================================================== ##
## 17. Write global outputs
## ================================================================== ##

print("Writing global architecture tables...")


write_tsv(
    nodes,
    OUT_NODES,
)


write_tsv(
    edges,
    OUT_EDGES,
)


write_tsv(
    summary,
    OUT_SUMMARY,
)


## ================================================================== ##
## 18. Write per-module architecture packages
## ================================================================== ##

print("Writing per-module architecture packages...")


for module in modules:

    outdir = (
        MODULE_DIR
        / module
    )


    outdir.mkdir(
        parents=True,
        exist_ok=True,
    )


    module_nodes = nodes[
        nodes[
            "module"
        ]
        ==
        module
    ].copy()


    module_edges = edges[
        edges[
            "module"
        ]
        ==
        module
    ].copy()


    module_summary = summary[
        summary[
            "module"
        ]
        ==
        module
    ].copy()


    write_tsv(
        module_nodes,
        outdir
        /
        "architecture_nodes.tsv",
    )


    write_tsv(
        module_edges,
        outdir
        /
        "architecture_edges_all_pairs.tsv",
    )


    write_tsv(
        module_edges[
            module_edges[
                "has_local_20kb_architecture_edge"
            ]
            ==
            1
        ],
        outdir
        /
        "architecture_edges_local_20kb.tsv",
    )


    write_tsv(
        module_summary,
        outdir
        /
        "architecture_summary.tsv",
    )


## ================================================================== ##
## 19. QC table
## ================================================================== ##

qc_rows = []


for label, (
    observed,
    expected,
) in checks.items():

    qc_rows.append(
        {
            "metric":
                label,

            "observed":
                observed,

            "expected":
                expected,

            "pass":
                int(
                    observed
                    ==
                    expected
                ),
        }
    )


qc_rows.extend(
    [
        {
            "metric":
                "architecture_nodes",

            "observed":
                len(
                    nodes
                ),

            "expected":
                EXPECTED_CLUSTERS,

            "pass":
                int(
                    len(
                        nodes
                    )
                    ==
                    EXPECTED_CLUSTERS
                ),
        },

        {
            "metric":
                "modules",

            "observed":
                len(
                    summary
                ),

            "expected":
                EXPECTED_MODULES,

            "pass":
                int(
                    len(
                        summary
                    )
                    ==
                    EXPECTED_MODULES
                ),
        },
    ]
)


qc = pd.DataFrame(
    qc_rows
)


write_tsv(
    qc,
    OUT_QC,
)


## ================================================================== ##
## 20. Terminal summary
## ================================================================== ##

print()
print("Module architecture summary")


print(
    f"  Modules:                         "
    f"{len(summary):,}"
)

print(
    f"  Member clusters / graph nodes:   "
    f"{len(nodes):,}"
)

print(
    f"  Possible within-module pairs:    "
    f"{len(edges):,}"
)

print(
    f"  Co-present pair families:        "
    f"{(edges['n_genomes_both'] > 0).sum():,}"
)

print(
    f"  Same-contig pair families:       "
    f"{(edges['n_genomes_same_contig'] > 0).sum():,}"
)

print(
    f"  <=20-kb architecture edges:      "
    f"{(edges['n_genomes_within_20kb'] > 0).sum():,}"
)

print(
    f"  <=10-kb pair families:           "
    f"{(edges['n_genomes_within_10kb'] > 0).sum():,}"
)

print(
    f"  <=5-kb pair families:            "
    f"{(edges['n_genomes_within_5kb'] > 0).sum():,}"
)

print(
    f"  Pair families with adjacency:    "
    f"{(edges['n_genomes_adjacent'] > 0).sum():,}"
)

print(
    f"  Genome-pair <=20-kb observations:"
    f" {edges['n_genomes_within_20kb'].sum():,}"
)


## ================================================================== ##
## 21. Module 20 / 35 spotlights
## ================================================================== ##

for spotlight_module in [
    "Module_20",
    "Module_35",
]:

    spot = edges[
        edges[
            "module"
        ]
        ==
        spotlight_module
    ].copy()


    print()
    print(
        f"{spotlight_module} architecture spotlight"
    )


    display = [
        "cluster_a",
        "cluster_b",

        "n_genomes_both",

        "within_20kb_support",
        "within_5kb_support",
        "adjacent_support",

        "dominant_local_orientation",
        "dominant_local_orientation_support",

        "dominant_transcriptional_order",
        "dominant_transcriptional_order_support",

        "consensus_geometry",

        "local_vs_consensus_orientation_agreement",
        "local_vs_consensus_order_agreement",

        "median_local_intergenic_gap_bp",
        "median_local_genes_between",
    ]


    print(
        spot[
            display
        ]
        .to_string(
            index=False
        )
    )


    spot_nodes = nodes[
        nodes[
            "module"
        ]
        ==
        spotlight_module
    ][
        [
            "cluster",
            "local_component",
            "local_component_size",
            "local_architecture_degree",
        ]
    ]


    print()
    print(
        spot_nodes.to_string(
            index=False
        )
    )


print()
print("=" * 80)
print("SUCCESS")
print("=" * 80)

print(
    f"Architecture nodes:   "
    f"{OUT_NODES}"
)

print(
    f"Architecture edges:   "
    f"{OUT_EDGES}"
)

print(
    f"Module summaries:     "
    f"{OUT_SUMMARY}"
)

print(
    f"QC:                   "
    f"{OUT_QC}"
)

print(
    f"Per-module packages:  "
    f"{MODULE_DIR}"
)
