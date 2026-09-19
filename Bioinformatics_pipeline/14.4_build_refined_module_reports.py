#!/usr/bin/env python3

from pathlib import Path
from collections import Counter
import sys

import numpy as np
import pandas as pd


## ================================================================== ##
## STAGE 14D - REFINED, GENOME-RESOLVED MODULE REPORTS
##
## PURPOSE
## -------
##
## Refine the Stage-14C reports by adding:
##
##   1. concise integrated module interpretation
##   2. genome-resolved module architecture
##   3. member-family annotation crosswalk
##   4. original Jaccard-network evidence
##   5. non-focal-module functional neighborhood context
##
##
## IMPORTANT ANALYTICAL PRINCIPLES
## -------------------------------
##
## * MMseqs2 families remain the primary sequence-family unit.
##
## * MCL modules remain genome-level co-occurrence communities and are
##   NOT assumed to be operons or single linear gene systems.
##
## * Physical architecture remains based on Stage-14B / Stage-13
##   evidence.
##
## * No new arbitrary conservation threshold is introduced.
##
## * Percentages remain coupled to raw numerator / denominator support.
##
## * GlobDB, FeGenie, FindMeHemes, SignalP and DeepTMHMM are independent
##   evidence layers.
##
## * A GlobDB COG attached to another MEMBER FAMILY is displayed in the
##   member-family annotation crosswalk, not misrepresented as an
##   ordinary neighboring gene.
##
## * Non-module COG context generated here is descriptive only:
##       +/-5 transcription-oriented genes
##       excluding genes belonging to the SAME focal MCL module.
##
##   It therefore does not replace the censoring-aware Stage-13C
##   annotation analysis.
##
##
## REPORT-LEVEL TOPOLOGY
## ---------------------
##
## DeepTMHMM BETA is interpreted as beta-barrel membrane topology.
##
## Upstream localization_class values are preserved and never rewritten.
## ================================================================== ##


## ================================================================== ##
## Paths
## ================================================================== ##

HERE = Path(__file__).resolve().parent
WORKFLOW = HERE.parent.parent


## ------------------------------------------------------------------ ##
## Stage 07
## ------------------------------------------------------------------ ##

NETWORK_DIR = (
    WORKFLOW
    / "07_cooccurrence_network"
)


MODULE_MEMBERSHIP = (
    NETWORK_DIR
    / "module_membership.tsv"
)


## ------------------------------------------------------------------ ##
## Stage 12
## ------------------------------------------------------------------ ##

PAIR_OBSERVATIONS = (
    WORKFLOW
    / "12_neighborhood_extraction"
    / "module_member_pair_observations.tsv"
)


OBSERVED_NEIGHBORHOODS = (
    WORKFLOW
    / "12_neighborhood_extraction"
    / "observed_neighborhood_genes.tsv"
)


## ------------------------------------------------------------------ ##
## Stage 14A
## ------------------------------------------------------------------ ##

REPORT_OVERVIEW = (
    WORKFLOW
    / "14_module_reports"
    / "14A_report_data"
    / "global"
    / "module_report_overview.tsv"
)


MODULE_PROTEINS = (
    WORKFLOW
    / "14_module_reports"
    / "14A_report_data"
    / "global"
    / "module_protein_evidence.tsv"
)


INTERVENING_EVIDENCE = (
    WORKFLOW
    / "14_module_reports"
    / "14A_report_data"
    / "global"
    / "local_intervening_gene_integrated_evidence.tsv"
)


## ------------------------------------------------------------------ ##
## Stage 14B
## ------------------------------------------------------------------ ##

ARCH_NODES = (
    WORKFLOW
    / "14_module_reports"
    / "14B_module_architecture"
    / "global"
    / "module_architecture_nodes.tsv"
)


ARCH_EDGES = (
    WORKFLOW
    / "14_module_reports"
    / "14B_module_architecture"
    / "global"
    / "module_architecture_edges.tsv"
)


ARCH_SUMMARY = (
    WORKFLOW
    / "14_module_reports"
    / "14B_module_architecture"
    / "global"
    / "module_architecture_summary.tsv"
)


## ------------------------------------------------------------------ ##
## Stage 14C
## ------------------------------------------------------------------ ##

CLUSTER_REPORT = (
    WORKFLOW
    / "14_module_reports"
    / "14C_rendered_reports"
    / "global"
    / "module_cluster_report_evidence.tsv"
)


## ================================================================== ##
## Output paths
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


OUT_REFINED_SUMMARY = (
    GLOBAL_DIR
    / "refined_module_summary.tsv"
)


OUT_GENOME_ARCH = (
    GLOBAL_DIR
    / "genome_resolved_module_architecture.tsv"
)


OUT_GENOME_PROTEINS = (
    GLOBAL_DIR
    / "genome_resolved_module_proteins.tsv"
)


OUT_NETWORK_EDGES = (
    GLOBAL_DIR
    / "module_network_edges.tsv"
)


OUT_NONMODULE_COG = (
    GLOBAL_DIR
    / "nonmodule_cog_context_5genes.tsv"
)


OUT_INDEX = (
    GLOBAL_DIR
    / "refined_module_report_index.tsv"
)


OUT_QC = (
    GLOBAL_DIR
    / "refined_module_report_qc.tsv"
)


## ================================================================== ##
## Locked upstream expectations
## ================================================================== ##

EXPECTED_MODULES = 35
EXPECTED_MODULE_CLUSTERS = 156
EXPECTED_MODULE_PROTEINS = 10_537

EXPECTED_GENOME_MODULE_PAIRS = 2_435
EXPECTED_FULL_GENOME_MODULE_PAIRS = 291

EXPECTED_PAIR_OBSERVATIONS = 34_037

EXPECTED_ARCH_PAIRS = 560
EXPECTED_LOCAL_20KB_EDGES = 192

EXPECTED_NETWORK_EDGES = 332

EXPECTED_INTERVENING_ROWS = 6_950
EXPECTED_UNIQUE_INTERVENING_PROTEINS = 3_291

EXPECTED_BETA_PROTEINS = 61


## ================================================================== ##
## Markdown display limits
##
## These affect display only. Full TSV files always retain all rows.
## ================================================================== ##

MAX_NETWORK_EDGES_DISPLAY = 30
MAX_LOCAL_EDGES_DISPLAY = 40

MAX_NONMODULE_COG_PER_FOCAL = 5

MAX_GENOME_ROWS_IN_MARKDOWN = 100

MAX_INTERVENING_ROWS_DISPLAY = 30


## ================================================================== ##
## Utility functions
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
            f"Required file does not exist:\n{path}"
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


def clean(value):

    if pd.isna(
        value
    ):

        return ""

    return str(
        value
    ).strip()


def numeric(value):

    if isinstance(
        value,
        pd.Series,
    ):

        return pd.to_numeric(
            value.replace(
                "",
                np.nan,
            ),
            errors="coerce",
        )


    value = clean(
        value
    )

    if value == "":

        return np.nan


    return pd.to_numeric(
        value,
        errors="coerce",
    )


def integer(value):

    value = numeric(
        value
    )

    if pd.isna(
        value
    ):

        return 0

    return int(
        value
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


def canonical_pair(
    a,
    b,
):

    if a <= b:

        return (
            a,
            b,
        )

    return (
        b,
        a,
    )


def compact_float(
    value,
    digits=1,
):

    value = numeric(
        value
    )

    if pd.isna(
        value
    ):

        return ""

    return f"{value:.{digits}f}"


def distribution_string(
    values,
):

    cleaned = [
        clean(
            value
        )

        for value
        in values

        if clean(
            value
        )
        not in [
            "",
            "NA",
            "nan",
        ]
    ]


    if len(
        cleaned
    ) == 0:

        return ""


    counter = Counter(
        cleaned
    )


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


def dominant_value(
    values,
):

    cleaned = [
        clean(
            value
        )

        for value
        in values

        if clean(
            value
        )
        not in [
            "",
            "NA",
            "nan",
        ]
    ]


    if len(
        cleaned
    ) == 0:

        return ""


    counter = Counter(
        cleaned
    )


    maximum = max(
        counter.values()
    )


    winners = sorted(
        value

        for value, count
        in counter.items()

        if count == maximum
    )


    if len(
        winners
    ) == 1:

        return winners[0]


    return "mixed_tie"


def md_escape(value):

    value = clean(
        value
    )

    return (
        value
        .replace(
            "|",
            "\\|",
        )
        .replace(
            "\n",
            " ",
        )
    )


def markdown_table(
    df,
    columns,
):

    columns = [
        column

        for column
        in columns

        if column
        in df.columns
    ]


    if len(
        df
    ) == 0:

        return "_No rows._\n"


    if len(
        columns
    ) == 0:

        return "_No displayable columns._\n"


    lines = []


    lines.append(
        "| "
        +
        " | ".join(
            md_escape(
                column
            )

            for column
            in columns
        )
        +
        " |"
    )


    lines.append(
        "| "
        +
        " | ".join(
            "---"

            for _ in columns
        )
        +
        " |"
    )


    for row in df[
        columns
    ].itertuples(
        index=False,
        name=None,
    ):

        lines.append(
            "| "
            +
            " | ".join(
                md_escape(
                    value
                )

                for value
                in row
            )
            +
            " |"
        )


    return (
        "\n".join(
            lines
        )
        +
        "\n"
    )


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
## Report-level topology
## ================================================================== ##

def report_topology_class(
    row,
):

    deeptmhmm = clean(
        row.get(
            "evidence_deeptmhmm_class",
            "",
        )
    )


    signalp = clean(
        row.get(
            "evidence_signalp_prediction",
            "",
        )
    )


    n_tm = numeric(
        row.get(
            "evidence_deeptmhmm_n_tm_helices",
            "",
        )
    )


    if deeptmhmm == "BETA":

        return "beta_barrel_membrane"


    if signalp in [
        "LIPO",
        "TATLIPO",
    ]:

        return "exported_lipoprotein"


    if deeptmhmm == "SP+TM":

        if not pd.isna(
            n_tm
        ):

            if n_tm == 1:

                return "exported_single_pass_membrane"


            if n_tm > 1:

                return "exported_multipass_membrane"


        return "exported_membrane_SP_plus_TM"


    if deeptmhmm == "TM":

        if not pd.isna(
            n_tm
        ):

            if n_tm == 1:

                return "single_pass_membrane"


            if n_tm > 1:

                return "multipass_membrane"


        return "membrane"


    if deeptmhmm == "SP":

        return "soluble_exported_periplasmic_like"


    if deeptmhmm == "GLOB":

        return "globular_nonmembrane_like"


    if signalp in [
        "SP",
        "TAT",
    ]:

        return "exported_topology_uncertain"


    return "uncertain"


## ================================================================== ##
## Per-genome transcriptional order
## ================================================================== ##

def transcriptional_order(
    left_cluster,
    right_cluster,
    left_strand,
    right_strand,
):

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


def broad_orientation(
    value,
):

    value = clean(
        value
    )


    if value.startswith(
        "codirectional"
    ):

        return "codirectional"


    if value == "convergent":

        return "convergent"


    if value == "divergent":

        return "divergent"


    return ""


## ================================================================== ##
## NETWORK EDGE DISCOVERY
##
## We deliberately do not hard-code a filename because the Stage-07
## final network filename may differ between workflow iterations.
##
## We search TSV/CSV files for:
##
##   * two plausible node columns
##   * one Jaccard column
##
## and select the candidate yielding exactly the locked 332 final
## network edges among the 156 MCL-assigned network nodes.
## ================================================================== ##

def detect_network_columns(
    columns,
):

    lower_to_original = {
        column.lower():
            column

        for column
        in columns
    }


    pair_candidates = [
        (
            "source",
            "target",
        ),
        (
            "cluster_a",
            "cluster_b",
        ),
        (
            "cluster1",
            "cluster2",
        ),
        (
            "cluster_1",
            "cluster_2",
        ),
        (
            "node1",
            "node2",
        ),
        (
            "node_1",
            "node_2",
        ),
        (
            "from",
            "to",
        ),
    ]


    pair = None


    for left, right in pair_candidates:

        if (
            left
            in lower_to_original
            and
            right
            in lower_to_original
        ):

            pair = (
                lower_to_original[
                    left
                ],
                lower_to_original[
                    right
                ],
            )

            break


    jaccard = None


    for column in columns:

        if "jaccard" in column.lower():

            jaccard = column

            break


    shared = None


    for column in columns:

        lower = column.lower()

        if (
            "shared"
            in lower
            and
            "genome"
            in lower
        ):

            shared = column

            break


    return (
        pair,
        jaccard,
        shared,
    )


def discover_network_edges(
    network_dir,
    module_clusters,
):

    print()
    print("Discovering authoritative final Jaccard network edge table...")


    candidates = []


    files = (
        list(
            network_dir.rglob(
                "*.tsv"
            )
        )
        +
        list(
            network_dir.rglob(
                "*.csv"
            )
        )
    )


    for path in sorted(
        files
    ):

        try:

            sep = (
                "\t"
                if path.suffix.lower() == ".tsv"
                else
                ","
            )


            header = pd.read_csv(
                path,
                sep=sep,
                nrows=0,
            )


            (
                pair,
                jaccard,
                shared,
            ) = detect_network_columns(
                header.columns
            )


            if (
                pair is None
                or
                jaccard is None
            ):

                continue


            df = pd.read_csv(
                path,
                sep=sep,
                dtype=str,
                keep_default_na=False,
            )


            left_col, right_col = pair


            subset = df[
                df[
                    left_col
                ].isin(
                    module_clusters
                )
                &
                df[
                    right_col
                ].isin(
                    module_clusters
                )
            ].copy()


            if len(
                subset
            ) == 0:

                continue


            canon = subset.apply(
                lambda row:
                    canonical_pair(
                        row[
                            left_col
                        ],
                        row[
                            right_col
                        ],
                    ),
                axis=1,
            )


            subset[
                "_a"
            ] = [
                item[0]

                for item
                in canon
            ]


            subset[
                "_b"
            ] = [
                item[1]

                for item
                in canon
            ]


            n_unique = (
                subset[
                    [
                        "_a",
                        "_b",
                    ]
                ]
                .drop_duplicates()
                .shape[0]
            )


            candidates.append(
                {
                    "path":
                        path,

                    "data":
                        df,

                    "left":
                        left_col,

                    "right":
                        right_col,

                    "jaccard":
                        jaccard,

                    "shared":
                        shared,

                    "n_network_node_edges":
                        n_unique,
                }
            )


        except Exception:

            continue


    exact = [
        candidate

        for candidate
        in candidates

        if candidate[
            "n_network_node_edges"
        ]
        ==
        EXPECTED_NETWORK_EDGES
    ]


    if len(
        exact
    ) == 0:

        print()
        print("Candidate network tables found:")


        for candidate in candidates:

            print(
                f"  {candidate['path']} "
                f"-> {candidate['n_network_node_edges']:,} "
                f"edges among MCL nodes"
            )


        fail(
            "Could not uniquely identify a Stage-07 network table "
            "with the expected 332 final edges."
        )


    ## Prefer filenames containing final/network/edge if multiple
    ## equivalent copies exist. ##

    def score(
        candidate,
    ):

        name = str(
            candidate[
                "path"
            ]
        ).lower()


        value = 0


        if "final" in name:

            value += 4


        if "network" in name:

            value += 3


        if "edge" in name:

            value += 2


        if "jaccard" in name:

            value += 1


        return value


    exact = sorted(
        exact,
        key=lambda candidate: (
            -score(
                candidate
            ),
            len(
                str(
                    candidate[
                        "path"
                    ]
                )
            ),
            str(
                candidate[
                    "path"
                ]
            ),
        ),
    )


    chosen = exact[0]


    print(
        f"  Selected: "
        f"{chosen['path']}"
    )


    df = chosen[
        "data"
    ].copy()


    left_col = chosen[
        "left"
    ]


    right_col = chosen[
        "right"
    ]


    df = df[
        df[
            left_col
        ].isin(
            module_clusters
        )
        &
        df[
            right_col
        ].isin(
            module_clusters
        )
    ].copy()


    canonical = df.apply(
        lambda row:
            canonical_pair(
                row[
                    left_col
                ],
                row[
                    right_col
                ],
            ),
        axis=1,
    )


    out = pd.DataFrame(
        {
            "cluster_a":
                [
                    item[0]

                    for item
                    in canonical
                ],

            "cluster_b":
                [
                    item[1]

                    for item
                    in canonical
                ],

            "jaccard":
                pd.to_numeric(
                    df[
                        chosen[
                            "jaccard"
                        ]
                    ],
                    errors="coerce",
                ),
        }
    )


    if chosen[
        "shared"
    ] is not None:

        out[
            "shared_genomes"
        ] = pd.to_numeric(
            df[
                chosen[
                    "shared"
                ]
            ],
            errors="coerce",
        )


    else:

        out[
            "shared_genomes"
        ] = np.nan


    out[
        "network_source_file"
    ] = str(
        chosen[
            "path"
        ]
    )


    out = (
        out
        .drop_duplicates(
            [
                "cluster_a",
                "cluster_b",
            ]
        )
        .reset_index(
            drop=True
        )
    )


    if len(
        out
    ) != EXPECTED_NETWORK_EDGES:

        fail(
            f"Expected {EXPECTED_NETWORK_EDGES} unique final "
            f"network edges, found {len(out)}."
        )


    return out


## ================================================================== ##
## Read core data
## ================================================================== ##

print("=" * 80)
print("STAGE 14D - REFINED, GENOME-RESOLVED MODULE REPORTS")
print("=" * 80)


print()
print("Reading core Stage-14 data...")


membership = read_tsv(
    MODULE_MEMBERSHIP
)


overview = read_tsv(
    REPORT_OVERVIEW
)


proteins = read_tsv(
    MODULE_PROTEINS
)


intervening = read_tsv(
    INTERVENING_EVIDENCE
)


nodes = read_tsv(
    ARCH_NODES
)


edges = read_tsv(
    ARCH_EDGES
)


arch_summary = read_tsv(
    ARCH_SUMMARY
)


cluster_report = read_tsv(
    CLUSTER_REPORT
)


pair_obs = read_tsv(
    PAIR_OBSERVATIONS
)


require_columns(
    membership,
    [
        "cluster",
        "module",
    ],
    MODULE_MEMBERSHIP.name,
)


require_columns(
    proteins,
    [
        "genome",
        "protein_id",
        "cluster",
        "module",
    ],
    MODULE_PROTEINS.name,
)


require_columns(
    edges,
    [
        "module",
        "cluster_a",
        "cluster_b",
        "n_genomes_both",
        "n_genomes_within_20kb",
        "within_20kb_support",
        "within_5kb_support",
        "adjacent_support",
    ],
    ARCH_EDGES.name,
)


require_columns(
    pair_obs,
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
    membership
) != EXPECTED_MODULE_CLUSTERS:

    fail(
        f"Expected {EXPECTED_MODULE_CLUSTERS} MCL clusters; "
        f"found {len(membership)}."
    )


if membership[
    "module"
].nunique() != EXPECTED_MODULES:

    fail(
        "Unexpected module count."
    )


if len(
    proteins
) != EXPECTED_MODULE_PROTEINS:

    fail(
        f"Expected {EXPECTED_MODULE_PROTEINS:,} module proteins; "
        f"found {len(proteins):,}."
    )


if len(
    edges
) != EXPECTED_ARCH_PAIRS:

    fail(
        f"Expected {EXPECTED_ARCH_PAIRS} architecture pairs; "
        f"found {len(edges)}."
    )


if len(
    pair_obs
) != EXPECTED_PAIR_OBSERVATIONS:

    fail(
        f"Expected {EXPECTED_PAIR_OBSERVATIONS:,} genome-pair "
        f"observations; found {len(pair_obs):,}."
    )


if len(
    intervening
) != EXPECTED_INTERVENING_ROWS:

    fail(
        f"Expected {EXPECTED_INTERVENING_ROWS:,} intervening rows; "
        f"found {len(intervening):,}."
    )


modules = sorted(
    membership[
        "module"
    ].unique()
)


module_clusters = set(
    membership[
        "cluster"
    ]
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


module_size = (
    membership
    .groupby(
        "module"
    )
    .size()
    .to_dict()
)


print(
    f"  Modules:                "
    f"{len(modules):,}"
)

print(
    f"  Member families:        "
    f"{len(membership):,}"
)

print(
    f"  Module proteins:        "
    f"{len(proteins):,}"
)

print(
    f"  Genome-pair rows:       "
    f"{len(pair_obs):,}"
)

print(
    f"  Architecture pairs:     "
    f"{len(edges):,}"
)

print(
    f"  Intervening rows:       "
    f"{len(intervening):,}"
)


## ================================================================== ##
## Report-level topology QC
## ================================================================== ##

print()
print("Deriving report-level protein topology...")


protein_records = proteins.to_dict(
    orient="records"
)


proteins[
    "report_topology_class"
] = [
    report_topology_class(
        record
    )

    for record
    in protein_records
]


beta = proteins[
    proteins[
        "evidence_deeptmhmm_class"
    ]
    ==
    "BETA"
]


if len(
    beta
) != EXPECTED_BETA_PROTEINS:

    fail(
        f"Expected {EXPECTED_BETA_PROTEINS} BETA proteins; "
        f"found {len(beta)}."
    )


if (
    beta[
        "report_topology_class"
    ]
    !=
    "beta_barrel_membrane"
).any():

    fail(
        "BETA topology QC failed."
    )


print(
    f"  BETA proteins:          "
    f"{len(beta):,}"
)

print(
    f"  BETA topology QC:       PASS"
)


## ================================================================== ##
## Discover original Stage-07 Jaccard network
## ================================================================== ##

network = discover_network_edges(
    NETWORK_DIR,
    module_clusters,
)


network[
    "module_a"
] = network[
    "cluster_a"
].map(
    cluster_to_module
)


network[
    "module_b"
] = network[
    "cluster_b"
].map(
    cluster_to_module
)


network[
    "same_mcl_module"
] = (
    network[
        "module_a"
    ]
    ==
    network[
        "module_b"
    ]
).astype(int)


module_network = network[
    network[
        "same_mcl_module"
    ]
    ==
    1
].copy()


module_network = module_network.rename(
    columns={
        "module_a":
            "module",
    }
)


module_network = module_network.drop(
    columns=[
        "module_b",
    ]
)


write_tsv(
    module_network,
    OUT_NETWORK_EDGES,
)


print(
    f"  Final network edges:    "
    f"{len(network):,}"
)

print(
    f"  Within-module edges:    "
    f"{len(module_network):,}"
)

print(
    f"  Cross-module edges:     "
    f"{len(network) - len(module_network):,}"
)


## ================================================================== ##
## Build genome-resolved protein table
## ================================================================== ##

print()
print("Building genome-resolved protein evidence...")


protein_output_columns = [
    "genome",
    "module",
    "cluster",
    "protein_id",

    "evidence_start",
    "evidence_end",
    "evidence_strand",
    "evidence_gene_rank",

    "evidence_number_of_hemes",

    "evidence_fegenie_positive",
    "evidence_fegenie_HMMs",
    "evidence_fegenie_categories",

    "evidence_findmehemes_positive",

    "evidence_signalp_prediction",
    "evidence_deeptmhmm_class",
    "evidence_deeptmhmm_n_tm_helices",

    "report_topology_class",

    "evidence_globdb_cog",
    "evidence_globdb_gene",
    "evidence_globdb_product",

    "evidence_annotation_match_type",
    "evidence_annotation_accepted",
]


protein_output_columns = [
    column

    for column
    in protein_output_columns

    if column
    in proteins.columns
]


genome_proteins = proteins[
    protein_output_columns
].copy()


write_tsv(
    genome_proteins,
    OUT_GENOME_PROTEINS,
)


## ================================================================== ##
## Prepare intervening-gene lookup by genome + module + pair
## ================================================================== ##

print("Preparing genome-level intervening-gene lookup...")


intervening[
    "_protein_key"
] = (
    intervening[
        "genome"
    ]
    +
    "\t"
    +
    intervening[
        "intervening_protein_id"
    ]
)


if (
    intervening[
        "_protein_key"
    ]
    .nunique()
    !=
    EXPECTED_UNIQUE_INTERVENING_PROTEINS
):

    fail(
        f"Expected {EXPECTED_UNIQUE_INTERVENING_PROTEINS:,} "
        f"unique intervening proteins; found "
        f"{intervening['_protein_key'].nunique():,}."
    )


intervening_pair_map = {}


for (
    genome,
    module,
    cluster_a,
    cluster_b,
), group in intervening.groupby(
    [
        "genome",
        "module",
        "cluster_a",
        "cluster_b",
    ],
    sort=False,
):

    ids = sorted(
        group[
            "intervening_protein_id"
        ].unique()
    )


    intervening_pair_map[
        (
            genome,
            module,
            cluster_a,
            cluster_b,
        )
    ] = ids


## ================================================================== ##
## Prepare pair observations for genome-resolved architecture
## ================================================================== ##

print("Building per-genome pair architecture...")


for column in [
    "any_same_contig",
    "any_adjacent",
    "any_within_5kb",
    "any_within_10kb",
    "any_within_20kb",
]:

    pair_obs[
        column
    ] = pd.to_numeric(
        pair_obs[
            column
        ],
        errors="raise",
    ).astype(int)


pair_obs[
    "_broad_orientation"
] = pair_obs[
    "closest_orientation_class"
].map(
    broad_orientation
)


pair_obs[
    "_transcriptional_order"
] = pair_obs.apply(
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
            "_broad_orientation"
        ]
        ==
        "codirectional"
        else
        "",
    axis=1,
)


def describe_local_pair(
    row,
):

    cluster_a = row[
        "cluster_a"
    ]


    cluster_b = row[
        "cluster_b"
    ]


    if row[
        "_transcriptional_order"
    ] != "":

        pair_label = row[
            "_transcriptional_order"
        ]


    else:

        pair_label = (
            f"{cluster_a}"
            f"<->{cluster_b}"
        )


    if row[
        "any_adjacent"
    ] == 1:

        scale = "adjacent"


    elif row[
        "any_within_5kb"
    ] == 1:

        scale = "<=5kb"


    elif row[
        "any_within_10kb"
    ] == 1:

        scale = "<=10kb"


    elif row[
        "any_within_20kb"
    ] == 1:

        scale = "<=20kb"


    else:

        scale = ""


    gap = clean(
        row[
            "closest_intergenic_gap_bp"
        ]
    )


    genes_between = clean(
        row[
            "closest_genes_between"
        ]
    )


    pair_key = canonical_pair(
        cluster_a,
        cluster_b,
    )


    intervening_ids = intervening_pair_map.get(
        (
            row[
                "genome"
            ],
            row[
                "module"
            ],
            pair_key[0],
            pair_key[1],
        ),
        [],
    )


    parts = [
        pair_label,
        scale,
    ]


    if gap != "":

        parts.append(
            f"gap={gap}bp"
        )


    if genes_between != "":

        parts.append(
            f"between={genes_between}"
        )


    if len(
        intervening_ids
    ) > 0:

        parts.append(
            "intervening="
            +
            ",".join(
                intervening_ids
            )
        )


    return " [" + "; ".join(
        part

        for part
        in parts[1:]

        if part != ""
    ) + "]" if pair_label == "" else (
        pair_label
        +
        " ["
        +
        "; ".join(
            part

            for part
            in parts[1:]

            if part != ""
        )
        +
        "]"
    )


pair_obs[
    "_local_pair_description"
] = pair_obs.apply(
    lambda row:
        describe_local_pair(
            row
        )
        if row[
            "any_within_20kb"
        ]
        ==
        1
        else
        "",
    axis=1,
)


## ================================================================== ##
## Build genome x module occurrence table directly from proteins
## ================================================================== ##

print("Building genome x module occurrence/completeness table...")


genome_arch_rows = []


for (
    genome,
    module,
), group in proteins.groupby(
    [
        "genome",
        "module",
    ],
    sort=True,
):

    clusters_present = sorted(
        group[
            "cluster"
        ].unique()
    )


    n_clusters_present = len(
        clusters_present
    )


    size = int(
        module_size[
            module
        ]
    )


    n_module_proteins = len(
        group
    )


    paralogue_excess = (
        n_module_proteins
        -
        n_clusters_present
    )


    completeness = (
        n_clusters_present
        /
        size
    )


    pair_group = pair_obs[
        (
            pair_obs[
                "genome"
            ]
            ==
            genome
        )
        &
        (
            pair_obs[
                "module"
            ]
            ==
            module
        )
    ]


    local_pairs = pair_group[
        pair_group[
            "any_within_20kb"
        ]
        ==
        1
    ]


    local_descriptions = sorted(
        description

        for description
        in local_pairs[
            "_local_pair_description"
        ]

        if description != ""
    )


    cluster_protein_strings = []


    for cluster, cluster_group in group.groupby(
        "cluster",
        sort=True,
    ):

        ids = sorted(
            cluster_group[
                "protein_id"
            ].unique()
        )


        cluster_protein_strings.append(
            f"{cluster}="
            +
            ",".join(
                ids
            )
        )


    fegenie_strings = []


    if (
        "evidence_fegenie_positive"
        in group.columns
        and
        "evidence_fegenie_HMMs"
        in group.columns
    ):

        fe_numeric = numeric(
            group[
                "evidence_fegenie_positive"
            ]
        )


        fe_group = group[
            fe_numeric
            ==
            1
        ]


        for row in fe_group.itertuples(
            index=False
        ):

            hmm = clean(
                getattr(
                    row,
                    "evidence_fegenie_HMMs",
                    "",
                )
            )


            fegenie_strings.append(
                f"{row.cluster}:"
                f"{row.protein_id}:"
                f"{hmm}"
            )


    genome_arch_rows.append(
        {
            "genome":
                genome,

            "module":
                module,

            "module_size":
                size,

            "n_clusters_present":
                n_clusters_present,

            "completeness":
                completeness,

            "is_full":
                int(
                    n_clusters_present
                    ==
                    size
                ),

            "n_module_proteins":
                n_module_proteins,

            "paralogue_excess":
                paralogue_excess,

            "clusters_present":
                "; ".join(
                    clusters_present
                ),

            "cluster_protein_ids":
                "; ".join(
                    cluster_protein_strings
                ),

            "fegenie_positive_module_proteins":
                "; ".join(
                    sorted(
                        fegenie_strings
                    )
                ),

            "n_copresent_cluster_pairs":
                len(
                    pair_group
                ),

            "n_same_contig_pairs":
                int(
                    pair_group[
                        "any_same_contig"
                    ].sum()
                )
                if len(
                    pair_group
                ) > 0
                else
                0,

            "n_local_20kb_pairs":
                int(
                    pair_group[
                        "any_within_20kb"
                    ].sum()
                )
                if len(
                    pair_group
                ) > 0
                else
                0,

            "n_local_5kb_pairs":
                int(
                    pair_group[
                        "any_within_5kb"
                    ].sum()
                )
                if len(
                    pair_group
                ) > 0
                else
                0,

            "n_adjacent_pairs":
                int(
                    pair_group[
                        "any_adjacent"
                    ].sum()
                )
                if len(
                    pair_group
                ) > 0
                else
                0,

            "local_pair_architecture":
                "; ".join(
                    local_descriptions
                ),
        }
    )


genome_arch = pd.DataFrame(
    genome_arch_rows
)


if len(
    genome_arch
) != EXPECTED_GENOME_MODULE_PAIRS:

    fail(
        f"Expected {EXPECTED_GENOME_MODULE_PAIRS:,} "
        f"genome-module occurrences; found "
        f"{len(genome_arch):,}."
    )


if int(
    genome_arch[
        "is_full"
    ].sum()
) != EXPECTED_FULL_GENOME_MODULE_PAIRS:

    fail(
        f"Expected {EXPECTED_FULL_GENOME_MODULE_PAIRS:,} "
        f"full genome-module occurrences; found "
        f"{int(genome_arch['is_full'].sum()):,}."
    )


write_tsv(
    genome_arch,
    OUT_GENOME_ARCH,
)


print(
    f"  Genome-module rows:     "
    f"{len(genome_arch):,}"
)

print(
    f"  Full module occurrences:"
    f" {int(genome_arch['is_full'].sum()):,}"
)


## ================================================================== ##
## Build non-focal-module COG context
##
## This solves the Module-20 ambiguity where Cluster_00035 itself may
## carry GlobDB COG3005 and would otherwise look like a separate NapC
## neighbor of Cluster_00048.
##
## Here we exclude:
##
##   * focal gene itself
##   * any gene belonging to the SAME MCL module as the focal protein
##
## and summarize only positions <=5 oriented genes away.
##
## This is descriptive support, not a censoring-adjusted conservation
## percentage.
## ================================================================== ##

print()
print("Building non-focal-module COG context within +/-5 genes...")


neigh = read_tsv(
    OBSERVED_NEIGHBORHOODS
)


require_columns(
    neigh,
    [
        "genome",
        "focal_protein_id",
        "focal_cluster",
        "focal_module",

        "neighbor_protein_id",
        "oriented_gene_offset",

        "is_focal",
        "neighbor_same_module_as_focal",

        "neighbor_globdb_cog",
        "neighbor_globdb_product",
    ],
    OBSERVED_NEIGHBORHOODS.name,
)


neigh[
    "_offset"
] = numeric(
    neigh[
        "oriented_gene_offset"
    ]
)


neigh[
    "_is_focal"
] = numeric(
    neigh[
        "is_focal"
    ]
)


neigh[
    "_same_module"
] = numeric(
    neigh[
        "neighbor_same_module_as_focal"
    ]
)


context = neigh[
    (
        neigh[
            "_is_focal"
        ]
        !=
        1
    )
    &
    (
        neigh[
            "_same_module"
        ]
        !=
        1
    )
    &
    (
        neigh[
            "_offset"
        ].abs()
        <=
        5
    )
    &
    (
        neigh[
            "neighbor_globdb_cog"
        ]
        !=
        ""
    )
].copy()


context[
    "direction"
] = np.where(
    context[
        "_offset"
    ]
    <
    0,
    "upstream",
    "downstream",
)


context[
    "_focal_occurrence_key"
] = (
    context[
        "genome"
    ]
    +
    "\t"
    +
    context[
        "focal_protein_id"
    ]
)


nonmodule_cog_rows = []


for (
    focal_module,
    focal_cluster,
    direction,
    cog,
), group in context.groupby(
    [
        "focal_module",
        "focal_cluster",
        "direction",
        "neighbor_globdb_cog",
    ],
    sort=True,
):

    products = [
        value

        for value
        in group[
            "neighbor_globdb_product"
        ]

        if clean(
            value
        )
        !=
        ""
    ]


    offsets = sorted(
        int(
            x
        )

        for x
        in group[
            "_offset"
        ].dropna()
    )


    nonmodule_cog_rows.append(
        {
            "module":
                focal_module,

            "focal_cluster":
                focal_cluster,

            "direction":
                direction,

            "neighbor_cog":
                cog,

            "n_neighbor_gene_rows":
                len(
                    group
                ),

            "n_focal_occurrences_with_cog":
                group[
                    "_focal_occurrence_key"
                ]
                .nunique(),

            "n_genomes_with_cog":
                group[
                    "genome"
                ]
                .nunique(),

            "minimum_absolute_gene_offset":
                min(
                    abs(
                        value
                    )

                    for value
                    in offsets
                )
                if len(
                    offsets
                ) > 0
                else
                np.nan,

            "oriented_offset_distribution":
                distribution_string(
                    [
                        str(
                            value
                        )

                        for value
                        in offsets
                    ]
                ),

            "dominant_product":
                dominant_value(
                    products
                ),

            "product_distribution":
                distribution_string(
                    products
                ),
        }
    )


nonmodule_cog = pd.DataFrame(
    nonmodule_cog_rows
)


write_tsv(
    nonmodule_cog,
    OUT_NONMODULE_COG,
)


print(
    f"  Non-module COG associations: "
    f"{len(nonmodule_cog):,}"
)


## ================================================================== ##
## Member-family annotation crosswalk
## ================================================================== ##

print()
print("Preparing member-family annotation crosswalk...")


crosswalk_columns = [
    "module",
    "cluster",

    "representative_protein",

    "n_proteins",
    "n_genomes",

    "mean_heme_count",
    "min_heme_count",
    "max_heme_count",

    "n_fegenie_positive",
    "n_fegenie_negative",

    "dominant_fegenie_HMM",
    "report_fegenie_HMM_distribution",

    "report_globdb_cog_distribution",
    "report_globdb_product_distribution",

    "report_signalp_distribution",
    "report_deeptmhmm_class_distribution",

    "dominant_report_topology",
    "dominant_report_topology_support",
    "report_topology_distribution",
]


crosswalk_columns = [
    column

    for column
    in crosswalk_columns

    if column
    in cluster_report.columns
]


crosswalk = cluster_report[
    crosswalk_columns
].copy()


## ================================================================== ##
## Build refined module summary
## ================================================================== ##

print("Building refined module summary...")


summary_rows = []


for module in modules:

    member = crosswalk[
        crosswalk[
            "module"
        ]
        ==
        module
    ]


    genomes = genome_arch[
        genome_arch[
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


    local_edges = module_edges[
        numeric(
            module_edges[
                "n_genomes_within_20kb"
            ]
        )
        >
        0
    ]


    network_edges = module_network[
        module_network[
            "module"
        ]
        ==
        module
    ]


    module_intervening = intervening[
        intervening[
            "module"
        ]
        ==
        module
    ]


    summary_rows.append(
        {
            "module":
                module,

            "n_member_clusters":
                len(
                    member
                ),

            "n_genomes":
                len(
                    genomes
                ),

            "n_full_genomes":
                int(
                    genomes[
                        "is_full"
                    ].sum()
                ),

            "pct_full_genomes":
                pct(
                    int(
                        genomes[
                            "is_full"
                        ].sum()
                    ),
                    len(
                        genomes
                    ),
                ),

            "n_possible_family_pairs":
                len(
                    module_edges
                ),

            "n_local_20kb_edges":
                len(
                    local_edges
                ),

            "n_direct_jaccard_edges_within_module":
                len(
                    network_edges
                ),

            "median_within_module_jaccard":
                (
                    network_edges[
                        "jaccard"
                    ].median()
                    if len(
                        network_edges
                    ) > 0
                    else
                    np.nan
                ),

            "maximum_within_module_jaccard":
                (
                    network_edges[
                        "jaccard"
                    ].max()
                    if len(
                        network_edges
                    ) > 0
                    else
                    np.nan
                ),

            "n_intervening_gene_rows":
                len(
                    module_intervening
                ),

            "n_unique_intervening_proteins":
                module_intervening[
                    "_protein_key"
                ]
                .nunique(),
        }
    )


refined_summary = pd.DataFrame(
    summary_rows
)


write_tsv(
    refined_summary,
    OUT_REFINED_SUMMARY,
)


## ================================================================== ##
## Integrated interpretation generator
## ================================================================== ##

def family_evidence_sentence(
    row,
):

    cluster = clean(
        row.get(
            "cluster",
            "",
        )
    )


    n_proteins = clean(
        row.get(
            "n_proteins",
            "",
        )
    )


    n_genomes = clean(
        row.get(
            "n_genomes",
            "",
        )
    )


    hemes = clean(
        row.get(
            "mean_heme_count",
            "",
        )
    )


    fegenie = clean(
        row.get(
            "report_fegenie_HMM_distribution",
            "",
        )
    )


    topology = clean(
        row.get(
            "dominant_report_topology",
            "",
        )
    )


    topology_support = clean(
        row.get(
            "dominant_report_topology_support",
            "",
        )
    )


    parts = []


    if (
        n_proteins != ""
        and
        n_genomes != ""
    ):

        parts.append(
            f"{n_proteins} proteins / {n_genomes} genomes"
        )


    if hemes != "":

        parts.append(
            f"mean heme count {compact_float(hemes, 1)}"
        )


    if fegenie != "":

        parts.append(
            f"FeGenie {fegenie}"
        )


    else:

        parts.append(
            "no FeGenie family assignment"
        )


    if topology != "":

        if topology_support != "":

            parts.append(
                f"{topology} ({topology_support})"
            )


        else:

            parts.append(
                topology
            )


    return (
        f"**{cluster}**: "
        +
        "; ".join(
            parts
        )
        +
        "."
    )


def architecture_edge_sentence(
    row,
):

    a = clean(
        row.get(
            "cluster_a",
            "",
        )
    )


    b = clean(
        row.get(
            "cluster_b",
            "",
        )
    )


    order = clean(
        row.get(
            "dominant_transcriptional_order",
            "",
        )
    )


    if order != "":

        label = order.replace(
            "->",
            " → ",
        )


    else:

        label = (
            f"{a} ↔ {b}"
        )


    within20 = clean(
        row.get(
            "within_20kb_support",
            "",
        )
    )


    within5 = clean(
        row.get(
            "within_5kb_support",
            "",
        )
    )


    adjacent = clean(
        row.get(
            "adjacent_support",
            "",
        )
    )


    orient = clean(
        row.get(
            "dominant_local_orientation",
            "",
        )
    )


    orient_support = clean(
        row.get(
            "dominant_local_orientation_support",
            "",
        )
    )


    order_support = clean(
        row.get(
            "dominant_transcriptional_order_support",
            "",
        )
    )


    parts = []


    if within20 != "":

        parts.append(
            f"≤20 kb {within20}"
        )


    if within5 != "":

        parts.append(
            f"≤5 kb {within5}"
        )


    if adjacent != "":

        parts.append(
            f"adjacent {adjacent}"
        )


    if (
        orient != ""
        and
        orient_support != ""
    ):

        parts.append(
            f"{orient} {orient_support}"
        )


    if (
        order != ""
        and
        order_support != ""
    ):

        parts.append(
            f"order {order_support}"
        )


    return (
        f"**{label}** — "
        +
        "; ".join(
            parts
        )
        +
        "."
    )


## ================================================================== ##
## Render refined reports
## ================================================================== ##

print()
print("Rendering refined module reports...")


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


    report_path = (
        outdir
        / "module_report.md"
    )


    member = crosswalk[
        crosswalk[
            "module"
        ]
        ==
        module
    ].copy()


    module_genomes = genome_arch[
        genome_arch[
            "module"
        ]
        ==
        module
    ].copy()


    module_protein_rows = genome_proteins[
        genome_proteins[
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


    local_edges = module_edges[
        numeric(
            module_edges[
                "n_genomes_within_20kb"
            ]
        )
        >
        0
    ].copy()


    module_nodes = nodes[
        nodes[
            "module"
        ]
        ==
        module
    ].copy()


    network_edges = module_network[
        module_network[
            "module"
        ]
        ==
        module
    ].copy()


    module_cog = nonmodule_cog[
        nonmodule_cog[
            "module"
        ]
        ==
        module
    ].copy()


    module_intervening = intervening[
        intervening[
            "module"
        ]
        ==
        module
    ].copy()


    summary = refined_summary[
        refined_summary[
            "module"
        ]
        ==
        module
    ].iloc[0]


    ## -------------------------------------------------------------- ##
    ## Write complete module-level supporting TSVs
    ## -------------------------------------------------------------- ##

    write_tsv(
        member,
        outdir
        /
        "member_family_annotation_crosswalk.tsv",
    )


    write_tsv(
        local_edges,
        outdir
        /
        "local_architecture.tsv",
    )


    write_tsv(
        network_edges,
        outdir
        /
        "network_edges.tsv",
    )


    write_tsv(
        module_genomes,
        outdir
        /
        "genome_resolved_architecture.tsv",
    )


    write_tsv(
        module_protein_rows,
        outdir
        /
        "genome_resolved_proteins.tsv",
    )


    write_tsv(
        module_cog,
        outdir
        /
        "nonmodule_cog_context_5genes.tsv",
    )


    write_tsv(
        module_intervening.drop(
            columns=[
                "_protein_key",
            ],
            errors="ignore",
        ),
        outdir
        /
        "intervening_gene_evidence.tsv",
    )


    ## -------------------------------------------------------------- ##
    ## Start Markdown report
    ## -------------------------------------------------------------- ##

    lines = []


    lines.append(
        f"# {module}"
    )


    lines.append(
        ""
    )


    ## -------------------------------------------------------------- ##
    ## Integrated interpretation
    ## -------------------------------------------------------------- ##

    lines.append(
        "## Integrated module interpretation"
    )


    lines.append(
        ""
    )


    n_genomes = int(
        summary[
            "n_genomes"
        ]
    )


    n_full = int(
        summary[
            "n_full_genomes"
        ]
    )


    size = int(
        summary[
            "n_member_clusters"
        ]
    )


    n_possible = int(
        summary[
            "n_possible_family_pairs"
        ]
    )


    n_local = int(
        summary[
            "n_local_20kb_edges"
        ]
    )


    lines.append(
        f"{module} contains **{size} MMseqs2 sequence families** "
        f"and occurs in **{n_genomes} genomes**. "
        f"All {size} member families occur together in "
        f"**{n_full}/{n_genomes} genomes "
        f"({pct(n_full, n_genomes):.1f}%)**."
    )


    lines.append(
        ""
    )


    lines.append(
        f"Of {n_possible} possible within-module family pairs, "
        f"**{n_local} have observed ≤20-kb physical association**. "
        f"The original Jaccard network contains "
        f"**{int(summary['n_direct_jaccard_edges_within_module'])} "
        f"direct within-module edges**."
    )


    lines.append(
        ""
    )


    ## Most recurrent local architecture edges. ##

    if len(
        local_edges
    ) > 0:

        lines.append(
            "**Most frequently observed local relationships:**"
        )


        lines.append(
            ""
        )


        ranked_local = local_edges.copy()


        ranked_local[
            "_rank"
        ] = numeric(
            ranked_local[
                "n_genomes_within_20kb"
            ]
        ).fillna(
            0
        )


        ranked_local = ranked_local.sort_values(
            [
                "_rank",
                "cluster_a",
                "cluster_b",
            ],
            ascending=[
                False,
                True,
                True,
            ],
            kind="stable",
        )


        for _, row in ranked_local.head(
            5
        ).iterrows():

            lines.append(
                "- "
                +
                architecture_edge_sentence(
                    row
                )
            )


        lines.append(
            ""
        )


    ## Family evidence. For small modules show all; for large modules
    ## show FeGenie-annotated families plus most prevalent families. ##

    if size <= 6:

        family_display = member.copy()


    else:

        annotated_mask = (
            member[
                "report_fegenie_HMM_distribution"
            ]
            !=
            ""
        ) if (
            "report_fegenie_HMM_distribution"
            in member.columns
        ) else pd.Series(
            False,
            index=member.index,
        )


        annotated = member[
            annotated_mask
        ].copy()


        remainder = member[
            ~annotated_mask
        ].copy()


        if "n_genomes" in remainder.columns:

            remainder[
                "_n_genomes"
            ] = numeric(
                remainder[
                    "n_genomes"
                ]
            ).fillna(
                0
            )


            remainder = remainder.sort_values(
                [
                    "_n_genomes",
                    "cluster",
                ],
                ascending=[
                    False,
                    True,
                ],
            )


        family_display = pd.concat(
            [
                annotated,
                remainder.head(
                    max(
                        0,
                        6
                        -
                        len(
                            annotated
                        ),
                    )
                ),
            ],
            ignore_index=True,
        ).drop_duplicates(
            "cluster"
        )


    lines.append(
        "**Key member-family evidence:**"
    )


    lines.append(
        ""
    )


    for _, row in family_display.iterrows():

        lines.append(
            "- "
            +
            family_evidence_sentence(
                row
            )
        )


    lines.append(
        ""
    )


    ## -------------------------------------------------------------- ##
    ## Occurrence / completeness
    ## -------------------------------------------------------------- ##

    lines.append(
        "## Occurrence and completeness"
    )


    lines.append(
        ""
    )


    lines.append(
        f"- Genomes with ≥1 module family: **{n_genomes}**"
    )


    lines.append(
        f"- Genomes containing all {size} families: "
        f"**{n_full}/{n_genomes} ({pct(n_full, n_genomes):.1f}%)**"
    )


    median_completeness = module_genomes[
        "completeness"
    ].median()


    lines.append(
        f"- Median completeness: "
        f"**{100 * median_completeness:.1f}%**"
    )


    lines.append(
        ""
    )


    ## -------------------------------------------------------------- ##
    ## Member annotation crosswalk
    ## -------------------------------------------------------------- ##

    lines.append(
        "## Member-family annotation crosswalk"
    )


    lines.append(
        ""
    )


    lines.append(
        "This table keeps sequence-family identity separate from "
        "independent annotation systems. A GlobDB NapC-like COG and a "
        "FeGenie MtoA assignment can therefore coexist for the same "
        "MMseqs2 family rather than being interpreted as separate genes."
    )


    lines.append(
        ""
    )


    crosswalk_display_columns = [
        "cluster",
        "n_proteins",
        "n_genomes",
        "mean_heme_count",

        "report_fegenie_HMM_distribution",

        "report_globdb_cog_distribution",
        "report_globdb_product_distribution",

        "dominant_report_topology",
        "dominant_report_topology_support",

        "report_signalp_distribution",
        "report_deeptmhmm_class_distribution",
    ]


    lines.append(
        markdown_table(
            member,
            crosswalk_display_columns,
        )
    )


    lines.append(
        ""
    )


    ## -------------------------------------------------------------- ##
    ## Physical architecture
    ## -------------------------------------------------------------- ##

    lines.append(
        "## Supported module architecture"
    )


    lines.append(
        ""
    )


    if len(
        local_edges
    ) == 0:

        lines.append(
            "No within-module family pair has observed ≤20-kb "
            "physical association."
        )


    else:

        local_display = local_edges.copy()


        local_display[
            "_rank"
        ] = numeric(
            local_display[
                "n_genomes_within_20kb"
            ]
        ).fillna(
            0
        )


        local_display = local_display.sort_values(
            [
                "_rank",
                "cluster_a",
                "cluster_b",
            ],
            ascending=[
                False,
                True,
                True,
            ],
        ).head(
            MAX_LOCAL_EDGES_DISPLAY
        )


        architecture_columns = [
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

            "median_local_intergenic_gap_bp",
            "median_local_genes_between",
        ]


        lines.append(
            markdown_table(
                local_display,
                architecture_columns,
            )
        )


        if len(
            local_edges
        ) > MAX_LOCAL_EDGES_DISPLAY:

            lines.append(
                ""
            )

            lines.append(
                f"_Showing {MAX_LOCAL_EDGES_DISPLAY} of "
                f"{len(local_edges)} local edges. "
                f"See `local_architecture.tsv` for all relationships._"
            )


    lines.append(
        ""
    )


    ## Local components only useful when >1. ##

    component_count = (
        module_nodes[
            "local_component"
        ]
        .nunique()
        if "local_component"
        in module_nodes.columns
        else
        0
    )


    if component_count > 1:

        lines.append(
            "### Local architecture components"
        )


        lines.append(
            ""
        )


        component_columns = [
            "cluster",
            "local_component",
            "local_component_size",
            "local_architecture_degree",
        ]


        lines.append(
            markdown_table(
                module_nodes,
                component_columns,
            )
        )


        lines.append(
            ""
        )


    ## -------------------------------------------------------------- ##
    ## Genome-resolved architecture
    ## -------------------------------------------------------------- ##

    lines.append(
        "## Genome-resolved architecture"
    )


    lines.append(
        ""
    )


    lines.append(
        "The complete table is stored in "
        "`genome_resolved_architecture.tsv`; exact protein-level "
        "annotations are stored in `genome_resolved_proteins.tsv`."
    )


    lines.append(
        ""
    )


    genome_display = module_genomes.copy()


    genome_display = genome_display.sort_values(
        [
            "is_full",
            "n_clusters_present",
            "n_local_20kb_pairs",
            "n_adjacent_pairs",
            "genome",
        ],
        ascending=[
            False,
            False,
            False,
            False,
            True,
        ],
        kind="stable",
    )


    if len(
        genome_display
    ) > MAX_GENOME_ROWS_IN_MARKDOWN:

        genome_display = genome_display.head(
            MAX_GENOME_ROWS_IN_MARKDOWN
        )


        lines.append(
            f"_Displaying the first {MAX_GENOME_ROWS_IN_MARKDOWN} "
            f"genomes after sorting by completeness and local "
            f"architecture. The TSV contains all {len(module_genomes)}._"
        )


        lines.append(
            ""
        )


    genome_columns = [
        "genome",

        "n_clusters_present",
        "module_size",
        "completeness",
        "is_full",

        "n_module_proteins",
        "paralogue_excess",

        "clusters_present",

        "n_local_20kb_pairs",
        "n_adjacent_pairs",

        "local_pair_architecture",
    ]


    lines.append(
        markdown_table(
            genome_display,
            genome_columns,
        )
    )


    lines.append(
        ""
    )


    ## -------------------------------------------------------------- ##
    ## Network evidence
    ## -------------------------------------------------------------- ##

    lines.append(
        "## Original Jaccard-network evidence"
    )


    lines.append(
        ""
    )


    lines.append(
        "These are the direct Stage-07 co-occurrence-network edges "
        "connecting member families within this MCL module. "
        "MCL membership does not require every pair of module families "
        "to have a direct network edge."
    )


    lines.append(
        ""
    )


    if len(
        network_edges
    ) == 0:

        lines.append(
            "_No direct within-module network edge recovered._"
        )


    else:

        network_display = network_edges.sort_values(
            [
                "jaccard",
                "cluster_a",
                "cluster_b",
            ],
            ascending=[
                False,
                True,
                True,
            ],
        ).head(
            MAX_NETWORK_EDGES_DISPLAY
        )


        network_columns = [
            "cluster_a",
            "cluster_b",
            "jaccard",
            "shared_genomes",
        ]


        lines.append(
            markdown_table(
                network_display,
                network_columns,
            )
        )


        if len(
            network_edges
        ) > MAX_NETWORK_EDGES_DISPLAY:

            lines.append(
                ""
            )


            lines.append(
                f"_Showing {MAX_NETWORK_EDGES_DISPLAY} of "
                f"{len(network_edges)} direct network edges. "
                f"See `network_edges.tsv` for the complete table._"
            )


    lines.append(
        ""
    )


    ## -------------------------------------------------------------- ##
    ## Non-module functional context
    ## -------------------------------------------------------------- ##

    lines.append(
        "## Non-module functional neighborhood context"
    )


    lines.append(
        ""
    )


    lines.append(
        "This section summarizes conventional COG/product annotations "
        "within ±5 transcription-oriented genes after excluding genes "
        "belonging to the same MCL module as the focal family. "
        "It is therefore intended to describe surrounding context "
        "without counting another module member's own GlobDB annotation "
        "as though it were a separate neighboring gene."
    )


    lines.append(
        ""
    )


    lines.append(
        "These are descriptive observed counts, not censoring-adjusted "
        "conservation percentages; the complete Stage-13C tables remain "
        "authoritative for annotation-conservation analysis."
    )


    lines.append(
        ""
    )


    cog_display_parts = []


    for focal_cluster, group in module_cog.groupby(
        "focal_cluster",
        sort=True,
    ):

        group = group.copy()


        group[
            "_rank"
        ] = numeric(
            group[
                "n_focal_occurrences_with_cog"
            ]
        ).fillna(
            0
        )


        group = group.sort_values(
            [
                "_rank",
                "n_genomes_with_cog",
                "neighbor_cog",
            ],
            ascending=[
                False,
                False,
                True,
            ],
        ).head(
            MAX_NONMODULE_COG_PER_FOCAL
        )


        cog_display_parts.append(
            group
        )


    if len(
        cog_display_parts
    ) == 0:

        lines.append(
            "_No non-module COG annotations observed within ±5 genes._"
        )


    else:

        cog_display = pd.concat(
            cog_display_parts,
            ignore_index=True,
        )


        cog_columns = [
            "focal_cluster",
            "direction",
            "neighbor_cog",

            "n_focal_occurrences_with_cog",
            "n_genomes_with_cog",

            "minimum_absolute_gene_offset",

            "dominant_product",
        ]


        lines.append(
            markdown_table(
                cog_display,
                cog_columns,
            )
        )


    lines.append(
        ""
    )


    ## -------------------------------------------------------------- ##
    ## Intervening genes
    ## -------------------------------------------------------------- ##

    lines.append(
        "## Genes between locally associated module families"
    )


    lines.append(
        ""
    )


    n_intervening_rows = len(
        module_intervening
    )


    n_unique_intervening = module_intervening[
        "_protein_key"
    ].nunique()


    lines.append(
        f"- Intervening-gene observations: "
        f"**{n_intervening_rows}**"
    )


    lines.append(
        f"- Unique intervening proteins: "
        f"**{n_unique_intervening}**"
    )


    if (
        "evidence_findmehemes_positive"
        in module_intervening.columns
    ):

        fmh = (
            numeric(
                module_intervening[
                    "evidence_findmehemes_positive"
                ]
            )
            ==
            1
        )


        lines.append(
            f"- FindMeHemes-positive observations: "
            f"**{int(fmh.sum())}**"
        )


    if (
        "evidence_fegenie_positive"
        in module_intervening.columns
    ):

        fe = (
            numeric(
                module_intervening[
                    "evidence_fegenie_positive"
                ]
            )
            ==
            1
        )


        lines.append(
            f"- FeGenie-positive observations: "
            f"**{int(fe.sum())}**"
        )


    lines.append(
        ""
    )


    if len(
        module_intervening
    ) > 0:

        intervening_display = module_intervening.head(
            MAX_INTERVENING_ROWS_DISPLAY
        )


        intervening_columns = [
            "genome",
            "cluster_a",
            "cluster_b",

            "intervening_protein_id",

            "evidence_fegenie_HMMs",
            "evidence_number_of_hemes",

            "evidence_signalp_prediction",
            "evidence_deeptmhmm_class",

            "evidence_globdb_cog",
            "evidence_globdb_product",
        ]


        lines.append(
            markdown_table(
                intervening_display,
                intervening_columns,
            )
        )


        if len(
            module_intervening
        ) > MAX_INTERVENING_ROWS_DISPLAY:

            lines.append(
                ""
            )


            lines.append(
                f"_Showing {MAX_INTERVENING_ROWS_DISPLAY} of "
                f"{len(module_intervening)} observations. "
                f"See `intervening_gene_evidence.tsv` for all rows._"
            )


    else:

        lines.append(
            "_No intervening genes occur between locally associated "
            "module families._"
        )


    lines.append(
        ""
    )


    ## -------------------------------------------------------------- ##
    ## Interpretation / methods caveats
    ## -------------------------------------------------------------- ##

    lines.append(
        "## Interpretation notes"
    )


    lines.append(
        ""
    )


    lines.append(
        "The MCL module is a co-occurrence community, not a predefined "
        "operon or biochemical system. Local architecture is inferred "
        "independently from physical genomic observations."
    )


    lines.append(
        ""
    )


    lines.append(
        "Sequence-family identity is defined by MMseqs2. FeGenie, "
        "GlobDB COG/product, FindMeHemes, SignalP and DeepTMHMM remain "
        "independent evidence layers and may disagree."
    )


    lines.append(
        ""
    )


    lines.append(
        "All conservation percentages should be interpreted together "
        "with their numerator and denominator. Contig censoring remains "
        "handled in the Stage-13 consensus analyses."
    )


    lines.append(
        ""
    )


    lines.append(
        "The report-level topology field is derived from raw predictors "
        "for interpretation; DeepTMHMM `BETA` is treated as beta-barrel "
        "membrane topology. No upstream localization field is modified."
    )


    lines.append(
        ""
    )


    report_path.write_text(
        "\n".join(
            lines
        )
    )


    index_rows.append(
        {
            "module":
                module,

            "report_path":
                str(
                    report_path.relative_to(
                        HERE
                    )
                ),

            "n_member_clusters":
                size,

            "n_genomes":
                n_genomes,

            "n_full_genomes":
                n_full,

            "n_local_20kb_edges":
                n_local,

            "n_direct_network_edges":
                len(
                    network_edges
                ),

            "n_nonmodule_cog_context_rows":
                len(
                    module_cog
                ),

            "n_intervening_gene_rows":
                n_intervening_rows,

            "n_unique_intervening_proteins":
                n_unique_intervening,
        }
    )


## ================================================================== ##
## Write index
## ================================================================== ##

index = pd.DataFrame(
    index_rows
)


write_tsv(
    index,
    OUT_INDEX,
)


## ================================================================== ##
## Final QC
## ================================================================== ##

print()
print("Running final Stage-14D QC...")


local_edge_count = int(
    (
        numeric(
            edges[
                "n_genomes_within_20kb"
            ]
        )
        >
        0
    )
    .sum()
)


checks = [
    (
        "modules",
        len(
            modules
        ),
        EXPECTED_MODULES,
    ),

    (
        "module_clusters",
        len(
            membership
        ),
        EXPECTED_MODULE_CLUSTERS,
    ),

    (
        "module_proteins",
        len(
            proteins
        ),
        EXPECTED_MODULE_PROTEINS,
    ),

    (
        "genome_module_pairs",
        len(
            genome_arch
        ),
        EXPECTED_GENOME_MODULE_PAIRS,
    ),

    (
        "full_genome_module_pairs",
        int(
            genome_arch[
                "is_full"
            ].sum()
        ),
        EXPECTED_FULL_GENOME_MODULE_PAIRS,
    ),

    (
        "architecture_pairs",
        len(
            edges
        ),
        EXPECTED_ARCH_PAIRS,
    ),

    (
        "local_20kb_edges",
        local_edge_count,
        EXPECTED_LOCAL_20KB_EDGES,
    ),

    (
        "network_edges",
        len(
            network
        ),
        EXPECTED_NETWORK_EDGES,
    ),

    (
        "intervening_gene_rows",
        len(
            intervening
        ),
        EXPECTED_INTERVENING_ROWS,
    ),

    (
        "unique_intervening_proteins",
        intervening[
            "_protein_key"
        ]
        .nunique(),
        EXPECTED_UNIQUE_INTERVENING_PROTEINS,
    ),

    (
        "beta_proteins",
        len(
            beta
        ),
        EXPECTED_BETA_PROTEINS,
    ),

    (
        "reports_rendered",
        len(
            index
        ),
        EXPECTED_MODULES,
    ),
]


qc_rows = []


for metric, observed, expected in checks:

    passed = int(
        observed
        ==
        expected
    )


    qc_rows.append(
        {
            "metric":
                metric,

            "observed":
                observed,

            "expected":
                expected,

            "pass":
                passed,
        }
    )


    if passed != 1:

        fail(
            f"QC failed for {metric}: "
            f"observed {observed}, expected {expected}."
        )


qc = pd.DataFrame(
    qc_rows
)


write_tsv(
    qc,
    OUT_QC,
)


## ================================================================== ##
## Terminal summary
## ================================================================== ##

print()
print("Refined module-report summary")


print(
    f"  Reports rendered:                  "
    f"{len(index):,}"
)

print(
    f"  Genome-module occurrences:         "
    f"{len(genome_arch):,}"
)

print(
    f"  Full genome-module occurrences:    "
    f"{int(genome_arch['is_full'].sum()):,}"
)

print(
    f"  Genome-resolved module proteins:   "
    f"{len(genome_proteins):,}"
)

print(
    f"  Final Jaccard network edges:        "
    f"{len(network):,}"
)

print(
    f"  Within-module Jaccard edges:        "
    f"{len(module_network):,}"
)

print(
    f"  Cross-module Jaccard edges:         "
    f"{len(network) - len(module_network):,}"
)

print(
    f"  <=20-kb architecture edges:         "
    f"{local_edge_count:,}"
)

print(
    f"  Non-module COG context rows:        "
    f"{len(nonmodule_cog):,}"
)

print(
    f"  Intervening-gene observations:      "
    f"{len(intervening):,}"
)

print(
    f"  Unique intervening proteins:        "
    f"{intervening['_protein_key'].nunique():,}"
)

print(
    f"  BETA -> beta-barrel topology:       "
    f"{len(beta):,}/{len(beta):,}"
)


## ================================================================== ##
## Module 20 / Module 35 spotlights
## ================================================================== ##

for spotlight in [
    "Module_20",
    "Module_35",
]:

    print()
    print(
        f"{spotlight} refined-report spotlight"
    )


    s = refined_summary[
        refined_summary[
            "module"
        ]
        ==
        spotlight
    ].iloc[0]


    print(
        f"  Member families:            "
        f"{int(s['n_member_clusters'])}"
    )

    print(
        f"  Genomes:                    "
        f"{int(s['n_genomes'])}"
    )

    print(
        f"  Full genomes:               "
        f"{int(s['n_full_genomes'])}"
    )

    print(
        f"  Local <=20-kb edges:        "
        f"{int(s['n_local_20kb_edges'])}"
    )

    print(
        f"  Direct Jaccard edges:       "
        f"{int(s['n_direct_jaccard_edges_within_module'])}"
    )

    print(
        f"  Intervening observations:   "
        f"{int(s['n_intervening_gene_rows'])}"
    )

    print(
        f"  Unique intervening proteins:"
        f" {int(s['n_unique_intervening_proteins'])}"
    )


    spot_edges = edges[
        edges[
            "module"
        ]
        ==
        spotlight
    ]


    spot_edges = spot_edges[
        numeric(
            spot_edges[
                "n_genomes_within_20kb"
            ]
        )
        >
        0
    ]


    if len(
        spot_edges
    ) > 0:

        display_columns = [
            "cluster_a",
            "cluster_b",

            "n_genomes_both",

            "within_20kb_support",
            "within_5kb_support",
            "adjacent_support",

            "dominant_transcriptional_order",
            "dominant_transcriptional_order_support",
        ]


        print()
        print(
            spot_edges[
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
    f"Refined summary:          "
    f"{OUT_REFINED_SUMMARY}"
)

print(
    f"Genome architecture:      "
    f"{OUT_GENOME_ARCH}"
)

print(
    f"Genome proteins:          "
    f"{OUT_GENOME_PROTEINS}"
)

print(
    f"Module network edges:     "
    f"{OUT_NETWORK_EDGES}"
)

print(
    f"Non-module COG context:   "
    f"{OUT_NONMODULE_COG}"
)

print(
    f"Report index:             "
    f"{OUT_INDEX}"
)

print(
    f"QC:                       "
    f"{OUT_QC}"
)

print(
    f"Refined reports:          "
    f"{MODULE_DIR}"
)
