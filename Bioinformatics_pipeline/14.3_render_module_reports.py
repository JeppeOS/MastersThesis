#!/usr/bin/env python3

from pathlib import Path
from collections import Counter
import sys

import numpy as np
import pandas as pd


## ================================================================== ##
## STAGE 14C - RENDER INTEGRATED MODULE REPORTS
##
## PURPOSE
## -------
##
## Convert the analytical outputs from Stages 09-14B into one
## human-readable evidence report for each MCL module.
##
##
## REPORT HIERARCHY
## ----------------
##
## MCL module
##   ↓
## member MMseqs2 sequence families
##   ↓
## genome occurrence/completeness
##   ↓
## supported local architecture
##   ↓
## focal-family consensus neighborhoods
##   ↓
## COG/product context
##   ↓
## intervening genes
##   ↓
## FeGenie / FindMeHemes / SignalP / DeepTMHMM evidence
##
##
## IMPORTANT
## ---------
##
## The report does NOT assume:
##
##   * an MCL module is an operon
##   * all module families are physically linked
##   * 100% (1/1) is equivalent to 100% (44/44)
##   * FeGenie-negative means functionally irrelevant
##
##
## DISPLAY LIMITS
## --------------
##
## Large COG/MMseqs tables are ranked for readability in the Markdown
## report only.
##
## No underlying analytical rows are discarded.
##
## The complete upstream TSV files remain authoritative.
##
##
## TOPOLOGY INTERPRETATION
## -----------------------
##
## A new report-level topology class is derived from the RAW predictors.
##
## Importantly:
##
##     DeepTMHMM BETA
##
## is interpreted as:
##
##     beta_barrel_membrane
##
## rather than inheriting the older
##
##     soluble_exported_periplasmic_like_candidate
##
## convenience label.
##
## The upstream localization_class is NOT modified.
## ================================================================== ##


## ================================================================== ##
## Paths
## ================================================================== ##

HERE = Path(__file__).resolve().parent
WORKFLOW = HERE.parent.parent


## Stage 14A ##

REPORT_OVERVIEW = (
    WORKFLOW
    / "14_module_reports"
    / "14A_report_data"
    / "global"
    / "module_report_overview.tsv"
)


CLUSTER_EVIDENCE = (
    WORKFLOW
    / "14_module_reports"
    / "14A_report_data"
    / "global"
    / "module_member_cluster_evidence.tsv"
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


## Stage 14B ##

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


## Stage 13C / 13D ##

INTERVENING_PATTERNS = (
    WORKFLOW
    / "13_neighborhood_conservation"
    / "13C_annotation_context"
    / "local_module_pair_intervening_pattern_summary.tsv"
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
## Outputs
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


OUT_SUMMARY = (
    GLOBAL_DIR
    / "module_report_summary.tsv"
)


OUT_CLUSTER_REPORT = (
    GLOBAL_DIR
    / "module_cluster_report_evidence.tsv"
)


OUT_INDEX = (
    GLOBAL_DIR
    / "module_report_index.tsv"
)


OUT_QC = (
    GLOBAL_DIR
    / "module_report_rendering_qc.tsv"
)


## ================================================================== ##
## Expected values
## ================================================================== ##

EXPECTED_MODULES = 35
EXPECTED_CLUSTERS = 156
EXPECTED_MODULE_PROTEINS = 10_537
EXPECTED_EDGES = 560
EXPECTED_LOCAL_EDGES = 192
EXPECTED_INTERVENING_ROWS = 6_950


## ================================================================== ##
## Display parameters
##
## These affect Markdown readability ONLY.
## They do NOT filter analytical outputs.
## ================================================================== ##

MAX_MMSEQ_RELATIONSHIPS_PER_FOCAL = 8
MAX_COG_RELATIONSHIPS_PER_FOCAL = 5
MAX_INTERVENING_PATTERNS = 20


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
            f"Required file not found:\n{path}"
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

    if pd.isna(value):

        return ""

    return str(
        value
    ).strip()


def number(value):

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

    value = number(
        value
    )


    if pd.isna(
        value
    ):

        return 0


    return int(
        value
    )


def is_positive(value):

    value = number(
        value
    )


    return (
        not pd.isna(
            value
        )
        and
        value == 1
    )


def escape_md(value):

    value = clean(
        value
    )


    value = value.replace(
        "|",
        "\\|",
    )


    value = value.replace(
        "\n",
        " ",
    )


    return value


def markdown_table(
    df,
    columns=None,
):

    if len(
        df
    ) == 0:

        return "_No rows._\n"


    if columns is None:

        columns = list(
            df.columns
        )


    columns = [
        column

        for column
        in columns

        if column
        in df.columns
    ]


    if len(
        columns
    ) == 0:

        return "_No displayable columns._\n"


    lines = []


    lines.append(
        "| "
        +
        " | ".join(
            escape_md(
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
                escape_md(
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
        !=
        ""
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
        key=lambda item: (
            -item[1],
            item[0],
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
        !=
        ""
    ]


    if len(
        cleaned
    ) == 0:

        return (
            "",
            0,
            0,
        )


    counter = Counter(
        cleaned
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
            len(
                cleaned
            ),
        )


    return (
        "mixed_tie",
        maximum,
        len(
            cleaned
        ),
    )


def support_string(
    numerator,
    denominator,
):

    return (
        f"{int(numerator)}"
        f"/"
        f"{int(denominator)}"
    )


## ================================================================== ##
## Report-level topology
##
## Raw predictor precedence:
##
## BETA
##   > SignalP LIPO
##   > SP+TM
##   > TM
##   > SP
##   > GLOB
##
## The older upstream localization_class is preserved separately.
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


    n_tm = number(
        row.get(
            "evidence_deeptmhmm_n_tm_helices",
            "",
        )
    )


    ## -------------------------------------------------------------- ##
    ## Beta-barrel membrane topology has highest precedence.
    ## -------------------------------------------------------------- ##

    if deeptmhmm == "BETA":

        return "beta_barrel_membrane"


    ## -------------------------------------------------------------- ##
    ## Lipoprotein signal
    ## -------------------------------------------------------------- ##

    if signalp in [
        "LIPO",
        "TATLIPO",
    ]:

        return "exported_lipoprotein"


    ## -------------------------------------------------------------- ##
    ## Signal peptide + alpha-helical TM topology
    ## -------------------------------------------------------------- ##

    if deeptmhmm == "SP+TM":

        if not pd.isna(
            n_tm
        ):

            if n_tm == 1:

                return (
                    "exported_single_pass_membrane"
                )


            if n_tm > 1:

                return (
                    "exported_multipass_membrane"
                )


        return (
            "exported_membrane_SP_plus_TM"
        )


    ## -------------------------------------------------------------- ##
    ## Membrane topology without DeepTMHMM SP class
    ## -------------------------------------------------------------- ##

    if deeptmhmm == "TM":

        if not pd.isna(
            n_tm
        ):

            if n_tm == 1:

                return "single_pass_membrane"


            if n_tm > 1:

                return "multipass_membrane"


        return "membrane"


    ## -------------------------------------------------------------- ##
    ## Exported but no detected mature membrane topology
    ## -------------------------------------------------------------- ##

    if deeptmhmm == "SP":

        return (
            "soluble_exported_periplasmic_like"
        )


    ## -------------------------------------------------------------- ##
    ## Globular / no membrane topology
    ## -------------------------------------------------------------- ##

    if deeptmhmm == "GLOB":

        return "globular_nonmembrane_like"


    ## -------------------------------------------------------------- ##
    ## SignalP-only fallback
    ## -------------------------------------------------------------- ##

    if signalp in [
        "SP",
        "TAT",
    ]:

        return (
            "exported_topology_uncertain"
        )


    return "uncertain"


## ================================================================== ##
## Read inputs
## ================================================================== ##

print("=" * 80)
print("STAGE 14C - RENDER INTEGRATED MODULE REPORTS")
print("=" * 80)


print()
print("Reading Stage-14 report data...")


overview = read_tsv(
    REPORT_OVERVIEW
)


clusters = read_tsv(
    CLUSTER_EVIDENCE
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


architecture = read_tsv(
    ARCH_SUMMARY
)


patterns = read_tsv(
    INTERVENING_PATTERNS
)


mmseq = read_tsv(
    MMSEQ_CONSENSUS
)


cogs = read_tsv(
    COG_CONSENSUS
)


require_columns(
    overview,
    [
        "module",
        "n_member_clusters",
        "n_module_proteins",
        "n_genomes_with_any_module_protein",
    ],
    REPORT_OVERVIEW.name,
)


require_columns(
    clusters,
    [
        "module",
        "cluster",
    ],
    CLUSTER_EVIDENCE.name,
)


require_columns(
    proteins,
    [
        "genome",
        "protein_id",
        "module",
        "cluster",
    ],
    MODULE_PROTEINS.name,
)


require_columns(
    intervening,
    [
        "genome",
        "module",
        "intervening_protein_id",
    ],
    INTERVENING_EVIDENCE.name,
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
        "has_local_20kb_architecture_edge",
    ],
    ARCH_EDGES.name,
)


require_columns(
    architecture,
    [
        "module",
        "n_member_clusters",
        "n_local_architecture_components",
        "largest_local_component_size",
    ],
    ARCH_SUMMARY.name,
)


modules = sorted(
    overview[
        "module"
    ].unique()
)


if len(
    modules
) != EXPECTED_MODULES:

    fail(
        f"Expected {EXPECTED_MODULES} modules; "
        f"found {len(modules)}."
    )


if len(
    clusters
) != EXPECTED_CLUSTERS:

    fail(
        f"Expected {EXPECTED_CLUSTERS} clusters; "
        f"found {len(clusters)}."
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
) != EXPECTED_EDGES:

    fail(
        f"Expected {EXPECTED_EDGES} architecture edges; "
        f"found {len(edges)}."
    )


if (
    pd.to_numeric(
        edges[
            "has_local_20kb_architecture_edge"
        ],
        errors="raise",
    )
    .sum()
    !=
    EXPECTED_LOCAL_EDGES
):

    fail(
        "Unexpected number of <=20-kb architecture edges."
    )


if len(
    intervening
) != EXPECTED_INTERVENING_ROWS:

    fail(
        f"Expected {EXPECTED_INTERVENING_ROWS:,} "
        f"intervening-gene rows; "
        f"found {len(intervening):,}."
    )


print(
    f"  Modules:                 "
    f"{len(modules):,}"
)

print(
    f"  Member clusters:         "
    f"{len(clusters):,}"
)

print(
    f"  Module proteins:         "
    f"{len(proteins):,}"
)

print(
    f"  Architecture pairs:      "
    f"{len(edges):,}"
)

print(
    f"  Local <=20-kb edges:     "
    f"{EXPECTED_LOCAL_EDGES:,}"
)

print(
    f"  Intervening-gene rows:   "
    f"{len(intervening):,}"
)


## ================================================================== ##
## Derive report-level topology for every module protein
## ================================================================== ##

print()
print("Deriving report-level topology from raw predictors...")


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
) > 0:

    if (
        beta[
            "report_topology_class"
        ]
        !=
        "beta_barrel_membrane"
    ).any():

        fail(
            "A DeepTMHMM BETA protein was not "
            "classified as beta_barrel_membrane."
        )


print(
    f"  DeepTMHMM BETA proteins: "
    f"{len(beta):,}"
)

print(
    f"  BETA -> beta barrel QC:  PASS"
)


## ================================================================== ##
## Summarize report topology per cluster
## ================================================================== ##

print()
print("Building report-level member-family evidence...")


topology_rows = []


for (
    module,
    cluster,
), group in proteins.groupby(
    [
        "module",
        "cluster",
    ],
    sort=True,
):

    (
        dominant_topology,
        n_dominant_topology,
        n_topology_informative,
    ) = dominant_value(
        group[
            "report_topology_class"
        ]
    )


    topology_rows.append(
        {
            "module":
                module,

            "cluster":
                cluster,

            "report_topology_distribution":
                distribution_string(
                    group[
                        "report_topology_class"
                    ]
                ),

            "dominant_report_topology":
                dominant_topology,

            "dominant_report_topology_support":
                support_string(
                    n_dominant_topology,
                    n_topology_informative,
                ),
        }
    )


topology_summary = pd.DataFrame(
    topology_rows
)


cluster_report = clusters.merge(
    topology_summary,
    on=[
        "module",
        "cluster",
    ],
    how="left",
    validate="one_to_one",
)


if len(
    cluster_report
) != EXPECTED_CLUSTERS:

    fail(
        "Cluster report evidence lost clusters."
    )


cluster_report.to_csv(
    OUT_CLUSTER_REPORT,
    sep="\t",
    index=False,
    na_rep="",
    float_format="%.6f",
)


## ================================================================== ##
## Intervening-gene module-level summary
##
## Keep BOTH:
##
##   row observations
##   unique genome + protein observations
##
## because one protein can potentially occur in more than one pair
## context.
## ================================================================== ##

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


summary_rows = []


for module in modules:

    ov = overview[
        overview[
            "module"
        ]
        ==
        module
    ].iloc[0]


    arch = architecture[
        architecture[
            "module"
        ]
        ==
        module
    ].iloc[0]


    mod_intervening = intervening[
        intervening[
            "module"
        ]
        ==
        module
    ].copy()


    n_intervening_rows = len(
        mod_intervening
    )


    n_unique_intervening = (
        mod_intervening[
            "_protein_key"
        ]
        .nunique()
    )


    if (
        "evidence_findmehemes_positive"
        in mod_intervening.columns
    ):

        fmh_mask = (
            pd.to_numeric(
                mod_intervening[
                    "evidence_findmehemes_positive"
                ].replace(
                    "",
                    np.nan,
                ),
                errors="coerce",
            )
            ==
            1
        )

    else:

        fmh_mask = pd.Series(
            False,
            index=mod_intervening.index,
        )


    if (
        "evidence_fegenie_positive"
        in mod_intervening.columns
    ):

        fegenie_mask = (
            pd.to_numeric(
                mod_intervening[
                    "evidence_fegenie_positive"
                ].replace(
                    "",
                    np.nan,
                ),
                errors="coerce",
            )
            ==
            1
        )

    else:

        fegenie_mask = pd.Series(
            False,
            index=mod_intervening.index,
        )


    summary_rows.append(
        {
            "module":
                module,

            "n_member_clusters":
                integer(
                    ov[
                        "n_member_clusters"
                    ]
                ),

            "n_module_proteins":
                integer(
                    ov[
                        "n_module_proteins"
                    ]
                ),

            "n_genomes_represented":
                integer(
                    ov[
                        "n_genomes_with_any_module_protein"
                    ]
                ),

            "n_possible_cluster_pairs":
                integer(
                    arch.get(
                        "n_possible_cluster_pairs",
                        0,
                    )
                ),

            "n_pairs_copresent":
                integer(
                    arch.get(
                        "n_pairs_copresent",
                        0,
                    )
                ),

            "n_local_20kb_edges":
                integer(
                    arch.get(
                        "n_pairs_within_20kb",
                        0,
                    )
                ),

            "n_local_10kb_edges":
                integer(
                    arch.get(
                        "n_pairs_within_10kb",
                        0,
                    )
                ),

            "n_local_5kb_edges":
                integer(
                    arch.get(
                        "n_pairs_within_5kb",
                        0,
                    )
                ),

            "n_adjacency_edges":
                integer(
                    arch.get(
                        "n_pairs_with_adjacency",
                        0,
                    )
                ),

            "n_local_architecture_components":
                integer(
                    arch.get(
                        "n_local_architecture_components",
                        0,
                    )
                ),

            "largest_local_component_size":
                integer(
                    arch.get(
                        "largest_local_component_size",
                        0,
                    )
                ),

            "n_intervening_gene_rows":
                n_intervening_rows,

            "n_unique_intervening_proteins":
                n_unique_intervening,

            "n_findmehemes_positive_intervening_rows":
                int(
                    fmh_mask.sum()
                ),

            "n_unique_findmehemes_positive_intervening_proteins":
                (
                    mod_intervening.loc[
                        fmh_mask,
                        "_protein_key",
                    ]
                    .nunique()
                ),

            "n_fegenie_positive_intervening_rows":
                int(
                    fegenie_mask.sum()
                ),

            "n_unique_fegenie_positive_intervening_proteins":
                (
                    mod_intervening.loc[
                        fegenie_mask,
                        "_protein_key",
                    ]
                    .nunique()
                ),
        }
    )


report_summary = pd.DataFrame(
    summary_rows
)


report_summary.to_csv(
    OUT_SUMMARY,
    sep="\t",
    index=False,
)


## ================================================================== ##
## Render reports
## ================================================================== ##

print()
print("Rendering module reports...")


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


    module_summary = report_summary[
        report_summary[
            "module"
        ]
        ==
        module
    ].iloc[0]


    member = cluster_report[
        cluster_report[
            "module"
        ]
        ==
        module
    ].copy()


    module_nodes = nodes[
        nodes[
            "module"
        ]
        ==
        module
    ].copy()


    local_edges = edges[
        (
            edges[
                "module"
            ]
            ==
            module
        )
        &
        (
            pd.to_numeric(
                edges[
                    "has_local_20kb_architecture_edge"
                ],
                errors="coerce",
            )
            ==
            1
        )
    ].copy()


    module_mmseq = mmseq[
        mmseq[
            "focal_module"
        ]
        ==
        module
    ].copy()


    module_cogs = cogs[
        cogs[
            "focal_module"
        ]
        ==
        module
    ].copy()


    module_patterns = patterns[
        patterns[
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


    ## -------------------------------------------------------------- ##
    ## Main report
    ## -------------------------------------------------------------- ##

    lines = []


    lines.append(
        f"# {module}"
    )


    lines.append(
        ""
    )


    lines.append(
        "## Interpretation framework"
    )


    lines.append(
        ""
    )


    lines.append(
        "This is an MCL-defined genome-level co-occurrence module. "
        "The module is not assumed to represent a single operon, "
        "protein complex, or linear conserved gene cassette."
    )


    lines.append(
        ""
    )


    lines.append(
        "MMseqs2 clusters are the primary sequence-family units. "
        "FeGenie, FindMeHemes, SignalP, DeepTMHMM and GlobDB "
        "annotations are retained as independent evidence layers."
    )


    lines.append(
        ""
    )


    lines.append(
        "Percentages should always be interpreted together with "
        "their numerator and denominator; for example, 100% (1/1) "
        "is not equivalent in evidential strength to 100% (44/44)."
    )


    lines.append(
        ""
    )


    ## -------------------------------------------------------------- ##
    ## Overview
    ## -------------------------------------------------------------- ##

    lines.append(
        "## Module overview"
    )


    lines.append(
        ""
    )


    lines.append(
        f"- Member MMseqs2 families: "
        f"{module_summary['n_member_clusters']}"
    )


    lines.append(
        f"- Module proteins: "
        f"{module_summary['n_module_proteins']}"
    )


    lines.append(
        f"- Genomes represented: "
        f"{module_summary['n_genomes_represented']}"
    )


    lines.append(
        f"- Possible within-module family pairs: "
        f"{module_summary['n_possible_cluster_pairs']}"
    )


    lines.append(
        f"- Co-present family pairs: "
        f"{module_summary['n_pairs_copresent']}"
    )


    lines.append(
        f"- Family pairs with ≤20-kb evidence: "
        f"{module_summary['n_local_20kb_edges']}"
    )


    lines.append(
        f"- Local architecture components: "
        f"{module_summary['n_local_architecture_components']}"
    )


    lines.append(
        f"- Largest local component: "
        f"{module_summary['largest_local_component_size']} families"
    )


    lines.append(
        ""
    )


    ## -------------------------------------------------------------- ##
    ## Member families
    ## -------------------------------------------------------------- ##

    lines.append(
        "## Member families"
    )


    lines.append(
        ""
    )


    member_columns = [
        "cluster",
        "n_proteins",
        "n_genomes",
        "mean_heme_count",
        "dominant_fegenie_HMM",
        "report_fegenie_HMM_distribution",
        "dominant_report_topology",
        "dominant_report_topology_support",
        "report_topology_distribution",
        "report_signalp_distribution",
        "report_deeptmhmm_class_distribution",
    ]


    lines.append(
        markdown_table(
            member,
            member_columns,
        )
    )


    lines.append(
        ""
    )


    lines.append(
        "The `dominant_report_topology` field is derived from the raw "
        "predictors for reporting. In particular, DeepTMHMM `BETA` "
        "is interpreted as beta-barrel membrane topology. The original "
        "upstream `localization_class` remains unchanged."
    )


    lines.append(
        ""
    )


    ## -------------------------------------------------------------- ##
    ## Local architecture components
    ## -------------------------------------------------------------- ##

    lines.append(
        "## Supported local architecture"
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
        "### Local components"
    )


    lines.append(
        ""
    )


    lines.append(
        markdown_table(
            module_nodes,
            component_columns,
        )
    )


    lines.append(
        ""
    )


    lines.append(
        "### Local family-pair relationships"
    )


    lines.append(
        ""
    )


    if len(
        local_edges
    ) == 0:

        lines.append(
            "_No within-module family pair was observed within ≤20 kb._"
        )


    else:

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
            "local_vs_consensus_orientation_agreement",
            "local_vs_consensus_order_agreement",
            "median_local_intergenic_gap_bp",
            "median_local_genes_between",
        ]


        local_edges = local_edges.sort_values(
            [
                "n_genomes_within_20kb",
                "n_genomes_both",
                "cluster_a",
                "cluster_b",
            ],
            ascending=[
                False,
                False,
                True,
                True,
            ],
            kind="stable",
        )


        lines.append(
            markdown_table(
                local_edges,
                architecture_columns,
            )
        )


    lines.append(
        ""
    )


    ## -------------------------------------------------------------- ##
    ## Focal-family MMseq consensus
    ## -------------------------------------------------------------- ##

    lines.append(
        "## Focal-family sequence-neighborhood consensus"
    )


    lines.append(
        ""
    )


    lines.append(
        "The following is a display-ranked subset of the full Stage-13D "
        "MMseqs2 consensus table. Ranking is for report readability only "
        "and is not an analytical filtering threshold."
    )


    lines.append(
        ""
    )


    same_module_mmseq = module_mmseq.copy()


    if (
        "neighbor_same_mcl_module_as_focal"
        in same_module_mmseq.columns
    ):

        same_module_mmseq = same_module_mmseq[
            pd.to_numeric(
                same_module_mmseq[
                    "neighbor_same_mcl_module_as_focal"
                ],
                errors="coerce",
            )
            ==
            1
        ]


    if (
        "direction_has_gene_window_hits"
        in same_module_mmseq.columns
    ):

        same_module_mmseq = same_module_mmseq[
            pd.to_numeric(
                same_module_mmseq[
                    "direction_has_gene_window_hits"
                ],
                errors="coerce",
            )
            ==
            1
        ]


    mmseq_display_parts = []


    for focal_cluster, group in same_module_mmseq.groupby(
        "focal_cluster",
        sort=True,
    ):

        group = group.copy()


        group[
            "_rank_genomes"
        ] = pd.to_numeric(
            group[
                "genome_saturation_n_positive"
            ],
            errors="coerce",
        ).fillna(
            0
        )


        group[
            "_rank_occurrences"
        ] = pd.to_numeric(
            group[
                "occurrence_saturation_n_positive"
            ],
            errors="coerce",
        ).fillna(
            0
        )


        group = group.sort_values(
            [
                "_rank_genomes",
                "_rank_occurrences",
                "neighbor_cluster",
            ],
            ascending=[
                False,
                False,
                True,
            ],
            kind="stable",
        ).head(
            MAX_MMSEQ_RELATIONSHIPS_PER_FOCAL
        )


        mmseq_display_parts.append(
            group
        )


    if len(
        mmseq_display_parts
    ) > 0:

        mmseq_display = pd.concat(
            mmseq_display_parts,
            ignore_index=True,
        )


        mmseq_columns = [
            "focal_cluster",
            "neighbor_cluster",
            "direction",
            "best_exact_gene_offset",
            "best_exact_occurrence_support",
            "occurrence_saturation_gene_radius",
            "occurrence_saturation_support",
            "genome_saturation_gene_radius",
            "genome_saturation_support",
            "physical_5kb_occurrence_support",
            "physical_10kb_occurrence_support",
            "physical_20kb_occurrence_support",
        ]


        lines.append(
            markdown_table(
                mmseq_display,
                mmseq_columns,
            )
        )


    else:

        lines.append(
            "_No same-module MMseqs2 neighborhood relationships "
            "with directional hits._"
        )


    lines.append(
        ""
    )


    ## -------------------------------------------------------------- ##
    ## COG context
    ## -------------------------------------------------------------- ##

    lines.append(
        "## Functional annotation context"
    )


    lines.append(
        ""
    )


    lines.append(
        "COG/product evidence is a secondary annotation layer. "
        "A missing COG is not interpreted as biological absence."
    )


    lines.append(
        ""
    )


    cog_hits = module_cogs.copy()


    if (
        "direction_has_gene_window_hits"
        in cog_hits.columns
    ):

        cog_hits = cog_hits[
            pd.to_numeric(
                cog_hits[
                    "direction_has_gene_window_hits"
                ],
                errors="coerce",
            )
            ==
            1
        ]


    cog_display_parts = []


    for focal_cluster, group in cog_hits.groupby(
        "focal_cluster",
        sort=True,
    ):

        group = group.copy()


        group[
            "_rank_genomes"
        ] = pd.to_numeric(
            group[
                "genome_saturation_n_positive"
            ],
            errors="coerce",
        ).fillna(
            0
        )


        group[
            "_rank_occurrences"
        ] = pd.to_numeric(
            group[
                "occurrence_saturation_n_positive"
            ],
            errors="coerce",
        ).fillna(
            0
        )


        group = group.sort_values(
            [
                "_rank_genomes",
                "_rank_occurrences",
                "neighbor_cog",
            ],
            ascending=[
                False,
                False,
                True,
            ],
            kind="stable",
        ).head(
            MAX_COG_RELATIONSHIPS_PER_FOCAL
        )


        cog_display_parts.append(
            group
        )


    if len(
        cog_display_parts
    ) > 0:

        cog_display = pd.concat(
            cog_display_parts,
            ignore_index=True,
        )


        cog_columns = [
            "focal_cluster",
            "neighbor_cog",
            "direction",
            "best_exact_gene_offset",
            "best_exact_occurrence_support",
            "occurrence_saturation_gene_radius",
            "occurrence_saturation_support",
            "genome_saturation_support",
            "dominant_product_at_occurrence_saturation_radius",
        ]


        lines.append(
            markdown_table(
                cog_display,
                cog_columns,
            )
        )


    else:

        lines.append(
            "_No COG neighborhood relationships with directional hits._"
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


    lines.append(
        f"- Intervening-gene observations: "
        f"{module_summary['n_intervening_gene_rows']}"
    )


    lines.append(
        f"- Unique intervening proteins: "
        f"{module_summary['n_unique_intervening_proteins']}"
    )


    lines.append(
        f"- FindMeHemes-positive observations: "
        f"{module_summary['n_findmehemes_positive_intervening_rows']}"
    )


    lines.append(
        f"- Unique FindMeHemes-positive proteins: "
        f"{module_summary['n_unique_findmehemes_positive_intervening_proteins']}"
    )


    lines.append(
        f"- FeGenie-positive observations: "
        f"{module_summary['n_fegenie_positive_intervening_rows']}"
    )


    lines.append(
        f"- Unique FeGenie-positive proteins: "
        f"{module_summary['n_unique_fegenie_positive_intervening_proteins']}"
    )


    lines.append(
        ""
    )


    if len(
        module_patterns
    ) > 0:

        module_patterns = module_patterns.copy()


        if (
            "n_local_observations_with_pattern"
            in module_patterns.columns
        ):

            module_patterns[
                "_pattern_rank"
            ] = pd.to_numeric(
                module_patterns[
                    "n_local_observations_with_pattern"
                ],
                errors="coerce",
            ).fillna(
                0
            )


            module_patterns = module_patterns.sort_values(
                [
                    "_pattern_rank",
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


        pattern_columns = [
            "cluster_a",
            "cluster_b",
            "n_intervening_genes",
            "intervening_functional_pattern",
            "intervening_cog_pattern",
            "n_local_observations_with_pattern",
            "support_among_local_nonadjacent",
            "support_among_all_local",
        ]


        lines.append(
            "### Most frequently observed local intervening patterns"
        )


        lines.append(
            ""
        )


        lines.append(
            markdown_table(
                module_patterns.head(
                    MAX_INTERVENING_PATTERNS
                ),
                pattern_columns,
            )
        )


    else:

        lines.append(
            "_No non-adjacent local pair occurrences were present._"
        )


    lines.append(
        ""
    )


    ## -------------------------------------------------------------- ##
    ## FeGenie-positive intervening genes
    ## -------------------------------------------------------------- ##

    if (
        len(
            module_intervening
        )
        >
        0
        and
        "evidence_fegenie_positive"
        in module_intervening.columns
    ):

        fe_mask = (
            pd.to_numeric(
                module_intervening[
                    "evidence_fegenie_positive"
                ].replace(
                    "",
                    np.nan,
                ),
                errors="coerce",
            )
            ==
            1
        )


        fe_intervening = module_intervening[
            fe_mask
        ].copy()


        if len(
            fe_intervening
        ) > 0:

            lines.append(
                "### FeGenie-positive intervening proteins"
            )


            lines.append(
                ""
            )


            fe_columns = [
                "genome",
                "cluster_a",
                "cluster_b",
                "intervening_protein_id",
                "evidence_fegenie_HMMs",
                "evidence_fegenie_categories",
                "evidence_number_of_hemes",
                "evidence_signalp_prediction",
                "evidence_deeptmhmm_class",
                "evidence_globdb_cog",
                "evidence_globdb_product",
            ]


            lines.append(
                markdown_table(
                    fe_intervening,
                    fe_columns,
                )
            )


            lines.append(
                ""
            )


    ## -------------------------------------------------------------- ##
    ## Report caveats
    ## -------------------------------------------------------------- ##

    lines.append(
        "## Report caveats"
    )


    lines.append(
        ""
    )


    lines.append(
        "- MCL membership describes genome-level co-occurrence, "
        "not necessarily physical linkage."
    )


    lines.append(
        "- Local architecture edges require observed ≤20-kb proximity "
        "but no arbitrary recurrence percentage threshold."
    )


    lines.append(
        "- Censoring from contig boundaries is retained in the Stage-13 "
        "informative denominators."
    )


    lines.append(
        "- Blank FeGenie/COG annotations do not imply that a protein "
        "lacks biological importance."
    )


    lines.append(
        "- The complete analytical TSVs should be used whenever the "
        "display-ranked Markdown tables are insufficient."
    )


    lines.append(
        ""
    )


    report_path.write_text(
        "\n".join(
            lines
        )
    )


    ## -------------------------------------------------------------- ##
    ## Write compact derived tables beside report
    ## -------------------------------------------------------------- ##

    member.to_csv(
        outdir
        /
        "member_family_report.tsv",
        sep="\t",
        index=False,
        na_rep="",
        float_format="%.6f",
    )


    local_edges.to_csv(
        outdir
        /
        "local_architecture_report.tsv",
        sep="\t",
        index=False,
        na_rep="",
        float_format="%.6f",
    )


    module_intervening.drop(
        columns=[
            "_protein_key",
        ],
        errors="ignore",
    ).to_csv(
        outdir
        /
        "intervening_gene_report.tsv",
        sep="\t",
        index=False,
        na_rep="",
        float_format="%.6f",
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
                module_summary[
                    "n_member_clusters"
                ],

            "n_genomes_represented":
                module_summary[
                    "n_genomes_represented"
                ],

            "n_local_20kb_edges":
                module_summary[
                    "n_local_20kb_edges"
                ],

            "n_local_architecture_components":
                module_summary[
                    "n_local_architecture_components"
                ],

            "largest_local_component_size":
                module_summary[
                    "largest_local_component_size"
                ],

            "n_unique_intervening_proteins":
                module_summary[
                    "n_unique_intervening_proteins"
                ],
        }
    )


## ================================================================== ##
## Index / QC
## ================================================================== ##

index = pd.DataFrame(
    index_rows
)


index.to_csv(
    OUT_INDEX,
    sep="\t",
    index=False,
)


if len(
    index
) != EXPECTED_MODULES:

    fail(
        "Did not render exactly 35 module reports."
    )


missing_reports = [
    module

    for module in modules

    if not (
        MODULE_DIR
        /
        module
        /
        "module_report.md"
    ).exists()
]


if len(
    missing_reports
) > 0:

    fail(
        "Missing rendered reports for: "
        +
        ", ".join(
            missing_reports
        )
    )


qc = pd.DataFrame(
    [
        [
            "modules_rendered",
            len(
                index
            ),
            EXPECTED_MODULES,
            int(
                len(
                    index
                )
                ==
                EXPECTED_MODULES
            ),
        ],

        [
            "member_clusters",
            len(
                cluster_report
            ),
            EXPECTED_CLUSTERS,
            int(
                len(
                    cluster_report
                )
                ==
                EXPECTED_CLUSTERS
            ),
        ],

        [
            "module_proteins",
            len(
                proteins
            ),
            EXPECTED_MODULE_PROTEINS,
            int(
                len(
                    proteins
                )
                ==
                EXPECTED_MODULE_PROTEINS
            ),
        ],

        [
            "architecture_pairs",
            len(
                edges
            ),
            EXPECTED_EDGES,
            int(
                len(
                    edges
                )
                ==
                EXPECTED_EDGES
            ),
        ],

        [
            "local_20kb_edges",
            int(
                pd.to_numeric(
                    edges[
                        "has_local_20kb_architecture_edge"
                    ],
                    errors="coerce",
                )
                .sum()
            ),
            EXPECTED_LOCAL_EDGES,
            int(
                pd.to_numeric(
                    edges[
                        "has_local_20kb_architecture_edge"
                    ],
                    errors="coerce",
                )
                .sum()
                ==
                EXPECTED_LOCAL_EDGES
            ),
        ],

        [
            "intervening_gene_rows",
            len(
                intervening
            ),
            EXPECTED_INTERVENING_ROWS,
            int(
                len(
                    intervening
                )
                ==
                EXPECTED_INTERVENING_ROWS
            ),
        ],

        [
            "beta_topology_qc_pass",
            int(
                (
                    beta[
                        "report_topology_class"
                    ]
                    ==
                    "beta_barrel_membrane"
                )
                .all()
            ),
            1,
            int(
                (
                    beta[
                        "report_topology_class"
                    ]
                    ==
                    "beta_barrel_membrane"
                )
                .all()
            ),
        ],
    ],
    columns=[
        "metric",
        "observed",
        "expected",
        "pass",
    ],
)


qc.to_csv(
    OUT_QC,
    sep="\t",
    index=False,
)


## ================================================================== ##
## Terminal summary
## ================================================================== ##

print()
print("Rendered module-report summary")


print(
    f"  Reports rendered:               "
    f"{len(index):,}"
)

print(
    f"  Member-family report rows:      "
    f"{len(cluster_report):,}"
)

print(
    f"  Local architecture edges:       "
    f"{EXPECTED_LOCAL_EDGES:,}"
)

print(
    f"  DeepTMHMM BETA proteins:        "
    f"{len(beta):,}"
)

print(
    f"  BETA called beta-barrel:        "
    f"{(beta['report_topology_class'] == 'beta_barrel_membrane').sum():,}"
)

print(
    f"  Intervening-gene observations:  "
    f"{len(intervening):,}"
)

print(
    f"  Unique intervening proteins:    "
    f"{intervening['_protein_key'].nunique():,}"
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
        f"{spotlight} rendered-report spotlight"
    )


    s = report_summary[
        report_summary[
            "module"
        ]
        ==
        spotlight
    ].iloc[0]


    print(
        f"  Member families:          "
        f"{s['n_member_clusters']}"
    )

    print(
        f"  Genomes represented:      "
        f"{s['n_genomes_represented']}"
    )

    print(
        f"  Local <=20-kb edges:      "
        f"{s['n_local_20kb_edges']}"
    )

    print(
        f"  Local components:         "
        f"{s['n_local_architecture_components']}"
    )

    print(
        f"  Intervening observations: "
        f"{s['n_intervening_gene_rows']}"
    )

    print(
        f"  Unique intervening genes: "
        f"{s['n_unique_intervening_proteins']}"
    )


    spot_member = cluster_report[
        cluster_report[
            "module"
        ]
        ==
        spotlight
    ]


    display_columns = [
        column

        for column in [
            "cluster",
            "n_proteins",
            "n_genomes",
            "dominant_fegenie_HMM",
            "report_fegenie_HMM_distribution",
            "dominant_report_topology",
            "dominant_report_topology_support",
            "report_topology_distribution",
        ]

        if column
        in spot_member.columns
    ]


    print()
    print(
        spot_member[
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
    f"Global report summary:  "
    f"{OUT_SUMMARY}"
)

print(
    f"Cluster report evidence:"
    f" {OUT_CLUSTER_REPORT}"
)

print(
    f"Report index:           "
    f"{OUT_INDEX}"
)

print(
    f"QC:                     "
    f"{OUT_QC}"
)

print(
    f"Rendered reports:       "
    f"{MODULE_DIR}"
)
