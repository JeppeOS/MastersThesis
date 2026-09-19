#!/usr/bin/env python3

from pathlib import Path
from collections import Counter
import sys

import numpy as np
import pandas as pd


## ================================================================== ##
## STAGE 15A - GENE-LEVEL RESOLUTION RECOVERY
##
## PURPOSE
## -------
##
## Recover the individual-gene resolution that was deliberately
## compressed during MMseqs2 clustering / network / MCL analysis.
##
##
## CORE PRINCIPLE
## --------------
##
## Across genomes:
##     MMseqs2 family / MCL module is the analytical unit.
##
## Within an individual genome:
##     the ACTUAL Prodigal gene/protein is the analytical unit.
##
##
## AUTOMATIC FOCAL GENES
## ---------------------
##
## Within the ORIGINAL Stage-12 module-associated regions:
##
##   1. every MCL/module protein
##   2. every FeGenie-positive gene
##   3. every FindMeHemes-positive gene
##
## becomes an individual focal gene.
##
##
## FOCAL ELIGIBILITY IS FIXED
## --------------------------
##
## A fresh focal neighborhood may extend outside the original
## module-associated region, but newly exposed genes outside that
## original region DO NOT automatically become new focal genes.
##
## This prevents recursive expansion across whole contigs.
##
##
## FRESH NEIGHBORHOOD EXTRACTION
## -----------------------------
##
## Every selected focal gene is re-centered against the complete
## Stage-11 Prodigal gene catalogue.
##
## A surrounding gene is retained if:
##
##   abs(gene-rank difference) <= 20
##
## OR
##
##   its CDS intersects:
##
##       focal CDS +/- 20,000 bp
##
##
## Every surrounding gene is retained, regardless of annotation.
##
##
## OUTPUT PHILOSOPHY
## -----------------
##
## Cluster/module identity becomes a PROPERTY of an individual gene.
##
## Thus two paralogues in the same MMseqs2 family remain two distinct
## focal genes.
##
## No FeGenie / FindMeHemes / SignalP / DeepTMHMM reruns occur here.
## Existing integrated evidence is simply recovered from the complete
## master gene catalogue.
## ================================================================== ##


## ================================================================== ##
## Paths
## ================================================================== ##

HERE = Path(__file__).resolve().parent
WORKFLOW = HERE.parent.parent


## Complete 2,002,656-gene Prodigal catalogue + annotations. ##

GENE_CATALOG = (
    WORKFLOW
    / "11_genome_annotations"
    / "annotated_gene_catalog.tsv"
)


## Original Stage-12 module-centered neighborhoods. ##

OBSERVED_NEIGHBORHOODS = (
    WORKFLOW
    / "12_neighborhood_extraction"
    / "observed_neighborhood_genes.tsv"
)


## Authoritative actual protein -> MMseq cluster -> MCL module mapping. ##

MODULE_PROTEINS = (
    WORKFLOW
    / "09_module_occurrence"
    / "protein_module_membership.tsv"
)


## ================================================================== ##
## Outputs
## ================================================================== ##

OUT_REGIONS = (
    HERE
    / "fixed_module_eligibility_regions.tsv"
)


OUT_REGION_GENES = (
    HERE
    / "module_region_genes.tsv"
)


OUT_FOCALS = (
    HERE
    / "focal_gene_catalog.tsv"
)


OUT_FOCALS_UNIQUE = (
    HERE
    / "focal_gene_catalog_unique.tsv"
)


OUT_NEIGHBORHOODS = (
    HERE
    / "focal_gene_neighborhoods.tsv"
)


OUT_REASON_SUMMARY = (
    HERE
    / "focal_gene_reason_summary.tsv"
)


OUT_QC = (
    HERE
    / "resolution_recovery_qc.tsv"
)


## ================================================================== ##
## Extraction parameters
## ================================================================== ##

GENE_RADIUS = 20
BP_RADIUS = 20_000


## ================================================================== ##
## Locked upstream expectations
## ================================================================== ##

EXPECTED_GENES = 2_002_656
EXPECTED_MODULE_PROTEINS = 10_537
EXPECTED_MODULES = 35

EXPECTED_STAGE12_ROWS = 253_261
EXPECTED_STAGE12_FOCALS = 10_537


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
            f"Required file does not exist:\n{path}"
        )


def clean(value):

    if pd.isna(value):

        return ""

    return str(value).strip()


def numeric(series):

    return pd.to_numeric(
        series.replace(
            "",
            np.nan,
        ),
        errors="coerce",
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


def resolve_column(
    columns,
    aliases,
    label,
    required=True,
):

    columns = list(columns)

    lower = {
        column.lower():
            column

        for column
        in columns
    }


    ## Exact alias first. ##

    for alias in aliases:

        if alias.lower() in lower:

            return lower[
                alias.lower()
            ]


    ## Case-insensitive normalized comparison. ##

    def normalize(value):

        return (
            value.lower()
            .replace("-", "_")
            .replace(" ", "_")
        )


    normalized = {
        normalize(column):
            column

        for column
        in columns
    }


    for alias in aliases:

        key = normalize(alias)

        if key in normalized:

            return normalized[key]


    if required:

        fail(
            f"Could not identify required column '{label}'.\n"
            f"Tried aliases: {aliases}\n"
            f"Available columns:\n"
            +
            "\n".join(
                columns
            )
        )


    return None


def flag_numeric(
    series,
):

    return (
        pd.to_numeric(
            series.replace(
                "",
                np.nan,
            ),
            errors="coerce",
        )
        ==
        1
    )


def distribution_string(values):

    cleaned = [
        clean(value)

        for value
        in values

        if clean(value) != ""
    ]


    if len(cleaned) == 0:

        return ""


    counts = Counter(cleaned)


    return "; ".join(
        f"{value}:{count}"

        for value, count
        in sorted(
            counts.items(),
            key=lambda item: (
                -item[1],
                item[0],
            ),
        )
    )


## ================================================================== ##
## Report-level topology interpretation
##
## This DOES NOT overwrite upstream localization_class.
## ================================================================== ##

def report_topology(
    deeptmhmm,
    signalp,
    n_tm,
):

    deeptmhmm = clean(deeptmhmm)
    signalp = clean(signalp)


    try:

        n_tm = float(n_tm)

    except Exception:

        n_tm = np.nan


    if deeptmhmm == "BETA":

        return "beta_barrel_membrane"


    if signalp in [
        "LIPO",
        "TATLIPO",
    ]:

        return "exported_lipoprotein"


    if deeptmhmm == "SP+TM":

        if not np.isnan(n_tm):

            if n_tm == 1:

                return "exported_single_pass_membrane"


            if n_tm > 1:

                return "exported_multipass_membrane"


        return "exported_membrane_SP_plus_TM"


    if deeptmhmm == "TM":

        if not np.isnan(n_tm):

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
## Merge overlapping original module-associated focal spans.
##
## Input spans are gene-rank intervals produced from the ORIGINAL
## Stage-12 module focal neighborhoods.
##
## Overlapping or directly touching rank intervals are merged.
## ================================================================== ##

def merge_region_spans(
    spans,
):

    region_rows = []


    region_counter = 0


    for (
        genome,
        module,
        contig,
    ), group in spans.groupby(
        [
            "genome",
            "module",
            "contig",
        ],
        sort=True,
    ):

        group = group.sort_values(
            [
                "region_rank_start",
                "region_rank_end",
            ],
            kind="stable",
        )


        current = None


        for row in group.itertuples(
            index=False
        ):

            rank_start = int(
                row.region_rank_start
            )

            rank_end = int(
                row.region_rank_end
            )


            bp_start = int(
                row.region_bp_start
            )

            bp_end = int(
                row.region_bp_end
            )


            focal_id = row.source_focal_protein_id


            if current is None:

                current = {
                    "genome":
                        genome,

                    "module":
                        module,

                    "contig":
                        contig,

                    "region_rank_start":
                        rank_start,

                    "region_rank_end":
                        rank_end,

                    "region_bp_start":
                        bp_start,

                    "region_bp_end":
                        bp_end,

                    "source_focals":
                        [
                            focal_id
                        ],
                }


                continue


            ## Merge overlapping / directly touching gene-rank spans. ##

            if (
                rank_start
                <=
                current[
                    "region_rank_end"
                ]
                +
                1
            ):

                current[
                    "region_rank_end"
                ] = max(
                    current[
                        "region_rank_end"
                    ],
                    rank_end,
                )


                current[
                    "region_bp_start"
                ] = min(
                    current[
                        "region_bp_start"
                    ],
                    bp_start,
                )


                current[
                    "region_bp_end"
                ] = max(
                    current[
                        "region_bp_end"
                    ],
                    bp_end,
                )


                current[
                    "source_focals"
                ].append(
                    focal_id
                )


            else:

                region_counter += 1


                current[
                    "region_id"
                ] = (
                    f"Region_{region_counter:06d}"
                )


                region_rows.append(
                    current
                )


                current = {
                    "genome":
                        genome,

                    "module":
                        module,

                    "contig":
                        contig,

                    "region_rank_start":
                        rank_start,

                    "region_rank_end":
                        rank_end,

                    "region_bp_start":
                        bp_start,

                    "region_bp_end":
                        bp_end,

                    "source_focals":
                        [
                            focal_id
                        ],
                }


        if current is not None:

            region_counter += 1


            current[
                "region_id"
            ] = (
                f"Region_{region_counter:06d}"
            )


            region_rows.append(
                current
            )


    out = pd.DataFrame(
        region_rows
    )


    out[
        "n_original_module_focals"
    ] = out[
        "source_focals"
    ].map(
        lambda values:
            len(
                set(values)
            )
    )


    out[
        "source_focal_protein_ids"
    ] = out[
        "source_focals"
    ].map(
        lambda values:
            ";".join(
                sorted(
                    set(values)
                )
            )
    )


    return out.drop(
        columns=[
            "source_focals",
        ]
    )


## ================================================================== ##
## Start
## ================================================================== ##

print("=" * 80)
print("STAGE 15A - GENE-LEVEL RESOLUTION RECOVERY")
print("=" * 80)


## ================================================================== ##
## 1. Read authoritative module proteins
## ================================================================== ##

print()
print("Reading authoritative module protein membership...")


module_proteins = pd.read_csv(
    MODULE_PROTEINS,
    sep="\t",
    dtype=str,
    keep_default_na=False,
)


required_module_columns = [
    "genome",
    "protein_id",
    "cluster",
    "module",
]


missing = (
    set(
        required_module_columns
    )
    -
    set(
        module_proteins.columns
    )
)


if missing:

    fail(
        "protein_module_membership.tsv missing: "
        +
        ", ".join(
            sorted(missing)
        )
    )


if len(
    module_proteins
) != EXPECTED_MODULE_PROTEINS:

    fail(
        f"Expected {EXPECTED_MODULE_PROTEINS:,} module proteins; "
        f"found {len(module_proteins):,}."
    )


if module_proteins[
    "module"
].nunique() != EXPECTED_MODULES:

    fail(
        f"Expected {EXPECTED_MODULES} modules; "
        f"found {module_proteins['module'].nunique()}."
    )


if module_proteins[
    [
        "genome",
        "protein_id",
    ]
].duplicated().any():

    fail(
        "Duplicate genome + protein_id keys "
        "in module protein membership."
    )


module_lookup = (
    module_proteins
    .set_index(
        [
            "genome",
            "protein_id",
        ]
    )[
        [
            "cluster",
            "module",
        ]
    ]
)


print(
    f"  Module proteins: "
    f"{len(module_proteins):,}"
)

print(
    f"  Modules:         "
    f"{module_proteins['module'].nunique():,}"
)


## ================================================================== ##
## 2. Read original Stage-12 neighborhoods
##
## These define the FIXED region in which automatic focal eligibility
## is assessed.
## ================================================================== ##

print()
print("Reading original Stage-12 module neighborhoods...")


stage12 = pd.read_csv(
    OBSERVED_NEIGHBORHOODS,
    sep="\t",
    dtype=str,
    keep_default_na=False,
)


required_stage12 = [
    "genome",

    "focal_protein_id",
    "focal_cluster",
    "focal_module",

    "focal_contig",

    "neighbor_protein_id",
    "neighbor_start",
    "neighbor_end",
    "neighbor_gene_rank",
]


missing = (
    set(
        required_stage12
    )
    -
    set(
        stage12.columns
    )
)


if missing:

    fail(
        "observed_neighborhood_genes.tsv missing: "
        +
        ", ".join(
            sorted(missing)
        )
    )


if len(
    stage12
) != EXPECTED_STAGE12_ROWS:

    fail(
        f"Expected {EXPECTED_STAGE12_ROWS:,} Stage-12 rows; "
        f"found {len(stage12):,}."
    )


if stage12[
    "focal_protein_id"
].nunique() != EXPECTED_STAGE12_FOCALS:

    fail(
        f"Expected {EXPECTED_STAGE12_FOCALS:,} "
        f"original module focal proteins; found "
        f"{stage12['focal_protein_id'].nunique():,}."
    )


stage12[
    "_neighbor_start"
] = pd.to_numeric(
    stage12[
        "neighbor_start"
    ],
    errors="raise",
).astype(int)


stage12[
    "_neighbor_end"
] = pd.to_numeric(
    stage12[
        "neighbor_end"
    ],
    errors="raise",
).astype(int)


stage12[
    "_neighbor_rank"
] = pd.to_numeric(
    stage12[
        "neighbor_gene_rank"
    ],
    errors="raise",
).astype(int)


print(
    f"  Stage-12 rows:       "
    f"{len(stage12):,}"
)

print(
    f"  Original focals:     "
    f"{stage12['focal_protein_id'].nunique():,}"
)

print(
    f"  Focal contigs:       "
    f"{stage12[['genome','focal_contig']].drop_duplicates().shape[0]:,}"
)


## ================================================================== ##
## 3. Build original focal spans, then merge into fixed module regions
## ================================================================== ##

print()
print("Constructing fixed module-associated eligibility regions...")


focal_spans = (
    stage12
    .groupby(
        [
            "genome",
            "focal_module",
            "focal_contig",
            "focal_protein_id",
        ],
        as_index=False,
    )
    .agg(
        region_rank_start=(
            "_neighbor_rank",
            "min",
        ),

        region_rank_end=(
            "_neighbor_rank",
            "max",
        ),

        region_bp_start=(
            "_neighbor_start",
            "min",
        ),

        region_bp_end=(
            "_neighbor_end",
            "max",
        ),
    )
    .rename(
        columns={
            "focal_module":
                "module",

            "focal_contig":
                "contig",

            "focal_protein_id":
                "source_focal_protein_id",
        }
    )
)


regions = merge_region_spans(
    focal_spans
)


regions = regions[
    [
        "region_id",

        "genome",
        "module",
        "contig",

        "region_rank_start",
        "region_rank_end",

        "region_bp_start",
        "region_bp_end",

        "n_original_module_focals",
        "source_focal_protein_ids",
    ]
]


write_tsv(
    regions,
    OUT_REGIONS,
)


print(
    f"  Fixed eligibility regions: "
    f"{len(regions):,}"
)


## ================================================================== ##
## 4. Build unique gene x module eligibility set directly from
##    original Stage-12 rows
##
## This is the exact union of genes actually contained in the original
## module-centered extraction.
## ================================================================== ##

print()
print("Building fixed module-region gene set...")


region_gene_keys = (
    stage12[
        [
            "genome",
            "focal_module",
            "focal_contig",
            "neighbor_protein_id",
            "_neighbor_rank",
        ]
    ]
    .rename(
        columns={
            "focal_module":
                "module",

            "focal_contig":
                "contig",

            "neighbor_protein_id":
                "protein_id",

            "_neighbor_rank":
                "gene_rank",
        }
    )
    .drop_duplicates(
        [
            "genome",
            "module",
            "protein_id",
        ]
    )
    .reset_index(
        drop=True
    )
)


## Assign each eligible gene to its merged fixed region. ##

region_assignment = []


region_dict = {}


for (
    genome,
    module,
    contig,
), group in regions.groupby(
    [
        "genome",
        "module",
        "contig",
    ],
    sort=False,
):

    region_dict[
        (
            genome,
            module,
            contig,
        )
    ] = group.sort_values(
        "region_rank_start"
    ).to_dict(
        orient="records"
    )


for row in region_gene_keys.itertuples(
    index=False
):

    candidates = region_dict.get(
        (
            row.genome,
            row.module,
            row.contig,
        ),
        [],
    )


    matches = [
        region

        for region
        in candidates

        if (
            int(
                row.gene_rank
            )
            >=
            int(
                region[
                    "region_rank_start"
                ]
            )
            and
            int(
                row.gene_rank
            )
            <=
            int(
                region[
                    "region_rank_end"
                ]
            )
        )
    ]


    if len(matches) != 1:

        fail(
            "Could not assign exactly one fixed eligibility region to "
            f"{row.genome} / {row.module} / {row.protein_id}. "
            f"Matches={len(matches)}"
        )


    region_assignment.append(
        matches[0][
            "region_id"
        ]
    )


region_gene_keys[
    "region_id"
] = region_assignment


print(
    f"  Gene x module region entries: "
    f"{len(region_gene_keys):,}"
)

print(
    f"  Unique actual genes in regions: "
    f"{region_gene_keys[['genome','protein_id']].drop_duplicates().shape[0]:,}"
)


## ================================================================== ##
## 5. Inspect complete gene catalogue header and resolve columns
## ================================================================== ##

print()
print("Resolving complete Stage-11 gene-catalogue schema...")


require_file(
    GENE_CATALOG
)


catalog_header = pd.read_csv(
    GENE_CATALOG,
    sep="\t",
    nrows=0,
)


columns = catalog_header.columns


COL_GENOME = resolve_column(
    columns,
    [
        "genome",
    ],
    "genome",
)


COL_PROTEIN = resolve_column(
    columns,
    [
        "protein_id",
        "prodigal_protein_id",
    ],
    "protein_id",
)


COL_CONTIG = resolve_column(
    columns,
    [
        "contig",
        "prodigal_contig",
    ],
    "contig",
)


COL_START = resolve_column(
    columns,
    [
        "start",
        "prodigal_start",
        "gene_start",
    ],
    "start",
)


COL_END = resolve_column(
    columns,
    [
        "end",
        "prodigal_end",
        "gene_end",
    ],
    "end",
)


COL_STRAND = resolve_column(
    columns,
    [
        "strand",
        "prodigal_strand",
    ],
    "strand",
)


COL_RANK = resolve_column(
    columns,
    [
        "contig_gene_rank",
        "gene_rank",
    ],
    "contig_gene_rank",
)


## Optional integrated evidence columns. ##

COL_CLUSTER = resolve_column(
    columns,
    [
        "cluster",
        "mmseq_cluster",
        "mmseqs_cluster",
    ],
    "MMseq cluster",
    required=False,
)


COL_MODULE = resolve_column(
    columns,
    [
        "module",
        "mcl_module",
    ],
    "MCL module",
    required=False,
)


COL_FE_POS = resolve_column(
    columns,
    [
        "fegenie_positive",
        "neighbor_fegenie_positive",
    ],
    "FeGenie positive",
    required=False,
)


COL_FE_HMM = resolve_column(
    columns,
    [
        "fegenie_HMMs",
        "fegenie_hmms",
        "fegenie_HMM",
        "fegenie_hmm",
    ],
    "FeGenie HMM",
    required=False,
)


COL_FE_CAT = resolve_column(
    columns,
    [
        "fegenie_categories",
        "fegenie_category",
    ],
    "FeGenie category",
    required=False,
)


COL_FMH = resolve_column(
    columns,
    [
        "findmehemes_positive",
    ],
    "FindMeHemes positive",
    required=False,
)


COL_HEMES = resolve_column(
    columns,
    [
        "number_of_hemes",
        "heme_count",
        "CXXCH_count",
    ],
    "heme count",
    required=False,
)


COL_SIGNALP = resolve_column(
    columns,
    [
        "signalp_prediction",
        "prediction",
    ],
    "SignalP prediction",
    required=False,
)


COL_DEEPTMHMM = resolve_column(
    columns,
    [
        "deeptmhmm_class",
    ],
    "DeepTMHMM class",
    required=False,
)


COL_NTM = resolve_column(
    columns,
    [
        "deeptmhmm_n_tm_helices",
        "n_tm_helices",
    ],
    "DeepTMHMM TM count",
    required=False,
)


COL_LOCALIZATION = resolve_column(
    columns,
    [
        "localization_class",
    ],
    "localization class",
    required=False,
)


COL_COG = resolve_column(
    columns,
    [
        "globdb_cog",
    ],
    "GlobDB COG",
    required=False,
)


COL_GENE = resolve_column(
    columns,
    [
        "globdb_gene",
    ],
    "GlobDB gene",
    required=False,
)


COL_PRODUCT = resolve_column(
    columns,
    [
        "globdb_product",
    ],
    "GlobDB product",
    required=False,
)


COL_MATCH = resolve_column(
    columns,
    [
        "annotation_match_type",
    ],
    "annotation match type",
    required=False,
)


COL_ACCEPT = resolve_column(
    columns,
    [
        "annotation_accepted",
    ],
    "annotation accepted",
    required=False,
)


COL_SOURCE = resolve_column(
    columns,
    [
        "candidate_source",
    ],
    "candidate source",
    required=False,
)


selected_columns = [
    value

    for value
    in [
        COL_GENOME,
        COL_PROTEIN,
        COL_CONTIG,
        COL_START,
        COL_END,
        COL_STRAND,
        COL_RANK,

        COL_CLUSTER,
        COL_MODULE,

        COL_FE_POS,
        COL_FE_HMM,
        COL_FE_CAT,

        COL_FMH,
        COL_HEMES,

        COL_SIGNALP,
        COL_DEEPTMHMM,
        COL_NTM,
        COL_LOCALIZATION,

        COL_COG,
        COL_GENE,
        COL_PRODUCT,

        COL_MATCH,
        COL_ACCEPT,

        COL_SOURCE,
    ]

    if value is not None
]


selected_columns = list(
    dict.fromkeys(
        selected_columns
    )
)


print(
    f"  Reading {len(selected_columns)} selected fields "
    f"from complete catalogue..."
)


catalog = pd.read_csv(
    GENE_CATALOG,
    sep="\t",
    dtype=str,
    keep_default_na=False,
    usecols=selected_columns,
)


if len(
    catalog
) != EXPECTED_GENES:

    fail(
        f"Expected {EXPECTED_GENES:,} complete-catalogue genes; "
        f"found {len(catalog):,}."
    )


## ================================================================== ##
## 6. Standardize complete catalogue column names
## ================================================================== ##

rename = {
    COL_GENOME:
        "genome",

    COL_PROTEIN:
        "protein_id",

    COL_CONTIG:
        "contig",

    COL_START:
        "start",

    COL_END:
        "end",

    COL_STRAND:
        "strand",

    COL_RANK:
        "gene_rank",
}


optional_rename = {
    COL_CLUSTER:
        "mmseq_cluster",

    COL_MODULE:
        "mcl_module",

    COL_FE_POS:
        "fegenie_positive",

    COL_FE_HMM:
        "fegenie_HMMs",

    COL_FE_CAT:
        "fegenie_categories",

    COL_FMH:
        "findmehemes_positive",

    COL_HEMES:
        "number_of_hemes",

    COL_SIGNALP:
        "signalp_prediction",

    COL_DEEPTMHMM:
        "deeptmhmm_class",

    COL_NTM:
        "deeptmhmm_n_tm_helices",

    COL_LOCALIZATION:
        "upstream_localization_class",

    COL_COG:
        "globdb_cog",

    COL_GENE:
        "globdb_gene",

    COL_PRODUCT:
        "globdb_product",

    COL_MATCH:
        "annotation_match_type",

    COL_ACCEPT:
        "annotation_accepted",

    COL_SOURCE:
        "candidate_source",
}


for old, new in optional_rename.items():

    if old is not None:

        rename[
            old
        ] = new


catalog = catalog.rename(
    columns=rename
)


## Add absent optional fields explicitly. ##

optional_standard_fields = [
    "mmseq_cluster",
    "mcl_module",

    "fegenie_positive",
    "fegenie_HMMs",
    "fegenie_categories",

    "findmehemes_positive",
    "number_of_hemes",

    "signalp_prediction",
    "deeptmhmm_class",
    "deeptmhmm_n_tm_helices",
    "upstream_localization_class",

    "globdb_cog",
    "globdb_gene",
    "globdb_product",

    "annotation_match_type",
    "annotation_accepted",

    "candidate_source",
]


for column in optional_standard_fields:

    if column not in catalog.columns:

        catalog[
            column
        ] = ""


catalog[
    "start"
] = pd.to_numeric(
    catalog[
        "start"
    ],
    errors="raise",
).astype(int)


catalog[
    "end"
] = pd.to_numeric(
    catalog[
        "end"
    ],
    errors="raise",
).astype(int)


catalog[
    "gene_rank"
] = pd.to_numeric(
    catalog[
        "gene_rank"
    ],
    errors="raise",
).astype(int)


if catalog[
    [
        "genome",
        "protein_id",
    ]
].duplicated().any():

    fail(
        "Complete gene catalogue has duplicate genome + protein_id keys."
    )


print(
    f"  Complete genes:       "
    f"{len(catalog):,}"
)

print(
    f"  Genomes:              "
    f"{catalog['genome'].nunique():,}"
)

print(
    f"  Contigs with genes:   "
    f"{catalog[['genome','contig']].drop_duplicates().shape[0]:,}"
)


## ================================================================== ##
## 7. Reassert authoritative MMseq cluster/module assignments
##
## Stage-09 mapping is authoritative for the 10,537 module proteins.
## ================================================================== ##

module_annotation = (
    module_proteins[
        [
            "genome",
            "protein_id",
            "cluster",
            "module",
        ]
    ]
    .rename(
        columns={
            "cluster":
                "_authoritative_module_cluster",

            "module":
                "_authoritative_module",
        }
    )
)


catalog = catalog.merge(
    module_annotation,
    on=[
        "genome",
        "protein_id",
    ],
    how="left",
    validate="one_to_one",
)


module_mask = (
    catalog[
        "_authoritative_module"
    ].notna()
)


catalog.loc[
    module_mask,
    "mmseq_cluster",
] = catalog.loc[
    module_mask,
    "_authoritative_module_cluster",
]


catalog.loc[
    module_mask,
    "mcl_module",
] = catalog.loc[
    module_mask,
    "_authoritative_module",
]


catalog = catalog.drop(
    columns=[
        "_authoritative_module_cluster",
        "_authoritative_module",
    ]
)


## ================================================================== ##
## 8. Add report-level topology to complete catalogue
## ================================================================== ##

print()
print("Deriving report-level topology for complete gene catalogue...")


catalog[
    "report_topology"
] = [
    report_topology(
        deeptmhmm,
        signalp,
        n_tm,
    )

    for deeptmhmm, signalp, n_tm
    in zip(
        catalog[
            "deeptmhmm_class"
        ],
        catalog[
            "signalp_prediction"
        ],
        catalog[
            "deeptmhmm_n_tm_helices"
        ],
    )
]


## ================================================================== ##
## 9. Join complete gene evidence to fixed module-region gene set
## ================================================================== ##

print()
print("Recovering complete evidence for genes in fixed module regions...")


region_genes = region_gene_keys.merge(
    catalog,
    on=[
        "genome",
        "protein_id",
    ],
    how="left",
    validate="many_to_one",
    indicator=True,
    suffixes=(
        "_region",
        "",
    ),
)


if (
    region_genes[
        "_merge"
    ]
    !=
    "both"
).any():

    n_missing = int(
        (
            region_genes[
                "_merge"
            ]
            !=
            "both"
        )
        .sum()
    )


    fail(
        f"{n_missing:,} fixed-region genes could not be recovered "
        f"from the complete Stage-11 catalogue."
    )


region_genes = region_genes.drop(
    columns=[
        "_merge",
    ]
)


## Validate region contig / rank. ##

if (
    region_genes[
        "contig_region"
    ]
    !=
    region_genes[
        "contig"
    ]
).any():

    fail(
        "Stage-12 and Stage-11 contig identities disagree "
        "for at least one fixed-region gene."
    )


if (
    region_genes[
        "gene_rank_region"
    ].astype(int)
    !=
    region_genes[
        "gene_rank"
    ].astype(int)
).any():

    fail(
        "Stage-12 and Stage-11 gene ranks disagree "
        "for at least one fixed-region gene."
    )


region_genes = region_genes.drop(
    columns=[
        "contig_region",
        "gene_rank_region",
    ]
)


## ================================================================== ##
## 10. Assign automatic focal reasons
## ================================================================== ##

print()
print("Assigning gene-level focal eligibility...")


module_pair_keys = set(
    zip(
        module_proteins[
            "genome"
        ],
        module_proteins[
            "protein_id"
        ],
        module_proteins[
            "module"
        ],
    )
)


region_genes[
    "focal_module_member"
] = [
    int(
        (
            genome,
            protein_id,
            module,
        )
        in
        module_pair_keys
    )

    for genome, protein_id, module
    in zip(
        region_genes[
            "genome"
        ],
        region_genes[
            "protein_id"
        ],
        region_genes[
            "module"
        ],
    )
]


region_genes[
    "focal_fegenie_positive"
] = flag_numeric(
    region_genes[
        "fegenie_positive"
    ]
).astype(int)


region_genes[
    "focal_findmehemes_positive"
] = flag_numeric(
    region_genes[
        "findmehemes_positive"
    ]
).astype(int)


## Reserved for later manual promotion without changing automatic rules. ##

region_genes[
    "focal_manual"
] = 0


region_genes[
    "is_automatic_focal"
] = (
    (
        region_genes[
            [
                "focal_module_member",
                "focal_fegenie_positive",
                "focal_findmehemes_positive",
            ]
        ]
        .max(
            axis=1
        )
    )
    ==
    1
).astype(int)


def focal_reason(row):

    reasons = []


    if row[
        "focal_module_member"
    ] == 1:

        reasons.append(
            "module_member"
        )


    if row[
        "focal_fegenie_positive"
    ] == 1:

        reasons.append(
            "FeGenie"
        )


    if row[
        "focal_findmehemes_positive"
    ] == 1:

        reasons.append(
            "FindMeHemes"
        )


    if row[
        "focal_manual"
    ] == 1:

        reasons.append(
            "manual"
        )


    return ";".join(
        reasons
    )


region_genes[
    "focal_reasons"
] = region_genes.apply(
    focal_reason,
    axis=1,
)


## ================================================================== ##
## 11. Save all genes in fixed module-associated regions
## ================================================================== ##

region_gene_output_columns = [
    "region_id",

    "genome",
    "module",
    "contig",

    "protein_id",
    "start",
    "end",
    "strand",
    "gene_rank",

    "focal_module_member",
    "focal_fegenie_positive",
    "focal_findmehemes_positive",
    "focal_manual",
    "is_automatic_focal",
    "focal_reasons",

    "mmseq_cluster",
    "mcl_module",

    "fegenie_positive",
    "fegenie_HMMs",
    "fegenie_categories",

    "findmehemes_positive",
    "number_of_hemes",

    "signalp_prediction",
    "deeptmhmm_class",
    "deeptmhmm_n_tm_helices",

    "report_topology",
    "upstream_localization_class",

    "globdb_cog",
    "globdb_gene",
    "globdb_product",

    "annotation_match_type",
    "annotation_accepted",

    "candidate_source",
]


write_tsv(
    region_genes[
        region_gene_output_columns
    ],
    OUT_REGION_GENES,
)


## ================================================================== ##
## 12. Build gene-level focal catalogue
##
## Unit:
##
##     genome x module x actual focal protein
##
## Same biological gene may legitimately appear once for more than one
## module if it lies in fixed eligibility regions of multiple modules.
## ================================================================== ##

focals = region_genes[
    region_genes[
        "is_automatic_focal"
    ]
    ==
    1
].copy()


if focals[
    [
        "genome",
        "module",
        "protein_id",
    ]
].duplicated().any():

    fail(
        "Duplicate genome + module + focal protein rows."
    )


focals[
    "focal_id"
] = (
    focals[
        "genome"
    ]
    +
    "|"
    +
    focals[
        "module"
    ]
    +
    "|"
    +
    focals[
        "protein_id"
    ]
)


focal_columns = [
    "focal_id",
    "region_id",

    "genome",
    "module",

    "protein_id",
    "contig",
    "start",
    "end",
    "strand",
    "gene_rank",

    "focal_module_member",
    "focal_fegenie_positive",
    "focal_findmehemes_positive",
    "focal_manual",
    "focal_reasons",

    "mmseq_cluster",
    "mcl_module",

    "fegenie_positive",
    "fegenie_HMMs",
    "fegenie_categories",

    "findmehemes_positive",
    "number_of_hemes",

    "signalp_prediction",
    "deeptmhmm_class",
    "deeptmhmm_n_tm_helices",

    "report_topology",
    "upstream_localization_class",

    "globdb_cog",
    "globdb_gene",
    "globdb_product",

    "annotation_match_type",
    "annotation_accepted",

    "candidate_source",
]


focals = focals[
    focal_columns
].sort_values(
    [
        "module",
        "genome",
        "contig",
        "gene_rank",
        "protein_id",
    ],
    kind="stable",
).reset_index(
    drop=True
)


write_tsv(
    focals,
    OUT_FOCALS,
)


print()
print("Gene-level focal catalogue")


print(
    f"  Focal genome-module genes: "
    f"{len(focals):,}"
)

print(
    f"  Unique actual focal genes:  "
    f"{focals[['genome','protein_id']].drop_duplicates().shape[0]:,}"
)

print(
    f"  Module-member focal rows:   "
    f"{focals['focal_module_member'].sum():,}"
)

print(
    f"  FeGenie-positive focal rows:"
    f" {focals['focal_fegenie_positive'].sum():,}"
)

print(
    f"  FindMeHemes-positive rows:  "
    f"{focals['focal_findmehemes_positive'].sum():,}"
)


## ================================================================== ##
## 13. Unique actual focal protein catalogue
##
## One row per genome + protein_id, with all eligible modules combined.
## ================================================================== ##

unique_rows = []


for (
    genome,
    protein_id,
), group in focals.groupby(
    [
        "genome",
        "protein_id",
    ],
    sort=True,
):

    first = group.iloc[0]


    unique_rows.append(
        {
            "genome":
                genome,

            "protein_id":
                protein_id,

            "eligible_modules":
                ";".join(
                    sorted(
                        group[
                            "module"
                        ].unique()
                    )
                ),

            "n_eligible_modules":
                group[
                    "module"
                ].nunique(),

            "focal_module_member":
                int(
                    group[
                        "focal_module_member"
                    ].max()
                ),

            "focal_fegenie_positive":
                int(
                    group[
                        "focal_fegenie_positive"
                    ].max()
                ),

            "focal_findmehemes_positive":
                int(
                    group[
                        "focal_findmehemes_positive"
                    ].max()
                ),

            "focal_reasons":
                ";".join(
                    reason

                    for reason
                    in [
                        "module_member"
                        if group[
                            "focal_module_member"
                        ].max()
                        ==
                        1
                        else
                        "",

                        "FeGenie"
                        if group[
                            "focal_fegenie_positive"
                        ].max()
                        ==
                        1
                        else
                        "",

                        "FindMeHemes"
                        if group[
                            "focal_findmehemes_positive"
                        ].max()
                        ==
                        1
                        else
                        "",
                    ]

                    if reason != ""
                ),

            "contig":
                first[
                    "contig"
                ],

            "start":
                first[
                    "start"
                ],

            "end":
                first[
                    "end"
                ],

            "strand":
                first[
                    "strand"
                ],

            "gene_rank":
                first[
                    "gene_rank"
                ],

            "mmseq_cluster":
                first[
                    "mmseq_cluster"
                ],

            "mcl_module":
                first[
                    "mcl_module"
                ],

            "fegenie_HMMs":
                first[
                    "fegenie_HMMs"
                ],

            "fegenie_categories":
                first[
                    "fegenie_categories"
                ],

            "findmehemes_positive":
                first[
                    "findmehemes_positive"
                ],

            "number_of_hemes":
                first[
                    "number_of_hemes"
                ],

            "signalp_prediction":
                first[
                    "signalp_prediction"
                ],

            "deeptmhmm_class":
                first[
                    "deeptmhmm_class"
                ],

            "report_topology":
                first[
                    "report_topology"
                ],

            "globdb_cog":
                first[
                    "globdb_cog"
                ],

            "globdb_product":
                first[
                    "globdb_product"
                ],
        }
    )


unique_focals = pd.DataFrame(
    unique_rows
)


write_tsv(
    unique_focals,
    OUT_FOCALS_UNIQUE,
)


## ================================================================== ##
## 14. Fresh extraction requires only contigs containing selected focals
## ================================================================== ##

print()
print("Preparing fresh gene-centered extraction...")


focal_contig_keys = set(
    (
        genome
        +
        "\t"
        +
        contig
    )

    for genome, contig
    in focals[
        [
            "genome",
            "contig",
        ]
    ]
    .drop_duplicates()
    .itertuples(
        index=False,
        name=None,
    )
)


catalog[
    "_contig_key"
] = (
    catalog[
        "genome"
    ]
    +
    "\t"
    +
    catalog[
        "contig"
    ]
)


focal_contig_catalog = catalog[
    catalog[
        "_contig_key"
    ].isin(
        focal_contig_keys
    )
].copy()


focal_contig_catalog = focal_contig_catalog.sort_values(
    [
        "genome",
        "contig",
        "gene_rank",
    ],
    kind="stable",
).reset_index(
    drop=True
)


print(
    f"  Focal-containing contigs: "
    f"{len(focal_contig_keys):,}"
)

print(
    f"  Genes on focal contigs:   "
    f"{len(focal_contig_catalog):,}"
)


## ================================================================== ##
## 15. Build efficient contig lookup
## ================================================================== ##

contig_lookup = {}


for (
    genome,
    contig,
), group in focal_contig_catalog.groupby(
    [
        "genome",
        "contig",
    ],
    sort=False,
):

    group = group.sort_values(
        "gene_rank",
        kind="stable",
    ).reset_index(
        drop=True
    )


    contig_lookup[
        (
            genome,
            contig,
        )
    ] = group


## ================================================================== ##
## 16. Fresh focal-centered extraction
## ================================================================== ##

print()
print(
    "Freshly extracting +/-20 genes OR +/-20 kb "
    "around every selected actual focal gene..."
)


neighborhood_parts = []


focal_annotation_columns = [
    "focal_id",
    "region_id",

    "genome",
    "module",

    "protein_id",
    "contig",
    "start",
    "end",
    "strand",
    "gene_rank",

    "focal_module_member",
    "focal_fegenie_positive",
    "focal_findmehemes_positive",
    "focal_reasons",

    "mmseq_cluster",
    "mcl_module",

    "fegenie_HMMs",
    "fegenie_categories",

    "number_of_hemes",

    "signalp_prediction",
    "deeptmhmm_class",
    "report_topology",

    "globdb_cog",
    "globdb_product",
]


for focal in focals[
    focal_annotation_columns
].itertuples(
    index=False
):

    contig_df = contig_lookup.get(
        (
            focal.genome,
            focal.contig,
        )
    )


    if contig_df is None:

        fail(
            f"Missing focal contig from catalogue: "
            f"{focal.genome} / {focal.contig}"
        )


    focal_rank = int(
        focal.gene_rank
    )


    focal_start = int(
        focal.start
    )


    focal_end = int(
        focal.end
    )


    rank_diff = (
        contig_df[
            "gene_rank"
        ]
        -
        focal_rank
    )


    within_gene = (
        rank_diff.abs()
        <=
        GENE_RADIUS
    )


    ## CDS intersects focal CDS +/- BP_RADIUS. ##

    low = (
        focal_start
        -
        BP_RADIUS
    )


    high = (
        focal_end
        +
        BP_RADIUS
    )


    within_bp = (
        (
            contig_df[
                "end"
            ]
            >=
            low
        )
        &
        (
            contig_df[
                "start"
            ]
            <=
            high
        )
    )


    selected = contig_df[
        within_gene
        |
        within_bp
    ].copy()


    selected[
        "focal_id"
    ] = focal.focal_id


    selected[
        "focal_region_id"
    ] = focal.region_id


    selected[
        "focal_context_module"
    ] = focal.module


    selected[
        "focal_protein_id"
    ] = focal.protein_id


    selected[
        "focal_contig"
    ] = focal.contig


    selected[
        "focal_start"
    ] = focal_start


    selected[
        "focal_end"
    ] = focal_end


    selected[
        "focal_strand"
    ] = focal.strand


    selected[
        "focal_gene_rank"
    ] = focal_rank


    selected[
        "focal_selection_reasons"
    ] = focal.focal_reasons


    selected[
        "focal_is_module_member"
    ] = focal.focal_module_member


    selected[
        "focal_is_fegenie_positive"
    ] = focal.focal_fegenie_positive


    selected[
        "focal_is_findmehemes_positive"
    ] = focal.focal_findmehemes_positive


    selected[
        "focal_mmseq_cluster"
    ] = focal.mmseq_cluster


    selected[
        "focal_mcl_module"
    ] = focal.mcl_module


    selected[
        "focal_fegenie_HMMs"
    ] = focal.fegenie_HMMs


    selected[
        "focal_fegenie_categories"
    ] = focal.fegenie_categories


    selected[
        "focal_number_of_hemes"
    ] = focal.number_of_hemes


    selected[
        "focal_signalp_prediction"
    ] = focal.signalp_prediction


    selected[
        "focal_deeptmhmm_class"
    ] = focal.deeptmhmm_class


    selected[
        "focal_report_topology"
    ] = focal.report_topology


    selected[
        "focal_globdb_cog"
    ] = focal.globdb_cog


    selected[
        "focal_globdb_product"
    ] = focal.globdb_product


    ## -------------------------------------------------------------- ##
    ## Relative gene / physical position
    ## -------------------------------------------------------------- ##

    selected[
        "genomic_gene_offset"
    ] = (
        selected[
            "gene_rank"
        ]
        -
        focal_rank
    )


    orientation_factor = (
        1
        if focal.strand == "+"
        else
        -1
    )


    selected[
        "oriented_gene_offset"
    ] = (
        selected[
            "genomic_gene_offset"
        ]
        *
        orientation_factor
    )


    focal_midpoint = (
        focal_start
        +
        focal_end
    ) / 2.0


    selected_midpoint = (
        selected[
            "start"
        ]
        +
        selected[
            "end"
        ]
    ) / 2.0


    selected[
        "genomic_midpoint_offset_bp"
    ] = (
        selected_midpoint
        -
        focal_midpoint
    )


    selected[
        "oriented_midpoint_offset_bp"
    ] = (
        selected[
            "genomic_midpoint_offset_bp"
        ]
        *
        orientation_factor
    )


    selected[
        "within_20_genes"
    ] = (
        selected[
            "genomic_gene_offset"
        ]
        .abs()
        <=
        GENE_RADIUS
    ).astype(int)


    selected[
        "within_20kb"
    ] = (
        (
            selected[
                "end"
            ]
            >=
            low
        )
        &
        (
            selected[
                "start"
            ]
            <=
            high
        )
    ).astype(int)


    selected[
        "extraction_basis"
    ] = np.select(
        [
            (
                selected[
                    "within_20_genes"
                ]
                ==
                1
            )
            &
            (
                selected[
                    "within_20kb"
                ]
                ==
                1
            ),

            selected[
                "within_20_genes"
            ]
            ==
            1,

            selected[
                "within_20kb"
            ]
            ==
            1,
        ],
        [
            "both",
            "gene_radius",
            "bp_radius",
        ],
        default="",
    )


    selected[
        "is_focal"
    ] = (
        selected[
            "protein_id"
        ]
        ==
        focal.protein_id
    ).astype(int)


    ## -------------------------------------------------------------- ##
    ## Is surrounding gene inside original fixed eligibility region?
    ##
    ## Important for distinguishing:
    ##
    ##   original eligibility territory
    ## vs
    ##   newly exposed fresh-neighborhood context.
    ## -------------------------------------------------------------- ##

    original_region_gene_set = set(
        region_genes.loc[
            (
                region_genes[
                    "genome"
                ]
                ==
                focal.genome
            )
            &
            (
                region_genes[
                    "module"
                ]
                ==
                focal.module
            ),
            "protein_id",
        ]
    )


    selected[
        "neighbor_inside_original_fixed_region"
    ] = selected[
        "protein_id"
    ].isin(
        original_region_gene_set
    ).astype(int)


    ## -------------------------------------------------------------- ##
    ## Neighbor focal eligibility ONLY if it belonged to the original
    ## fixed module-associated eligibility region.
    ## -------------------------------------------------------------- ##

    selected[
        "neighbor_is_selected_focal_for_this_module"
    ] = selected[
        "protein_id"
    ].isin(
        set(
            focals.loc[
                (
                    focals[
                        "genome"
                    ]
                    ==
                    focal.genome
                )
                &
                (
                    focals[
                        "module"
                    ]
                    ==
                    focal.module
                ),
                "protein_id",
            ]
        )
    ).astype(int)


    neighborhood_parts.append(
        selected
    )


if len(
    neighborhood_parts
) == 0:

    fail(
        "No focal neighborhoods were produced."
    )


neighborhoods = pd.concat(
    neighborhood_parts,
    ignore_index=True,
)


## ================================================================== ##
## 17. Rename surrounding-gene fields explicitly
## ================================================================== ##

neighbor_rename = {
    "protein_id":
        "neighbor_protein_id",

    "contig":
        "neighbor_contig",

    "start":
        "neighbor_start",

    "end":
        "neighbor_end",

    "strand":
        "neighbor_strand",

    "gene_rank":
        "neighbor_gene_rank",

    "mmseq_cluster":
        "neighbor_mmseq_cluster",

    "mcl_module":
        "neighbor_mcl_module",

    "fegenie_positive":
        "neighbor_fegenie_positive",

    "fegenie_HMMs":
        "neighbor_fegenie_HMMs",

    "fegenie_categories":
        "neighbor_fegenie_categories",

    "findmehemes_positive":
        "neighbor_findmehemes_positive",

    "number_of_hemes":
        "neighbor_number_of_hemes",

    "signalp_prediction":
        "neighbor_signalp_prediction",

    "deeptmhmm_class":
        "neighbor_deeptmhmm_class",

    "deeptmhmm_n_tm_helices":
        "neighbor_deeptmhmm_n_tm_helices",

    "report_topology":
        "neighbor_report_topology",

    "upstream_localization_class":
        "neighbor_upstream_localization_class",

    "globdb_cog":
        "neighbor_globdb_cog",

    "globdb_gene":
        "neighbor_globdb_gene",

    "globdb_product":
        "neighbor_globdb_product",

    "annotation_match_type":
        "neighbor_annotation_match_type",

    "annotation_accepted":
        "neighbor_annotation_accepted",

    "candidate_source":
        "neighbor_candidate_source",
}


neighborhoods = neighborhoods.rename(
    columns=neighbor_rename
)


neighborhoods = neighborhoods.drop(
    columns=[
        "_contig_key",
    ],
    errors="ignore",
)


## ================================================================== ##
## 18. Sort and write fresh focal neighborhoods
## ================================================================== ##

neighborhoods = neighborhoods.sort_values(
    [
        "focal_context_module",
        "genome",
        "focal_contig",
        "focal_gene_rank",
        "neighbor_gene_rank",
    ],
    kind="stable",
).reset_index(
    drop=True
)


write_tsv(
    neighborhoods,
    OUT_NEIGHBORHOODS,
)


## ================================================================== ##
## 19. Focal-reason summary
## ================================================================== ##

reason_summary = (
    focals
    .groupby(
        [
            "focal_module_member",
            "focal_fegenie_positive",
            "focal_findmehemes_positive",
            "focal_reasons",
        ],
        as_index=False,
    )
    .agg(
        n_focal_rows=(
            "focal_id",
            "size",
        ),

        n_unique_actual_genes=(
            "protein_id",
            "nunique",
        ),

        n_genomes=(
            "genome",
            "nunique",
        ),

        n_modules=(
            "module",
            "nunique",
        ),
    )
    .sort_values(
        [
            "n_focal_rows",
            "focal_reasons",
        ],
        ascending=[
            False,
            True,
        ],
    )
)


write_tsv(
    reason_summary,
    OUT_REASON_SUMMARY,
)


## ================================================================== ##
## 20. QC
## ================================================================== ##

print()
print("Running resolution-recovery QC...")


## Every module protein must remain a focal gene in its own module. ##

module_focal_keys = set(
    zip(
        focals.loc[
            focals[
                "focal_module_member"
            ]
            ==
            1,
            "genome",
        ],
        focals.loc[
            focals[
                "focal_module_member"
            ]
            ==
            1,
            "protein_id",
        ],
        focals.loc[
            focals[
                "focal_module_member"
            ]
            ==
            1,
            "module",
        ],
    )
)


expected_module_keys = set(
    zip(
        module_proteins[
            "genome"
        ],
        module_proteins[
            "protein_id"
        ],
        module_proteins[
            "module"
        ],
    )
)


if module_focal_keys != expected_module_keys:

    missing = (
        expected_module_keys
        -
        module_focal_keys
    )


    extra = (
        module_focal_keys
        -
        expected_module_keys
    )


    fail(
        f"Module-focal reconstruction mismatch. "
        f"Missing={len(missing):,}; extra={len(extra):,}"
    )


## Every focal must have exactly one self row. ##

self_counts = (
    neighborhoods[
        neighborhoods[
            "is_focal"
        ]
        ==
        1
    ]
    .groupby(
        "focal_id"
    )
    .size()
)


if len(
    self_counts
) != len(
    focals
):

    fail(
        "Not every focal gene has a self row."
    )


if (
    self_counts
    !=
    1
).any():

    fail(
        "At least one focal gene has !=1 self row."
    )


## No fresh neighborhood row may be on another contig. ##

if (
    neighborhoods[
        "focal_contig"
    ]
    !=
    neighborhoods[
        "neighbor_contig"
    ]
).any():

    fail(
        "Fresh extraction contains a cross-contig neighbor."
    )


## Extraction basis must be valid. ##

if not neighborhoods[
    "extraction_basis"
].isin(
    [
        "both",
        "gene_radius",
        "bp_radius",
    ]
).all():

    fail(
        "Unexpected extraction_basis value."
    )


## Confirm automatic focal conditions. ##

bad_focal = focals[
    (
        focals[
            [
                "focal_module_member",
                "focal_fegenie_positive",
                "focal_findmehemes_positive",
            ]
        ]
        .max(
            axis=1
        )
        !=
        1
    )
]


if len(
    bad_focal
) > 0:

    fail(
        "A selected automatic focal gene has no automatic focal reason."
    )


## Confirm newly exposed genes outside fixed region did not acquire
## automatic focal status merely by being exposed by fresh extraction. ##

n_outside_context = int(
    (
        neighborhoods[
            "neighbor_inside_original_fixed_region"
        ]
        ==
        0
    )
    .sum()
)


n_outside_marked_focal = int(
    (
        (
            neighborhoods[
                "neighbor_inside_original_fixed_region"
            ]
            ==
            0
        )
        &
        (
            neighborhoods[
                "neighbor_is_selected_focal_for_this_module"
            ]
            ==
            1
        )
    )
    .sum()
)


if n_outside_marked_focal != 0:

    fail(
        "Recursive focal expansion detected: a newly exposed gene "
        "outside the fixed eligibility region became focal."
    )


qc = pd.DataFrame(
    [
        [
            "complete_gene_catalog_rows",
            len(
                catalog
            ),
            EXPECTED_GENES,
            int(
                len(
                    catalog
                )
                ==
                EXPECTED_GENES
            ),
        ],

        [
            "original_stage12_rows",
            len(
                stage12
            ),
            EXPECTED_STAGE12_ROWS,
            int(
                len(
                    stage12
                )
                ==
                EXPECTED_STAGE12_ROWS
            ),
        ],

        [
            "original_module_focals",
            stage12[
                "focal_protein_id"
            ].nunique(),
            EXPECTED_STAGE12_FOCALS,
            int(
                stage12[
                    "focal_protein_id"
                ].nunique()
                ==
                EXPECTED_STAGE12_FOCALS
            ),
        ],

        [
            "module_proteins_recovered_as_focals",
            len(
                module_focal_keys
            ),
            EXPECTED_MODULE_PROTEINS,
            int(
                len(
                    module_focal_keys
                )
                ==
                EXPECTED_MODULE_PROTEINS
            ),
        ],

        [
            "modules",
            focals[
                "module"
            ].nunique(),
            EXPECTED_MODULES,
            int(
                focals[
                    "module"
                ].nunique()
                ==
                EXPECTED_MODULES
            ),
        ],

        [
            "fixed_eligibility_regions",
            len(
                regions
            ),
            "",
            1,
        ],

        [
            "fixed_region_gene_module_rows",
            len(
                region_genes
            ),
            "",
            1,
        ],

        [
            "unique_fixed_region_genes",
            region_genes[
                [
                    "genome",
                    "protein_id",
                ]
            ]
            .drop_duplicates()
            .shape[0],
            "",
            1,
        ],

        [
            "gene_level_focal_rows",
            len(
                focals
            ),
            "",
            1,
        ],

        [
            "unique_actual_focal_genes",
            len(
                unique_focals
            ),
            "",
            1,
        ],

        [
            "fegenie_positive_focal_rows",
            int(
                focals[
                    "focal_fegenie_positive"
                ].sum()
            ),
            "",
            1,
        ],

        [
            "findmehemes_positive_focal_rows",
            int(
                focals[
                    "focal_findmehemes_positive"
                ].sum()
            ),
            "",
            1,
        ],

        [
            "fresh_focal_neighborhood_rows",
            len(
                neighborhoods
            ),
            "",
            1,
        ],

        [
            "fresh_focal_self_rows",
            int(
                (
                    neighborhoods[
                        "is_focal"
                    ]
                    ==
                    1
                )
                .sum()
            ),
            len(
                focals
            ),
            int(
                (
                    neighborhoods[
                        "is_focal"
                    ]
                    ==
                    1
                )
                .sum()
                ==
                len(
                    focals
                )
            ),
        ],

        [
            "fresh_context_rows_outside_original_fixed_region",
            n_outside_context,
            "",
            1,
        ],

        [
            "outside_region_rows_marked_as_selected_focal",
            n_outside_marked_focal,
            0,
            int(
                n_outside_marked_focal
                ==
                0
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


write_tsv(
    qc,
    OUT_QC,
)


## ================================================================== ##
## 21. Terminal summaries
## ================================================================== ##

print()
print("Resolution-recovery summary")


print(
    f"  Complete Prodigal genes:              "
    f"{len(catalog):,}"
)

print(
    f"  Fixed module eligibility regions:     "
    f"{len(regions):,}"
)

print(
    f"  Gene x module region entries:         "
    f"{len(region_genes):,}"
)

print(
    f"  Unique genes in fixed regions:        "
    f"{region_genes[['genome','protein_id']].drop_duplicates().shape[0]:,}"
)

print(
    f"  Gene-level focal rows:                "
    f"{len(focals):,}"
)

print(
    f"  Unique actual focal genes:            "
    f"{len(unique_focals):,}"
)

print(
    f"  Module-member focal rows:             "
    f"{int(focals['focal_module_member'].sum()):,}"
)

print(
    f"  FeGenie-positive focal rows:          "
    f"{int(focals['focal_fegenie_positive'].sum()):,}"
)

print(
    f"  FindMeHemes-positive focal rows:      "
    f"{int(focals['focal_findmehemes_positive'].sum()):,}"
)

print(
    f"  Fresh focal-neighborhood rows:        "
    f"{len(neighborhoods):,}"
)

print(
    f"  Newly exposed context rows outside "
    f"fixed regions: {n_outside_context:,}"
)

print(
    f"  Recursive outside-region focals:      "
    f"{n_outside_marked_focal:,}"
)


print()
print("Focal-reason combinations")


print(
    reason_summary.to_string(
        index=False
    )
)


## ================================================================== ##
## 22. Module 20 / Module 35 gene-resolution spotlights
## ================================================================== ##

for spotlight in [
    "Module_20",
    "Module_35",
]:

    spot = focals[
        focals[
            "module"
        ]
        ==
        spotlight
    ].copy()


    print()
    print(
        f"{spotlight} gene-level focal spotlight"
    )


    print(
        f"  Focal rows:          "
        f"{len(spot):,}"
    )

    print(
        f"  Unique genes:        "
        f"{spot[['genome','protein_id']].drop_duplicates().shape[0]:,}"
    )

    print(
        f"  Module members:      "
        f"{int(spot['focal_module_member'].sum()):,}"
    )

    print(
        f"  FeGenie-positive:    "
        f"{int(spot['focal_fegenie_positive'].sum()):,}"
    )

    print(
        f"  FindMeHemes-positive:"
        f" {int(spot['focal_findmehemes_positive'].sum()):,}"
    )


    display = [
        "genome",
        "protein_id",

        "focal_reasons",

        "mmseq_cluster",
        "mcl_module",

        "fegenie_HMMs",
        "number_of_hemes",

        "signalp_prediction",
        "deeptmhmm_class",
        "report_topology",

        "globdb_cog",
        "globdb_product",
    ]


    ## Prefer additional biologically interesting genes first. ##

    spot[
        "_priority"
    ] = (
        spot[
            "focal_fegenie_positive"
        ]
        +
        spot[
            "focal_findmehemes_positive"
        ]
        -
        spot[
            "focal_module_member"
        ]
    )


    print()
    print(
        spot.sort_values(
            [
                "_priority",
                "focal_module_member",
                "genome",
                "gene_rank",
            ],
            ascending=[
                False,
                False,
                True,
                True,
            ],
        )[
            display
        ]
        .head(
            30
        )
        .to_string(
            index=False
        )
    )


print()
print("=" * 80)
print("SUCCESS")
print("=" * 80)


print(
    f"Fixed regions:            "
    f"{OUT_REGIONS}"
)

print(
    f"All genes in regions:     "
    f"{OUT_REGION_GENES}"
)

print(
    f"Gene-level focals:        "
    f"{OUT_FOCALS}"
)

print(
    f"Unique actual focals:     "
    f"{OUT_FOCALS_UNIQUE}"
)

print(
    f"Fresh focal neighborhoods:"
    f" {OUT_NEIGHBORHOODS}"
)

print(
    f"Focal reason summary:     "
    f"{OUT_REASON_SUMMARY}"
)

print(
    f"QC:                       "
    f"{OUT_QC}"
)
