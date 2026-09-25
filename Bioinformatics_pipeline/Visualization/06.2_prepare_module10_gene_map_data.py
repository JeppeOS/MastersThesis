#!/usr/bin/env python3


## ================================================================== ##
## STAGE 16A - MODULE 10 / Cluster_00050 Cyc2 GENE-MAP DATA
##
## Purpose
## -------
##
## Convert the deduplicated Stage-15A Cluster_00050 neighborhoods into
## the focal-oriented, plot-ready format expected by the generalized
## Stage-16A2 local architecture pipeline.
##
##
## Authoritative biological focal unit
## -----------------------------------
##
##     genome + focal_protein_id
##
## The preceding QC established:
##
##     53 unique biological Cluster_00050 focal loci
##     53 genomes
##
##
## Focal orientation
## -----------------
##
## Every locus is normalized so that Cluster_00050 / Cyc2 points
## left-to-right.
##
## If focal strand is "+":
##
##     plot coordinate = genomic coordinate - focal midpoint
##
## If focal strand is "-":
##
##     plot coordinate = focal midpoint - genomic coordinate
##
## and gene strands are reversed accordingly.
##
##
## Plotting / local-clustering window
## ----------------------------------
##
## Every actual CDS intersecting:
##
##     -10,000 bp ... +10,000 bp
##
## relative to the focal Cyc2 midpoint is retained.
##
## This matches the principle used for the Module 20 visualization:
##
##     filter visual prominence, not biological context
##
## i.e. ordinary genes remain in the local-family analysis rather
## than restricting analysis only to cytochromes or module members.
##
##
## IMPORTANT
## ---------
##
## display_class and display_label are visualization metadata only.
##
## They DO NOT modify:
##
##     global MMseq family identity
##     MCL module identity
##     FeGenie
##     FindMeHemes
##     SignalP
##     DeepTMHMM
##     GlobDB annotations
##
## Those source columns remain authoritative.
## ================================================================== ##


from pathlib import Path

import numpy as np
import pandas as pd


## ================================================================== ##
## Constants
## ================================================================== ##

WORKFLOW = Path(
    "~/methanotrophs/methanotroph_project/jeppe/"
    "genome_analysis_workflow"
).expanduser()


FOCAL_CLUSTER = "Cluster_00050"

MCL_MODULE = "Module_10"

WINDOW_BP = 10_000


MODULE_CLUSTERS = {
    "Cluster_00050",
    "Cluster_00064",
    "Cluster_00069",
    "Cluster_00091",
    "Cluster_00222",
}


EXPECTED_FOCAL_REGIONS = 53


## ================================================================== ##
## Inputs
## ================================================================== ##

INPUT_FILE = (
    WORKFLOW
    / "16_visualization"
    / "16A_cluster00050_cyc2_gene_map_framework"
    / "Cluster_00050_unique_focal_neighborhoods.tsv"
)


TAXONOMY_FILE = (
    WORKFLOW
    / "16_visualization"
    / "globdb_r226_taxonomy.tsv"
)


## ================================================================== ##
## Outputs
## ================================================================== ##

OUTPUT_DIR = (
    WORKFLOW
    / "16_visualization"
    / "16A_cluster00050_cyc2_gene_map_framework"
)


OUTPUT_FILE = (
    OUTPUT_DIR
    / "Module_10_gene_map_data.tsv"
)


FOCAL_REGIONS_FILE = (
    OUTPUT_DIR
    / "Module_10_focal_regions.tsv"
)


QC_FILE = (
    OUTPUT_DIR
    / "Module_10_gene_map_qc.tsv"
)


## ================================================================== ##
## Helpers
## ================================================================== ##

def fail(message):

    raise RuntimeError(
        message
    )


def clean(series):

    return (
        series
        .fillna("")
        .astype(str)
        .str.strip()
    )


def numeric(series, name):

    out = pd.to_numeric(
        series,
        errors="coerce"
    )

    if out.isna().any():

        bad = series.loc[
            out.isna()
        ].head(
            20
        )

        fail(
            f"Could not parse numeric column {name}. "
            f"Example problematic values:\n{bad}"
        )

    return out


def parse_taxonomy(
    text,
    prefix
):

    text = str(
        text
    )

    for token in text.split(";"):

        token = token.strip()

        if token.startswith(
            prefix
        ):

            value = token[
                len(prefix):
            ].strip()

            if value:

                return value

    return ""


def flip_strand(value):

    value = str(
        value
    ).strip()

    if value == "+":
        return "-"

    if value == "-":
        return "+"

    return value


## ================================================================== ##
## Read input
## ================================================================== ##

print(
    "=" * 100
)

print(
    "READ DEDUPLICATED Cluster_00050 NEIGHBORHOODS"
)

print(
    "=" * 100
)


if not INPUT_FILE.exists():

    fail(
        f"Input file does not exist:\n{INPUT_FILE}"
    )


x = pd.read_csv(
    INPUT_FILE,
    sep="\t",
    dtype=str,
    keep_default_na=False
)


print(
    f"Input rows: {len(x):,}"
)


## ================================================================== ##
## Required input columns
## ================================================================== ##

required = {

    "genome",

    "neighbor_contig",
    "neighbor_protein_id",
    "neighbor_start",
    "neighbor_end",
    "neighbor_strand",

    "neighbor_mmseq_cluster",
    "neighbor_mcl_module",

    "neighbor_fegenie_positive",
    "neighbor_fegenie_HMMs",

    "neighbor_findmehemes_positive",
    "neighbor_number_of_hemes",

    "neighbor_signalp_prediction",
    "neighbor_deeptmhmm_class",
    "neighbor_report_topology",

    "neighbor_globdb_cog",
    "neighbor_globdb_product",

    "focal_protein_id",
    "focal_contig",
    "focal_start",
    "focal_end",
    "focal_strand",

    "focal_mmseq_cluster",
    "focal_mcl_module",

    "is_focal",

    "plot_region_id",

}


missing = (
    required
    -
    set(
        x.columns
    )
)


if missing:

    fail(
        "Missing required columns:\n"
        +
        "\n".join(
            sorted(
                missing
            )
        )
    )


## ================================================================== ##
## Clean identifiers
## ================================================================== ##

for col in [

    "genome",
    "neighbor_protein_id",
    "neighbor_contig",
    "neighbor_strand",

    "neighbor_mmseq_cluster",
    "neighbor_mcl_module",

    "focal_protein_id",
    "focal_contig",
    "focal_strand",

    "focal_mmseq_cluster",
    "focal_mcl_module",

    "plot_region_id",

]:

    x[
        col
    ] = clean(
        x[
            col
        ]
    )


## ================================================================== ##
## Numerical coordinates
## ================================================================== ##

x[
    "neighbor_start_num"
] = numeric(
    x[
        "neighbor_start"
    ],
    "neighbor_start"
)


x[
    "neighbor_end_num"
] = numeric(
    x[
        "neighbor_end"
    ],
    "neighbor_end"
)


x[
    "focal_start_num"
] = numeric(
    x[
        "focal_start"
    ],
    "focal_start"
)


x[
    "focal_end_num"
] = numeric(
    x[
        "focal_end"
    ],
    "focal_end"
)


## ================================================================== ##
## Validate focal metadata
## ================================================================== ##

print()

print(
    "=" * 100
)

print(
    "FOCAL METADATA QC"
)

print(
    "=" * 100
)


region_meta = (

    x.groupby(
        "plot_region_id",
        dropna=False
    )

    .agg(

        genome=(
            "genome",
            "first"
        ),

        n_genomes=(
            "genome",
            "nunique"
        ),

        focal_protein_id=(
            "focal_protein_id",
            "first"
        ),

        n_focal_proteins=(
            "focal_protein_id",
            "nunique"
        ),

        focal_contig=(
            "focal_contig",
            "first"
        ),

        n_focal_contigs=(
            "focal_contig",
            "nunique"
        ),

        focal_start=(
            "focal_start_num",
            "first"
        ),

        n_focal_starts=(
            "focal_start_num",
            "nunique"
        ),

        focal_end=(
            "focal_end_num",
            "first"
        ),

        n_focal_ends=(
            "focal_end_num",
            "nunique"
        ),

        focal_strand=(
            "focal_strand",
            "first"
        ),

        n_focal_strands=(
            "focal_strand",
            "nunique"
        ),

        focal_mmseq_cluster=(
            "focal_mmseq_cluster",
            "first"
        ),

        n_focal_clusters=(
            "focal_mmseq_cluster",
            "nunique"
        ),

    )

    .reset_index()

)


bad_meta = region_meta.loc[

    (
        region_meta[
            "n_genomes"
        ]
        !=
        1
    )

    |

    (
        region_meta[
            "n_focal_proteins"
        ]
        !=
        1
    )

    |

    (
        region_meta[
            "n_focal_contigs"
        ]
        !=
        1
    )

    |

    (
        region_meta[
            "n_focal_starts"
        ]
        !=
        1
    )

    |

    (
        region_meta[
            "n_focal_ends"
        ]
        !=
        1
    )

    |

    (
        region_meta[
            "n_focal_strands"
        ]
        !=
        1
    )

    |

    (
        region_meta[
            "n_focal_clusters"
        ]
        !=
        1
    )

]


if not bad_meta.empty:

    print(
        bad_meta.to_string(
            index=False
        )
    )

    fail(
        "Inconsistent focal metadata within one or more "
        "biological focal regions."
    )


if len(
    region_meta
) != EXPECTED_FOCAL_REGIONS:

    fail(
        f"Expected {EXPECTED_FOCAL_REGIONS} biological focal regions, "
        f"found {len(region_meta)}."
    )


if set(
    region_meta[
        "focal_mmseq_cluster"
    ]
) != {
    FOCAL_CLUSTER
}:

    fail(
        "One or more focal regions do not use Cluster_00050."
    )


print(
    f"Biological focal regions: {len(region_meta)}"
)

print(
    f"Genomes: {region_meta['genome'].nunique()}"
)


## ================================================================== ##
## Focal midpoint
## ================================================================== ##

x[
    "focal_midpoint"
] = (
    x[
        "focal_start_num"
    ]
    +
    x[
        "focal_end_num"
    ]
) / 2.0


## ================================================================== ##
## Normalize coordinates relative to focal Cyc2 orientation
## ================================================================== ##

plus = (
    x[
        "focal_strand"
    ]
    ==
    "+"
)


minus = (
    x[
        "focal_strand"
    ]
    ==
    "-"
)


if not (
    plus
    |
    minus
).all():

    bad = sorted(
        x.loc[
            ~(
                plus
                |
                minus
            ),
            "focal_strand"
        ].unique()
    )

    fail(
        f"Unexpected focal strand values: {bad}"
    )


## ------------------------------------------------------------------ ##
## Focal points right:
##
## '+' focal:
##
##     ordinary genomic coordinate subtraction
##
## '-' focal:
##
##     reverse coordinate system around focal midpoint
## ------------------------------------------------------------------ ##

x[
    "plot_start_bp"
] = np.where(

    plus,

    x[
        "neighbor_start_num"
    ]
    -
    x[
        "focal_midpoint"
    ],

    x[
        "focal_midpoint"
    ]
    -
    x[
        "neighbor_end_num"
    ]

)


x[
    "plot_end_bp"
] = np.where(

    plus,

    x[
        "neighbor_end_num"
    ]
    -
    x[
        "focal_midpoint"
    ],

    x[
        "focal_midpoint"
    ]
    -
    x[
        "neighbor_start_num"
    ]

)


x[
    "plot_strand"
] = np.where(

    plus,

    x[
        "neighbor_strand"
    ],

    x[
        "neighbor_strand"
    ].map(
        flip_strand
    )

)


## ================================================================== ##
## Basic coordinate QC
## ================================================================== ##

if (
    x[
        "plot_start_bp"
    ]
    >
    x[
        "plot_end_bp"
    ]
).any():

    fail(
        "At least one normalized gene has plot_start_bp > plot_end_bp."
    )


x[
    "gene_midpoint_bp"
] = (
    x[
        "plot_start_bp"
    ]
    +
    x[
        "plot_end_bp"
    ]
) / 2.0


## ================================================================== ##
## Keep every actual CDS intersecting +/-10 kb
## ================================================================== ##

x = x.loc[

    (
        x[
            "plot_end_bp"
        ]
        >=
        -WINDOW_BP
    )

    &

    (
        x[
            "plot_start_bp"
        ]
        <=
        WINDOW_BP
    )

].copy()


print()

print(
    "=" * 100
)

print(
    "PLOTTING WINDOW"
)

print(
    "=" * 100
)


print(
    f"Genes intersecting +/-{WINDOW_BP:,} bp: {len(x):,}"
)


## ================================================================== ##
## Standardized protein ID
## ================================================================== ##

x[
    "protein_id"
] = x[
    "neighbor_protein_id"
]


## ================================================================== ##
## Focal marker
## ================================================================== ##

is_focal_flag = pd.to_numeric(
    x[
        "is_focal"
    ],
    errors="coerce"
).fillna(
    0
).astype(
    int
)


x[
    "is_focal_plot"
] = (

    (
        is_focal_flag
        ==
        1
    )

    |

    (
        (
            x[
                "neighbor_protein_id"
            ]
            ==
            x[
                "focal_protein_id"
            ]
        )

        &

        (
            x[
                "neighbor_mmseq_cluster"
            ]
            ==
            FOCAL_CLUSTER
        )
    )

)


## ================================================================== ##
## Exactly one focal gene per biological region
## ================================================================== ##

focal_counts = (

    x.loc[
        x[
            "is_focal_plot"
        ]
    ]

    .groupby(
        "plot_region_id"
    )

    .size()

)


missing_focals = (
    set(
        region_meta[
            "plot_region_id"
        ]
    )
    -
    set(
        focal_counts.index
    )
)


if missing_focals:

    fail(
        "Some focal loci disappeared from the +/-10 kb window:\n"
        +
        "\n".join(
            sorted(
                missing_focals
            )
        )
    )


bad_focal_counts = focal_counts.loc[
    focal_counts
    !=
    1
]


if not bad_focal_counts.empty:

    print(
        bad_focal_counts
    )

    fail(
        "Expected exactly one focal Cyc2 gene per plot region."
    )


## ================================================================== ##
## Every focal must point right after normalization
## ================================================================== ##

bad_focal_strand = x.loc[

    x[
        "is_focal_plot"
    ]

    &

    (
        x[
            "plot_strand"
        ]
        !=
        "+"
    )

]


if not bad_focal_strand.empty:

    fail(
        "At least one focal Cyc2 does not point right after "
        "coordinate normalization."
    )


## ================================================================== ##
## Taxonomy
## ================================================================== ##

print()

print(
    "=" * 100
)

print(
    "ADD GLOBDB r226 TAXONOMY"
)

print(
    "=" * 100
)


if not TAXONOMY_FILE.exists():

    fail(
        f"Taxonomy file missing:\n{TAXONOMY_FILE}"
    )


taxonomy = pd.read_csv(
    TAXONOMY_FILE,
    sep="\t",
    header=None,
    names=[
        "genome",
        "taxonomy"
    ],
    usecols=[
        0,
        1
    ],
    dtype=str,
    keep_default_na=False
)


taxonomy = taxonomy.drop_duplicates(
    "genome"
)


taxonomy[
    "taxonomy_family"
] = taxonomy[
    "taxonomy"
].map(
    lambda value:
    parse_taxonomy(
        value,
        "f__"
    )
)


taxonomy[
    "taxonomy_genus"
] = taxonomy[
    "taxonomy"
].map(
    lambda value:
    parse_taxonomy(
        value,
        "g__"
    )
)


taxonomy[
    "taxonomy_species"
] = taxonomy[
    "taxonomy"
].map(
    lambda value:
    parse_taxonomy(
        value,
        "s__"
    )
)


taxonomy[
    "taxonomy_display"
] = np.where(

    taxonomy[
        "taxonomy_species"
    ]
    !=
    "",

    taxonomy[
        "taxonomy_species"
    ],

    np.where(

        taxonomy[
            "taxonomy_genus"
        ]
        !=
        "",

        taxonomy[
            "taxonomy_genus"
        ],

        taxonomy[
            "genome"
        ]

    )

)


taxonomy[
    "priority_taxon"
] = np.select(

    [

        taxonomy[
            "taxonomy_genus"
        ].str.startswith(
            "Methylobacter",
            na=False
        ),

        taxonomy[
            "taxonomy_genus"
        ].str.startswith(
            "Crenothrix",
            na=False
        ),

    ],

    [
        "Methylobacter",
        "Crenothrix",
    ],

    default="Other"

)


x = x.merge(

    taxonomy[
        [
            "genome",
            "taxonomy",

            "taxonomy_family",
            "taxonomy_genus",
            "taxonomy_species",
            "taxonomy_display",

            "priority_taxon",
        ]
    ],

    on="genome",

    how="left",

    validate="many_to_one"

)


missing_taxonomy = x.loc[
    x[
        "taxonomy_display"
    ].isna()
    |
    (
        x[
            "taxonomy_display"
        ]
        ==
        ""
    ),
    "genome"
].drop_duplicates()


if not missing_taxonomy.empty:

    fail(
        "Taxonomy missing for:\n"
        +
        "\n".join(
            missing_taxonomy
        )
    )


## ================================================================== ##
## Visualization-level classes
##
## These are deliberately broad at this stage.
##
## Local-family architecture analysis will later determine which
## recurrent families deserve explicit colors / biological labels.
## ================================================================== ##

fegenie_positive = pd.to_numeric(
    x[
        "neighbor_fegenie_positive"
    ],
    errors="coerce"
).fillna(
    0
).astype(
    int
) == 1


findmehemes_positive = pd.to_numeric(
    x[
        "neighbor_findmehemes_positive"
    ],
    errors="coerce"
).fillna(
    0
).astype(
    int
) == 1


same_module = (

    x[
        "neighbor_mmseq_cluster"
    ].isin(
        MODULE_CLUSTERS
    )

    &

    ~x[
        "is_focal_plot"
    ]

)


x[
    "display_class"
] = np.select(

    [

        x[
            "is_focal_plot"
        ],

        same_module,

        fegenie_positive,

        findmehemes_positive,

    ],

    [

        "focal_cluster00050",

        "same_module_family",

        "fegenie_other",

        "findmehemes_other",

    ],

    default="background_gene"

)


## ================================================================== ##
## Initial labels
##
## These are intentionally conservative.
##
## We will refine module-partner labels after inspecting the exhaustive
## local-family architecture rather than guessing functions now.
## ================================================================== ##

module_label_map = {

    "Cluster_00064":
        "Cluster_00064",

    "Cluster_00069":
        "Cluster_00069",

    "Cluster_00091":
        "Cluster_00091",

    "Cluster_00222":
        "Cluster_00222",

}


x[
    "display_label"
] = ""


x.loc[
    x[
        "is_focal_plot"
    ],
    "display_label"
] = "Cyc2"


for cluster, label in module_label_map.items():

    x.loc[
        (
            x[
                "neighbor_mmseq_cluster"
            ]
            ==
            cluster
        )
        &
        ~x[
            "is_focal_plot"
        ],
        "display_label"
    ] = label


## ------------------------------------------------------------------ ##
## Other FeGenie hits:
## retain their authoritative HMM name as the initial label.
## ------------------------------------------------------------------ ##

mask = (

    x[
        "display_class"
    ]
    ==
    "fegenie_other"
)


x.loc[
    mask,
    "display_label"
] = clean(
    x.loc[
        mask,
        "neighbor_fegenie_HMMs"
    ]
)


x.loc[
    mask
    &
    (
        x[
            "display_label"
        ]
        ==
        ""
    ),
    "display_label"
] = "Other FeGenie hit"


## ------------------------------------------------------------------ ##
## Other heme proteins:
## generic label only.
## ------------------------------------------------------------------ ##

x.loc[
    x[
        "display_class"
    ]
    ==
    "findmehemes_other",
    "display_label"
] = "Other heme protein"


## ================================================================== ##
## Duplicate gene observations within region QC
## ================================================================== ##

duplicates = x.loc[
    x.duplicated(
        [
            "plot_region_id",
            "protein_id"
        ],
        keep=False
    )
]


if not duplicates.empty:

    print(
        duplicates[
            [
                "plot_region_id",
                "genome",
                "protein_id",
            ]
        ].to_string(
            index=False
        )
    )

    fail(
        "Duplicate plot_region_id + protein_id rows remain."
    )


## ================================================================== ##
## Focal region output
## ================================================================== ##

focal_regions = (

    x.loc[
        x[
            "is_focal_plot"
        ],

        [
            "plot_region_id",

            "genome",
            "focal_protein_id",
            "focal_contig",

            "taxonomy_family",
            "taxonomy_genus",
            "taxonomy_species",
            "taxonomy_display",
            "priority_taxon",

        ]

    ]

    .drop_duplicates()

    .sort_values(
        [
            "taxonomy_family",
            "taxonomy_genus",
            "taxonomy_species",
            "genome"
        ]
    )

)


if len(
    focal_regions
) != EXPECTED_FOCAL_REGIONS:

    fail(
        f"Expected {EXPECTED_FOCAL_REGIONS} final focal-region rows, "
        f"found {len(focal_regions)}."
    )


## ================================================================== ##
## Final ordering
## ================================================================== ##

x = x.sort_values(

    [
        "plot_region_id",
        "plot_start_bp",
        "plot_end_bp",
        "protein_id",
    ]

).reset_index(
    drop=True
)


## ================================================================== ##
## Write output
## ================================================================== ##

x.to_csv(
    OUTPUT_FILE,
    sep="\t",
    index=False
)


focal_regions.to_csv(
    FOCAL_REGIONS_FILE,
    sep="\t",
    index=False
)


## ================================================================== ##
## QC table
## ================================================================== ##

qc = pd.DataFrame(
    [
        {
            "metric":
                "n_regions",

            "value":
                x[
                    "plot_region_id"
                ].nunique()
        },
        {
            "metric":
                "n_genomes",

            "value":
                x[
                    "genome"
                ].nunique()
        },
        {
            "metric":
                "n_rows",

            "value":
                len(
                    x
                )
        },
        {
            "metric":
                "n_unique_proteins",

            "value":
                x[
                    [
                        "genome",
                        "protein_id"
                    ]
                ].drop_duplicates().shape[0]
        },
        {
            "metric":
                "n_focal_rows",

            "value":
                int(
                    x[
                        "is_focal_plot"
                    ].sum()
                )
        },
        {
            "metric":
                "n_module10_partner_rows",

            "value":
                int(
                    (
                        x[
                            "display_class"
                        ]
                        ==
                        "same_module_family"
                    ).sum()
                )
        },
        {
            "metric":
                "n_fegenie_other_rows",

            "value":
                int(
                    (
                        x[
                            "display_class"
                        ]
                        ==
                        "fegenie_other"
                    ).sum()
                )
        },
        {
            "metric":
                "n_findmehemes_other_rows",

            "value":
                int(
                    (
                        x[
                            "display_class"
                        ]
                        ==
                        "findmehemes_other"
                    ).sum()
                )
        },
        {
            "metric":
                "n_background_rows",

            "value":
                int(
                    (
                        x[
                            "display_class"
                        ]
                        ==
                        "background_gene"
                    ).sum()
                )
        },
    ]
)


qc.to_csv(
    QC_FILE,
    sep="\t",
    index=False
)


## ================================================================== ##
## Console summary
## ================================================================== ##

print()

print(
    "=" * 100
)

print(
    "MODULE 10 GENE-MAP DATA COMPLETE"
)

print(
    "=" * 100
)


print(
    qc.to_string(
        index=False
    )
)


print()

print(
    "Display classes:"
)


print(
    x[
        "display_class"
    ].value_counts().to_string()
)


print()

print(
    "Module-10 family occurrence in +/-10 kb plot data:"
)


print(

    x.loc[
        x[
            "neighbor_mmseq_cluster"
        ].isin(
            MODULE_CLUSTERS
        )
    ]

    .groupby(
        "neighbor_mmseq_cluster"
    )

    .agg(

        n_proteins=(
            "protein_id",
            "nunique"
        ),

        n_regions=(
            "plot_region_id",
            "nunique"
        ),

        n_genomes=(
            "genome",
            "nunique"
        ),

    )

    .reset_index()

    .sort_values(
        "neighbor_mmseq_cluster"
    )

    .to_string(
        index=False
    )

)


print()

print(
    f"Gene-map data:\n{OUTPUT_FILE}"
)

print()

print(
    f"Focal regions:\n{FOCAL_REGIONS_FILE}"
)

print()

print(
    f"QC:\n{QC_FILE}"
)
