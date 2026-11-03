#!/usr/bin/env python3

from pathlib import Path
import pandas as pd


## ================================================================== ##
## PURPOSE
##
## Build the authoritative manifest for the standardized AlphaFold3
## figures.
##
## Input:
##   Stage-34 source audit
##
## Output:
##   One row per plotted structure containing:
##
##   - figure set
##   - panel
##   - genome/taxon display label
##   - source type
##   - protein-family annotation
##   - protein ID
##   - heme count
##   - sequence length
##   - true AlphaFold3 pTM
##   - selected structure path
##   - source-workflow provenance
##   - reference-panel designation
##
## Structural similarity metrics are NOT added here.
##
## They will be calculated independently with US-align in Stage 36,
## then joined back onto this manifest.
## ================================================================== ##


BASE = (
    Path.home()
    / "methanotrophs"
    / "methanotroph_project"
    / "jeppe"
)


AUDIT = (
    BASE
    / "structure_standardization"
    / "00_audit"
    / "alphafold_structure_source_audit.tsv"
)


OUTDIR = (
    BASE
    / "structure_standardization"
    / "01_manifest"
)


OUTDIR.mkdir(
    parents=True,
    exist_ok=True
)


OUTFILE = (
    OUTDIR
    / "standard_structure_manifest.tsv"
)


## ================================================================== ##
## Read source audit
## ================================================================== ##


audit = pd.read_csv(
    AUDIT,
    sep="\t"
)


if len(audit) != 15:

    raise RuntimeError(
        f"Expected 15 audited structures; "
        f"found {len(audit)}."
    )


if not audit[
    "ptm_matches_thesis"
].all():

    raise RuntimeError(
        "Not all audited pTM values match the thesis values. "
        "Resolve Stage 34 before continuing."
    )


## ================================================================== ##
## Standardized biological labels
##
## Design principle:
##
## Line 1 in final figure:
##     panel + genome/taxon
##
## Line 2:
##     protein family / annotation
##
## Line 3:
##     pTM + structural-comparison metric
##
## Protein IDs, lengths and heme counts remain available in the
## manifest but do not need to dominate the visual panel heading.
## ================================================================== ##


LABELS = {


    ## ============================================================= ##
    ## Cluster_00035
    ## Original genome_analysis_workflow
    ## ============================================================= ##

    ("cluster00035", "A"): {

        "figure_title":
            "Cluster_00035",

        "source_type":
            "GlobDB genome",

        "genome_label":
            "GCA_903905335",

        "annotation_label":
            "MtoA-annotated Cluster_00035 homologue",

        "short_annotation":
            "MtoA",

        "reference_panel":
            True,
    },


    ("cluster00035", "B"): {

        "figure_title":
            "Cluster_00035",

        "source_type":
            "GlobDB genome",

        "genome_label":
            "Methylobacter_C sp003584895",

        "annotation_label":
            "FeGenie-negative Cluster_00035 homologue",

        "short_annotation":
            "FeGenie-negative homologue",

        "reference_panel":
            False,
    },


    ## ============================================================= ##
    ## Cyc2
    ## ============================================================= ##

    ("cyc2", "A"): {

        "figure_title":
            "Cyc2 representatives",

        "source_type":
            "MAG",

        "genome_label":
            "Methylobacter_A",

        "annotation_label":
            "Cluster_00024 — FeGenie: Cyc2_repCluster2",

        "short_annotation":
            "Cluster_00024",

        "reference_panel":
            True,
    },


    ("cyc2", "B"): {

        "figure_title":
            "Cyc2 representatives",

        "source_type":
            "MAG",

        "genome_label":
            "Methylumidiphilus",

        "annotation_label":
            "Cluster_00050 — FeGenie: Cyc2_repCluster2",

        "short_annotation":
            "Cluster_00050",

        "reference_panel":
            False,
    },


    ("cyc2", "C"): {

        "figure_title":
            "Cyc2 representatives",

        "source_type":
            "GlobDB genome",

        "genome_label":
            "Methylobacter_A sp023229605",

        "annotation_label":
            "Cluster_00050 — FeGenie: Cyc2_repCluster2",

        "short_annotation":
            "Cluster_00050",

        "reference_panel":
            False,
    },


    ("cyc2", "D"): {

        "figure_title":
            "Cyc2 representatives",

        "source_type":
            "GlobDB genome",

        "genome_label":
            "PGZD01 sp002843155",

        "annotation_label":
            "Cluster_00206 — FeGenie: Cyc2_repCluster2",

        "short_annotation":
            "Cluster_00206",

        "reference_panel":
            False,
    },


    ## ============================================================= ##
    ## MCA0421 / Cluster_00063
    ## Original genome_analysis_workflow
    ## ============================================================= ##

    ("mca0421", "A"): {

        "figure_title":
            "MCA0421 / Cluster_00063",

        "source_type":
            "Reference genome",

        "genome_label":
            "Methylococcus capsulatus Bath",

        "annotation_label":
            "MCA0421 — Cluster_00063 9-heme cytochrome",

        "short_annotation":
            "MCA0421",

        "reference_panel":
            True,
    },


    ("mca0421", "B"): {

        "figure_title":
            "MCA0421 / Cluster_00063",

        "source_type":
            "GlobDB genome",

        "genome_label":
            "Methylococcus capsulatus",

        "annotation_label":
            "Cluster_00063 — MCA0421 homologue",

        "short_annotation":
            "Cluster_00063",

        "reference_panel":
            False,
    },


    ("mca0421", "C"): {

        "figure_title":
            "MCA0421 / Cluster_00063",

        "source_type":
            "GlobDB genome",

        "genome_label":
            "Methylococcus sp. 016106025",

        "annotation_label":
            "Cluster_00063 — MCA0421 homologue",

        "short_annotation":
            "Cluster_00063",

        "reference_panel":
            False,
    },


    ("mca0421", "D"): {

        "figure_title":
            "MCA0421 / Cluster_00063",

        "source_type":
            "GlobDB genome",

        "genome_label":
            "Methylobacter_A sp. 026399025",

        "annotation_label":
            "Cluster_00063 — MCA0421 homologue",

        "short_annotation":
            "Cluster_00063",

        "reference_panel":
            False,
    },


    ## ============================================================= ##
    ## Cluster_00121
    ## ============================================================= ##

    ("cluster00121", "A"): {

        "figure_title":
            "Cluster_00121",

        "source_type":
            "GlobDB genome",

        "genome_label":
            "SPIREOTU_00519314 BCRBG_22465",

        "annotation_label":
            "Cluster_00121 — 9-heme cytochrome",

        "short_annotation":
            "Cluster_00121",

        "reference_panel":
            True,
    },


    ("cluster00121", "B"): {

        "figure_title":
            "Cluster_00121",

        "source_type":
            "MAG",

        "genome_label":
            "Methylumidiphilus",

        "annotation_label":
            "Cluster_00121 — 9-heme cytochrome",

        "short_annotation":
            "Cluster_00121",

        "reference_panel":
            False,
    },


    ("cluster00121", "C"): {

        "figure_title":
            "Cluster_00121",

        "source_type":
            "MAG",

        "genome_label":
            "Methylumidiphilus",

        "annotation_label":
            "Cluster_00121 — 9-heme cytochrome",

        "short_annotation":
            "Cluster_00121",

        "reference_panel":
            False,
    },


    ("cluster00121", "D"): {

        "figure_title":
            "Cluster_00121",

        "source_type":
            "MAG",

        "genome_label":
            "Methylovulum",

        "annotation_label":
            "Cluster_00121 — 9-heme cytochrome",

        "short_annotation":
            "Cluster_00121",

        "reference_panel":
            False,
    },


    ("cluster00121", "E"): {

        "figure_title":
            "Cluster_00121",

        "source_type":
            "MAG",

        "genome_label":
            "Methylomonas",

        "annotation_label":
            "Cluster_00121 — 9-heme cytochrome",

        "short_annotation":
            "Cluster_00121",

        "reference_panel":
            False,
    },
}


## ================================================================== ##
## Optional Foldseek context
##
## This is retained in the manifest but deliberately separated from
## the primary annotation because Foldseek similarity is not direct
## functional validation.
## ================================================================== ##


FOLDSEEK_CONTEXT = {

    ("cluster00035", "A"):
        "MtrA-like structural hits",

    ("cluster00035", "B"):
        "MtrA-like structural hits",


    ("cluster00121", "A"):
        "IPT/TIG-like; MtrC/MtrF-like structural hits",

    ("cluster00121", "B"):
        "IPT/TIG-like; MtrC/MtrF-like structural hits",

    ("cluster00121", "C"):
        "IPT/TIG-like; MtrC/MtrF-like structural hits",

    ("cluster00121", "D"):
        "IPT/TIG-like; MtrC/MtrF-like structural hits",

    ("cluster00121", "E"):
        "IPT/TIG-like; MtrC/MtrF-like structural hits",
}


## ================================================================== ##
## Build standardized manifest
## ================================================================== ##


rows = []


for _, row in audit.iterrows():

    figure_set = str(
        row[
            "figure_set"
        ]
    )


    panel = str(
        row[
            "panel"
        ]
    )


    key = (
        figure_set,
        panel
    )


    if key not in LABELS:

        raise RuntimeError(
            f"No standardized labels defined for "
            f"{figure_set} panel {panel}."
        )


    label = LABELS[
        key
    ]


    source_suffix = (
        f" — {label['source_type']}"
        if label[
            "source_type"
        ]
        else ""
    )


    panel_heading = (
        f"{panel}. "
        f"{label['genome_label']}"
        f"{source_suffix}"
    )


    rows.append({

        ## ---------------------------------------------------------- ##
        ## Plot ordering
        ## ---------------------------------------------------------- ##

        "figure_set":
            figure_set,

        "figure_title":
            label[
                "figure_title"
            ],

        "panel":
            panel,

        "reference_panel":
            bool(
                label[
                    "reference_panel"
                ]
            ),


        ## ---------------------------------------------------------- ##
        ## Standardized displayed information
        ## ---------------------------------------------------------- ##

        "panel_heading":
            panel_heading,

        "genome_label":
            label[
                "genome_label"
            ],

        "source_type":
            label[
                "source_type"
            ],

        "annotation_label":
            label[
                "annotation_label"
            ],

        "short_annotation":
            label[
                "short_annotation"
            ],

        "foldseek_context":
            FOLDSEEK_CONTEXT.get(
                key,
                ""
            ),


        ## ---------------------------------------------------------- ##
        ## Protein information
        ## ---------------------------------------------------------- ##

        "protein_id":
            row[
                "protein_id"
            ],

        "sequence_length":
            row[
                "sequence_length"
            ],

        "expected_hemes":
            row[
                "expected_hemes"
            ],


        ## ---------------------------------------------------------- ##
        ## AlphaFold model information
        ## ---------------------------------------------------------- ##

        "model_number":
            int(
                row[
                    "model_number"
                ]
            ),

        "ptm":
            float(
                row[
                    "ptm"
                ]
            ),

        "structure_cif":
            row[
                "structure_cif"
            ],

        "summary_json":
            row[
                "summary_json"
            ],


        ## ---------------------------------------------------------- ##
        ## Provenance
        ## ---------------------------------------------------------- ##

        "source_workflow":
            row[
                "source_workflow"
            ],

        "source_root":
            row[
                "source_root"
            ],

        "genome_id":
            row[
                "genome_id"
            ],
    })


manifest = pd.DataFrame(
    rows
)


## ================================================================== ##
## Figure ordering
## ================================================================== ##


figure_order = {

    "cluster00035":
        1,

    "cyc2":
        2,

    "mca0421":
        3,

    "cluster00121":
        4,
}


manifest[
    "_figure_order"
] = manifest[
    "figure_set"
].map(
    figure_order
)


manifest = manifest.sort_values(
    [
        "_figure_order",
        "panel"
    ]
).drop(
    columns=[
        "_figure_order"
    ]
).reset_index(
    drop=True
)


## ================================================================== ##
## Validation
## ================================================================== ##


if len(
    manifest
) != 15:

    raise RuntimeError(
        f"Expected 15 standardized structures; "
        f"found {len(manifest)}."
    )


## Exactly one reference per figure set

reference_counts = (
    manifest
    .groupby(
        "figure_set"
    )
    ["reference_panel"]
    .sum()
)


bad_reference = (
    reference_counts[
        reference_counts != 1
    ]
)


if len(
    bad_reference
) > 0:

    raise RuntimeError(
        "Each figure set must contain exactly one "
        "reference panel:\n"
        +
        bad_reference.to_string()
    )


## Every reference must be panel A

bad_ref_panel = manifest[
    manifest[
        "reference_panel"
    ]
    &
    (
        manifest[
            "panel"
        ] != "A"
    )
]


if len(
    bad_ref_panel
) > 0:

    raise RuntimeError(
        "Reference panel is not A for:\n"
        +
        bad_ref_panel[
            [
                "figure_set",
                "panel"
            ]
        ].to_string(
            index=False
        )
    )


## No missing pTM

if manifest[
    "ptm"
].isna().any():

    raise RuntimeError(
        "Missing pTM values in standardized manifest."
    )


## No duplicate protein IDs within a figure set

duplicates = manifest[
    manifest.duplicated(
        subset=[
            "figure_set",
            "protein_id"
        ],
        keep=False
    )
]


if len(
    duplicates
) > 0:

    raise RuntimeError(
        "Duplicate proteins within a figure set:\n"
        +
        duplicates[
            [
                "figure_set",
                "panel",
                "protein_id"
            ]
        ].to_string(
            index=False
        )
    )


## ================================================================== ##
## Write manifest
## ================================================================== ##


manifest.to_csv(
    OUTFILE,
    sep="\t",
    index=False
)


## ================================================================== ##
## Console report
## ================================================================== ##


print()

print(
    "=" * 120
)

print(
    "STANDARDIZED STRUCTURE MANIFEST"
)

print(
    "=" * 120
)

print()


display = manifest[
    [
        "figure_set",
        "panel",
        "genome_label",
        "source_type",
        "short_annotation",
        "protein_id",
        "expected_hemes",
        "ptm",
        "reference_panel"
    ]
]


print(
    display.to_string(
        index=False
    )
)


print()

print(
    "=" * 120
)


for figure_set, group in manifest.groupby(
    "figure_set",
    sort=False
):

    print()

    print(
        figure_set
    )

    print(
        "-" * len(
            figure_set
        )
    )


    for _, row in group.iterrows():

        print(
            row[
                "panel_heading"
            ]
        )

        print(
            "    "
            + row[
                "annotation_label"
            ]
        )

        print(
            f"    pTM {row['ptm']:.2f}"
            +
            (
                " | reference"
                if row[
                    "reference_panel"
                ]
                else
                " | comparison to A pending US-align"
            )
        )


print()

print(
    "=" * 120
)

print(
    f"Manifest written:\n{OUTFILE}"
)

print()

print(
    "Stage 36 can now calculate structural similarity "
    "against panel A for each figure set."
)
