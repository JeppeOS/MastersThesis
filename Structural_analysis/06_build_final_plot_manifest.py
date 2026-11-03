#!/usr/bin/env python3

from pathlib import Path
import pandas as pd


## ================================================================== ##
## PURPOSE
##
## Convert the standardized structure manifest + US-align results into
## the exact manifest used for final figure rendering.
##
## No structural calculations occur here.
##
## Each panel gets three standardized information lines:
##
##   1. Genome / taxon + source
##   2. Protein-family annotation
##   3. pTM + TM-score/RMSD relative to panel A
##
## Panel A is the reference and therefore has no pairwise metric.
##
## Original holo CIFs remain the rendering source so modeled Heme C
## ligands are retained.
## ================================================================== ##


BASE = (
    Path.home()
    / "methanotrophs"
    / "methanotroph_project"
    / "jeppe"
)


INFILE = (
    BASE
    / "structure_standardization"
    / "02_structural_comparison"
    / "tables"
    / "standard_structure_manifest_with_USalign.tsv"
)


OUTDIR = (
    BASE
    / "structure_standardization"
    / "03_plot_manifest"
)


OUTDIR.mkdir(
    parents=True,
    exist_ok=True
)


OUTFILE = (
    OUTDIR
    / "final_structure_plot_manifest.tsv"
)


## ================================================================== ##
## Read
## ================================================================== ##


df = pd.read_csv(
    INFILE,
    sep="\t"
)


if len(df) != 15:

    raise RuntimeError(
        f"Expected 15 structures; found {len(df)}."
    )


## ================================================================== ##
## Compact source labels for figure headings
## ================================================================== ##


SOURCE_LABEL = {

    "MAG":
        "MAG",

    "GlobDB genome":
        "GlobDB",

    "Reference genome":
        "reference",
}


df[
    "source_display"
] = df[
    "source_type"
].map(
    SOURCE_LABEL
)


if df[
    "source_display"
].isna().any():

    bad = df.loc[
        df[
            "source_display"
        ].isna(),
        "source_type"
    ].unique()

    raise RuntimeError(
        "Unknown source_type values: "
        f"{bad}"
    )


## ================================================================== ##
## Standardized displayed text
## ================================================================== ##


df[
    "display_line_1"
] = (

    df[
        "panel"
    ].astype(
        str
    )

    + ". "

    + df[
        "genome_label"
    ].astype(
        str
    )

    + " — "

    + df[
        "source_display"
    ]
)


df[
    "display_line_2"
] = df[
    "annotation_label"
]


def metric_line(
    row
):

    ptm = float(
        row[
            "ptm"
        ]
    )


    if bool(
        row[
            "reference_panel"
        ]
    ):

        return (
            f"pTM {ptm:.2f}  |  reference"
        )


    tm = float(
        row[
            "tm_score_vs_A"
        ]
    )


    rmsd = float(
        row[
            "rmsd_vs_A_A"
        ]
    )


    return (
        f"pTM {ptm:.2f}"
        f"  |  vs A: TM-score {tm:.3f}"
        f"  |  RMSD {rmsd:.2f} Å"
    )


df[
    "display_line_3"
] = df.apply(
    metric_line,
    axis=1
)


## ================================================================== ##
## Figure titles and subtitles
## ================================================================== ##


FIGURE_TEXT = {

    "cluster00035": {

        "title":
            "Cluster_00035 structural comparison",

        "subtitle":
            (
                "MtoA-annotated and FeGenie-negative "
                "homologues"
            ),
    },


    "cyc2": {

        "title":
            "Cyc2_repCluster2 structural comparison",

        "subtitle":
            (
                "Representative proteins from three "
                "MMseqs2 protein families"
            ),
    },


    "mca0421": {

        "title":
            "MCA0421 / Cluster_00063 structural comparison",

        "subtitle":
            (
                "MCA0421 and representative homologues "
                "across Methylococcales"
            ),
    },


    "cluster00121": {

        "title":
            "Cluster_00121 structural comparison",

        "subtitle":
            (
                "Representative nine-heme cytochromes "
                "from GlobDB and enrichment MAGs"
            ),
    },
}


df[
    "figure_display_title"
] = df[
    "figure_set"
].map(
    lambda x:
        FIGURE_TEXT[
            x
        ][
            "title"
        ]
)


df[
    "figure_display_subtitle"
] = df[
    "figure_set"
].map(
    lambda x:
        FIGURE_TEXT[
            x
        ][
            "subtitle"
        ]
)


## ================================================================== ##
## Standardized pLDDT palette
##
## Stored in each row for explicit provenance.
##
## Final renderer will use:
##
## >90     #0053D6
## 70-90   #65CBF3
## 50-70   #FFDB13
## <50     #FF7D45
## ================================================================== ##


df[
    "plddt_very_high"
] = "#0053D6"


df[
    "plddt_confident"
] = "#65CBF3"


df[
    "plddt_low"
] = "#FFDB13"


df[
    "plddt_very_low"
] = "#FF7D45"


## ================================================================== ##
## Explicit holo-rendering settings
##
## These make it impossible to accidentally switch to the protein-only
## alignment PDB during figure rendering.
## ================================================================== ##


df[
    "render_structure"
] = df[
    "structure_cif"
]


df[
    "alignment_structure"
] = df[
    "protein_pdb"
]


df[
    "render_heme_c"
] = True


df[
    "render_protein_by_plddt"
] = True


## ================================================================== ##
## Validation
## ================================================================== ##


if (
    df[
        "render_structure"
    ]
    ==
    df[
        "alignment_structure"
    ]
).any():

    raise RuntimeError(
        "A rendering structure is identical to an alignment PDB. "
        "Final rendering must use holo CIFs."
    )


for figure, group in df.groupby(
    "figure_set"
):

    if (
        group[
            "reference_panel"
        ].sum()
        != 1
    ):

        raise RuntimeError(
            f"{figure}: expected exactly one reference panel."
        )


    reference = group.loc[
        group[
            "reference_panel"
        ],
        "panel"
    ].iloc[
        0
    ]


    if reference != "A":

        raise RuntimeError(
            f"{figure}: reference panel is not A."
        )


## ================================================================== ##
## Keep useful provenance and output columns
## ================================================================== ##


columns = [

    "figure_set",
    "figure_display_title",
    "figure_display_subtitle",

    "panel",
    "reference_panel",

    "genome_label",
    "source_type",
    "source_display",

    "annotation_label",
    "short_annotation",
    "foldseek_context",

    "protein_id",
    "sequence_length",
    "expected_hemes",

    "ptm",
    "tm_score_vs_A",
    "rmsd_vs_A_A",
    "aligned_length_vs_A",
    "seq_identity_vs_A",

    "display_line_1",
    "display_line_2",
    "display_line_3",

    "render_structure",
    "alignment_structure",
    "render_heme_c",
    "render_protein_by_plddt",

    "plddt_very_high",
    "plddt_confident",
    "plddt_low",
    "plddt_very_low",

    "source_workflow",
    "source_root",
    "summary_json",
]


df = df[
    columns
]


df.to_csv(
    OUTFILE,
    sep="\t",
    index=False
)


## ================================================================== ##
## Console preview
## ================================================================== ##


print()

print(
    "=" * 110
)

print(
    "FINAL STRUCTURE PLOT MANIFEST"
)

print(
    "=" * 110
)


for figure, group in df.groupby(
    "figure_set",
    sort=False
):

    print()

    print(
        group[
            "figure_display_title"
        ].iloc[
            0
        ]
    )

    print(
        group[
            "figure_display_subtitle"
        ].iloc[
            0
        ]
    )

    print(
        "-" * 90
    )


    for _, row in group.iterrows():

        print(
            row[
                "display_line_1"
            ]
        )

        print(
            "  "
            + row[
                "display_line_2"
            ]
        )

        print(
            "  "
            + row[
                "display_line_3"
            ]
        )

        print()


print(
    "=" * 110
)

print()

print(
    "Canonical pLDDT colours:"
)

print(
    "  Very high (>90):    #0053D6"
)

print(
    "  Confident (70-90):  #65CBF3"
)

print(
    "  Low (50-70):        #FFDB13"
)

print(
    "  Very low (<50):     #FF7D45"
)

print()

print(
    "Final rendering source:"
)

print(
    "  original mature holo AlphaFold3 CIF"
)

print()

print(
    "Structural-alignment source:"
)

print(
    "  protein-only PDB"
)

print()

print(
    "Heme C ligands WILL be retained in final figures."
)

print()

print(
    f"Written:\n{OUTFILE}"
)
