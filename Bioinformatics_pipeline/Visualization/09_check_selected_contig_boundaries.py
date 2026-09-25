#!/usr/bin/env python3


## ================================================================== ##
## CROSS-FIGURE CONTIG-BOUNDARY QC
##
## Purpose
## -------
##
## Determine whether apparently short gene-map neighborhoods are short
## because:
##
##   A. the focal gene lies close to a contig boundary, or
##
##   B. the contig continues but no protein-coding CDS is present in
##      part of the displayed window.
##
##
## The script examines the currently plotted:
##
##   Module 20
##   Module 10
##   Module 35
##
## neighborhoods.
##
##
## Coordinates are normalized so that the focal gene points right.
##
## Output:
##
##   selected_gene_map_contig_boundary_qc.tsv
##
## This output can later be joined to the R plotting scripts to draw a
## thin sequence-availability baseline underneath each gene track.
## ================================================================== ##


from pathlib import Path
import gzip

import pandas as pd


## ================================================================== ##
## Paths
## ================================================================== ##

WORKFLOW = Path(
    "~/methanotrophs/methanotroph_project/jeppe/"
    "genome_analysis_workflow"
).expanduser()


FASTA_DIR = Path(
    "~/methanotrophs/methanotroph_project/jeppe/"
    "methylococcales_fastas"
).expanduser()


OUTPUT_DIR = (
    WORKFLOW
    / "16_visualization"
    / "16C_cross_figure_qc"
)


OUTPUT_FILE = (
    OUTPUT_DIR
    / "selected_gene_map_contig_boundary_qc.tsv"
)


DATASETS = [

    {
        "figure": "Module_20",

        "file":
            WORKFLOW
            / "16_visualization"
            / "16A3_module20_representative_selection"
            / "Module_20_selected_gene_map_data.tsv",

        "focal_markers": [
            (
                "display_class",
                "focal_cluster00035"
            ),
        ],
    },

    {
        "figure": "Module_10",

        "file":
            WORKFLOW
            / "16_visualization"
            / "16A3_module10_representative_selection"
            / "Module_10_selected_gene_map_data.tsv",

        "focal_markers": [
            (
                "display_class",
                "focal_cluster00050"
            ),
        ],
    },

    {
        "figure": "Module_35",

        "file":
            WORKFLOW
            / "16_visualization"
            / "16A2_local_architecture_pipeline"
            / "runs"
            / "Module_35"
            / "Module_35_gene_map_data_with_local_families.tsv",

        "focal_markers": [
            (
                "local_family_id",
                "M35_local_001"
            ),
            (
                "display_class",
                "focal_cyc2"
            ),
        ],
    },

]


## ================================================================== ##
## Helpers
## ================================================================== ##

def fail(message):

    raise RuntimeError(
        message
    )


def clean(value):

    if pd.isna(
        value
    ):

        return ""

    return str(
        value
    ).strip()


def open_fasta(path):

    if str(
        path
    ).endswith(
        ".gz"
    ):

        return gzip.open(
            path,
            "rt"
        )

    return open(
        path,
        "rt"
    )


def find_genome_fasta(genome):

    candidates = [

        FASTA_DIR
        / f"{genome}.fa.gz",

        FASTA_DIR
        / f"{genome}.fa",

        FASTA_DIR
        / f"{genome}.fna.gz",

        FASTA_DIR
        / f"{genome}.fna",

        FASTA_DIR
        / f"{genome}.fasta.gz",

        FASTA_DIR
        / f"{genome}.fasta",

    ]


    existing = [

        path

        for path in candidates

        if path.exists()

    ]


    if len(
        existing
    ) == 0:

        fail(
            f"No genome FASTA found for {genome}."
        )


    if len(
        existing
    ) > 1:

        fail(
            "More than one possible FASTA found for "
            f"{genome}:\n"
            +
            "\n".join(
                str(
                    path
                )
                for path in existing
            )
        )


    return existing[
        0
    ]


def read_fasta_lengths(path):

    lengths = {}

    current_id = None

    current_length = 0


    with open_fasta(
        path
    ) as handle:

        for raw in handle:

            line = raw.strip()


            if not line:

                continue


            if line.startswith(
                ">"
            ):

                if current_id is not None:

                    lengths[
                        current_id
                    ] = current_length


                current_id = (
                    line[
                        1:
                    ]
                    .split()[0]
                )

                current_length = 0


            else:

                if current_id is None:

                    fail(
                        f"FASTA sequence before header in {path}"
                    )


                current_length += len(
                    line
                )


    if current_id is not None:

        lengths[
            current_id
        ] = current_length


    return lengths


def choose_region_column(df):

    for col in [

        "plot_region_id",
        "focal_region_id",

    ]:

        if col in df.columns:

            return col


    ## Module 35 currently has one plotted focal region per genome.
    return "genome"


def get_focal_mask(
    df,
    marker_candidates
):

    for column, value in marker_candidates:

        if column not in df.columns:

            continue


        mask = (
            df[
                column
            ]
            .fillna("")
            .astype(str)
            .str.strip()
            ==
            value
        )


        if mask.any():

            return (
                mask,
                f"{column}={value}"
            )


    ## -------------------------------------------------------------- ##
    ## Generic fallback: explicit focal marker.
    ## -------------------------------------------------------------- ##

    if "is_focal_plot" in df.columns:

        values = pd.to_numeric(
            df[
                "is_focal_plot"
            ],
            errors="coerce"
        ).fillna(
            0
        )


        if (
            values == 1
        ).any():

            return (
                values == 1,
                "is_focal_plot=1"
            )


    if "is_focal" in df.columns:

        values = pd.to_numeric(
            df[
                "is_focal"
            ],
            errors="coerce"
        ).fillna(
            0
        )


        if (
            values == 1
        ).any():

            return (
                values == 1,
                "is_focal=1"
            )


    ## -------------------------------------------------------------- ##
    ## Final fallback:
    ## protein_id equals focal_protein_id.
    ## -------------------------------------------------------------- ##

    if {
        "protein_id",
        "focal_protein_id"
    }.issubset(
        df.columns
    ):

        mask = (
            df[
                "protein_id"
            ].astype(str)
            ==
            df[
                "focal_protein_id"
            ].astype(str)
        )


        if mask.any():

            return (
                mask,
                "protein_id=focal_protein_id"
            )


    fail(
        "Could not identify focal rows."
    )


def focal_raw_coordinates(row):

    ## -------------------------------------------------------------- ##
    ## Prefer the neighbor/raw gene coordinates from the focal row.
    ## -------------------------------------------------------------- ##

    if all(
        col in row.index
        for col in [
            "neighbor_contig",
            "neighbor_start",
            "neighbor_end",
            "neighbor_strand",
        ]
    ):

        contig = clean(
            row[
                "neighbor_contig"
            ]
        )

        start = pd.to_numeric(
            row[
                "neighbor_start"
            ],
            errors="coerce"
        )

        end = pd.to_numeric(
            row[
                "neighbor_end"
            ],
            errors="coerce"
        )

        strand = clean(
            row[
                "neighbor_strand"
            ]
        )


        if (
            contig != ""
            and
            not pd.isna(
                start
            )
            and
            not pd.isna(
                end
            )
            and
            strand in {
                "+",
                "-"
            }
        ):

            return (
                contig,
                float(
                    start
                ),
                float(
                    end
                ),
                strand,
                "neighbor_coordinates"
            )


    ## -------------------------------------------------------------- ##
    ## Fallback to focal_* metadata.
    ## -------------------------------------------------------------- ##

    required = [

        "focal_contig",
        "focal_start",
        "focal_end",
        "focal_strand",

    ]


    if not all(
        col in row.index
        for col in required
    ):

        fail(
            "Neither focal-row neighbor coordinates nor focal_* "
            "coordinates were available."
        )


    contig = clean(
        row[
            "focal_contig"
        ]
    )

    start = pd.to_numeric(
        row[
            "focal_start"
        ],
        errors="coerce"
    )

    end = pd.to_numeric(
        row[
            "focal_end"
        ],
        errors="coerce"
    )

    strand = clean(
        row[
            "focal_strand"
        ]
    )


    if (
        contig == ""
        or
        pd.isna(
            start
        )
        or
        pd.isna(
            end
        )
        or
        strand not in {
            "+",
            "-"
        }
    ):

        fail(
            "Invalid focal coordinate metadata."
        )


    return (
        contig,
        float(
            start
        ),
        float(
            end
        ),
        strand,
        "focal_coordinates"
    )


## ================================================================== ##
## Main
## ================================================================== ##

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)


results = []


for config in DATASETS:

    figure = config[
        "figure"
    ]

    path = config[
        "file"
    ]


    print()

    print(
        "=" * 110
    )

    print(
        figure
    )

    print(
        "=" * 110
    )


    if not path.exists():

        fail(
            f"Input missing:\n{path}"
        )


    df = pd.read_csv(
        path,
        sep="\t",
        dtype=str,
        keep_default_na=False
    )


    region_col = choose_region_column(
        df
    )


    focal_mask, focal_method = get_focal_mask(

        df,

        config[
            "focal_markers"
        ]

    )


    focals = df.loc[
        focal_mask
    ].copy()


    ## -------------------------------------------------------------- ##
    ## Remove accidental exact repeat focal rows before checking.
    ## -------------------------------------------------------------- ##

    focal_id_col = (

        "protein_id"

        if "protein_id" in focals.columns

        else "neighbor_protein_id"

    )


    focals = focals.drop_duplicates(

        [
            region_col,
            focal_id_col
        ]

    )


    focal_counts = (

        focals

        .groupby(
            region_col
        )

        .size()

    )


    bad = focal_counts.loc[
        focal_counts != 1
    ]


    if not bad.empty:

        print(
            bad
        )

        fail(
            f"{figure}: expected exactly one focal gene per region."
        )


    fasta_cache = {}


    for _, focal in focals.iterrows():

        region = clean(
            focal[
                region_col
            ]
        )

        genome = clean(
            focal[
                "genome"
            ]
        )


        region_genes = df.loc[
            df[
                region_col
            ].astype(str)
            ==
            region
        ].copy()


        if region_genes.empty:

            fail(
                f"{figure}: no genes found for region {region}"
            )


        (
            contig,
            focal_start,
            focal_end,
            focal_strand,
            coordinate_source,

        ) = focal_raw_coordinates(
            focal
        )


        ## ---------------------------------------------------------- ##
        ## Load contig lengths from authoritative genome FASTA.
        ## ---------------------------------------------------------- ##

        if genome not in fasta_cache:

            fasta_path = find_genome_fasta(
                genome
            )

            fasta_cache[
                genome
            ] = (
                fasta_path,
                read_fasta_lengths(
                    fasta_path
                )
            )


        fasta_path, lengths = fasta_cache[
            genome
        ]


        if contig not in lengths:

            print()

            print(
                f"Genome: {genome}"
            )

            print(
                f"Requested contig: {contig}"
            )

            print(
                "First FASTA contig IDs:"
            )

            print(
                "\n".join(
                    list(
                        lengths.keys()
                    )[
                        :20
                    ]
                )
            )


            fail(
                f"Contig {contig!r} was not found in FASTA "
                f"{fasta_path}"
            )


        contig_length = lengths[
            contig
        ]


        ## ---------------------------------------------------------- ##
        ## Focal midpoint in original genomic coordinates.
        ## ---------------------------------------------------------- ##

        focal_midpoint = (
            focal_start
            +
            focal_end
        ) / 2.0


        genomic_left_available = (
            focal_midpoint
            -
            1.0
        )


        genomic_right_available = (
            contig_length
            -
            focal_midpoint
        )


        ## ---------------------------------------------------------- ##
        ## Convert availability to the orientation used in gene maps.
        ##
        ## Every focal gene points to the right in the plots.
        ## ---------------------------------------------------------- ##

        if focal_strand == "+":

            plot_left_available = (
                genomic_left_available
            )

            plot_right_available = (
                genomic_right_available
            )


        else:

            plot_left_available = (
                genomic_right_available
            )

            plot_right_available = (
                genomic_left_available
            )


        contig_plot_start_bp = (
            -plot_left_available
        )


        contig_plot_end_bp = (
            plot_right_available
        )


        ## ---------------------------------------------------------- ##
        ## Actual extent of CDS arrows currently present in input.
        ## ---------------------------------------------------------- ##

        starts = pd.to_numeric(
            region_genes[
                "plot_start_bp"
            ],
            errors="coerce"
        )


        ends = pd.to_numeric(
            region_genes[
                "plot_end_bp"
            ],
            errors="coerce"
        )


        if (
            starts.isna().any()
            or
            ends.isna().any()
        ):

            fail(
                f"{figure} / {region}: non-numeric plot coordinates."
            )


        leftmost_gene_bp = starts.min()

        rightmost_gene_bp = ends.max()


        ## ---------------------------------------------------------- ##
        ## Taxon / track label
        ## ---------------------------------------------------------- ##

        if "track_label" in focal.index:

            track_label = clean(
                focal[
                    "track_label"
                ]
            )

        elif "taxonomy_display" in focal.index:

            track_label = clean(
                focal[
                    "taxonomy_display"
                ]
            )

        else:

            track_label = genome


        plot_order = ""


        if "plot_order" in focal.index:

            plot_order = clean(
                focal[
                    "plot_order"
                ]
            )


        ## ---------------------------------------------------------- ##
        ## Censoring at useful reference windows.
        ##
        ## We report both 10 kb and 15 kb because the current figures
        ## do not all necessarily use identical displayed spans.
        ## ---------------------------------------------------------- ##

        left_censored_10kb = (
            plot_left_available
            <
            10_000
        )

        right_censored_10kb = (
            plot_right_available
            <
            10_000
        )


        left_censored_15kb = (
            plot_left_available
            <
            15_000
        )

        right_censored_15kb = (
            plot_right_available
            <
            15_000
        )


        ## ---------------------------------------------------------- ##
        ## If a contig does reach +/-10 kb, quantify how much of the
        ## outer 10-kb plotting window contains no plotted CDS.
        ##
        ## These are NOT necessarily truly gene-free sequence if some
        ## other feature type was not part of the protein catalogue.
        ## They mean only "no plotted protein-coding CDS".
        ## ---------------------------------------------------------- ##

        if not left_censored_10kb:

            left_no_plotted_CDS_tail_10kb = max(
                0.0,
                leftmost_gene_bp
                -
                (
                    -10_000
                )
            )

        else:

            left_no_plotted_CDS_tail_10kb = None


        if not right_censored_10kb:

            right_no_plotted_CDS_tail_10kb = max(
                0.0,
                10_000
                -
                rightmost_gene_bp
            )

        else:

            right_no_plotted_CDS_tail_10kb = None


        results.append(
            {

                "figure":
                    figure,

                "plot_order":
                    plot_order,

                "genome":
                    genome,

                "track_label":
                    track_label,

                "region_id":
                    region,

                "focal_marker_method":
                    focal_method,

                "coordinate_source":
                    coordinate_source,

                "focal_contig":
                    contig,

                "contig_length_bp":
                    contig_length,

                "focal_start":
                    focal_start,

                "focal_end":
                    focal_end,

                "focal_strand":
                    focal_strand,

                "focal_midpoint":
                    focal_midpoint,

                "plot_left_available_bp":
                    plot_left_available,

                "plot_right_available_bp":
                    plot_right_available,

                "contig_plot_start_bp":
                    contig_plot_start_bp,

                "contig_plot_end_bp":
                    contig_plot_end_bp,

                "leftmost_plotted_gene_bp":
                    leftmost_gene_bp,

                "rightmost_plotted_gene_bp":
                    rightmost_gene_bp,

                "left_censored_at_10kb":
                    left_censored_10kb,

                "right_censored_at_10kb":
                    right_censored_10kb,

                "left_censored_at_15kb":
                    left_censored_15kb,

                "right_censored_at_15kb":
                    right_censored_15kb,

                "left_no_plotted_CDS_tail_within_10kb_bp":
                    left_no_plotted_CDS_tail_10kb,

                "right_no_plotted_CDS_tail_within_10kb_bp":
                    right_no_plotted_CDS_tail_10kb,

                "genome_fasta":
                    str(
                        fasta_path
                    ),

            }
        )


## ================================================================== ##
## Output
## ================================================================== ##

out = pd.DataFrame(
    results
)


out[
    "plot_order_numeric"
] = pd.to_numeric(
    out[
        "plot_order"
    ],
    errors="coerce"
)


out = out.sort_values(

    [
        "figure",
        "plot_order_numeric",
        "track_label",
    ],

    na_position="last"

).drop(
    columns=[
        "plot_order_numeric"
    ]
)


out.to_csv(
    OUTPUT_FILE,
    sep="\t",
    index=False
)


## ================================================================== ##
## Console report
## ================================================================== ##

report = out[
    [
        "figure",
        "plot_order",
        "track_label",

        "focal_contig",
        "contig_length_bp",

        "plot_left_available_bp",
        "plot_right_available_bp",

        "leftmost_plotted_gene_bp",
        "rightmost_plotted_gene_bp",

        "left_censored_at_10kb",
        "right_censored_at_10kb",

        "left_censored_at_15kb",
        "right_censored_at_15kb",
    ]
].copy()


for col in [

    "plot_left_available_bp",
    "plot_right_available_bp",

    "leftmost_plotted_gene_bp",
    "rightmost_plotted_gene_bp",

]:

    report[
        col
    ] = pd.to_numeric(
        report[
            col
        ],
        errors="coerce"
    ).round(
        1
    )


print()

print(
    "=" * 150
)

print(
    "CONTIG-BOUNDARY QC"
)

print(
    "=" * 150
)


print(
    report.to_string(
        index=False
    )
)


print()

print(
    "=" * 150
)

print(
    "REGIONS CENSORED WITHIN +/-10 kb"
)

print(
    "=" * 150
)


censored = out.loc[

    out[
        "left_censored_at_10kb"
    ]

    |

    out[
        "right_censored_at_10kb"
    ]

]


if censored.empty:

    print(
        "[none]"
    )

else:

    print(

        censored[
            [
                "figure",
                "plot_order",
                "track_label",

                "plot_left_available_bp",
                "plot_right_available_bp",

                "left_censored_at_10kb",
                "right_censored_at_10kb",
            ]
        ].to_string(
            index=False
        )

    )


print()

print(
    f"Full QC table:\n{OUTPUT_FILE}"
)
