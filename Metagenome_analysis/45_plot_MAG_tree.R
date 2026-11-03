#!/usr/bin/env Rscript


## ============================================================
## MAG + ALLOCRENOTHRIX PHYLOGENY
##
## Tips:
##   16 methanotroph MAGs
##   2 Allocrenothrix methanica genomes
##   1 Methylophaga outgroup
##
## Tree:
##   pruned from the validated rooted 650-genome bac120 tree
##
## Tracks:
##   1. Methane monooxygenase repertoire
##   2. Completeness
##   3. Contamination
## ============================================================


suppressPackageStartupMessages({

    library(tidyverse)
    library(ape)
    library(ggtree)
    library(ggnewscale)

})


## ============================================================
## PATHS
## ============================================================

PROJECT <- file.path(
    Sys.getenv("HOME"),
    "methanotrophs",
    "methanotroph_project",
    "jeppe",
    "metagenome"
)


PHYLO <- file.path(
    PROJECT,
    "comparative_analysis",
    "final_species_phylogeny"
)


TREE_FILE <- file.path(
    PHYLO,
    "05_rooted",
    "final_methylococcales_GlobDB_MAGs_Umezawa_bac120_ROOTED.tree"
)


META_FILE <- file.path(
    PHYLO,
    "08_visualization",
    "final_650_tip_metadata.tsv"
)


MMO_FILE <- file.path(
    PHYLO,
    "08_visualization",
    "extra_tracks",
    "all_650_pmo_pxm_smmo_QC.tsv"
)


OUTDIR <- file.path(
    PHYLO,
    "08_visualization",
    "MAG_subtree"
)


dir.create(
    OUTDIR,
    recursive = TRUE,
    showWarnings = FALSE
)


## ============================================================
## OUTPUT FILES
## ============================================================

OUT_PDF <- file.path(
    OUTDIR,
    "MAGs_Allocrenothrix_outgroup_phylogeny.pdf"
)


OUT_PNG <- file.path(
    OUTDIR,
    "MAGs_Allocrenothrix_outgroup_phylogeny.png"
)


OUT_TREE <- file.path(
    OUTDIR,
    "MAGs_Allocrenothrix_outgroup_subtree.tree"
)


OUT_TRACKS <- file.path(
    OUTDIR,
    "MAGs_Allocrenothrix_outgroup_tracks.tsv"
)


## ============================================================
## COLORS
##
## Same basic logic as the final 650-tip tree.
## ============================================================

MMO_COLORS <- c(

    "Particulate MMO" =
        "#7CAE00",      # lime green

    "sMMO" =
        "#C77CFF",      # lavender

    "Both" =
        "#111111",      # near-black

    "Not detected" =
        "#D9D9D9"       # light grey
)


COMP_LOW <- "#F7FBFF"
COMP_HIGH <- "#08519C"


CONT_LOW <- "#FFF5F0"
CONT_HIGH <- "#A50F15"


TRACK_OUTLINE <- "grey35"


## ============================================================
## READ TREE
## ============================================================

full_tree <- read.tree(
    TREE_FILE
)


cat(
    "\nFull tree tips:",
    length(
        full_tree$tip.label
    ),
    "\n"
)


## ============================================================
## READ METADATA
## ============================================================

meta <- read_tsv(
    META_FILE,
    show_col_types = FALSE
)


if (
    nrow(meta) != 650
) {

    stop(
        "Expected 650 metadata rows; found ",
        nrow(meta)
    )
}


## ============================================================
## SELECT THE 19 TIPS
## ============================================================

mag_ids <- meta %>%
    filter(
        str_detect(
            genome,
            "^MAG_"
        )
    ) %>%
    pull(
        genome
    )


if (
    length(
        mag_ids
    ) != 16
) {

    stop(
        "Expected exactly 16 MAG IDs; found ",
        length(
            mag_ids
        )
    )
}


umezawa_ids <- c(

    "UME_AF98_Allocrenothrix_methanica",

    "UME_MI19235_Allocrenothrix_methanica"
)


outgroup_id <-
    "OUT_Methylophaga_nitratireducenticrescens"


keep_ids <- c(

    mag_ids,

    umezawa_ids,

    outgroup_id
)


if (
    length(
        keep_ids
    ) != 19
) {

    stop(
        "Expected 19 selected genomes; found ",
        length(
            keep_ids
        )
    )
}


missing_from_tree <- setdiff(
    keep_ids,
    full_tree$tip.label
)


if (
    length(
        missing_from_tree
    ) > 0
) {

    stop(
        paste0(
            "The following requested tips are missing from the tree:\n",
            paste(
                missing_from_tree,
                collapse = "\n"
            )
        )
    )
}


## ============================================================
## PRUNE THE VALIDATED TREE
##
## This preserves the topology/branch lengths from the full
## bac120 phylogeny rather than estimating a second tree.
## ============================================================

subtree <- drop.tip(

    full_tree,

    setdiff(
        full_tree$tip.label,
        keep_ids
    )
)


subtree <- ladderize(
    subtree,
    right = FALSE
)


if (
    length(
        subtree$tip.label
    ) != 19
) {

    stop(
        "Pruned tree does not contain exactly 19 tips."
    )
}


write.tree(
    subtree,
    file = OUT_TREE
)


## ============================================================
## READ MMO DATA
## ============================================================

mmo <- read_tsv(
    MMO_FILE,
    show_col_types = FALSE
)


## ============================================================
## FINAL MMO DEFINITION
##
## Particulate MMO:
##     complete canonical pMMO OR complete pXMO
##
## sMMO:
##     complete soluble system
##
## Partial systems do not count as complete detection.
## ============================================================

mmo <- mmo %>%
    mutate(

        particulate_MMO =
            canonical_pMMO == 1 |
            pXMO == 1,

        soluble_MMO =
            sMMO == 1,

        MMO_status =
            case_when(

                particulate_MMO &
                    soluble_MMO ~

                    "Both",

                particulate_MMO ~

                    "Particulate MMO",

                soluble_MMO ~

                    "sMMO",

                TRUE ~

                    "Not detected"
            ),

        MMO_status =
            factor(

                MMO_status,

                levels = c(
                    "Particulate MMO",
                    "sMMO",
                    "Both",
                    "Not detected"
                )
            )
    )


## ============================================================
## BUILD TRACK TABLE
## ============================================================

tracks <- meta %>%
    filter(
        genome %in%
            keep_ids
    ) %>%
    left_join(

        mmo %>%
            select(
                genome,
                MMO_status
            ),

        by = "genome"
    )


if (
    nrow(
        tracks
    ) != 19
) {

    stop(
        "Expected 19 rows in plotting table; found ",
        nrow(
            tracks
        )
    )
}


if (
    any(
        is.na(
            tracks$MMO_status
        )
    )
) {

    stop(
        "At least one selected genome lacks an MMO classification."
    )
}


## ============================================================
## LABELS
##
## Keep the MAG ID so individual MAGs remain distinguishable.
## Use the existing species label when available.
## ============================================================

if (
    "species_label" %in%
        names(tracks)
) {

    tracks <- tracks %>%
        mutate(

            display_label =
                case_when(

                    genome ==
                        "UME_AF98_Allocrenothrix_methanica" ~

                        "AF98 | Allocrenothrix methanica",

                    genome ==
                        "UME_MI19235_Allocrenothrix_methanica" ~

                        "MI19235 | Allocrenothrix methanica",

                    genome ==
                        outgroup_id ~

                        "Methylophaga nitratireducenticrescens | outgroup",

                    str_detect(
                        genome,
                        "^MAG_"
                    ) &
                        !is.na(
                            species_label
                        ) &
                        species_label != "" ~

                        paste0(
                            genome,
                            " | ",
                            species_label
                        ),

                    TRUE ~

                        genome
                )
        )

} else {

    tracks <- tracks %>%
        mutate(

            display_label =
                case_when(

                    genome ==
                        "UME_AF98_Allocrenothrix_methanica" ~

                        "AF98 | Allocrenothrix methanica",

                    genome ==
                        "UME_MI19235_Allocrenothrix_methanica" ~

                        "MI19235 | Allocrenothrix methanica",

                    genome ==
                        outgroup_id ~

                        "Methylophaga nitratireducenticrescens | outgroup",

                    TRUE ~

                        genome
                )
        )
}


tracks <- tracks %>%
    mutate(

        label_fontface =
            if_else(
                str_detect(
                    genome,
                    "^MAG_"
                ),
                "bold",
                "plain"
            )
    )



## ============================================================
## FINAL TIP LABELS
##
## MAG:
##
##     genus + species where available | experimental sample
##
## Published Allocrenothrix:
##
##     Allocrenothrix methanica AF98 | Umezawa et al. 2026
##
## Outgroup:
##
##     Methylophaga nitratireducenticrescens | outgroup
##
## GTDB placeholder suffixes such as Methylobacter_A and
## sp002256465 are retained because they carry taxonomic
## information.
##
## SemiBin IDs are shown ONLY when two otherwise identical
## genus-level labels occur in the same sample.
## ============================================================


## ============================================================
## Experimental sample names
## ============================================================

sample_labels <- c(

    "barcode10" =
        "BRt5_2",

    "barcode11" =
        "ULt4_2",

    "barcode12" =
        "FHt_2",

    "barcode13" =
        "UL_3 (Anoxic)",

    "barcode14" =
        "UL_SAt_2",

    "barcode15" =
        "UL_SAt_3",

    "barcode16" =
        "UL_SA_3"
)


## ============================================================
## Taxonomic labels for the 16 MAGs
##
## These correspond to the final GTDB-Tk classifications used
## in the project.
## ============================================================

mag_taxonomy <- tribble(

    ~genome,        ~taxon_label,                           ~barcode,

    "MAG_b10_SB14", "Methylomonas koyamae_A",               "barcode10",
    "MAG_b10_SB1",  "Methylumidiphilus",                    "barcode10",
    "MAG_b10_SB8",  "Methylumidiphilus",                    "barcode10",

    "MAG_b12_SB5",  "Methylomonas",                         "barcode12",
    "MAG_b12_SB0",  "Methylovulum",                         "barcode12",
    "MAG_b12_SB6",  "Methylovulum",                         "barcode12",

    "MAG_b13_SB1",  "Methylobacter_A tundripaludum_A",      "barcode13",

    "MAG_b14_SB1",  "Methylobacter_A",                      "barcode14",
    "MAG_b14_SB0",  "Methylobacter_A",                      "barcode14",

    "MAG_b15_SB0",  "Methylobacter_A tundripaludum_A",      "barcode15",
    "MAG_b15_SB11", "Methylobacter_A",                      "barcode15",
    "MAG_b15_SB3",  "Methylobacter_C sp002256465",          "barcode15",
    "MAG_b15_SB10", "Methylomonas sp018335995",             "barcode15",

    "MAG_b16_SB1",  "Methylobacter_C sp002256465",          "barcode16",
    "MAG_b16_SB12", "Methylomonas",                         "barcode16",
    "MAG_b16_SB11", "UBA4132",                              "barcode16"
)


## ============================================================
## Add MAG taxonomy and sample information
## ============================================================

tracks <- tracks %>%
    left_join(
        mag_taxonomy,
        by =
            "genome"
    ) %>%
    mutate(

        sample_label =
            if_else(

                !is.na(
                    barcode
                ),

                unname(
                    sample_labels[
                        barcode
                    ]
                ),

                NA_character_
            )
    )


## ============================================================
## Identify MAGs requiring bin-number disambiguation
##
## Example:
##
##     Methylumidiphilus [SB1] | BRt5_2
##     Methylumidiphilus [SB8] | BRt5_2
##
## rather than two identical:
##
##     Methylumidiphilus | BRt5_2
## ============================================================

duplicate_mag_taxa <- tracks %>%
    filter(
        str_detect(
            genome,
            "^MAG_"
        )
    ) %>%
    count(
        taxon_label,
        sample_label,
        name =
            "n"
    ) %>%
    filter(
        n >
            1
    )


tracks <- tracks %>%
    left_join(

        duplicate_mag_taxa %>%
            mutate(
                duplicate_taxon =
                    TRUE
            ) %>%
            select(
                taxon_label,
                sample_label,
                duplicate_taxon
            ),

        by =
            c(
                "taxon_label",
                "sample_label"
            )
    ) %>%
    mutate(

        duplicate_taxon =
            replace_na(
                duplicate_taxon,
                FALSE
            ),

        bin_label =
            str_extract(
                genome,
                "SB[0-9]+$"
            )
    )


## ============================================================
## Construct final displayed tip labels
## ============================================================

tracks <- tracks %>%
    mutate(

        display_label =
            case_when(

                genome ==
                    "UME_AF98_Allocrenothrix_methanica" ~

                    paste0(
                        "Allocrenothrix methanica AF98",
                        " | Umezawa et al. 2026"
                    ),


                genome ==
                    "UME_MI19235_Allocrenothrix_methanica" ~

                    paste0(
                        "Allocrenothrix methanica MI19235",
                        " | Umezawa et al. 2026"
                    ),


                genome ==
                    outgroup_id ~

                    paste0(
                        "Methylophaga nitratireducenticrescens",
                        " | outgroup"
                    ),


                str_detect(
                    genome,
                    "^MAG_"
                ) &
                    duplicate_taxon ~

                    paste0(
                        taxon_label,
                        " [",
                        bin_label,
                        "] | ",
                        sample_label
                    ),


                str_detect(
                    genome,
                    "^MAG_"
                ) ~

                    paste0(
                        taxon_label,
                        " | ",
                        sample_label
                    ),


                TRUE ~

                    genome
            ),


        label_fontface =
            if_else(

                str_detect(
                    genome,
                    "^MAG_"
                ),

                "bold",

                "plain"
            )
    )


## ============================================================
## QC
## ============================================================

if (
    any(
        is.na(
            tracks$display_label
        )
    )
) {

    stop(
        "One or more final tip labels is NA."
    )
}


if (
    anyDuplicated(
        tracks$display_label
    ) >
        0
) {

    stop(
        "Duplicate final display labels remain."
    )
}


cat(
    "\nFinal subtree tip labels:\n\n"
)


print(
    tracks %>%
        select(
            genome,
            display_label
        ),
    n =
        Inf
)


write_tsv(
    tracks,
    OUT_TRACKS
)


## ============================================================
## BUILD RECTANGULAR TREE
## ============================================================

p_tree <- ggtree(

    subtree,

    layout = "rectangular",

    linewidth = 0.55
)


tip_positions <- p_tree$data %>%
    filter(
        isTip
    ) %>%
    transmute(

        genome =
            label,

        x_tree =
            x,

        y =
            y
    )


plot_dat <- tip_positions %>%
    left_join(
        tracks,
        by = "genome"
    )


if (
    nrow(
        plot_dat
    ) != 19
) {

    stop(
        "Failed to map all 19 plotting rows onto the tree."
    )
}


## ============================================================
## TRACK GEOMETRY
## ============================================================

TREE_MAX <- max(
    p_tree$data$x,
    na.rm = TRUE
)


TREE_MIN <- min(
    p_tree$data$x,
    na.rm = TRUE
)


TREE_RANGE <- TREE_MAX - TREE_MIN


if (
    TREE_RANGE <= 0
) {

    TREE_RANGE <- 1
}


TRACK_WIDTH <-
    TREE_RANGE *
    0.035


TRACK_GAP <-
    TREE_RANGE *
    0.025


MMO_X <-
    TREE_MAX +
    TREE_RANGE *
    0.07


COMP_X <-
    MMO_X +
    TRACK_WIDTH +
    TRACK_GAP


CONT_X <-
    COMP_X +
    TRACK_WIDTH +
    TRACK_GAP


LABEL_X <-
    CONT_X +
    TRACK_WIDTH +
    TREE_RANGE *
    0.05


PLOT_XMAX <-
    LABEL_X +
    TREE_RANGE *
    0.90


## Track-heading position

TRACK_TITLE_Y <-
    max(
        plot_dat$y
    ) +
    0.75


## ============================================================
## PLOT
## ============================================================

p <- p_tree +


    ## --------------------------------------------------------
    ## MMO TRACK
    ## --------------------------------------------------------

    geom_tile(

        data =
            plot_dat,

        aes(
            x =
                MMO_X,
            y =
                y,
            fill =
                MMO_status
        ),

        inherit.aes =
            FALSE,

        width =
            TRACK_WIDTH,

        height =
            0.98,

        color =
            TRACK_OUTLINE,

        linewidth =
            0.20
    ) +


    ## FORCE_SMMO_LEGEND_KEY
    ##
    ## No sMMO-only genome occurs in this subtree.
    ## This invisible off-plot tile forces ggplot to render the
    ## sMMO legend key using the same color as the big tree.

    geom_tile(

        data =
            tibble(

                x =
                    MMO_X,

                y =
                    -1,

                MMO_status =
                    factor(

                        "sMMO",

                        levels =
                            c(
                                "Particulate MMO",
                                "sMMO",
                                "Both",
                                "Not detected"
                            )
                    )
            ),

        aes(
            x =
                x,
            y =
                y,
            fill =
                MMO_status
        ),

        inherit.aes =
            FALSE,

        width =
            TRACK_WIDTH,

        height =
            0.98,

        color =
            TRACK_OUTLINE,

        linewidth =
            0.20,

        show.legend =
            TRUE
    ) +


    scale_fill_manual(

        values =
            MMO_COLORS,

        limits =
            c(
                "Particulate MMO",
                "sMMO",
                "Both",
                "Not detected"
            ),

        breaks =
            c(
                "Particulate MMO",
                "sMMO",
                "Both",
                "Not detected"
            ),

        labels =
            c(
                "Particulate MMO (pmo/pxm)",
                "sMMO",
                "Both",
                "Not detected"
            ),

        drop =
            FALSE,

        name =
            "Methane monooxygenase"
    ) +


    ggnewscale::new_scale_fill() +


    ## --------------------------------------------------------
    ## COMPLETENESS
    ## --------------------------------------------------------

    geom_tile(

        data =
            plot_dat,

        aes(
            x =
                COMP_X,
            y =
                y,
            fill =
                completeness
        ),

        inherit.aes =
            FALSE,

        width =
            TRACK_WIDTH,

        height =
            0.98,

        color =
            TRACK_OUTLINE,

        linewidth =
            0.20
    ) +


    scale_fill_gradient(

        low =
            COMP_LOW,

        high =
            COMP_HIGH,

        limits =
            c(
                0,
                100
            ),

        na.value =
            "grey95",

        name =
            "Completeness (%)"
    ) +


    ggnewscale::new_scale_fill() +


    ## --------------------------------------------------------
    ## CONTAMINATION
    ## --------------------------------------------------------

    geom_tile(

        data =
            plot_dat,

        aes(
            x =
                CONT_X,
            y =
                y,
            fill =
                contamination
        ),

        inherit.aes =
            FALSE,

        width =
            TRACK_WIDTH,

        height =
            0.98,

        color =
            TRACK_OUTLINE,

        linewidth =
            0.20
    ) +


    scale_fill_gradient(

        low =
            CONT_LOW,

        high =
            CONT_HIGH,

        na.value =
            "grey95",

        name =
            "Contamination (%)"
    ) +


    ## --------------------------------------------------------
    ## TIP LABELS
    ## --------------------------------------------------------

    geom_text(

        data =
            plot_dat,

        aes(
            x =
                LABEL_X,
            y =
                y,
            label =
                display_label,
            fontface =
                label_fontface
        ),

        inherit.aes =
            FALSE,

        hjust =
            0,

        size =
            3.3
    ) +


    ## --------------------------------------------------------
    ## TRACK HEADINGS
    ##
    ## Tilted and placed directly above their respective tracks.
    ## --------------------------------------------------------

    annotate(

        "text",

        x =
            MMO_X,

        y =
            TRACK_TITLE_Y,

        label =
            "MMO",

        angle =
            45,

        hjust =
            0,

        vjust =
            0.5,

        size =
            3.2,

        fontface =
            "bold"
    ) +


    annotate(

        "text",

        x =
            COMP_X,

        y =
            TRACK_TITLE_Y,

        label =
            "Completeness",

        angle =
            45,

        hjust =
            0,

        vjust =
            0.5,

        size =
            3.2,

        fontface =
            "bold"
    ) +


    annotate(

        "text",

        x =
            CONT_X,

        y =
            TRACK_TITLE_Y,

        label =
            "Contamination",

        angle =
            45,

        hjust =
            0,

        vjust =
            0.5,

        size =
            3.2,

        fontface =
            "bold"
    ) +


    ## --------------------------------------------------------
    ## SCALE BAR
    ## --------------------------------------------------------

    geom_treescale(

        x =
            TREE_MIN,

        y =
            0,

        fontsize =
            3
    ) +


    coord_cartesian(

        xlim =
            c(
                TREE_MIN,
                PLOT_XMAX
            ),

        ylim =
            c(
                0,
                TRACK_TITLE_Y +
                    1.0
            ),

        clip =
            "off"
    ) +


    theme_tree2() +


    theme(

        plot.margin =
            margin(
                15,
                35,
                15,
                15
            ),

        legend.position =
            "right",

        legend.title =
            element_text(
                size =
                    10,
                face =
                    "bold"
            ),

        legend.text =
            element_text(
                size =
                    9
            ),

        axis.text.x =
            element_text(
                size =
                    8
            ),

        axis.title.x =
            element_blank()
    )


## ============================================================
## SAVE
## ============================================================

ggsave(

    filename =
        OUT_PDF,

    plot =
        p,

    width =
        13,

    height =
        8.5,

    units =
        "in",

    device =
        cairo_pdf
)


ggsave(

    filename =
        OUT_PNG,

    plot =
        p,

    width =
        13,

    height =
        8.5,

    units =
        "in",

    dpi =
        400
)


## ============================================================
## REPORT
## ============================================================

cat(
    "\n============================================================\n"
)

cat(
    "MAG SUBTREE COMPLETE\n"
)

cat(
    "============================================================\n\n"
)


cat(
    "Tips:",
    length(
        subtree$tip.label
    ),
    "\n\n"
)


cat(
    "MMO status:\n"
)


print(
    table(
        tracks$MMO_status,
        useNA = "ifany"
    )
)


cat(
    "\nCompleteness range:\n"
)


print(
    range(
        tracks$completeness,
        na.rm = TRUE
    )
)


cat(
    "\nContamination range:\n"
)


print(
    range(
        tracks$contamination,
        na.rm = TRUE
    )
)


cat(
    "\nOutputs:\n"
)

cat(
    "  ",
    OUT_PDF,
    "\n",
    sep = ""
)

cat(
    "  ",
    OUT_PNG,
    "\n",
    sep = ""
)

cat(
    "  ",
    OUT_TREE,
    "\n",
    sep = ""
)

cat(
    "  ",
    OUT_TRACKS,
    "\n",
    sep = ""
)
