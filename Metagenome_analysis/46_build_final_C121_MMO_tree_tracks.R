#!/usr/bin/env Rscript


## ================================================================== ##
## STAGE 93H
##
## FINAL EXTRA TRACKS FOR THE 650-TIP SPECIES TREE
##
## Track 1:
##
##   Cluster_00121
##       Present
##       absence implicit
##
## Track 2:
##
##   Methane monooxygenase repertoire
##
##       Particulate MMO
##           complete canonical pMMO OR complete pXMO
##
##       sMMO
##           complete sMMO, no complete particulate system
##
##       Both
##           complete particulate system AND complete sMMO
##
##       Not detected
##           neither complete particulate nor complete soluble system
##
##
## IMPORTANT:
##
## Partial pMMO / pXMO / sMMO component sets are conservatively
## classified as "Not detected" because this track represents
## COMPLETE detected systems.
##
## The detailed partial-component information remains available in
## the Stage 93F QC tables.
## ================================================================== ##


suppressPackageStartupMessages({

    library(dplyr)
    library(readr)

})


## ================================================================== ##
## Paths
## ================================================================== ##

ROOT <- file.path(
    Sys.getenv("HOME"),
    "methanotrophs",
    "methanotroph_project",
    "jeppe"
)


PROJECT <- file.path(
    ROOT,
    "metagenome"
)


PHYLO <- file.path(
    PROJECT,
    "comparative_analysis",
    "final_species_phylogeny"
)


EXTRA <- file.path(
    PHYLO,
    "08_visualization",
    "extra_tracks"
)


META_FILE <- file.path(
    PHYLO,
    "08_visualization",
    "final_650_tip_metadata.tsv"
)


MMO_QC_FILE <- file.path(
    EXTRA,
    "all_650_pmo_pxm_smmo_QC.tsv"
)


C121_FOCAL <- file.path(
    PROJECT,
    "comparative_analysis",
    "Cluster_00121",
    "Cluster_00121_focal_regions.tsv"
)


OUT_FILE <- file.path(
    EXTRA,
    "final_C121_MMO_tree_tracks.tsv"
)


## ================================================================== ##
## MAG tree-ID crosswalk
## ================================================================== ##

mag_crosswalk <- tribble(

    ~original_genome,                              ~genome,

    "semibin2__flye__barcode10_seqs__SemiBin_14", "MAG_b10_SB14",
    "semibin2__flye__barcode10_seqs__SemiBin_1",  "MAG_b10_SB1",
    "semibin2__flye__barcode10_seqs__SemiBin_8",  "MAG_b10_SB8",

    "semibin2__flye__barcode12_seqs__SemiBin_5",  "MAG_b12_SB5",
    "semibin2__flye__barcode12_seqs__SemiBin_0",  "MAG_b12_SB0",
    "semibin2__flye__barcode12_seqs__SemiBin_6",  "MAG_b12_SB6",

    "semibin2__flye__barcode13_seqs__SemiBin_1",  "MAG_b13_SB1",

    "semibin2__flye__barcode14_seqs__SemiBin_1",  "MAG_b14_SB1",
    "semibin2__flye__barcode14_seqs__SemiBin_0",  "MAG_b14_SB0",

    "semibin2__flye__barcode15_seqs__SemiBin_0",  "MAG_b15_SB0",
    "semibin2__flye__barcode15_seqs__SemiBin_11", "MAG_b15_SB11",
    "semibin2__flye__barcode15_seqs__SemiBin_3",  "MAG_b15_SB3",
    "semibin2__flye__barcode15_seqs__SemiBin_10", "MAG_b15_SB10",

    "semibin2__flye__barcode16_seqs__SemiBin_1",  "MAG_b16_SB1",
    "semibin2__flye__barcode16_seqs__SemiBin_12", "MAG_b16_SB12",
    "semibin2__flye__barcode16_seqs__SemiBin_11", "MAG_b16_SB11"
)


## ================================================================== ##
## Input QC
## ================================================================== ##

for (
    path in
    c(
        META_FILE,
        MMO_QC_FILE,
        C121_FOCAL
    )
) {

    if (
        !file.exists(
            path
        )
    ) {

        stop(
            "Missing required file: ",
            path
        )
    }
}


## ================================================================== ##
## Read final 650-tip metadata
## ================================================================== ##

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


if (
    n_distinct(
        meta$genome
    ) != 650
) {

    stop(
        "Metadata does not contain 650 unique genome IDs."
    )
}


## ================================================================== ##
## MMO calls from Stage 93F
## ================================================================== ##

mmo <- read_tsv(
    MMO_QC_FILE,
    show_col_types = FALSE
)


if (
    nrow(mmo) != 650
) {

    stop(
        "Expected 650 Stage-93F MMO rows; found ",
        nrow(mmo)
    )
}


if (
    n_distinct(
        mmo$genome
    ) != 650
) {

    stop(
        "Stage-93F MMO table does not contain 650 unique genomes."
    )
}


## ================================================================== ##
## Define FINAL methane-monooxygenase track
##
## particulate =
##
##     canonical pMMO
##          OR
##     pXMO
##
## Complete systems only.
## ================================================================== ##

mmo_final <- mmo %>%
    mutate(

        particulate_MMO =
            canonical_pMMO ==
                1 |
            pXMO ==
                1,

        soluble_MMO =
            sMMO ==
                1,

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
            )
    ) %>%
    select(

        genome,

        MMO_status,

        canonical_pMMO,

        pXMO,

        sMMO,

        partial_pMMO,

        partial_pXMO,

        partial_sMMO
    )


## ================================================================== ##
## Cluster_00121
## ================================================================== ##

c121 <- read_tsv(
    C121_FOCAL,
    show_col_types = FALSE
)


if (
    !"genome" %in%
        names(c121)
) {

    stop(
        "Cluster_00121 focal table lacks a genome column."
    )
}


c121_positive <- c121 %>%
    transmute(

        genome =
            recode(
                as.character(
                    genome
                ),

                !!!setNames(
                    mag_crosswalk$genome,
                    mag_crosswalk$original_genome
                )
            )
    ) %>%
    distinct(
        genome
    )


cat(
    "\nCluster_00121 unique focal genomes: ",
    nrow(
        c121_positive
    ),
    "\n",
    sep = ""
)


if (
    nrow(
        c121_positive
    ) !=
        23
) {

    warning(
        "Expected 23 Cluster_00121-positive genomes; found ",
        nrow(
            c121_positive
        ),
        ". Check before final interpretation."
    )
}


not_in_tree <- setdiff(
    c121_positive$genome,
    meta$genome
)


if (
    length(
        not_in_tree
    ) >
        0
) {

    cat(
        "\nCluster_00121 IDs absent from tree:\n"
    )


    print(
        not_in_tree
    )


    stop(
        paste0(
            "One or more Cluster_00121 focal genomes could not be mapped ",
            "to the final tree."
        )
    )
}


## ================================================================== ##
## Build final 650-tip track table
## ================================================================== ##

final_tracks <- meta %>%
    select(
        genome
    ) %>%
    left_join(
        mmo_final,
        by = "genome"
    ) %>%
    mutate(

        Cluster00121_plot =
            if_else(

                genome %in%
                    c121_positive$genome,

                "Present",

                NA_character_
            ),

        MMO_status =
            factor(

                MMO_status,

                levels =
                    c(
                        "Particulate MMO",
                        "sMMO",
                        "Both",
                        "Not detected"
                    )
            )
    )


## ================================================================== ##
## Final QC
## ================================================================== ##

if (
    any(
        is.na(
            final_tracks$MMO_status
        )
    )
) {

    stop(
        "One or more tree tips lacks a final MMO classification."
    )
}


if (
    sum(
        final_tracks$Cluster00121_plot ==
            "Present",
        na.rm = TRUE
    ) !=
        23
) {

    stop(
        paste0(
            "Final Cluster_00121 track does not contain exactly 23 ",
            "positive genomes."
        )
    )
}


write_tsv(
    final_tracks,
    OUT_FILE
)


## ================================================================== ##
## Report
## ================================================================== ##

cat(
    "\n============================================================\n"
)

cat(
    "STAGE 93H COMPLETE\n"
)

cat(
    "============================================================\n\n"
)


cat(
    "Final MMO repertoire:\n"
)


print(
    table(
        final_tracks$MMO_status,
        useNA = "ifany"
    )
)


cat(
    "\nCluster_00121:\n"
)


print(
    table(
        final_tracks$Cluster00121_plot,
        useNA = "ifany"
    )
)


cat(
    "\nDetailed particulate-system composition:\n"
)


cat(
    "Canonical pMMO only: ",
    sum(
        final_tracks$canonical_pMMO == 1 &
        final_tracks$pXMO == 0,
        na.rm = TRUE
    ),
    "\n",
    sep = ""
)


cat(
    "pXMO only: ",
    sum(
        final_tracks$canonical_pMMO == 0 &
        final_tracks$pXMO == 1,
        na.rm = TRUE
    ),
    "\n",
    sep = ""
)


cat(
    "Canonical pMMO + pXMO: ",
    sum(
        final_tracks$canonical_pMMO == 1 &
        final_tracks$pXMO == 1,
        na.rm = TRUE
    ),
    "\n",
    sep = ""
)


cat(
    "\nOutput:\n  ",
    OUT_FILE,
    "\n",
    sep = ""
)
