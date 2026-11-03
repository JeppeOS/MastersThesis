#!/usr/bin/env Rscript


## ============================================================
## APPENDIX TABLE - FOCAL CYTOCHROME CLUSTERS
##
## Clusters:
##   Cluster_00024  No MCL module
##   Cluster_00035  Module_20
##   Cluster_00042  No MCL module
##   Cluster_00050  Module_10
##   Cluster_00063  Module_06
##   Cluster_00121  No MCL module
##   Cluster_00206  Module_35
##   Cluster_00288  No MCL module
##
## Primary input:
##   genome_analysis_workflow/06_clustering/cluster_membership.tsv
##
## Outputs:
##   TSV
##   landscape A4 DOCX
## ============================================================


suppressPackageStartupMessages({

    library(tidyverse)
    library(officer)
    library(flextable)

})


## ============================================================
## PATHS
## ============================================================

ROOT <- file.path(
    Sys.getenv("HOME"),
    "methanotrophs",
    "methanotroph_project",
    "jeppe"
)


WORKFLOW <- file.path(
    ROOT,
    "genome_analysis_workflow"
)


METAGENOME <- file.path(
    ROOT,
    "metagenome"
)


INPUT_FILE <- file.path(
    WORKFLOW,
    "06_clustering",
    "cluster_membership.tsv"
)


MAG_FUNCTIONAL_MASTER <- file.path(
    METAGENOME,
    "functional_analysis",
    "integrated_annotations",
    "MASTER_PROTEIN_ANNOTATIONS.tsv"
)


OUTDIR <- file.path(
    ROOT,
    "metagenome",
    "appendix_tables"
)


dir.create(
    OUTDIR,
    recursive = TRUE,
    showWarnings = FALSE
)


OUT_TSV <- file.path(
    OUTDIR,
    "Appendix_Focal_Cluster_Characteristics.tsv"
)


OUT_DOCX <- file.path(
    OUTDIR,
    "Appendix_Focal_Cluster_Characteristics.docx"
)


## ============================================================
## FOCAL CLUSTERS
## ============================================================

cluster_info <- tribble(

    ~cluster,
    ~module,
    ~functional_interpretation,

    "Cluster_00024",
    "No MCL module",
    "Large exported Cyc2 family",

    "Cluster_00035",
    "Module_20",
    "MtoA/MtrA-like multiheme cytochromes",

    "Cluster_00042",
    "No MCL module",
    "Multiheme cytochrome family",

    "Cluster_00050",
    "Module_10",
    "Cyc2 family",

    "Cluster_00063",
    "Module_06",
    "MCA0421/c553O-like multiheme cytochromes",

    "Cluster_00121",
    "No MCL module",
    "Conserved exported 9-heme cytochrome family",

    "Cluster_00206",
    "Module_35",
    "Cyc2 family",

    "Cluster_00288",
    "No MCL module",
    "Multiheme cytochrome family"
)


FOCAL_CLUSTERS <- cluster_info$cluster


## ============================================================
## Additional MAG homologues assigned to original GlobDB
## protein clusters
##
## Cluster_00050 contains two homologues recovered from the
## SemiBin2 Methylumidiphilus MAGs.
## ============================================================

mag_cluster_membership <- tribble(

    ~cluster,
    ~genome,
    ~protein_id,

    "Cluster_00050",
    "semibin2__flye__barcode10_seqs__SemiBin_1",
    "contig_25_317",

    "Cluster_00050",
    "semibin2__flye__barcode10_seqs__SemiBin_8",
    "contig_26_118"
)


## ============================================================
## READ DATA
## ============================================================

if (
    !file.exists(
        INPUT_FILE
    )
) {

    stop(
        "Input file not found: ",
        INPUT_FILE
    )
}


dat <- read_tsv(
    INPUT_FILE,
    show_col_types = FALSE
)


cat(
    "\nOriginal GlobDB cluster-membership proteins: ",
    nrow(dat),
    "\n",
    sep = ""
)


if (
    !file.exists(
        MAG_FUNCTIONAL_MASTER
    )
) {

    stop(
        "MAG functional annotation table not found: ",
        MAG_FUNCTIONAL_MASTER
    )
}


mag_annotations <- read_tsv(
    MAG_FUNCTIONAL_MASTER,
    show_col_types = FALSE
)


cat(
    "MAG functional annotation rows: ",
    nrow(mag_annotations),
    "\n",
    sep = ""
)


## ============================================================
## REQUIRED COLUMNS
## ============================================================

required <- c(

    "cluster",
    "genome",
    "protein_id",
    "length",

    "fegenie_positive",
    "fegenie_HMMs",

    "findmehemes_positive",
    "number_of_hemes",

    "signalp_prediction",

    "deeptmhmm_class",
    "deeptmhmm_n_tm_helices",

    "localization_class"
)


missing <- setdiff(
    required,
    names(dat)
)


if (
    length(
        missing
    ) > 0
) {

    stop(
        paste0(
            "Missing required columns:\n",
            paste(
                missing,
                collapse = "\n"
            )
        )
    )
}



## ============================================================
## Validate MAG functional annotation columns
##
## "cluster" is supplied by mag_cluster_membership rather than
## by MASTER_PROTEIN_ANNOTATIONS.tsv.
## ============================================================

mag_required <- setdiff(
    required,
    "cluster"
)


missing_mag_columns <- setdiff(
    mag_required,
    names(
        mag_annotations
    )
)


if (
    length(
        missing_mag_columns
    ) > 0
) {

    stop(
        paste0(
            "MAG functional annotation table is missing columns:\n",
            paste(
                missing_mag_columns,
                collapse = "\n"
            )
        )
    )
}


## ============================================================
## KEEP ONLY FOCAL CLUSTERS
## ============================================================

## ------------------------------------------------------------
## Original GlobDB proteins
## ------------------------------------------------------------

focal_globdb <- dat %>%
    filter(
        cluster %in%
            FOCAL_CLUSTERS
    ) %>%
    select(
        all_of(
            required
        )
    ) %>%
    distinct(
        cluster,
        genome,
        protein_id,
        .keep_all = TRUE
    )


## ------------------------------------------------------------
## Retrieve complete annotation records for the two MAG
## Cluster_00050 homologues.
## ------------------------------------------------------------

mag_focal <- mag_cluster_membership %>%

    inner_join(
        mag_annotations,
        by =
            c(
                "genome",
                "protein_id"
            )
    ) %>%

    select(
        all_of(
            required
        )
    ) %>%

    distinct(
        cluster,
        genome,
        protein_id,
        .keep_all = TRUE
    )


## ------------------------------------------------------------
## QC: both intended MAG proteins must be recovered.
## ------------------------------------------------------------

if (
    nrow(
        mag_focal
    ) !=
        2
) {

    stop(
        paste0(
            "Expected exactly 2 MAG proteins for Cluster_00050, ",
            "but recovered ",
            nrow(
                mag_focal
            ),
            "."
        )
    )
}


expected_mag_pairs <- mag_cluster_membership %>%
    transmute(
        key =
            paste(
                genome,
                protein_id,
                sep = " | "
            )
    ) %>%
    pull(
        key
    )


observed_mag_pairs <- mag_focal %>%
    transmute(
        key =
            paste(
                genome,
                protein_id,
                sep = " | "
            )
    ) %>%
    pull(
        key
    )


if (
    !setequal(
        expected_mag_pairs,
        observed_mag_pairs
    )
) {

    stop(
        "Recovered MAG Cluster_00050 proteins do not match expected proteins."
    )
}


cat(
    "\nMAG proteins added to Cluster_00050:\n"
)


print(
    mag_focal %>%
        select(
            cluster,
            genome,
            protein_id,
            length,
            fegenie_HMMs,
            number_of_hemes,
            signalp_prediction,
            deeptmhmm_class,
            deeptmhmm_n_tm_helices,
            localization_class
        )
)


## ------------------------------------------------------------
## Combined dataset used for every cluster-level statistic.
##
## Therefore the MAG proteins contribute to:
##
##   - protein count
##   - genome count
##   - copy-number statistics
##   - FeGenie composition
##   - CXXCH mean/range
##   - protein-length mean/range
##   - DeepTMHMM composition
##   - TM-helix mean/range
##   - SignalP composition
##   - localization composition
## ------------------------------------------------------------

focal <- bind_rows(
    focal_globdb,
    mag_focal
) %>%
    distinct(
        cluster,
        genome,
        protein_id,
        .keep_all = TRUE
    )


cat(
    "\nCombined focal proteins by cluster:\n"
)


print(
    focal %>%
        count(
            cluster
        )
)


## ============================================================
## HELPERS
## ============================================================


## ------------------------------------------------------------
## Detect TRUE-like values robustly.
## ------------------------------------------------------------

is_true <- function(x) {

    str_to_lower(
        str_trim(
            as.character(x)
        )
    ) %in%
        c(
            "true",
            "t",
            "1",
            "yes",
            "y"
        )
}


## ------------------------------------------------------------
## Clean annotation fields.
## ------------------------------------------------------------

clean_annotation <- function(x) {

    x <- as.character(x)

    x[
        is.na(x) |
        x == "" |
        x == "NA" |
        x == "None"
    ] <- "Not predicted"

    x
}


## ------------------------------------------------------------
## Compact composition:
##
## MtoA: 61; MtrA: 3; FeGenie-negative: 10
## ------------------------------------------------------------

composition_string <- function(x) {

    x <- as.character(x)

    counts <- sort(
        table(x),
        decreasing = TRUE
    )


    paste0(
        names(counts),
        ": ",
        as.integer(counts),
        collapse = "; "
    )
}


## ------------------------------------------------------------
## Numeric mean + range.
##
## Example:
##
## 8.7 (7-9)
## ------------------------------------------------------------

mean_range <- function(
    x,
    digits = 1
) {

    x <- as.numeric(x)

    x <- x[
        !is.na(x)
    ]


    if (
        length(x) == 0
    ) {

        return("—")
    }


    mean_value <- round(
        mean(x),
        digits
    )


    min_value <- round(
        min(x),
        digits
    )


    max_value <- round(
        max(x),
        digits
    )


    paste0(
        mean_value,
        " (",
        min_value,
        "–",
        max_value,
        ")"
    )
}


## ------------------------------------------------------------
## Integer-like mean + integer range.
## ------------------------------------------------------------

mean_integer_range <- function(x) {

    x <- as.numeric(x)

    x <- x[
        !is.na(x)
    ]


    if (
        length(x) == 0
    ) {

        return("—")
    }


    paste0(

        round(
            mean(x),
            1
        ),

        " (",

        min(x),

        "–",

        max(x),

        ")"
    )
}


## ============================================================
## STANDARDIZE ANNOTATION FIELDS
## ============================================================

focal <- focal %>%
    mutate(

        length =
            as.numeric(
                length
            ),


        ## ----------------------------------------------------
        ## FeGenie annotation.
        ##
        ## Negative proteins are explicitly represented rather
        ## than disappearing as NA.
        ## ----------------------------------------------------

        FeGenie_annotation =
            case_when(

                is_true(
                    fegenie_positive
                ) &
                    !is.na(
                        fegenie_HMMs
                    ) &
                    fegenie_HMMs != "" ~

                    as.character(
                        fegenie_HMMs
                    ),

                TRUE ~

                    "FeGenie-negative"
            ),


        ## ----------------------------------------------------
        ## CXXCH / heme-motif count.
        ##
        ## number_of_hemes derives from FindMeHemes.
        ##
        ## Proteins explicitly negative by FindMeHemes receive
        ## zero motifs rather than NA.
        ## ----------------------------------------------------

        CXXCH =
            case_when(

                !is.na(
                    suppressWarnings(
                        as.numeric(
                            number_of_hemes
                        )
                    )
                ) ~

                    suppressWarnings(
                        as.numeric(
                            number_of_hemes
                        )
                    ),

                !is_true(
                    findmehemes_positive
                ) ~

                    0,

                TRUE ~

                    NA_real_
            ),


        DeepTMHMM =
            clean_annotation(
                deeptmhmm_class
            ),


        TM_helices =
            suppressWarnings(
                as.numeric(
                    deeptmhmm_n_tm_helices
                )
            ),


        SignalP =
            clean_annotation(
                signalp_prediction
            ),


        Localization =
            clean_annotation(
                localization_class
            )
    )


## ============================================================
## COPY-NUMBER INFORMATION
## ============================================================

copy_summary <- focal %>%
    count(
        cluster,
        genome,
        name =
            "copies"
    ) %>%
    group_by(
        cluster
    ) %>%
    summarise(

        max_copies_per_genome =
            max(
                copies
            ),

        genomes_multicopy =
            sum(
                copies >
                    1
            ),

        .groups =
            "drop"
    )


## ============================================================
## CLUSTER SUMMARY
## ============================================================

summary <- focal %>%
    group_by(
        cluster
    ) %>%
    summarise(

        n_proteins =
            n(),

        n_genomes =
            n_distinct(
                genome
            ),

        FeGenie_composition =
            composition_string(
                FeGenie_annotation
            ),

        CXXCH_mean_range =
            mean_range(
                CXXCH,
                digits =
                    1
            ),

        length_mean_range =
            mean_integer_range(
                length
            ),

        DeepTMHMM_composition =
            composition_string(
                DeepTMHMM
            ),

        TM_helices_mean_range =
            mean_range(
                TM_helices,
                digits =
                    1
            ),

        SignalP_composition =
            composition_string(
                SignalP
            ),

        Localization_composition =
            composition_string(
                Localization
            ),

        .groups =
            "drop"
    ) %>%
    left_join(
        copy_summary,
        by =
            "cluster"
    ) %>%
    left_join(
        cluster_info,
        by =
            "cluster"
    )


## ============================================================
## ORDER
## ============================================================

summary <- summary %>%
    mutate(

        cluster =
            factor(
                cluster,
                levels =
                    FOCAL_CLUSTERS
            )
    ) %>%
    arrange(
        cluster
    )


## ============================================================
## FINAL DISPLAY TABLE
## ============================================================

appendix <- summary %>%
    transmute(

        `Cluster` =
            as.character(
                cluster
            ),

        `Module` =
            module,

        `Functional interpretation` =
            functional_interpretation,

        `Proteins` =
            n_proteins,

        `Genomes` =
            n_genomes,

        `Max copies/genome` =
            max_copies_per_genome,

        `Multicopy genomes` =
            genomes_multicopy,

        `FeGenie annotation composition` =
            FeGenie_composition,

        `CXXCH motifs: mean (range)` =
            CXXCH_mean_range,

        `Protein length, aa: mean (range)` =
            length_mean_range,

        `DeepTMHMM class composition` =
            DeepTMHMM_composition,

        `TM helices: mean (range)` =
            TM_helices_mean_range,

        `SignalP composition` =
            SignalP_composition,

        `Localization composition` =
            Localization_composition
    )


## ============================================================
## QC
## ============================================================

if (
    nrow(
        appendix
    ) !=
        8
) {

    stop(
        "Expected exactly eight focal cluster rows."
    )
}


cat(
    "\n============================================================\n"
)

cat(
    "FOCAL CLUSTER SUMMARY\n"
)

cat(
    "============================================================\n\n"
)


print(
    appendix,
    width =
        Inf
)


## Known workflow sizes - warnings only.
##
## These are useful guards against accidentally using a
## different cluster-membership file.

expected_sizes <- c(

    "Cluster_00024" =
        116,

    "Cluster_00035" =
        74,

    "Cluster_00042" =
        63,

    "Cluster_00050" =
        55,

    "Cluster_00063" =
        46,

    "Cluster_00121" =
        19,

    "Cluster_00206" =
        8,

    "Cluster_00288" =
        5
)


observed_sizes <- setNames(
    appendix$Proteins,
    appendix$Cluster
)


for (
    cl in
    names(
        expected_sizes
    )
) {

    if (
        observed_sizes[
            cl
        ] !=
            expected_sizes[
                cl
            ]
    ) {

        warning(
            cl,
            ": expected ",
            expected_sizes[
                cl
            ],
            " proteins from the validated workflow, but found ",
            observed_sizes[
                cl
            ],
            "."
        )
    }
}


## ============================================================
## WRITE TSV
## ============================================================

write_tsv(
    appendix,
    OUT_TSV
)


## ============================================================
## WORD TABLE
## ============================================================

ft <- flextable(
    appendix
)


ft <- set_header_labels(

    ft,

    `Max copies/genome` =
        "Max copies/\ngenome",

    `Multicopy genomes` =
        "Multicopy\ngenomes",

    `FeGenie annotation composition` =
        "FeGenie annotation\ncomposition",

    `CXXCH motifs: mean (range)` =
        "CXXCH motifs\nmean (range)",

    `Protein length, aa: mean (range)` =
        "Protein length (aa)\nmean (range)",

    `DeepTMHMM class composition` =
        "DeepTMHMM class\ncomposition",

    `TM helices: mean (range)` =
        "TM helices\nmean (range)",

    `SignalP composition` =
        "SignalP\ncomposition",

    `Localization composition` =
        "Localization\ncomposition"
)


## ------------------------------------------------------------
## Typography
## ------------------------------------------------------------

ft <- font(
    ft,
    fontname =
        "Arial",
    part =
        "all"
)


ft <- fontsize(
    ft,
    size =
        7,
    part =
        "body"
)


ft <- fontsize(
    ft,
    size =
        7.5,
    part =
        "header"
)


ft <- bold(
    ft,
    part =
        "header"
)


## ------------------------------------------------------------
## Alignment
## ------------------------------------------------------------

ft <- align(
    ft,
    align =
        "center",
    part =
        "header"
)


ft <- align(

    ft,

    j =
        c(
            "Cluster",
            "Module",
            "Functional interpretation",
            "FeGenie annotation composition",
            "DeepTMHMM class composition",
            "SignalP composition",
            "Localization composition"
        ),

    align =
        "left",

    part =
        "body"
)


ft <- align(

    ft,

    j =
        c(
            "Proteins",
            "Genomes",
            "Max copies/genome",
            "Multicopy genomes",
            "CXXCH motifs: mean (range)",
            "Protein length, aa: mean (range)",
            "TM helices: mean (range)"
        ),

    align =
        "center",

    part =
        "body"
)


## ------------------------------------------------------------
## Header / row shading
## ------------------------------------------------------------

ft <- bg(
    ft,
    bg =
        "#D9D9D9",
    part =
        "header"
)


ft <- bg(
    ft,
    i =
        c(
            2,
            4,
            6,
            8
        ),
    bg =
        "#F5F5F5",
    part =
        "body"
)


## ------------------------------------------------------------
## Borders
## ------------------------------------------------------------

thin_border <- fp_border(
    color =
        "#BFBFBF",
    width =
        0.4
)


ft <- border(
    ft,
    border =
        thin_border,
    part =
        "all"
)


## ------------------------------------------------------------
## Padding
## ------------------------------------------------------------

ft <- padding(

    ft,

    padding.top =
        2,

    padding.bottom =
        2,

    padding.left =
        2,

    padding.right =
        2,

    part =
        "all"
)


## ------------------------------------------------------------
## Column widths
## ------------------------------------------------------------

ft <- width(
    ft,
    j =
        "Cluster",
    width =
        0.75
)


ft <- width(
    ft,
    j =
        "Module",
    width =
        0.60
)


ft <- width(
    ft,
    j =
        "Functional interpretation",
    width =
        1.25
)


ft <- width(
    ft,
    j =
        c(
            "Proteins",
            "Genomes"
        ),
    width =
        0.45
)


ft <- width(
    ft,
    j =
        c(
            "Max copies/genome",
            "Multicopy genomes"
        ),
    width =
        0.55
)


ft <- width(
    ft,
    j =
        "FeGenie annotation composition",
    width =
        1.45
)


ft <- width(
    ft,
    j =
        "CXXCH motifs: mean (range)",
    width =
        0.85
)


ft <- width(
    ft,
    j =
        "Protein length, aa: mean (range)",
    width =
        0.95
)


ft <- width(
    ft,
    j =
        "DeepTMHMM class composition",
    width =
        1.15
)


ft <- width(
    ft,
    j =
        "TM helices: mean (range)",
    width =
        0.75
)


ft <- width(
    ft,
    j =
        "SignalP composition",
    width =
        1.10
)


ft <- width(
    ft,
    j =
        "Localization composition",
    width =
        1.15
)


ft <- set_table_properties(

    ft,

    layout =
        "fixed",

    opts_word =
        list(
            split =
                FALSE,
            repeat_headers =
                TRUE
        )
)


## ============================================================
## WORD DOCUMENT
## ============================================================

doc <- read_docx()


landscape_section <- prop_section(

    page_size =
        page_size(
            orient =
                "landscape",
            width =
                8.27,
            height =
                11.69
        ),

    page_margins =
        page_mar(
            top =
                0.35,
            bottom =
                0.35,
            left =
                0.30,
            right =
                0.30
        )
)


doc <- body_set_default_section(
    doc,
    landscape_section
)


doc <- body_add_par(

    doc,

    paste0(
        "Appendix Table. Characteristics of the eight focal ",
        "cytochrome protein clusters."
    ),

    style =
        "Normal"
)


doc <- body_add_par(
    doc,
    ""
)


doc <- body_add_flextable(
    doc,
    ft
)


print(
    doc,
    target =
        OUT_DOCX
)


cat(
    "\nOutputs:\n"
)

cat(
    "  ",
    OUT_TSV,
    "\n",
    sep = ""
)

cat(
    "  ",
    OUT_DOCX,
    "\n",
    sep = ""
)
