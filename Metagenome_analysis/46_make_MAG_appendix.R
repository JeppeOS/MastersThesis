#!/usr/bin/env Rscript


## ============================================================
## APPENDIX TABLE - ALL 187 SEMIBIN2 MAGs
##
## Output:
##   1. Landscape A4 Word table
##   2. TSV copy of the exact displayed data
##
## Columns:
##   MAG ID
##   Sample
##   Quality
##   Completeness (%)
##   Contamination (%)
##   Domain
##   Phylum
##   Class
##   Order
##   Family
##   Genus
##   Species
##   Genome size (Mbp)
##   GC (%)
##   Total contigs
## ============================================================


suppressPackageStartupMessages({

    library(tidyverse)
    library(officer)
    library(flextable)

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


INPUT_FILE <- file.path(
    PROJECT,
    "semibin2_analysis",
    "MAG_master_table_SemiBin2.tsv"
)


OUTDIR <- file.path(
    PROJECT,
    "appendix_tables"
)


dir.create(
    OUTDIR,
    recursive = TRUE,
    showWarnings = FALSE
)


OUT_DOCX <- file.path(
    OUTDIR,
    "Appendix_All_187_SemiBin2_MAGs.docx"
)


OUT_TSV <- file.path(
    OUTDIR,
    "Appendix_All_187_SemiBin2_MAGs.tsv"
)


## ============================================================
## SAMPLE LABELS
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


sample_order <- names(
    sample_labels
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


mag_raw <- read_tsv(
    INPUT_FILE,
    show_col_types = FALSE
)


cat(
    "\nInput rows: ",
    nrow(mag_raw),
    "\n",
    sep = ""
)


## ============================================================
## COLUMN QC
## ============================================================

required_columns <- c(

    "Genome_ID",
    "Sample_short",
    "Quality",
    "Completeness",
    "Contamination",

    "Domain",
    "Phylum",
    "Class",
    "Order",
    "Family",
    "Genus",
    "Species",

    "Genome_Size",
    "GC_Content",
    "Total_Contigs",
    "GTDB_domain_source"
)


missing_columns <- setdiff(
    required_columns,
    names(mag_raw)
)


if (
    length(
        missing_columns
    ) >
        0
) {

    stop(
        paste0(
            "Missing required columns:\n",
            paste(
                missing_columns,
                collapse = "\n"
            )
        )
    )
}


## ============================================================
## KEEP SEMIBIN2 MAGs
##
## The master table should already contain these, but filtering
## explicitly makes the appendix definition reproducible.
## ============================================================

mag <- mag_raw


if (
    "Binner" %in%
        names(mag)
) {

    mag <- mag %>%
        filter(
            Binner ==
                "SemiBin2"
        )
}


cat(
    "SemiBin2 MAGs: ",
    nrow(mag),
    "\n",
    sep = ""
)


if (
    nrow(mag) !=
        187
) {

    stop(
        paste0(
            "Expected exactly 187 SemiBin2 MAGs, but found ",
            nrow(mag),
            ". Check the master table before producing the appendix."
        )
    )
}


## ============================================================
## CHECK SAMPLE MAPPING
## ============================================================

unmapped_samples <- setdiff(
    unique(
        mag$Sample_short
    ),
    names(
        sample_labels
    )
)


if (
    length(
        unmapped_samples
    ) >
        0
) {

    stop(
        paste0(
            "The following Sample_short values do not have a ",
            "sample label:\n",
            paste(
                unmapped_samples,
                collapse = "\n"
            )
        )
    )
}


## ============================================================
## BUILD COMPACT MAG ID
##
## Example:
##
## semibin2__flye__barcode10_seqs__SemiBin_12
##
## becomes:
##
## MAG_b10_SB12
## ============================================================

mag <- mag %>%
    mutate(

        barcode_number =
            str_match(
                Genome_ID,
                "barcode([0-9]+)"
            )[, 2],

        bin_number =
            str_match(
                Genome_ID,
                "SemiBin_([0-9]+)"
            )[, 2],

        MAG_ID =
            paste0(
                "MAG_b",
                barcode_number,
                "_SB",
                bin_number
            )
    )


if (
    any(
        is.na(
            mag$barcode_number
        )
    ) ||
    any(
        is.na(
            mag$bin_number
        )
    )
) {

    stop(
        "Could not derive barcode/bin numbers from one or more Genome_ID values."
    )
}


if (
    n_distinct(
        mag$MAG_ID
    ) !=
        187
) {

    stop(
        "Compact MAG IDs are not unique."
    )
}


## ============================================================
## CLEAN MISSING TAXONOMY
## ============================================================

clean_taxonomy <- function(x) {

    x <- as.character(x)

    x[
        is.na(x) |
        x == "" |
        x == "NA"
    ] <- "—"

    x
}


## ============================================================
## FORMAT APPENDIX DATA
## ============================================================

appendix <- mag %>%
    mutate(

        Sample =
            unname(
                sample_labels[
                    Sample_short
                ]
            ),

        Sample_short =
            factor(
                Sample_short,
                levels =
                    sample_order
            ),

        bin_number_numeric =
            as.integer(
                bin_number
            ),

        Completeness =
            round(
                as.numeric(
                    Completeness
                ),
                2
            ),

        Contamination =
            round(
                as.numeric(
                    Contamination
                ),
                2
            ),

        Genome_Size_Mbp =
            round(
                as.numeric(
                    Genome_Size
                ) /
                    1e6,
                2
            ),

        ## GC_Content is stored as a fraction in the master table,
        ## e.g. 0.56 = 56%.
        GC_percent =
            round(
                as.numeric(
                    GC_Content
                ) *
                    100,
                1
            ),

        Total_Contigs =
            as.integer(
                Total_Contigs
            ),

        ## ------------------------------------------------------------
        ## Resolve GTDB domain.
        ##
        ## Completely unclassified genomes may not have a usable
        ## classification string, but GTDB_domain_source still stores
        ## whether GTDB-Tk treated the genome as Bacteria or Archaea.
        ## ------------------------------------------------------------

        Domain =
            case_when(

                !is.na(classification) &
                    str_starts(
                        classification,
                        "d__Bacteria"
                    ) ~

                    "Bacteria",


                !is.na(classification) &
                    str_starts(
                        classification,
                        "d__Archaea"
                    ) ~

                    "Archaea",


                str_to_lower(
                    str_trim(
                        as.character(
                            GTDB_domain_source
                        )
                    )
                ) ==
                    "bacteria" ~

                    "Bacteria",


                str_to_lower(
                    str_trim(
                        as.character(
                            GTDB_domain_source
                        )
                    )
                ) %in%
                    c(
                        "archaea",
                        "archea"
                    ) ~

                    "Archaea",


                str_to_lower(
                    str_trim(
                        as.character(
                            Domain
                        )
                    )
                ) ==
                    "bacteria" ~

                    "Bacteria",


                str_to_lower(
                    str_trim(
                        as.character(
                            Domain
                        )
                    )
                ) %in%
                    c(
                        "archaea",
                        "archea"
                    ) ~

                    "Archaea",


                TRUE ~

                    Domain
            ),

        across(
            c(
                Domain,
                Phylum,
                Class,
                Order,
                Family,
                Genus,
                Species
            ),
            clean_taxonomy
        )
    ) %>%
    arrange(
        Sample_short,
        bin_number_numeric
    ) %>%
    transmute(

        `MAG ID` =
            MAG_ID,

        `Sample` =
            Sample,

        `Quality` =
            Quality,

        `Completeness (%)` =
            Completeness,

        `Contamination (%)` =
            Contamination,

        `Domain` =
            Domain,

        `Phylum` =
            Phylum,

        `Class` =
            Class,

        `Order` =
            Order,

        `Family` =
            Family,

        `Genus` =
            Genus,

        `Species` =
            Species,

        `Genome size (Mbp)` =
            Genome_Size_Mbp,

        `GC (%)` =
            GC_percent,

        `Total contigs` =
            Total_Contigs
    )


## ============================================================
## FINAL TABLE QC
## ============================================================

stopifnot(
    nrow(
        appendix
    ) ==
        187
)


cat(
    "\nFinal appendix rows: ",
    nrow(appendix),
    "\n",
    sep = ""
)



cat(
    "\nGTDB_domain_source in raw SemiBin2 table:\n"
)


print(
    table(
        mag$GTDB_domain_source,
        useNA =
            "ifany"
    )
)


cat(
    "\nFinal Domain values in appendix:\n"
)


print(
    table(
        appendix$Domain,
        useNA =
            "ifany"
    )
)


raw_archaea <- sum(

    str_to_lower(
        str_trim(
            as.character(
                mag$GTDB_domain_source
            )
        )
    ) %in%
        c(
            "archaea",
            "archea"
        ),

    na.rm =
        TRUE
)


final_archaea <- sum(
    appendix$Domain ==
        "Archaea",
    na.rm =
        TRUE
)


cat(
    "\nArchaeal MAGs according to GTDB_domain_source: ",
    raw_archaea,
    "\n",
    sep = ""
)


cat(
    "Archaeal MAGs in final appendix table: ",
    final_archaea,
    "\n",
    sep = ""
)


if (
    raw_archaea !=
        final_archaea
) {

    stop(
        paste0(
            "Domain QC failed: raw table contains ",
            raw_archaea,
            " archaeal MAGs but appendix contains ",
            final_archaea,
            "."
        )
    )
}


cat(
    "\nMAGs per sample:\n"
)


print(
    table(
        appendix$Sample
    )
)


cat(
    "\nQuality classes:\n"
)


print(
    table(
        appendix$Quality
    )
)


## ============================================================
## WRITE TSV
## ============================================================

write_tsv(
    appendix,
    OUT_TSV,
    na = ""
)


## ============================================================
## WORD TABLE
## ============================================================

ft <- flextable(
    appendix
)


## ------------------------------------------------------------
## Header labels with line breaks to conserve width
## ------------------------------------------------------------

ft <- set_header_labels(

    ft,

    `Completeness (%)` =
        "Completeness\n(%)",

    `Contamination (%)` =
        "Contamination\n(%)",

    `Genome size (Mbp)` =
        "Genome size\n(Mbp)",

    `GC (%)` =
        "GC\n(%)",

    `Total contigs` =
        "Total\ncontigs"
)


## ------------------------------------------------------------
## General formatting
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
        6.5,
    part =
        "body"
)


ft <- fontsize(
    ft,
    size =
        7,
    part =
        "header"
)


ft <- bold(
    ft,
    part =
        "header"
)


ft <- align(
    ft,
    align =
        "center",
    part =
        "header"
)


## Left-align taxonomy / text fields

ft <- align(

    ft,

    j =
        c(
            "MAG ID",
            "Sample",
            "Quality",
            "Domain",
            "Phylum",
            "Class",
            "Order",
            "Family",
            "Genus",
            "Species"
        ),

    align =
        "left",

    part =
        "body"
)


## Numeric columns centered

ft <- align(

    ft,

    j =
        c(
            "Completeness (%)",
            "Contamination (%)",
            "Genome size (Mbp)",
            "GC (%)",
            "Total contigs"
        ),

    align =
        "center",

    part =
        "body"
)


## ------------------------------------------------------------
## Header background
## ------------------------------------------------------------

ft <- bg(
    ft,
    bg =
        "#D9D9D9",
    part =
        "header"
)


## ------------------------------------------------------------
## Light alternating rows
## ------------------------------------------------------------

even_rows <- seq(
    2,
    nrow(
        appendix
    ),
    by =
        2
)


ft <- bg(
    ft,
    i =
        even_rows,
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
        0.35
)


header_border <- fp_border(
    color =
        "#666666",
    width =
        0.7
)


ft <- border(
    ft,
    border =
        thin_border,
    part =
        "body"
)


ft <- border(
    ft,
    border =
        header_border,
    part =
        "header"
)


## ------------------------------------------------------------
## Compact padding
## ------------------------------------------------------------

ft <- padding(
    ft,
    padding.top =
        1.2,
    padding.bottom =
        1.2,
    padding.left =
        1.5,
    padding.right =
        1.5,
    part =
        "all"
)


## ------------------------------------------------------------
## Fixed widths for landscape A4
##
## Total ~10.95 inches.
## ------------------------------------------------------------

ft <- width(
    ft,
    j =
        "MAG ID",
    width =
        0.75
)


ft <- width(
    ft,
    j =
        "Sample",
    width =
        0.85
)


ft <- width(
    ft,
    j =
        "Quality",
    width =
        0.85
)


ft <- width(
    ft,
    j =
        "Completeness (%)",
    width =
        0.55
)


ft <- width(
    ft,
    j =
        "Contamination (%)",
    width =
        0.55
)


ft <- width(
    ft,
    j =
        "Domain",
    width =
        0.55
)


ft <- width(
    ft,
    j =
        "Phylum",
    width =
        0.80
)


ft <- width(
    ft,
    j =
        "Class",
    width =
        0.80
)


ft <- width(
    ft,
    j =
        "Order",
    width =
        0.85
)


ft <- width(
    ft,
    j =
        "Family",
    width =
        0.85
)


ft <- width(
    ft,
    j =
        "Genus",
    width =
        0.85
)


ft <- width(
    ft,
    j =
        "Species",
    width =
        1.15
)


ft <- width(
    ft,
    j =
        "Genome size (Mbp)",
    width =
        0.65
)


ft <- width(
    ft,
    j =
        "GC (%)",
    width =
        0.45
)


ft <- width(
    ft,
    j =
        "Total contigs",
    width =
        0.45
)


## ------------------------------------------------------------
## Repeat header across pages and prevent individual rows
## splitting across two pages.
## ------------------------------------------------------------

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
## CREATE LANDSCAPE A4 WORD DOCUMENT
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
                0.28,
            bottom =
                0.28,
            left =
                0.28,
            right =
                0.28,
            header =
                0.15,
            footer =
                0.15
        )
)


doc <- body_set_default_section(
    doc,
    landscape_section
)


## ============================================================
## CAPTION
## ============================================================

doc <- body_add_par(

    doc,

    value =
        paste0(
            "Appendix Table. Taxonomic classification and genome ",
            "quality statistics for all 187 SemiBin2 metagenome-",
            "assembled genomes (MAGs)."
        ),

    style =
        "Normal"
)


doc <- body_add_par(
    doc,
    value =
        ""
)


## ============================================================
## ADD TABLE
## ============================================================

doc <- body_add_flextable(
    doc,
    value =
        ft
)


## ============================================================
## SAVE
## ============================================================

print(
    doc,
    target =
        OUT_DOCX
)


cat(
    "\n============================================================\n"
)

cat(
    "APPENDIX TABLE COMPLETE\n"
)

cat(
    "============================================================\n\n"
)


cat(
    "Word document:\n  ",
    OUT_DOCX,
    "\n\n",
    sep = ""
)


cat(
    "TSV table:\n  ",
    OUT_TSV,
    "\n",
    sep = ""
)
