#!/usr/bin/env Rscript


## ================================================================== ##
## STAGE 93G
##
## QC THE 209 ORIGINAL MMO-"NONE" GENOMES
##
## Questions:
##
##   1. Are MMO-none genomes lower completeness?
##   2. Are they concentrated in particular genera?
##   3. Are partial-pMMO genomes lower completeness than genomes with
##      complete pMMO?
##   4. Which high-quality genomes genuinely have no detectable
##      pmo / pxm / mmo components?
##
## No sequence analysis is performed here.
## ================================================================== ##


suppressPackageStartupMessages({

    library(dplyr)
    library(readr)
    library(stringr)

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


WF <- file.path(
    ROOT,
    "genome_analysis_workflow"
)


PHYLO <- file.path(
    PROJECT,
    "comparative_analysis",
    "final_species_phylogeny"
)


QC_FILE <- file.path(
    PHYLO,
    "08_visualization",
    "extra_tracks",
    "all_650_pmo_pxm_smmo_QC.tsv"
)


META_FILE <- file.path(
    PHYLO,
    "08_visualization",
    "final_650_tip_metadata.tsv"
)


GLOBDB_TAX <- file.path(
    WF,
    "16_visualization",
    "globdb_r226_taxonomy.tsv"
)


MAG_MASTER <- file.path(
    PROJECT,
    "semibin2_analysis",
    "MAG_master_table_SemiBin2.tsv"
)


OUTDIR <- file.path(
    PHYLO,
    "08_visualization",
    "extra_tracks",
    "MMO_none_QC"
)


dir.create(
    OUTDIR,
    recursive = TRUE,
    showWarnings = FALSE
)


## ================================================================== ##
## Helpers
## ================================================================== ##

extract_rank <- function(
    x,
    rank_prefix
) {

    pattern <- paste0(
        "(?:^|;)",
        rank_prefix,
        "([^;]*)"
    )


    hit <- regexec(
        pattern,
        x,
        perl = TRUE
    )


    matches <- regmatches(
        x,
        hit
    )


    vapply(

        matches,

        function(z) {

            if (
                length(z) >= 2 &&
                nzchar(z[2])
            ) {

                z[2]

            } else {

                NA_character_
            }
        },

        character(1)
    )
}


clean_taxon <- function(x) {

    x <- as.character(x)

    x[
        is.na(x) |
        x %in%
            c(
                "",
                "NA",
                "N/A"
            )
    ] <- NA_character_


    x <- str_remove(
        x,
        "^g__"
    )


    x
}


## ================================================================== ##
## Read main data
## ================================================================== ##

qc <- read_tsv(
    QC_FILE,
    show_col_types = FALSE
)


meta <- read_tsv(
    META_FILE,
    show_col_types = FALSE
) %>%
    mutate(

        completeness =
            as.numeric(
                completeness
            ),

        contamination =
            as.numeric(
                contamination
            )
    )


## ================================================================== ##
## GlobDB genus
## ================================================================== ##

glob_tax <- read.delim(

    GLOBDB_TAX,

    header = FALSE,

    sep = "\t",

    quote = "",

    stringsAsFactors = FALSE
)


colnames(
    glob_tax
)[1:2] <- c(
    "genome",
    "taxonomy"
)


glob_genus <- glob_tax %>%
    transmute(

        genome,

        genus =
            clean_taxon(
                extract_rank(
                    taxonomy,
                    "g__"
                )
            )
    )


## ================================================================== ##
## MAG genus
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


mag <- read_tsv(
    MAG_MASTER,
    show_col_types = FALSE
)


genome_col <- names(mag)[
    sapply(
        mag,
        function(x) {
            sum(
                as.character(x) %in%
                    mag_crosswalk$original_genome,
                na.rm = TRUE
            )
        }
    ) ==
        16
][1]


tax_col <- names(mag)[
    str_detect(
        names(mag),
        regex(
            "taxonomy|classification",
            ignore_case = TRUE
        )
    )
][1]


mag_genus <- mag %>%
    transmute(

        original_genome =
            as.character(
                .data[[genome_col]]
            ),

        taxonomy =
            as.character(
                .data[[tax_col]]
            )
    ) %>%
    inner_join(
        mag_crosswalk,
        by = "original_genome"
    ) %>%
    transmute(

        genome,

        genus =
            clean_taxon(
                extract_rank(
                    taxonomy,
                    "g__"
                )
            )
    )


## ================================================================== ##
## Special tips
## ================================================================== ##

special_genus <- tribble(

    ~genome,                                     ~genus,

    "UME_AF98_Allocrenothrix_methanica",         "Allocrenothrix",

    "UME_MI19235_Allocrenothrix_methanica",      "Allocrenothrix",

    "OUT_Methylophaga_nitratireducenticrescens", "Methylophaga"
)


taxonomy <- bind_rows(

    glob_genus,

    mag_genus,

    special_genus

) %>%
    distinct(
        genome,
        .keep_all = TRUE
    )


## ================================================================== ##
## Combined QC table
## ================================================================== ##

dat <- qc %>%
    left_join(
        meta %>%
            select(
                genome,
                source,
                completeness,
                contamination
            ),
        by = "genome"
    ) %>%
    left_join(
        taxonomy,
        by = "genome"
    ) %>%
    mutate(

        genus =
            coalesce(
                genus,
                "Unclassified"
            ),

        quality_bin =
            case_when(

                completeness >= 90 &
                contamination <= 5 ~

                    ">=90% complete, <=5% contam.",

                completeness >= 70 ~

                    "70-<90% complete",

                completeness >= 50 ~

                    "50-<70% complete",

                completeness < 50 ~

                    "<50% complete",

                TRUE ~

                    "Quality unavailable"
            ),

        detection_group =
            case_when(

                canonical_pMMO == 1 &
                sMMO == 1 ~

                    "pMMO + sMMO",

                canonical_pMMO == 1 ~

                    "pMMO",

                sMMO == 1 ~

                    "sMMO",

                pXMO == 1 ~

                    "pXMO only",

                partial_pMMO == 1 ~

                    "partial pMMO",

                partial_pXMO == 1 ~

                    "partial pXMO",

                partial_sMMO == 1 ~

                    "partial sMMO",

                TRUE ~

                    "no MMO component detected"
            )
    )


write_tsv(
    dat,
    file.path(
        OUTDIR,
        "all_650_MMO_taxonomy_quality_QC.tsv"
    )
)


## ================================================================== ##
## Original None genomes
## ================================================================== ##

none <- dat %>%
    filter(
        original_MMO_status ==
            "None"
    )


write_tsv(
    none,
    file.path(
        OUTDIR,
        "original_None_with_taxonomy_quality.tsv"
    )
)


## ================================================================== ##
## Genus summary
## ================================================================== ##

genus_summary <- none %>%
    count(
        genus,
        detection_group,
        name = "n"
    ) %>%
    arrange(
        desc(n),
        genus
    )


write_tsv(
    genus_summary,
    file.path(
        OUTDIR,
        "original_None_by_genus.tsv"
    )
)


## ================================================================== ##
## Quality summary
## ================================================================== ##

quality_summary <- none %>%
    count(
        quality_bin,
        detection_group,
        name = "n"
    ) %>%
    arrange(
        quality_bin,
        detection_group
    )


write_tsv(
    quality_summary,
    file.path(
        OUTDIR,
        "original_None_by_quality.tsv"
    )
)


## ================================================================== ##
## High-quality genomes with absolutely no MMO component
## ================================================================== ##

hq_none <- none %>%
    filter(

        detection_group ==
            "no MMO component detected",

        completeness >=
            90,

        contamination <=
            5
    ) %>%
    arrange(
        genus,
        genome
    )


write_tsv(
    hq_none,
    file.path(
        OUTDIR,
        "HIGH_QUALITY_no_MMO_component.tsv"
    )
)


## ================================================================== ##
## Report
## ================================================================== ##

cat(
    "\n============================================================\n"
)

cat(
    "STAGE 93G - MMO NONE QC\n"
)

cat(
    "============================================================\n\n"
)


cat(
    "Original None genomes: ",
    nrow(none),
    "\n\n",
    sep = ""
)


cat(
    "Detection classes among original None:\n"
)


print(
    table(
        none$detection_group
    )
)


cat(
    "\nCompleteness summary for original None:\n"
)


print(
    summary(
        none$completeness
    )
)


cat(
    "\nQuality bins:\n"
)


print(
    table(
        none$quality_bin,
        useNA = "ifany"
    )
)


cat(
    "\nHigh-quality genomes with NO detectable MMO component: ",
    nrow(hq_none),
    "\n",
    sep = ""
)


cat(
    "\nTop genera among original None:\n"
)


print(
    none %>%
        count(
            genus,
            sort = TRUE
        ) %>%
        head(
            25
        ),
    n = 25
)


cat(
    "\nHigh-quality no-component genomes by genus:\n"
)


print(
    hq_none %>%
        count(
            genus,
            sort = TRUE
        ),
    n = Inf
)


cat(
    "\nOutputs:\n"
)

cat(
    "  ",
    OUTDIR,
    "\n",
    sep = ""
)
