#!/usr/bin/env Rscript


## ================================================================== ##
## STAGE 91
##
## BUILD UNIFIED METADATA FOR FINAL 650-TIP SPECIES TREE
##
## Sources:
##
##   631 GlobDB Methylococcales
##    16 SemiBin2 MAGs
##     2 Umezawa Allocrenothrix genomes
##     1 Methylophaga outgroup
##
## Tracks prepared here:
##
##   1. MCA0421-family homolog detected / not detected
##   2. CheckM2 completeness
##   3. CheckM2 contamination
##
## MCA0421 family:
##
##   old GlobDB family Cluster_00063
##
## For the incomplete MAGs, "Not detected" should be interpreted as:
##
##   no Cluster_00063 homolog detected in the assembled MAG
##
## and not necessarily biological absence.
## ================================================================== ##


suppressPackageStartupMessages({

    library(ape)
    library(readr)
    library(dplyr)
    library(stringr)
    library(tidyr)

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


WF <- file.path(
    ROOT,
    "genome_analysis_workflow"
)


MG <- file.path(
    ROOT,
    "metagenome"
)


PHYLO <- file.path(
    MG,
    "comparative_analysis",
    "final_species_phylogeny"
)


TREE_FILE <- file.path(
    PHYLO,
    "05_rooted",
    "final_methylococcales_GlobDB_MAGs_Umezawa_bac120_ROOTED.tree"
)


OUTDIR <- file.path(
    PHYLO,
    "08_visualization"
)


dir.create(
    OUTDIR,
    recursive = TRUE,
    showWarnings = FALSE
)


OUTPUT <- file.path(
    OUTDIR,
    "final_650_tip_metadata.tsv"
)


## ================================================================== ##
## Existing project files
## ================================================================== ##

OLD_MANIFEST <- file.path(
    WF,
    "17_phylogeny",
    "03_bac120_fasttree",
    "00_inputs",
    "genomes_632.tsv"
)


GLOBDB_TAX <- file.path(
    WF,
    "16_visualization",
    "globdb_r226_taxonomy.tsv"
)


GLOBDB_MODULES <- file.path(
    WF,
    "09_module_occurrence",
    "protein_module_membership.tsv"
)


GLOBDB_CHECKM2 <- file.path(
    WF,
    "17_phylogeny",
    "05_checkm2_quality",
    "01_results",
    "quality_report.tsv"
)


MAG_MASTER <- file.path(
    MG,
    "semibin2_analysis",
    "MAG_master_table_SemiBin2.tsv"
)


MAG_CLUSTER_MAP <- file.path(
    MG,
    "functional_analysis",
    "globdb_cluster_mapping",
    "Methylococcales_GlobDB_cluster_assignments.tsv"
)


UME_CHECKM2 <- file.path(
    PHYLO,
    "06_checkm2_umezawa",
    "01_results",
    "quality_report.tsv"
)


UME_MCA <- file.path(
    PHYLO,
    "07_MCA0421_umezawa",
    "01_search",
    "MCA0421_presence_absence.tsv"
)


OUTGROUP <- "OUT_Methylophaga_nitratireducenticrescens"

MCA_CLUSTER <- "Cluster_00063"


## ================================================================== ##
## Helpers
## ================================================================== ##

strip_fasta_extension <- function(x) {

    x %>%
        basename() %>%
        str_remove("\\.gz$") %>%
        str_remove("\\.(fna|fa|fasta)$")
}


extract_rank <- function(x, rank_prefix) {

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


find_column <- function(df, patterns, description) {

    nm <- names(df)

    hit <- nm[
        Reduce(
            `|`,
            lapply(
                patterns,
                function(p) {
                    str_detect(
                        nm,
                        regex(
                            p,
                            ignore_case = TRUE
                        )
                    )
                }
            )
        )
    ]


    if (
        length(hit) == 0
    ) {

        stop(
            "Could not identify ",
            description,
            " column. Available columns:\n",
            paste(
                nm,
                collapse = "\n"
            )
        )
    }


    hit[1]
}


## ================================================================== ##
## Tree
## ================================================================== ##

tree <- read.tree(
    TREE_FILE
)


if (
    Ntip(tree) != 650
) {

    stop(
        "Expected 650 tree tips; found ",
        Ntip(tree)
    )
}


if (
    !OUTGROUP %in%
        tree$tip.label
) {

    stop(
        "Outgroup is missing from tree."
    )
}


tree_tips <- tibble(
    genome =
        tree$tip.label
)


cat(
    "\nTree tips: ",
    nrow(tree_tips),
    "\n",
    sep = ""
)


## ================================================================== ##
## Explicit MAG tip-name crosswalk
##
## Tree IDs were deliberately shortened in Stage 84.
## ================================================================== ##

mag_crosswalk <- tribble(

    ~original_genome,                                      ~genome,

    "semibin2__flye__barcode10_seqs__SemiBin_14",         "MAG_b10_SB14",
    "semibin2__flye__barcode10_seqs__SemiBin_1",          "MAG_b10_SB1",
    "semibin2__flye__barcode10_seqs__SemiBin_8",          "MAG_b10_SB8",

    "semibin2__flye__barcode12_seqs__SemiBin_5",          "MAG_b12_SB5",
    "semibin2__flye__barcode12_seqs__SemiBin_0",          "MAG_b12_SB0",
    "semibin2__flye__barcode12_seqs__SemiBin_6",          "MAG_b12_SB6",

    "semibin2__flye__barcode13_seqs__SemiBin_1",          "MAG_b13_SB1",

    "semibin2__flye__barcode14_seqs__SemiBin_1",          "MAG_b14_SB1",
    "semibin2__flye__barcode14_seqs__SemiBin_0",          "MAG_b14_SB0",

    "semibin2__flye__barcode15_seqs__SemiBin_0",          "MAG_b15_SB0",
    "semibin2__flye__barcode15_seqs__SemiBin_11",         "MAG_b15_SB11",
    "semibin2__flye__barcode15_seqs__SemiBin_3",          "MAG_b15_SB3",
    "semibin2__flye__barcode15_seqs__SemiBin_10",         "MAG_b15_SB10",

    "semibin2__flye__barcode16_seqs__SemiBin_1",          "MAG_b16_SB1",
    "semibin2__flye__barcode16_seqs__SemiBin_12",         "MAG_b16_SB12",
    "semibin2__flye__barcode16_seqs__SemiBin_11",         "MAG_b16_SB11"
)


if (
    nrow(mag_crosswalk) != 16
) {

    stop(
        "MAG crosswalk should contain exactly 16 genomes."
    )
}


## ================================================================== ##
## GlobDB taxonomy
## ================================================================== ##

glob_tax <- read.delim(
    GLOBDB_TAX,
    header = FALSE,
    sep = "\t",
    quote = "",
    stringsAsFactors = FALSE
)


colnames(glob_tax)[1:2] <- c(
    "genome",
    "taxonomy"
)


glob_tax <- glob_tax %>%
    transmute(

        genome,

        genus =
            extract_rank(
                taxonomy,
                "g__"
            ),

        source =
            "GlobDB"
    )


## ================================================================== ##
## GlobDB MCA0421 / Cluster_00063
## ================================================================== ##

pm <- read_tsv(
    GLOBDB_MODULES,
    show_col_types = FALSE
)


glob_mca_present <- pm %>%
    filter(
        cluster ==
            MCA_CLUSTER
    ) %>%
    distinct(
        genome
    ) %>%
    mutate(
        MCA0421_status =
            "Present"
    )


## ================================================================== ##
## GlobDB CheckM2
##
## Use old manifest to connect CheckM2 file basenames to tree IDs.
## ================================================================== ##

old_manifest <- read_tsv(
    OLD_MANIFEST,
    col_names =
        c(
            "fasta",
            "genome"
        ),
    show_col_types =
        FALSE
) %>%
    filter(
        genome !=
            OUTGROUP
    ) %>%
    mutate(
        checkm_name =
            strip_fasta_extension(
                fasta
            )
    )


glob_quality_raw <- read_tsv(
    GLOBDB_CHECKM2,
    show_col_types =
        FALSE
)


glob_name_col <- find_column(
    glob_quality_raw,
    c("^Name$", "^Genome$"),
    "GlobDB CheckM2 genome-name"
)


glob_comp_col <- find_column(
    glob_quality_raw,
    c("Completeness"),
    "GlobDB completeness"
)


glob_cont_col <- find_column(
    glob_quality_raw,
    c("Contamination"),
    "GlobDB contamination"
)


glob_quality <- glob_quality_raw %>%
    transmute(

        checkm_name =
            strip_fasta_extension(
                .data[[glob_name_col]]
            ),

        completeness =
            as.numeric(
                .data[[glob_comp_col]]
            ),

        contamination =
            as.numeric(
                .data[[glob_cont_col]]
            )
    ) %>%
    right_join(
        old_manifest %>%
            select(
                genome,
                checkm_name
            ),
        by =
            "checkm_name"
    ) %>%
    select(
        genome,
        completeness,
        contamination
    )


if (
    any(
        is.na(
            glob_quality$completeness
        )
    )
) {

    bad <- glob_quality %>%
        filter(
            is.na(
                completeness
            )
        )


    print(bad)


    stop(
        "One or more GlobDB genomes did not match the old CheckM2 report."
    )
}


## ================================================================== ##
## MAG master table:
##
## automatically identify genome, taxonomy, completeness and
## contamination columns.
## ================================================================== ##

mag_master <- read_tsv(
    MAG_MASTER,
    show_col_types =
        FALSE
)


## Find the column containing the full SemiBin genome names.

mag_genome_scores <- sapply(
    mag_master,
    function(x) {

        sum(
            as.character(x) %in%
                mag_crosswalk$original_genome,
            na.rm = TRUE
        )
    }
)


mag_genome_col <- names(
    which.max(
        mag_genome_scores
    )
)


if (
    max(
        mag_genome_scores
    ) <
        16
) {

    stop(
        "Could not identify all 16 MAG genome IDs in MAG_master_table_SemiBin2.tsv.\n",
        "Best matching column: ",
        mag_genome_col,
        " with ",
        max(mag_genome_scores),
        " matches."
    )
}


mag_comp_col <- find_column(
    mag_master,
    c("Completeness"),
    "MAG completeness"
)


mag_cont_col <- find_column(
    mag_master,
    c("Contamination"),
    "MAG contamination"
)


mag_tax_col <- find_column(
    mag_master,
    c("taxonomy", "classification"),
    "MAG taxonomy"
)


mag_meta <- mag_master %>%
    transmute(

        original_genome =
            as.character(
                .data[[mag_genome_col]]
            ),

        taxonomy =
            as.character(
                .data[[mag_tax_col]]
            ),

        completeness =
            as.numeric(
                .data[[mag_comp_col]]
            ),

        contamination =
            as.numeric(
                .data[[mag_cont_col]]
            )
    ) %>%
    inner_join(
        mag_crosswalk,
        by =
            "original_genome"
    ) %>%
    mutate(

        genus =
            extract_rank(
                taxonomy,
                "g__"
            ),

        source =
            "Metagenome MAG"
    )


if (
    nrow(
        mag_meta
    ) !=
        16
) {

    stop(
        "Expected 16 MAG metadata rows; found ",
        nrow(mag_meta)
    )
}


## ================================================================== ##
## MAG MCA0421 mapping
##
## A MAG counts as detected when any passing old-family assignment
## contains Cluster_00063.
## ================================================================== ##

mag_map <- read_tsv(
    MAG_CLUSTER_MAP,
    show_col_types =
        FALSE
)


if (
    !"genome" %in%
        names(
            mag_map
        )
) {

    stop(
        "MAG cluster mapping lacks a 'genome' column."
    )
}


cluster_text_cols <- intersect(
    c(
        "top_scoring_cluster",
        "all_clusters_hit"
    ),
    names(
        mag_map
    )
)


if (
    length(
        cluster_text_cols
    ) ==
        0
) {

    stop(
        "Could not find top_scoring_cluster/all_clusters_hit in MAG mapping table."
    )
}


mag_mca_present <- mag_map %>%
    mutate(

        MCA_hit =
            if_any(
                all_of(
                    cluster_text_cols
                ),
                ~ str_detect(
                    coalesce(
                        as.character(
                            .x
                        ),
                        ""
                    ),
                    fixed(
                        MCA_CLUSTER
                    )
                )
            )
    ) %>%
    filter(
        MCA_hit
    ) %>%
    distinct(
        original_genome =
            genome
    ) %>%
    inner_join(
        mag_crosswalk,
        by =
            "original_genome"
    ) %>%
    transmute(
        genome,
        MCA0421_status =
            "Present"
    )


## ================================================================== ##
## Umezawa quality
## ================================================================== ##

ume_q_raw <- read_tsv(
    UME_CHECKM2,
    show_col_types =
        FALSE
)


ume_name_col <- find_column(
    ume_q_raw,
    c("^Name$", "^Genome$"),
    "Umezawa CheckM2 genome-name"
)


ume_comp_col <- find_column(
    ume_q_raw,
    c("Completeness"),
    "Umezawa completeness"
)


ume_cont_col <- find_column(
    ume_q_raw,
    c("Contamination"),
    "Umezawa contamination"
)


ume_quality <- ume_q_raw %>%
    transmute(

        genome =
            as.character(
                .data[[ume_name_col]]
            ),

        completeness =
            as.numeric(
                .data[[ume_comp_col]]
            ),

        contamination =
            as.numeric(
                .data[[ume_cont_col]]
            )
    )


## ================================================================== ##
## Umezawa MCA0421
## ================================================================== ##

ume_mca <- read_tsv(
    UME_MCA,
    show_col_types =
        FALSE
) %>%
    transmute(

        genome,

        MCA0421_status =
            if_else(
                MCA0421_present ==
                    1,
                "Present",
                "Not detected"
            )
    )


## ================================================================== ##
## Assemble GlobDB metadata
## ================================================================== ##

globdb_tips <- old_manifest %>%
    select(
        genome
    )


globdb_meta <- globdb_tips %>%
    left_join(
        glob_tax,
        by =
            "genome"
    ) %>%
    left_join(
        glob_quality,
        by =
            "genome"
    ) %>%
    left_join(
        glob_mca_present,
        by =
            "genome"
    ) %>%
    mutate(

        MCA0421_status =
            replace_na(
                MCA0421_status,
                "Not detected"
            ),

        source =
            "GlobDB"
    )


## ================================================================== ##
## Assemble MAG metadata
## ================================================================== ##

mag_meta <- mag_meta %>%
    left_join(
        mag_mca_present,
        by =
            "genome"
    ) %>%
    mutate(

        MCA0421_status =
            replace_na(
                MCA0421_status,
                "Not detected"
            )
    )


## ================================================================== ##
## Assemble Umezawa metadata
## ================================================================== ##

ume_meta <- tibble(

    genome =
        c(
            "UME_AF98_Allocrenothrix_methanica",
            "UME_MI19235_Allocrenothrix_methanica"
        ),

    genus =
        c(
            "Allocrenothrix",
            "Allocrenothrix"
        ),

    source =
        c(
            "Umezawa isolate",
            "Umezawa isolate"
        )
) %>%
    left_join(
        ume_quality,
        by =
            "genome"
    ) %>%
    left_join(
        ume_mca,
        by =
            "genome"
    )


## ================================================================== ##
## Outgroup
##
## No CheckM2 or MCA0421 search has been performed for the outgroup.
## Keep these as NA / Not assessed rather than calling absence.
## ================================================================== ##

outgroup_meta <- tibble(

    genome =
        OUTGROUP,

    genus =
        "Methylophaga",

    source =
        "Outgroup",

    completeness =
        NA_real_,

    contamination =
        NA_real_,

    MCA0421_status =
        "Not assessed"
)


## ================================================================== ##
## Final 650-genome metadata
## ================================================================== ##

meta <- bind_rows(

    globdb_meta %>%
        select(
            genome,
            genus,
            source,
            completeness,
            contamination,
            MCA0421_status
        ),

    mag_meta %>%
        select(
            genome,
            genus,
            source,
            completeness,
            contamination,
            MCA0421_status
        ),

    ume_meta %>%
        select(
            genome,
            genus,
            source,
            completeness,
            contamination,
            MCA0421_status
        ),

    outgroup_meta

) %>%
    right_join(
        tree_tips,
        by =
            "genome"
    ) %>%
    mutate(

        display_label =
            paste0(
                genus,
                " | ",
                genome
            ),

        MCA0421_status =
            factor(
                MCA0421_status,
                levels =
                    c(
                        "Present",
                        "Not detected",
                        "Not assessed"
                    )
            )
    )


## ================================================================== ##
## QC
## ================================================================== ##

cat(
    "\n============================================================\n"
)

cat(
    "STAGE 91 - FINAL TREE METADATA QC\n"
)

cat(
    "============================================================\n\n"
)


cat(
    "Rows: ",
    nrow(meta),
    "\n",
    sep = ""
)


cat(
    "Unique genomes: ",
    n_distinct(meta$genome),
    "\n\n",
    sep = ""
)


cat(
    "Sources:\n"
)

print(
    table(
        meta$source,
        useNA =
            "ifany"
    )
)


cat(
    "\nMCA0421:\n"
)

print(
    table(
        meta$MCA0421_status,
        useNA =
            "ifany"
    )
)


cat(
    "\nCompleteness:\n"
)

print(
    summary(
        meta$completeness
    )
)


cat(
    "\nContamination:\n"
)

print(
    summary(
        meta$contamination
    )
)


if (
    nrow(meta) != 650 ||
    n_distinct(meta$genome) != 650
) {

    stop(
        "Final metadata does not contain exactly 650 unique genomes."
    )
}


missing_genus <- meta %>%
    filter(
        is.na(genus) |
        genus ==
            ""
    )


if (
    nrow(
        missing_genus
    ) >
        0
) {

    cat(
        "\nGenomes lacking genus:\n"
    )

    print(
        missing_genus %>%
            select(
                genome,
                source
            ),
        n =
            Inf
    )


    stop(
        "One or more displayed tips lack genus labels."
    )
}


quality_missing <- meta %>%
    filter(
        source !=
            "Outgroup",
        is.na(
            completeness
        ) |
        is.na(
            contamination
        )
    )


if (
    nrow(
        quality_missing
    ) >
        0
) {

    cat(
        "\nNon-outgroup genomes lacking quality values:\n"
    )

    print(
        quality_missing,
        n =
            Inf
    )


    stop(
        "Missing CheckM2 values outside the outgroup."
    )
}


## ================================================================== ##
## Save
## ================================================================== ##

write_tsv(
    meta,
    OUTPUT
)


cat(
    "\nUmezawa genomes:\n"
)


print(
    meta %>%
        filter(
            source ==
                "Umezawa isolate"
        ) %>%
        select(
            genome,
            completeness,
            contamination,
            MCA0421_status
        )
)


cat(
    "\nOutput:\n  ",
    OUTPUT,
    "\n",
    sep = ""
)


cat(
    "\nStage 91 complete.\n"
)
