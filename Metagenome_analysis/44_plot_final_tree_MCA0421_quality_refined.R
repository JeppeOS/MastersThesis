#!/usr/bin/env Rscript

## ================================================================== ##
## STAGE 92
##
## FINAL bac120 TREE - REFINED VERSION
##
## 650 tips:
##
##   631 GlobDB Methylococcales
##    16 SemiBin2 MAGs
##     2 Umezawa Allocrenothrix isolates
##     1 Methylophaga outgroup
##
## Rings, tree -> outward:
##
##   1. Cyc2 focal-module origin
##   2. Module 20 / Cluster_00035
##   3. MCA0421-family homolog
##   4. CheckM2 completeness
##   5. CheckM2 contamination
##   6. Species labels
##
## Refinements:
##
##   - distinct categorical colors across tracks
##   - only completeness uses blue
##   - MCA0421 absence is implicit white / empty
##   - MCA0421 legend shows only "Present"
##   - tip labels use GTDB species rather than accession
##   - MAG labels are bold
##   - larger tip labels
##   - no floating in-gap track labels
##   - alternating grey genus-level shading extends from tree centre
##     to the inner edge of the first track
## ================================================================== ##

suppressPackageStartupMessages({
    library(ape)
    library(ggtree)
    library(ggplot2)
    library(dplyr)
    library(readr)
    library(stringr)
    library(tidyr)
    library(ggnewscale)
    library(scales)
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

PROJECT <- file.path(
    ROOT,
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

OUTDIR <- file.path(
    PHYLO,
    "08_visualization"
)

META_FILE <- file.path(
    OUTDIR,
    "final_650_tip_metadata.tsv"
)

GLOBDB_TAX <- file.path(
    WF,
    "16_visualization",
    "globdb_r226_taxonomy.tsv"
)

GLOBDB_MASTER <- file.path(
    WF,
    "05_integrated_annotations",
    "MASTER_PROTEIN_ANNOTATIONS.tsv"
)

MODULE_FILE <- file.path(
    WF,
    "09_module_occurrence",
    "protein_module_membership.tsv"
)

MAG_MASTER_FILE <- file.path(
    PROJECT,
    "semibin2_analysis",
    "MAG_master_table_SemiBin2.tsv"
)

MAG_FUNCTIONAL_MASTER <- file.path(
    PROJECT,
    "functional_analysis",
    "integrated_annotations",
    "MASTER_PROTEIN_ANNOTATIONS.tsv"
)

MAG_CLUSTER_MAP <- file.path(
    PROJECT,
    "functional_analysis",
    "globdb_cluster_mapping",
    "Methylococcales_GlobDB_cluster_assignments.tsv"
)

PDF_OUT <- file.path(
    OUTDIR,
    "final_bac120_Cyc2_Module20_MCA0421_quality_refined.pdf"
)

PNG_OUT <- file.path(
    OUTDIR,
    "final_bac120_Cyc2_Module20_MCA0421_quality_refined.png"
)

TRACK_OUT <- file.path(
    OUTDIR,
    "final_650_tip_functional_tracks_refined.tsv"
)


FINAL_C121_MMO_TRACK_FILE <- file.path(
    OUTDIR,
    "extra_tracks",
    "final_C121_MMO_tree_tracks.tsv"
)

## ================================================================== ##
## Constants
## ================================================================== ##

OUTGROUP <- "OUT_Methylophaga_nitratireducenticrescens"

FOCAL_M20_CLUSTER <- "Cluster_00035"

## ================================================================== ##
## Figure controls
## ================================================================== ##

OPEN_ANGLE <- 24

TREE_LINEWIDTH <- 0.18

RING_OFFSET_FRACTION <- 0.020
RING_WIDTH_FRACTION  <- 0.018
RING_GAP_FRACTION    <- 0.006

LABEL_PAD_FRACTION <- 0.007

TIP_LABEL_SIZE_NONMAG <- 1.10
TIP_LABEL_SIZE_MAG    <- 1.20

OUTLINE_COLOR <- "grey85"
OUTLINE_LINEWIDTH <- 0.10

## ================================================================== ##
## Colors
##
## Completeness is the ONLY blue track.
## ================================================================== ##

CYC2_COLORS <- c(
    "Module 10"      = "#E66101",
    "Module 35"      = "#5E3C99",
    "Module 10 + 35" = "#E7298A"
)

MODULE20_COLORS <- c(
    "MtoA"                    = "#A6761D",
    "MtrA"                    = "#E6AB02",
    "FeGenie-negative homolog" = "#666666",
    "Multiple types"          = "#FB8072"
)

MCA_COLOR <- c(
    "Present" = "#1B7837"
)



## Cluster_00121

CLUSTER00121_COLOR <- c(

    "Present" =
        "#00A6A6"       # teal
)


## Methane monooxygenase repertoire

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

COMPLETENESS_LOW  <- "#F7FBFF"
COMPLETENESS_HIGH <- "#08519C"

CONTAMINATION_LOW  <- "#FFF5F0"
CONTAMINATION_HIGH <- "#A50F15"

## Alternating genus shading
GENUS_SHADE_LIGHT <- "grey96"
GENUS_SHADE_DARK  <- "grey88"

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
## Helpers
## ================================================================== ##

extract_rank <- function(x, rank_prefix) {
    pattern <- paste0("(?:^|;)", rank_prefix, "([^;]*)")
    hit <- regexec(pattern, x, perl = TRUE)
    matches <- regmatches(x, hit)

    vapply(
        matches,
        function(z) {
            if (length(z) >= 2 && nzchar(z[2])) {
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

    x[is.na(x) | x %in% c("", "NA", "N/A")] <- NA_character_

    x <- str_remove(x, "^s__")
    x <- str_remove(x, "^g__")
    x <- str_replace_all(x, "_", " ")
    x <- str_squish(x)

    x
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
                        regex(p, ignore_case = TRUE)
                    )
                }
            )
        )
    ]

    if (length(hit) == 0) {
        stop(
            "Could not identify ", description, " column.\nAvailable columns:\n",
            paste(nm, collapse = "\n")
        )
    }

    hit[1]
}

make_outline <- function(cells, xmin_value, xmax_value) {
    cells %>%
        mutate(
            xmin = xmin_value,
            xmax = xmax_value
        )
}

## ================================================================== ##
## Tree
## ================================================================== ##

tree <- read.tree(TREE_FILE)

if (Ntip(tree) != 650) {
    stop("Expected 650 displayed tips; found ", Ntip(tree))
}

if (!OUTGROUP %in% tree$tip.label) {
    stop("Outgroup not found in tree.")
}

cat("Displayed tips: ", Ntip(tree), "\n", sep = "")

## ================================================================== ##
## General metadata
## ================================================================== ##

meta <- read_tsv(
    META_FILE,
    show_col_types = FALSE
)

if (nrow(meta) != 650) {
    stop("Metadata should contain exactly 650 rows.")
}

if (!setequal(meta$genome, tree$tip.label)) {
    stop("Tree and metadata genome IDs do not match.")
}

meta <- meta %>%
    mutate(
        completeness = as.numeric(completeness),
        contamination = as.numeric(contamination),
        is_MAG =
            str_detect(source, regex("MAG|Metagenome", ignore_case = TRUE)) |
            str_detect(genome, "^MAG_")
    )

## ================================================================== ##
## Species names - GlobDB
##
## Use GTDB species rank s__ directly.
## ================================================================== ##

glob_tax_raw <- read.delim(
    GLOBDB_TAX,
    header = FALSE,
    sep = "\t",
    quote = "",
    stringsAsFactors = FALSE
)

if (ncol(glob_tax_raw) < 2) {
    stop("GlobDB taxonomy file has fewer than two columns.")
}

colnames(glob_tax_raw)[1:2] <- c("genome", "taxonomy")

glob_species <- glob_tax_raw %>%
    transmute(
        genome,
        genus = clean_taxon(extract_rank(taxonomy, "g__")),
        species = clean_taxon(extract_rank(taxonomy, "s__"))
    ) %>%
    mutate(
        species_label = case_when(
            !is.na(species) & species != "" ~ species,
            !is.na(genus)   & genus   != "" ~ genus,
            TRUE ~ genome
        )
    ) %>%
    select(genome, genus, species_label)

## ================================================================== ##
## Species names - MAGs
## ================================================================== ##

mag_master_tax <- read_tsv(
    MAG_MASTER_FILE,
    show_col_types = FALSE
)

mag_genome_scores <- sapply(
    mag_master_tax,
    function(x) {
        sum(as.character(x) %in% mag_crosswalk$original_genome, na.rm = TRUE)
    }
)

mag_genome_col <- names(which.max(mag_genome_scores))

if (max(mag_genome_scores) < 16) {
    stop("Could not identify all 16 MAG genome IDs in MAG master table.")
}

mag_tax_col <- find_column(
    mag_master_tax,
    c("taxonomy", "classification"),
    "MAG taxonomy"
)

mag_species <- mag_master_tax %>%
    transmute(
        original_genome = as.character(.data[[mag_genome_col]]),
        taxonomy = as.character(.data[[mag_tax_col]])
    ) %>%
    inner_join(
        mag_crosswalk,
        by = "original_genome"
    ) %>%
    mutate(
        genus = clean_taxon(extract_rank(taxonomy, "g__")),
        species = clean_taxon(extract_rank(taxonomy, "s__")),
        species_label = case_when(
            !is.na(species) & species != "" ~ species,
            !is.na(genus)   & genus   != "" ~ genus,
            TRUE ~ genome
        )
    ) %>%
    select(genome, genus, species_label)

if (nrow(mag_species) != 16) {
    stop("Expected 16 MAG species-label rows; found ", nrow(mag_species))
}

## ================================================================== ##
## Species labels - Umezawa + outgroup
## ================================================================== ##

special_species <- tribble(
    ~genome,                                      ~genus,              ~species_label,
    "UME_AF98_Allocrenothrix_methanica",          "Allocrenothrix",   "Allocrenothrix methanica AF98",
    "UME_MI19235_Allocrenothrix_methanica",       "Allocrenothrix",   "Allocrenothrix methanica MI19235",
    "OUT_Methylophaga_nitratireducenticrescens",  "Methylophaga",     "Methylophaga nitratireducenticrescens"
)

## ================================================================== ##
## Combine all tip species labels
## ================================================================== ##

species_labels <- bind_rows(
    glob_species,
    mag_species,
    special_species
) %>%
    distinct(
        genome,
        .keep_all = TRUE
    )

meta <- meta %>%
    select(-any_of(c("species_label", "genus"))) %>%
    left_join(
        species_labels,
        by = "genome"
    )

missing_species <- meta %>%
    filter(is.na(species_label) | species_label == "")

if (nrow(missing_species) > 0) {
    cat("\nTips lacking species labels:\n")
    print(missing_species %>% select(genome, source), n = Inf)
    stop("One or more tips lack species labels.")
}

## ================================================================== ##
## GlobDB integrated annotations
## ================================================================== ##

glob_master <- read_tsv(
    GLOBDB_MASTER,
    show_col_types = FALSE
)

glob_master <- glob_master %>%
    mutate(
        fegenie_HMMs = coalesce(as.character(fegenie_HMMs), ""),
        fegenie_positive_norm =
            fegenie_positive %in% c(1, "1", TRUE, "TRUE", "True")
    ) %>%
    select(
        genome,
        protein_id,
        fegenie_positive_norm,
        fegenie_HMMs
    )

## ================================================================== ##
## GlobDB protein-module membership
## ================================================================== ##

pm <- read_tsv(
    MODULE_FILE,
    show_col_types = FALSE
) %>%
    select(
        genome,
        protein_id,
        cluster,
        module
    )

cluster_module_lookup <- pm %>%
    filter(!is.na(cluster), !is.na(module)) %>%
    distinct(cluster, module)

## ================================================================== ##
## GlobDB Cyc2 focal-module origin
## ================================================================== ##

glob_cyc2_proteins <- glob_master %>%
    filter(
        fegenie_positive_norm,
        str_detect(fegenie_HMMs, fixed("Cyc2", ignore_case = TRUE))
    ) %>%
    distinct(genome, protein_id)

glob_cyc2_joined <- glob_cyc2_proteins %>%
    left_join(
        pm %>% select(genome, protein_id, module),
        by = c("genome", "protein_id")
    )

glob_cyc2_by_genome <- glob_cyc2_joined %>%
    group_by(genome) %>%
    summarise(
        has_Module10 = any(module == "Module_10", na.rm = TRUE),
        has_Module35 = any(module == "Module_35", na.rm = TRUE),
        .groups = "drop"
    ) %>%
    mutate(
        Cyc2_module = case_when(
            has_Module10 & has_Module35 ~ "Module 10 + 35",
            has_Module10 ~ "Module 10",
            has_Module35 ~ "Module 35",
            TRUE ~ NA_character_
        )
    ) %>%
    select(genome, Cyc2_module)

## ================================================================== ##
## GlobDB Module 20 / Cluster_00035
## ================================================================== ##

glob_m20_focal <- pm %>%
    filter(
        module == "Module_20",
        cluster == FOCAL_M20_CLUSTER
    ) %>%
    left_join(
        glob_master %>%
            select(
                genome,
                protein_id,
                fegenie_positive_norm,
                fegenie_HMMs
            ),
        by = c("genome", "protein_id")
    ) %>%
    mutate(
        protein_type = case_when(
            str_detect(fegenie_HMMs, fixed("MtoA", ignore_case = TRUE)) ~ "MtoA",
            str_detect(fegenie_HMMs, fixed("MtrA", ignore_case = TRUE)) ~ "MtrA",
            !fegenie_positive_norm ~ "FeGenie-negative homolog",
            TRUE ~ "Other"
        )
    )

glob_module20_by_genome <- glob_m20_focal %>%
    filter(
        protein_type %in% c("MtoA", "MtrA", "FeGenie-negative homolog")
    ) %>%
    group_by(genome) %>%
    summarise(
        types = list(unique(na.omit(protein_type))),
        .groups = "drop"
    ) %>%
    mutate(
        Module20_type = vapply(
            types,
            function(x) {
                if (length(x) > 1) {
                    "Multiple types"
                } else if (length(x) == 1) {
                    x[[1]]
                } else {
                    NA_character_
                }
            },
            character(1)
        )
    ) %>%
    select(genome, Module20_type)

## ================================================================== ##
## MAG functional annotations
## ================================================================== ##

mag_master <- read_tsv(
    MAG_FUNCTIONAL_MASTER,
    show_col_types = FALSE
) %>%
    mutate(
        fegenie_HMMs = coalesce(as.character(fegenie_HMMs), ""),
        fegenie_positive_norm =
            fegenie_positive %in% c(1, "1", TRUE, "TRUE", "True")
    )

mag_map <- read_tsv(
    MAG_CLUSTER_MAP,
    show_col_types = FALSE
) %>%
    mutate(
        top_scoring_cluster = as.character(top_scoring_cluster)
    )

## ================================================================== ##
## MAG Cyc2 focal-module origin
## ================================================================== ##

mag_cyc2 <- mag_master %>%
    filter(
        fegenie_positive_norm,
        str_detect(fegenie_HMMs, fixed("Cyc2", ignore_case = TRUE))
    ) %>%
    select(
        original_genome = genome,
        protein_id
    ) %>%
    left_join(
        mag_map %>%
            select(
                original_genome = genome,
                protein_id,
                top_scoring_cluster
            ),
        by = c("original_genome", "protein_id")
    ) %>%
    left_join(
        cluster_module_lookup,
        by = c("top_scoring_cluster" = "cluster")
    ) %>%
    inner_join(
        mag_crosswalk,
        by = "original_genome"
    )

mag_cyc2_by_genome <- mag_cyc2 %>%
    group_by(genome) %>%
    summarise(
        has_Module10 = any(module == "Module_10", na.rm = TRUE),
        has_Module35 = any(module == "Module_35", na.rm = TRUE),
        .groups = "drop"
    ) %>%
    mutate(
        Cyc2_module = case_when(
            has_Module10 & has_Module35 ~ "Module 10 + 35",
            has_Module10 ~ "Module 10",
            has_Module35 ~ "Module 35",
            TRUE ~ NA_character_
        )
    ) %>%
    select(genome, Cyc2_module)

## ================================================================== ##
## MAG Module 20 / Cluster_00035
## ================================================================== ##

mag_m20 <- mag_map %>%
    filter(top_scoring_cluster == FOCAL_M20_CLUSTER) %>%
    select(
        original_genome = genome,
        protein_id
    ) %>%
    left_join(
        mag_master %>%
            transmute(
                original_genome = genome,
                protein_id,
                fegenie_positive_norm,
                fegenie_HMMs
            ),
        by = c("original_genome", "protein_id")
    ) %>%
    inner_join(
        mag_crosswalk,
        by = "original_genome"
    ) %>%
    mutate(
        protein_type = case_when(
            str_detect(fegenie_HMMs, fixed("MtoA", ignore_case = TRUE)) ~ "MtoA",
            str_detect(fegenie_HMMs, fixed("MtrA", ignore_case = TRUE)) ~ "MtrA",
            !fegenie_positive_norm ~ "FeGenie-negative homolog",
            TRUE ~ "Other"
        )
    )

mag_module20_by_genome <- mag_m20 %>%
    filter(
        protein_type %in% c("MtoA", "MtrA", "FeGenie-negative homolog")
    ) %>%
    group_by(genome) %>%
    summarise(
        types = list(unique(na.omit(protein_type))),
        .groups = "drop"
    ) %>%
    mutate(
        Module20_type = vapply(
            types,
            function(x) {
                if (length(x) > 1) {
                    "Multiple types"
                } else if (length(x) == 1) {
                    x[[1]]
                } else {
                    NA_character_
                }
            },
            character(1)
        )
    ) %>%
    select(genome, Module20_type)

## ================================================================== ##
## Merge functional tracks
## ================================================================== ##

cyc2_by_genome <- bind_rows(
    glob_cyc2_by_genome,
    mag_cyc2_by_genome
) %>%
    group_by(genome) %>%
    summarise(
        has_Module10 = any(Cyc2_module %in% c("Module 10", "Module 10 + 35")),
        has_Module35 = any(Cyc2_module %in% c("Module 35", "Module 10 + 35")),
        .groups = "drop"
    ) %>%
    mutate(
        Cyc2_module = case_when(
            has_Module10 & has_Module35 ~ "Module 10 + 35",
            has_Module10 ~ "Module 10",
            has_Module35 ~ "Module 35",
            TRUE ~ NA_character_
        )
    ) %>%
    select(genome, Cyc2_module)

module20_by_genome <- bind_rows(
    glob_module20_by_genome,
    mag_module20_by_genome
) %>%
    group_by(genome) %>%
    summarise(
        types = list(unique(na.omit(Module20_type))),
        .groups = "drop"
    ) %>%
    mutate(
        Module20_type = vapply(
            types,
            function(x) {
                if (length(x) > 1) {
                    "Multiple types"
                } else if (length(x) == 1) {
                    x[[1]]
                } else {
                    NA_character_
                }
            },
            character(1)
        )
    ) %>%
    select(genome, Module20_type)

## ================================================================== ##
## Combined 650-tip annotation table
## ================================================================== ##

ann <- meta %>%
    left_join(
        cyc2_by_genome,
        by = "genome"
    ) %>%
    left_join(
        module20_by_genome,
        by = "genome"
    ) %>%
    mutate(
        Cyc2_module = factor(
            Cyc2_module,
            levels = c("Module 10", "Module 35", "Module 10 + 35")
        ),
        Module20_type = factor(
            Module20_type,
            levels = c("MtoA", "MtrA", "FeGenie-negative homolog", "Multiple types")
        ),
        MCA0421_plot = if_else(
            as.character(MCA0421_status) == "Present",
            "Present",
            NA_character_
        )
    )

## ================================================================== ##
## Final Cluster_00121 + methane-monooxygenase tracks
## ================================================================== ##

if (
    !file.exists(
        FINAL_C121_MMO_TRACK_FILE
    )
) {

    stop(
        "Missing final Cluster_00121/MMO track file: ",
        FINAL_C121_MMO_TRACK_FILE,
        "\nRun Stage 93H first."
    )
}


final_extra_tracks <- read_tsv(
    FINAL_C121_MMO_TRACK_FILE,
    show_col_types = FALSE
)


if (
    nrow(
        final_extra_tracks
    ) != 650
) {

    stop(
        "Expected 650 rows in final extra-track table; found ",
        nrow(
            final_extra_tracks
        )
    )
}


if (
    n_distinct(
        final_extra_tracks$genome
    ) != 650
) {

    stop(
        "Final extra-track table does not contain 650 unique genomes."
    )
}


ann <- ann %>%
    left_join(

        final_extra_tracks %>%
            select(
                genome,
                Cluster00121_plot,
                MMO_status
            ),

        by = "genome"
    ) %>%
    mutate(

        Cluster00121_plot =
            if_else(
                Cluster00121_plot ==
                    "Present",
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


if (
    any(
        is.na(
            ann$MMO_status
        )
    )
) {

    stop(
        "One or more final tree tips lacks an MMO status."
    )
}


write_tsv(
    ann,
    TRACK_OUT
)

## ================================================================== ##
## Base fan tree
## ================================================================== ##

p_base <- ggtree(
    tree,
    layout = "fan",
    open.angle = OPEN_ANGLE,
    linewidth = TREE_LINEWIDTH
)

tree_data <- p_base$data

tree_radius <- max(
    tree_data$x,
    na.rm = TRUE
)

tip_geom <- tree_data %>%
    filter(isTip) %>%
    transmute(
        genome = label,
        y = y,
        angle_raw = angle
    ) %>%
    left_join(
        ann %>%
            select(
                genome,
                species_label,
                genus,
                source,
                is_MAG,
                Cyc2_module,
                Module20_type,
                MCA0421_plot,

                Cluster00121_plot,

                MMO_status,
                completeness,
                contamination
            ),
        by = "genome"
    )

## ================================================================== ##
## Ring geometry
## ================================================================== ##

RING_OFFSET <- tree_radius * RING_OFFSET_FRACTION
RING_WIDTH  <- tree_radius * RING_WIDTH_FRACTION
RING_GAP    <- tree_radius * RING_GAP_FRACTION
LABEL_PAD   <- tree_radius * LABEL_PAD_FRACTION

CYC2_XMIN <- tree_radius + RING_OFFSET
CYC2_XMAX <- CYC2_XMIN + RING_WIDTH

MODULE20_XMIN <- CYC2_XMAX + RING_GAP
MODULE20_XMAX <- MODULE20_XMIN + RING_WIDTH

MCA_XMIN <- MODULE20_XMAX + RING_GAP
MCA_XMAX <- MCA_XMIN + RING_WIDTH

C121_XMIN <-
    MCA_XMAX +
    RING_GAP


C121_XMAX <-
    C121_XMIN +
    RING_WIDTH


MMO_XMIN <-
    C121_XMAX +
    RING_GAP


MMO_XMAX <-
    MMO_XMIN +
    RING_WIDTH


COMP_XMIN <-
    MMO_XMAX +
    RING_GAP


COMP_XMAX <-
    COMP_XMIN +
    RING_WIDTH

CONT_XMIN <- COMP_XMAX + RING_GAP
CONT_XMAX <- CONT_XMIN + RING_WIDTH

LABEL_X <- CONT_XMAX + LABEL_PAD

## ================================================================== ##
## Tip-label orientation
## ================================================================== ##

tip_geom <- tip_geom %>%
    mutate(
        text_angle = angle_raw,
        text_hjust = 0
    )

flip <- tip_geom$text_angle > 90 & tip_geom$text_angle < 270

tip_geom$text_angle[flip] <-
    (tip_geom$text_angle[flip] + 180) %% 360

tip_geom$text_hjust[flip] <- 1

## ================================================================== ##
## Genus-level alternating background sectors
##
## Taxonomy / tip labels remain untouched.
##
## "shade_group" is ONLY used for the grey background. This allows
## visually patching small holes where a placeholder genus sits inside
## an otherwise continuous clade.
## ================================================================== ##

tip_shading <- tip_geom %>%
    arrange(y) %>%
    mutate(

        ## Start from actual GTDB genus when available.
        shade_group = case_when(
            !is.na(genus) & genus != "" ~ genus,
            TRUE ~ species_label
        ),

        ## Manual visual-group corrections only for shading.
        shade_group = case_when(

            ## Patch SPIREOTU_00336771 into the surrounding
            ## Methylovulum shading block.
            genus == "SPIREOTU_00336771" ~ "Methylovulum",

            ## Treat UBA4132 and Allocrenothrix as one continuous
            ## shading block.
            genus %in% c("UBA4132", "Allocrenothrix") ~
                "UBA4132_Allocrenothrix",

            TRUE ~ shade_group
        )
    )

## ================================================================== ##
## Patch isolated one-tip holes automatically
##
## If one tip is flanked on both sides by the same shading group,
## assign it to that surrounding shading group. This only affects
## background shading, not taxonomy or labels.
## ================================================================== ##

tip_shading <- tip_shading %>%
    mutate(
        previous_group = lag(shade_group),
        next_group = lead(shade_group)
    ) %>%
    mutate(
        shade_group = if_else(
            !is.na(previous_group) &
                !is.na(next_group) &
                previous_group == next_group &
                shade_group != previous_group,
            previous_group,
            shade_group
        )
    ) %>%
    select(
        -previous_group,
        -next_group
    )

## ================================================================== ##
## Collapse contiguous tips belonging to the same shading group
## ================================================================== ##

genus_blocks <- tip_shading %>%
    mutate(
        new_block = row_number() == 1L |
            shade_group != lag(shade_group),

        genus_block = cumsum(new_block)
    ) %>%
    group_by(genus_block, shade_group) %>%
    summarise(
        ymin = min(y) - 0.5,
        ymax = max(y) + 0.5,
        .groups = "drop"
    ) %>%
    mutate(
        xmin = 0,
        xmax = CYC2_XMIN,
        shade_group_number = genus_block %% 2L
    )

## ================================================================== ##
## Grey shading rectangles
## ================================================================== ##

shade_light <- geom_rect(
    data = genus_blocks %>%
        filter(shade_group_number == 0L),
    aes(
        xmin = xmin,
        xmax = xmax,
        ymin = ymin,
        ymax = ymax
    ),
    inherit.aes = FALSE,
    fill = GENUS_SHADE_LIGHT,
    color = NA
)

shade_dark <- geom_rect(
    data = genus_blocks %>%
        filter(shade_group_number == 1L),
    aes(
        xmin = xmin,
        xmax = xmax,
        ymin = ymin,
        ymax = ymax
    ),
    inherit.aes = FALSE,
    fill = GENUS_SHADE_DARK,
    color = NA
)

## Put shading behind the tree.
p_base$layers <- c(
    list(shade_light, shade_dark),
    p_base$layers
)

cat(
    "\nVisual genus-shading blocks: ",
    nrow(genus_blocks),
    "\n",
    sep = ""
)

## ================================================================== ##
## Ring outline cells
## ================================================================== ##

outline_cells <- tip_geom %>%
    transmute(
        genome,
        y,
        ymin = y - 0.48,
        ymax = y + 0.48
    )

cyc2_outline <- make_outline(outline_cells, CYC2_XMIN, CYC2_XMAX)
module20_outline <- make_outline(outline_cells, MODULE20_XMIN, MODULE20_XMAX)
mca_outline <- make_outline(outline_cells, MCA_XMIN, MCA_XMAX)
c121_outline <- make_outline(
    outline_cells,
    C121_XMIN,
    C121_XMAX
)


mmo_outline <- make_outline(
    outline_cells,
    MMO_XMIN,
    MMO_XMAX
)


comp_outline <- make_outline(outline_cells, COMP_XMIN, COMP_XMAX)
cont_outline <- make_outline(outline_cells, CONT_XMIN, CONT_XMAX)

## ================================================================== ##
## Ring data
## ================================================================== ##

cyc2_ring <- tip_geom %>%
    filter(!is.na(Cyc2_module)) %>%
    select(genome, Cyc2_module) %>%
    left_join(outline_cells, by = "genome") %>%
    mutate(
        xmin = CYC2_XMIN,
        xmax = CYC2_XMAX
    )

module20_ring <- tip_geom %>%
    filter(!is.na(Module20_type)) %>%
    select(genome, Module20_type) %>%
    left_join(outline_cells, by = "genome") %>%
    mutate(
        xmin = MODULE20_XMIN,
        xmax = MODULE20_XMAX
    )

mca_ring <- tip_geom %>%
    filter(MCA0421_plot == "Present") %>%
    select(genome, MCA0421_plot) %>%
    left_join(outline_cells, by = "genome") %>%
    mutate(
        xmin = MCA_XMIN,
        xmax = MCA_XMAX
    )

## Cluster_00121 - positive cells only

c121_ring <- tip_geom %>%
    filter(
        Cluster00121_plot ==
            "Present"
    ) %>%
    select(
        genome,
        Cluster00121_plot
    ) %>%
    left_join(
        outline_cells,
        by = "genome"
    ) %>%
    mutate(
        xmin =
            C121_XMIN,
        xmax =
            C121_XMAX
    )


## Methane monooxygenase repertoire - all 650 genomes

mmo_ring <- tip_geom %>%
    select(
        genome,
        MMO_status
    ) %>%
    left_join(
        outline_cells,
        by = "genome"
    ) %>%
    mutate(
        xmin =
            MMO_XMIN,
        xmax =
            MMO_XMAX
    )


comp_ring <- tip_geom %>%
    select(genome, completeness) %>%
    left_join(outline_cells, by = "genome") %>%
    mutate(
        xmin = COMP_XMIN,
        xmax = COMP_XMAX
    )

cont_ring <- tip_geom %>%
    select(genome, contamination) %>%
    left_join(outline_cells, by = "genome") %>%
    mutate(
        xmin = CONT_XMIN,
        xmax = CONT_XMAX
    )

## ================================================================== ##
## Plot
## ================================================================== ##

p <- p_base +

    ## -------------------------------------------------------------- ##
    ## Ring outlines
    ## -------------------------------------------------------------- ##

    geom_rect(
        data = cyc2_outline,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        inherit.aes = FALSE,
        fill = NA,
        color = OUTLINE_COLOR,
        linewidth = OUTLINE_LINEWIDTH
    ) +

    geom_rect(
        data = module20_outline,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        inherit.aes = FALSE,
        fill = NA,
        color = OUTLINE_COLOR,
        linewidth = OUTLINE_LINEWIDTH
    ) +

    geom_rect(
        data = mca_outline,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        inherit.aes = FALSE,
        fill = NA,
        color = OUTLINE_COLOR,
        linewidth = OUTLINE_LINEWIDTH
    ) +

    geom_rect(
        data =
            c121_outline,
        aes(
            xmin =
                xmin,
            xmax =
                xmax,
            ymin =
                ymin,
            ymax =
                ymax
        ),
        inherit.aes =
            FALSE,
        fill =
            NA,
        color =
            OUTLINE_COLOR,
        linewidth =
            OUTLINE_LINEWIDTH
    ) +


    geom_rect(
        data =
            mmo_outline,
        aes(
            xmin =
                xmin,
            xmax =
                xmax,
            ymin =
                ymin,
            ymax =
                ymax
        ),
        inherit.aes =
            FALSE,
        fill =
            NA,
        color =
            OUTLINE_COLOR,
        linewidth =
            OUTLINE_LINEWIDTH
    ) +



    geom_rect(
        data = comp_outline,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        inherit.aes = FALSE,
        fill = NA,
        color = OUTLINE_COLOR,
        linewidth = OUTLINE_LINEWIDTH
    ) +

    geom_rect(
        data = cont_outline,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        inherit.aes = FALSE,
        fill = NA,
        color = OUTLINE_COLOR,
        linewidth = OUTLINE_LINEWIDTH
    ) +

    ## -------------------------------------------------------------- ##
    ## 1. Cyc2 focal-module origin
    ## -------------------------------------------------------------- ##

    geom_rect(
        data = cyc2_ring,
        aes(
            xmin = xmin,
            xmax = xmax,
            ymin = ymin,
            ymax = ymax,
            fill = Cyc2_module
        ),
        inherit.aes = FALSE,
        color = NA
    ) +

    scale_fill_manual(
        values = CYC2_COLORS,
        breaks = c("Module 10", "Module 35", "Module 10 + 35"),
        drop = TRUE,
        na.value = "white",
        name = "Cyc2 focal-module\norigin"
    ) +

    ggnewscale::new_scale_fill() +

    ## -------------------------------------------------------------- ##
    ## 2. Module 20 / Cluster_00035
    ## -------------------------------------------------------------- ##

    geom_rect(
        data = module20_ring,
        aes(
            xmin = xmin,
            xmax = xmax,
            ymin = ymin,
            ymax = ymax,
            fill = Module20_type
        ),
        inherit.aes = FALSE,
        color = NA
    ) +

    scale_fill_manual(
        values = MODULE20_COLORS,
        breaks = c("MtoA", "MtrA", "FeGenie-negative homolog", "Multiple types"),
        drop = TRUE,
        na.value = "white",
        name = "Module 20\nCluster_00035"
    ) +

    ggnewscale::new_scale_fill() +

    ## -------------------------------------------------------------- ##
    ## 3. MCA0421
    ## -------------------------------------------------------------- ##

    geom_rect(
        data = mca_ring,
        aes(
            xmin = xmin,
            xmax = xmax,
            ymin = ymin,
            ymax = ymax,
            fill = MCA0421_plot
        ),
        inherit.aes = FALSE,
        color = NA
    ) +

    scale_fill_manual(
        values = MCA_COLOR,
        breaks = "Present",
        labels = "Present",
        drop = FALSE,
        name = "MCA0421-family\nhomolog"
    ) +

    ggnewscale::new_scale_fill() +

    ## -------------------------------------------------------------- ##
    ## 4. Cluster_00121
    ##
    ## Presence only. Absence remains white / empty.
    ## -------------------------------------------------------------- ##

    geom_rect(
        data =
            c121_ring,
        aes(
            xmin =
                xmin,
            xmax =
                xmax,
            ymin =
                ymin,
            ymax =
                ymax,
            fill =
                Cluster00121_plot
        ),
        inherit.aes =
            FALSE,
        color =
            NA
    ) +


    scale_fill_manual(

        values =
            CLUSTER00121_COLOR,

        breaks =
            "Present",

        labels =
            "Present",

        drop =
            FALSE,

        name =
            "Cluster_00121"
    ) +


    ggnewscale::new_scale_fill() +


    ## -------------------------------------------------------------- ##
    ## 5. Methane monooxygenase repertoire
    ##
    ## All 650 genomes.
    ##
    ## Particulate MMO = complete canonical pMMO OR complete pXMO.
    ## -------------------------------------------------------------- ##

    geom_rect(
        data =
            mmo_ring,
        aes(
            xmin =
                xmin,
            xmax =
                xmax,
            ymin =
                ymin,
            ymax =
                ymax,
            fill =
                MMO_status
        ),
        inherit.aes =
            FALSE,
        color =
            NA
    ) +


    scale_fill_manual(

        values =
            MMO_COLORS,

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
            "Methane\nmonooxygenase"
    ) +


    ggnewscale::new_scale_fill() +


    ## -------------------------------------------------------------- ##
    ## 6. Completeness
    ## -------------------------------------------------------------- ##

    geom_rect(
        data = comp_ring,
        aes(
            xmin = xmin,
            xmax = xmax,
            ymin = ymin,
            ymax = ymax,
            fill = completeness
        ),
        inherit.aes = FALSE,
        color = NA
    ) +

    scale_fill_gradient(
        low = COMPLETENESS_LOW,
        high = COMPLETENESS_HIGH,
        limits = c(0, 100),
        breaks = c(0, 25, 50, 75, 100),
        oob = squish,
        na.value = "grey85",
        name = "CheckM2\ncompleteness (%)"
    ) +

    ggnewscale::new_scale_fill() +

    ## -------------------------------------------------------------- ##
    ## 7. Contamination
    ## -------------------------------------------------------------- ##

    geom_rect(
        data = cont_ring,
        aes(
            xmin = xmin,
            xmax = xmax,
            ymin = ymin,
            ymax = ymax,
            fill = contamination
        ),
        inherit.aes = FALSE,
        color = NA
    ) +

    scale_fill_gradient(
        low = CONTAMINATION_LOW,
        high = CONTAMINATION_HIGH,
        na.value = "grey85",
        name = "CheckM2\ncontamination (%)"
    ) +

    ## -------------------------------------------------------------- ##
    ## Species labels - ordinary genomes
    ## -------------------------------------------------------------- ##

    geom_text(
        data = tip_geom %>% filter(!is_MAG),
        aes(
            x = LABEL_X,
            y = y,
            label = species_label,
            angle = text_angle,
            hjust = text_hjust
        ),
        inherit.aes = FALSE,
        size = TIP_LABEL_SIZE_NONMAG,
        fontface = "plain"
    ) +

    ## -------------------------------------------------------------- ##
    ## Species labels - MAGs
    ## -------------------------------------------------------------- ##

    geom_text(
        data = tip_geom %>% filter(is_MAG),
        aes(
            x = LABEL_X,
            y = y,
            label = species_label,
            angle = text_angle,
            hjust = text_hjust
        ),
        inherit.aes = FALSE,
        size = TIP_LABEL_SIZE_MAG,
        fontface = "bold"
    ) +

    ## -------------------------------------------------------------- ##
    ## Figure extent
    ## -------------------------------------------------------------- ##

    xlim(
        0,
        LABEL_X + tree_radius * 0.47
    ) +

    labs(
        title = "Methylococcales bac120 phylogeny"
    ) +

    theme(
        plot.title = element_text(
            hjust = 0.5,
            face = "bold",
            size = 15
        ),
        legend.position = "right",
        legend.box = "vertical",
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8),
        plot.margin = margin(20, 35, 20, 20)
    )

## ================================================================== ##
## Save
## ================================================================== ##

ggsave(
    PDF_OUT,
    p,
    width = 30,
    height = 30,
    limitsize = FALSE
)

ggsave(
    PNG_OUT,
    p,
    width = 30,
    height = 30,
    dpi = 250,
    limitsize = FALSE
)

## ================================================================== ##
## Report
## ================================================================== ##

cat("\n============================================================\n")
cat("STAGE 92 REFINED COMPLETE\n")
cat("============================================================\n\n")

cat("Displayed tips: ", Ntip(tree), "\n", sep = "")
cat("\nMAG tips: ", sum(tip_geom$is_MAG, na.rm = TRUE), "\n", sep = "")
cat("\nGenus shading blocks: ", nrow(genus_blocks), "\n", sep = "")

cat("\nCyc2 focal-module origin:\n")
print(table(ann$Cyc2_module, useNA = "ifany"))

cat("\nModule 20 / Cluster_00035:\n")
print(table(ann$Module20_type, useNA = "ifany"))

cat("\nMCA0421 positive genomes:\n")
print(sum(ann$MCA0421_plot == "Present", na.rm = TRUE))

cat("\nExample MAG species labels:\n")
print(
    tip_geom %>%
        filter(is_MAG) %>%
        select(genome, species_label),
    n = Inf
)

cat("\nUmezawa labels:\n")
print(
    tip_geom %>%
        filter(
            genome %in% c(
                "UME_AF98_Allocrenothrix_methanica",
                "UME_MI19235_Allocrenothrix_methanica"
            )
        ) %>%
        select(genome, species_label)
)

cat("\nTrack table:\n  ", TRACK_OUT, "\n", sep = "")
cat("\nPNG:\n  ", PNG_OUT, "\n", sep = "")
cat("\nPDF:\n  ", PDF_OUT, "\n", sep = "")
