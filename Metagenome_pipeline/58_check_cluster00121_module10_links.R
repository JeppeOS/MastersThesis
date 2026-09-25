#!/usr/bin/env Rscript


## ================================================================== ##
## STAGE 82
##
## TARGETED Cluster_00121 <-> MODULE-10 PROXIMITY CHECK
##
## No new clustering.
##
## Search the complete existing Cluster_00121 neighborhood extraction
## for proteins assigned to Cyc2 / Module-10-associated old GlobDB
## cytochrome families.
##
## Families:
##
##   Cluster_00050  Cyc2
##   Cluster_00091  CytC553-like
##   Cluster_00064  CytC553-like
##   Cluster_00069  1-heme
##   Cluster_00222  CytC553-like
##   Cluster_00141  accessory CytC553-like
##
## Question:
##
## Is the MOTU40_034781 linkage unique, or do other Cluster_00121
## loci contain the same families elsewhere in the extracted window?
## ================================================================== ##


suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
})


project <- file.path(
  Sys.getenv("HOME"),
  "methanotrophs",
  "methanotroph_project",
  "jeppe",
  "metagenome"
)


indir <- file.path(
  project,
  "comparative_analysis",
  "Cluster_00121"
)


input_file <- file.path(
  indir,
  "Cluster_00121_combined_neighborhood_genes_local_families.tsv"
)


outdir <- file.path(
  indir,
  "module10_proximity_check"
)


dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)


out_hits <- file.path(
  outdir,
  "Cluster_00121_Module10_family_proximity_hits.tsv"
)


out_regions <- file.path(
  outdir,
  "Cluster_00121_Module10_family_proximity_by_region.tsv"
)


## ================================================================== ##
## Families of interest
## ================================================================== ##

module10_families <- c(

  "Cluster_00050" =
    "Cyc2",

  "Cluster_00091" =
    "CytC553-like C00091",

  "Cluster_00064" =
    "CytC553-like C00064",

  "Cluster_00069" =
    "1-heme C00069",

  "Cluster_00222" =
    "CytC553-like C00222",

  "Cluster_00141" =
    "Accessory CytC553-like"
)


## ================================================================== ##
## Read
## ================================================================== ##

x <- read_tsv(
  input_file,
  show_col_types = FALSE
)


required <- c(

  "source_dataset",

  "focal_id",

  "genome",

  "focal_protein_id",

  "neighbor_protein_id",

  "oriented_gene_offset",

  "oriented_midpoint_offset_bp",

  "oriented_neighbor_strand",

  "is_focal",

  "local_family_id",

  "old_global_mmseq_cluster",

  "number_of_hemes",

  "comparison_cog",

  "comparison_product"
)


missing <- setdiff(
  required,
  names(x)
)


if (
  length(missing) >
    0
) {

  stop(
    "Missing columns: ",
    paste(
      missing,
      collapse = ", "
    )
  )
}


if (
  n_distinct(
    x$focal_id
  ) !=
    23
) {

  stop(
    "Expected 23 Cluster_00121 regions."
  )
}


## ================================================================== ##
## Clean
## ================================================================== ##

x <- x %>%
  mutate(

    oriented_gene_offset =
      as.numeric(
        oriented_gene_offset
      ),

    oriented_midpoint_offset_bp =
      as.numeric(
        oriented_midpoint_offset_bp
      ),

    number_of_hemes =
      replace_na(
        as.numeric(
          number_of_hemes
        ),
        0
      ),

    old_global_mmseq_cluster =
      replace_na(
        as.character(
          old_global_mmseq_cluster
        ),
        ""
      )
  )


## ================================================================== ##
## Find Module-10-related proteins
## ================================================================== ##

hits <- x %>%
  filter(

    is_focal == 0,

    old_global_mmseq_cluster %in%
      names(
        module10_families
      )
  ) %>%
  mutate(

    module10_family_label =
      unname(
        module10_families[
          old_global_mmseq_cluster
        ]
      ),

    absolute_gene_distance =
      abs(
        oriented_gene_offset
      ),

    absolute_bp_distance =
      abs(
        oriented_midpoint_offset_bp
      )
  ) %>%
  arrange(
    focal_id,
    absolute_gene_distance,
    oriented_gene_offset
  )


write_tsv(
  hits,
  out_hits
)


## ================================================================== ##
## Per-region summary
## ================================================================== ##

region_meta <- x %>%
  filter(
    is_focal ==
      1
  ) %>%
  distinct(
    focal_id,
    source_dataset,
    genome,
    focal_protein_id
  )


region_summary <- hits %>%
  group_by(
    focal_id
  ) %>%
  summarise(

    n_module10_proteins =
      n(),

    n_module10_families =
      n_distinct(
        old_global_mmseq_cluster
      ),

    has_Cyc2 =
      any(
        old_global_mmseq_cluster ==
          "Cluster_00050"
      ),

    closest_module10_gene_offset =
      oriented_gene_offset[
        which.min(
          absolute_gene_distance
        )
      ],

    closest_module10_bp =
      oriented_midpoint_offset_bp[
        which.min(
          absolute_bp_distance
        )
      ],

    families =
      paste(
        unique(
          paste0(
            old_global_mmseq_cluster,
            " [",
            module10_family_label,
            "]"
          )
        ),
        collapse =
          "; "
      ),

    ordered_hits =
      paste(
        paste0(
          old_global_mmseq_cluster,
          "(",
          if_else(
            oriented_gene_offset >
              0,
            paste0(
              "+",
              oriented_gene_offset
            ),
            as.character(
              oriented_gene_offset
            )
          ),
          ",",
          oriented_neighbor_strand,
          ")"
        ),
        collapse =
          " -> "
      ),

    .groups =
      "drop"
  )


region_summary <- region_meta %>%
  left_join(
    region_summary,
    by =
      "focal_id"
  ) %>%
  mutate(

    n_module10_proteins =
      replace_na(
        n_module10_proteins,
        0L
      ),

    n_module10_families =
      replace_na(
        n_module10_families,
        0L
      ),

    has_Cyc2 =
      replace_na(
        has_Cyc2,
        FALSE
      ),

    families =
      replace_na(
        families,
        ""
      ),

    ordered_hits =
      replace_na(
        ordered_hits,
        ""
      )
  ) %>%
  arrange(
    desc(
      has_Cyc2
    ),
    desc(
      n_module10_families
    ),
    desc(
      n_module10_proteins
    ),
    source_dataset,
    genome
  )


write_tsv(
  region_summary,
  out_regions
)


## ================================================================== ##
## Report
## ================================================================== ##

cat(
  "\n",
  paste(
    rep(
      "=",
      90
    ),
    collapse = ""
  ),
  "\n",
  sep = ""
)


cat(
  "STAGE 82 - Cluster_00121 / MODULE-10 PROXIMITY\n"
)


cat(
  paste(
    rep(
      "=",
      90
    ),
    collapse = ""
  ),
  "\n\n",
  sep = ""
)


cat(
  "Cluster_00121 regions:                  ",
  n_distinct(
    x$focal_id
  ),
  "\n",
  sep = ""
)


cat(
  "Regions with Module-10-family protein:  ",
  sum(
    region_summary$n_module10_proteins >
      0
  ),
  "\n",
  sep = ""
)


cat(
  "Regions with nearby Cluster_00050 Cyc2: ",
  sum(
    region_summary$has_Cyc2
  ),
  "\n",
  sep = ""
)


cat(
  "Total Module-10-family proteins:        ",
  nrow(
    hits
  ),
  "\n\n",
  sep = ""
)


cat(
  "REGIONS WITH HITS\n\n"
)


print(
  region_summary %>%
    filter(
      n_module10_proteins >
        0
    ) %>%
    select(
      source_dataset,
      focal_id,
      genome,
      focal_protein_id,
      n_module10_proteins,
      n_module10_families,
      has_Cyc2,
      closest_module10_gene_offset,
      closest_module10_bp,
      ordered_hits
    ),
  n = Inf,
  width = Inf
)


cat(
  "\n",
  paste(
    rep(
      "-",
      90
    ),
    collapse = ""
  ),
  "\n"
)


cat(
  "INDIVIDUAL HITS\n"
)


cat(
  paste(
    rep(
      "-",
      90
    ),
    collapse = ""
  ),
  "\n\n"
)


print(
  hits %>%
    select(

      source_dataset,

      focal_id,

      genome,

      focal_protein_id,

      neighbor_protein_id,

      oriented_gene_offset,

      oriented_midpoint_offset_bp,

      oriented_neighbor_strand,

      old_global_mmseq_cluster,

      module10_family_label,

      number_of_hemes,

      comparison_product
    ),
  n = Inf,
  width = Inf
)


cat(
  "\nOutputs:\n"
)


cat(
  "  ",
  out_hits,
  "\n",
  sep = ""
)


cat(
  "  ",
  out_regions,
  "\n",
  sep = ""
)
