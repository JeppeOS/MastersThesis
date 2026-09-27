#!/usr/bin/env Rscript


## ================================================================== ##
## STAGE 81
##
## TARGETED Cluster_00121 LOCAL-CONTEXT CHECK
##
## No new clustering is performed.
##
## Questions:
##
##   1. What occurs within ±5 genes of Cluster_00121 across all
##      23 focal regions?
##
##   2. Do recurrent neighboring local families form a conserved
##      ordered subarchitecture?
##
##   3. Are there any non-focal FindMeHemes / FeGenie proteins near
##      Cluster_00121, and if so do they already map to old
##      Cluster_XXXXX families?
##
## Outputs:
##
##   Cluster_00121_plusminus5_gene_context.tsv
##   Cluster_00121_recurrent_subarchitecture_by_region.tsv
##
## This is a diagnostic / interpretation stage only.
## ================================================================== ##


suppressPackageStartupMessages({

  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)

})


## ================================================================== ##
## Paths
## ================================================================== ##

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


local_dir <- file.path(
  indir,
  "local_clustering"
)


outdir <- file.path(
  indir,
  "local_context_check"
)


dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)


neighborhood_file <- file.path(
  indir,
  "Cluster_00121_combined_neighborhood_genes_local_families.tsv"
)


family_summary_file <- file.path(
  local_dir,
  "C121_local_family_summary.tsv"
)


out_context <- file.path(
  outdir,
  "Cluster_00121_plusminus5_gene_context.tsv"
)


out_architecture <- file.path(
  outdir,
  "Cluster_00121_recurrent_subarchitecture_by_region.tsv"
)


## ================================================================== ##
## Constants
## ================================================================== ##

FOCAL_FAMILY <- "C121_local_001"

GENE_RADIUS <- 5


## ================================================================== ##
## Helpers
## ================================================================== ##

require_columns <- function(
  df,
  required,
  table_name
) {

  missing <- setdiff(
    required,
    names(df)
  )


  if (
    length(
      missing
    ) >
      0
  ) {

    stop(
      table_name,
      " is missing required columns: ",
      paste(
        missing,
        collapse = ", "
      )
    )
  }
}


clean_character <- function(x) {

  x <- as.character(
    x
  )


  x[
    is.na(
      x
    )
  ] <- ""


  x
}


## ================================================================== ##
## Read data
## ================================================================== ##

genes <- read_tsv(
  neighborhood_file,
  show_col_types = FALSE
)


families <- read_tsv(
  family_summary_file,
  show_col_types = FALSE
)


## ================================================================== ##
## Validate schema
## ================================================================== ##

required_gene_columns <- c(

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

  "comparison_cog",

  "comparison_product",

  "old_global_mmseq_cluster",

  "fegenie_positive",

  "fegenie_HMMs",

  "findmehemes_positive",

  "number_of_hemes"
)


require_columns(
  genes,
  required_gene_columns,
  "Cluster_00121 neighborhood table"
)


require_columns(
  families,
  c(
    "local_family_id",
    "n_focal_regions",
    "prevalence_all_23",
    "dominant_cog",
    "dominant_product",
    "recurrent_family"
  ),
  "Cluster_00121 family summary"
)


if (
  n_distinct(
    genes$focal_id
  ) !=
    23
) {

  stop(
    "Expected 23 Cluster_00121 focal regions; found ",
    n_distinct(
      genes$focal_id
    ),
    "."
  )
}


## ================================================================== ##
## Numeric cleanup
## ================================================================== ##

genes <- genes %>%
  mutate(

    oriented_gene_offset =
      as.numeric(
        oriented_gene_offset
      ),

    oriented_midpoint_offset_bp =
      as.numeric(
        oriented_midpoint_offset_bp
      ),

    fegenie_positive =
      replace_na(
        as.numeric(
          fegenie_positive
        ),
        0
      ),

    findmehemes_positive =
      replace_na(
        as.numeric(
          findmehemes_positive
        ),
        0
      ),

    number_of_hemes =
      replace_na(
        as.numeric(
          number_of_hemes
        ),
        0
      ),

    comparison_cog =
      clean_character(
        comparison_cog
      ),

    comparison_product =
      clean_character(
        comparison_product
      ),

    old_global_mmseq_cluster =
      clean_character(
        old_global_mmseq_cluster
      ),

    fegenie_HMMs =
      clean_character(
        fegenie_HMMs
      )
  )


## ================================================================== ##
## Identify recurrent non-focal families
##
## Current analysis showed C121_local_002–006 in five regions.
##
## Keep this data-driven:
##
##   recurrent
##   non-focal
##   >= 5 / 23 focal regions
## ================================================================== ##

recurrent_families <- families %>%
  filter(

    local_family_id !=
      FOCAL_FAMILY,

    recurrent_family ==
      1,

    n_focal_regions >=
      5
  ) %>%
  arrange(
    desc(
      n_focal_regions
    ),
    local_family_id
  )


recurrent_ids <- recurrent_families$local_family_id


if (
  length(
    recurrent_ids
  ) ==
    0
) {

  stop(
    "No recurrent non-focal local families were recovered."
  )
}


## ================================================================== ##
## ±5 gene context
## ================================================================== ##

context <- genes %>%
  filter(

    abs(
      oriented_gene_offset
    ) <=
      GENE_RADIUS

  ) %>%
  mutate(

    relationship = case_when(

      is_focal ==
        1 ~

        "focal",

      oriented_gene_offset <
        0 ~

        "upstream",

      oriented_gene_offset >
        0 ~

        "downstream",

      TRUE ~

        "other"
    ),


    recurrent_context_family =
      local_family_id %in%
      recurrent_ids,


    nonfocal_heme_or_fegenie =

      is_focal ==
        0 &

      (
        findmehemes_positive ==
          1 |

        fegenie_positive ==
          1
      )
  ) %>%
  arrange(
    source_dataset,
    focal_id,
    oriented_gene_offset
  )


write_tsv(
  context,
  out_context
)


## ================================================================== ##
## Build recurrent-family architecture for every focal region
##
## Only focal + recurrent neighboring families are included in the
## architecture signature.
##
## Examples:
##
##   C121_local_006(-3,-) ->
##   C121_local_005(-2,+) ->
##   C121_local_004(-1,+) ->
##   C121_local_001(0,+)  ->
##   C121_local_003(+1,+)
##
## The exact observed order is retained rather than assuming such a
## pattern in advance.
## ================================================================== ##

architecture_members <- genes %>%
  filter(

    local_family_id %in%
      c(
        FOCAL_FAMILY,
        recurrent_ids
      )

  ) %>%
  arrange(
    focal_id,
    oriented_gene_offset
  ) %>%
  mutate(

    architecture_token =
      paste0(
        local_family_id,
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
      )
  )


## ================================================================== ##
## Presence matrix
## ================================================================== ##

presence <- architecture_members %>%
  distinct(
    focal_id,
    local_family_id
  ) %>%
  mutate(
    present =
      1L
  ) %>%
  complete(

    focal_id =
      unique(
        genes$focal_id
      ),

    local_family_id =
      c(
        FOCAL_FAMILY,
        recurrent_ids
      ),

    fill =
      list(
        present =
          0L
      )
  ) %>%
  pivot_wider(

    names_from =
      local_family_id,

    values_from =
      present,

    names_prefix =
      "has_"
  )


## ================================================================== ##
## Region metadata
## ================================================================== ##

region_metadata <- genes %>%
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


## ================================================================== ##
## Ordered architecture signatures
## ================================================================== ##

signatures <- architecture_members %>%
  group_by(
    focal_id
  ) %>%
  summarise(

    recurrent_family_count =
      n_distinct(
        local_family_id[
          local_family_id !=
            FOCAL_FAMILY
        ]
      ),

    ordered_recurrent_architecture =
      paste(
        architecture_token,
        collapse =
          " -> "
      ),

    .groups =
      "drop"
  )


## ================================================================== ##
## Determine whether all recurrent non-focal families coexist
## ================================================================== ##

all_recurrent_ids <- recurrent_ids


complete_subarchitecture_regions <- architecture_members %>%
  filter(
    local_family_id %in%
      recurrent_ids
  ) %>%
  distinct(
    focal_id,
    local_family_id
  ) %>%
  count(
    focal_id,
    name =
      "n_recurrent_families_present"
  ) %>%
  mutate(

    has_complete_recurrent_set =
      n_recurrent_families_present ==
      length(
        all_recurrent_ids
      )
  )


## ================================================================== ##
## Final per-region table
## ================================================================== ##

architecture_by_region <- region_metadata %>%
  left_join(
    signatures,
    by =
      "focal_id"
  ) %>%
  left_join(
    complete_subarchitecture_regions,
    by =
      "focal_id"
  ) %>%
  left_join(
    presence,
    by =
      "focal_id"
  ) %>%
  mutate(

    recurrent_family_count =
      replace_na(
        recurrent_family_count,
        0L
      ),

    n_recurrent_families_present =
      replace_na(
        n_recurrent_families_present,
        0L
      ),

    has_complete_recurrent_set =
      replace_na(
        has_complete_recurrent_set,
        FALSE
      ),

    ordered_recurrent_architecture =
      replace_na(
        ordered_recurrent_architecture,
        FOCAL_FAMILY
      )
  ) %>%
  arrange(

    desc(
      has_complete_recurrent_set
    ),

    desc(
      recurrent_family_count
    ),

    source_dataset,

    genome
  )


write_tsv(
  architecture_by_region,
  out_architecture
)


## ================================================================== ##
## Positional summary for recurrent families
## ================================================================== ##

position_summary <- genes %>%
  filter(
    local_family_id %in%
      recurrent_ids
  ) %>%
  group_by(
    local_family_id
  ) %>%
  summarise(

    n_regions =
      n_distinct(
        focal_id
      ),

    median_gene_offset =
      median(
        oriented_gene_offset,
        na.rm = TRUE
      ),

    min_gene_offset =
      min(
        oriented_gene_offset,
        na.rm = TRUE
      ),

    max_gene_offset =
      max(
        oriented_gene_offset,
        na.rm = TRUE
      ),

    median_midpoint_bp =
      median(
        oriented_midpoint_offset_bp,
        na.rm = TRUE
      ),

    dominant_strand =
      names(
        sort(
          table(
            oriented_neighbor_strand
          ),
          decreasing =
            TRUE
        )
      )[1],

    .groups =
      "drop"
  ) %>%
  left_join(

    recurrent_families %>%
      select(
        local_family_id,
        dominant_cog,
        dominant_product
      ),

    by =
      "local_family_id"
  )


## ================================================================== ##
## Nearby non-focal heme / FeGenie proteins
## ================================================================== ##

nearby_candidates <- context %>%
  filter(
    nonfocal_heme_or_fegenie
  ) %>%
  select(

    source_dataset,

    focal_id,

    genome,

    focal_protein_id,

    neighbor_protein_id,

    oriented_gene_offset,

    oriented_neighbor_strand,

    local_family_id,

    fegenie_positive,

    fegenie_HMMs,

    findmehemes_positive,

    number_of_hemes,

    old_global_mmseq_cluster,

    comparison_cog,

    comparison_product
  )


## ================================================================== ##
## Recurrent architecture pattern frequencies
## ================================================================== ##

pattern_summary <- architecture_by_region %>%
  count(

    ordered_recurrent_architecture,

    name =
      "n_regions",

    sort =
      TRUE
  )


## ================================================================== ##
## Terminal report
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
  "STAGE 81 - Cluster_00121 LOCAL-CONTEXT CHECK\n"
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
  "Focal regions:                         ",
  n_distinct(
    genes$focal_id
  ),
  "\n",
  sep = ""
)


cat(
  "Genes inspected within ±5 positions:   ",
  nrow(
    context
  ),
  "\n",
  sep = ""
)


cat(
  "Recurrent non-focal families:          ",
  length(
    recurrent_ids
  ),
  "\n",
  sep = ""
)


cat(
  "Regions with complete recurrent set:   ",
  sum(
    architecture_by_region$has_complete_recurrent_set
  ),
  "\n\n",
  sep = ""
)


cat(
  "Recurrent families:\n\n"
)


print(
  position_summary,
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
  "\n",
  sep = ""
)


cat(
  "MOST COMMON RECURRENT ARCHITECTURE PATTERNS\n"
)


cat(
  paste(
    rep(
      "-",
      90
    ),
    collapse = ""
  ),
  "\n\n",
  sep = ""
)


print(
  pattern_summary,
  n = 15,
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
  "\n",
  sep = ""
)


cat(
  "NON-FOCAL FeGenie / FindMeHemes PROTEINS WITHIN ±5 GENES\n"
)


cat(
  paste(
    rep(
      "-",
      90
    ),
    collapse = ""
  ),
  "\n\n",
  sep = ""
)


if (
  nrow(
    nearby_candidates
  ) ==
    0
) {

  cat(
    "None detected.\n"
  )

} else {

  print(
    nearby_candidates,
    n = Inf,
    width = Inf
  )
}


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
  "OUTPUTS\n"
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
  "±5-gene context:\n  ",
  out_context,
  "\n\n",
  sep = ""
)


cat(
  "Recurrent architecture by region:\n  ",
  out_architecture,
  "\n",
  sep = ""
)
