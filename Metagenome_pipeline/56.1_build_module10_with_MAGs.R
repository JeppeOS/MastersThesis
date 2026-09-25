#!/usr/bin/env Rscript

## ================================================================== ##
## STAGE 80A - ADD TWO Cluster_00050 MAG NEIGHBORHOODS TO MODULE 10
##
## Preserves the original 12 selected GlobDB Module-10 regions and
## appends:
##
##   barcode10 SemiBin_1 | Methylumidiphilus | contig_25_317
##   barcode10 SemiBin_8 | Methylumidiphilus | contig_26_118
##
## The output uses exactly the compact plotting schema required by the
## existing Stage-16B Module-10 plotting script.
## ================================================================== ##

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
})


## ================================================================== ##
## Paths
## ================================================================== ##

project <- path.expand(
  "~/methanotrophs/methanotroph_project/jeppe/metagenome"
)

old_workflow <- path.expand(
  "~/methanotrophs/methanotroph_project/jeppe/genome_analysis_workflow"
)

old_gene_map <- file.path(
  old_workflow,
  "16_visualization",
  "16A3_module10_representative_selection",
  "Module_10_selected_gene_map_data.tsv"
)

old_contig_qc <- file.path(
  old_workflow,
  "16_visualization",
  "16C_cross_figure_qc",
  "selected_gene_map_contig_boundary_qc.tsv"
)

mg_neighborhood_file <- file.path(
  project,
  "functional_analysis",
  "methylococcales_neighborhoods",
  "priority_observed_neighborhood_genes_local_families.tsv"
)

cluster_assignment_file <- file.path(
  project,
  "functional_analysis",
  "globdb_cluster_mapping",
  "Methylococcales_GlobDB_cluster_assignments.tsv"
)

outdir <- file.path(
  project,
  "comparative_analysis",
  "Module_10_with_MAGs"
)

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

output_table <- file.path(
  outdir,
  "Module_10_selected_gene_map_data_with_MAGs.tsv"
)

output_contig_qc <- file.path(
  outdir,
  "Module_10_with_MAGs_contig_boundary_qc.tsv"
)


## ================================================================== ##
## Helpers
## ================================================================== ##

require_columns <- function(df, required, table_name) {

  missing <- setdiff(
    required,
    names(df)
  )

  if (length(missing) > 0) {

    stop(
      table_name,
      " is missing required columns: ",
      paste(missing, collapse = ", ")
    )
  }
}


first_existing_column <- function(
  df,
  candidates,
  regex = NULL,
  required = FALSE
) {

  hit <- candidates[
    candidates %in% names(df)
  ]

  if (length(hit) > 0) {
    return(hit[1])
  }


  if (!is.null(regex)) {

    hit <- names(df)[
      str_detect(
        names(df),
        regex(regex, ignore_case = TRUE)
      )
    ]

    if (length(hit) > 0) {
      return(hit[1])
    }
  }


  if (required) {

    stop(
      "Could not identify required column. Candidates: ",
      paste(candidates, collapse = ", ")
    )
  }


  NA_character_
}


first_nonempty <- function(
  df,
  candidates
) {

  present <- candidates[
    candidates %in% names(df)
  ]


  out <- rep(
    "",
    nrow(df)
  )


  for (column in present) {

    value <- as.character(
      df[[column]]
    )

    use <- (
      is.na(out) |
      out == ""
    ) &
      !is.na(value) &
      value != ""


    out[use] <- value[use]
  }


  out
}


numeric_column <- function(
  df,
  candidates,
  default = 0
) {

  column <- first_existing_column(
    df,
    candidates
  )


  if (is.na(column)) {
    return(rep(default, nrow(df)))
  }


  suppressWarnings(
    as.numeric(
      df[[column]]
    )
  ) %>%
    replace_na(default)
}


character_column <- function(
  df,
  candidates,
  default = ""
) {

  column <- first_existing_column(
    df,
    candidates
  )


  if (is.na(column)) {
    return(rep(default, nrow(df)))
  }


  out <- as.character(
    df[[column]]
  )

  out[is.na(out)] <- default

  out
}


## ================================================================== ##
## Read old validated Module-10 table
## ================================================================== ##

message("Reading original 12-region Module-10 figure table...")

old <- read_tsv(
  old_gene_map,
  show_col_types = FALSE
)


required_plot_columns <- c(
  "genome",
  "protein_id",
  "plot_region_id",
  "plot_order",
  "track_label",
  "taxonomy_display",
  "plot_start_bp",
  "plot_end_bp",
  "plot_strand",
  "display_class",
  "local_family_id",
  "neighbor_mmseq_cluster",
  "neighbor_fegenie_HMMs",
  "neighbor_number_of_hemes",
  "neighbor_globdb_cog",
  "neighbor_globdb_product"
)


require_columns(
  old,
  required_plot_columns,
  "Original Module-10 table"
)


if (n_distinct(old$plot_region_id) != 12) {
  stop("Original Module-10 table does not contain 12 regions.")
}


old_plot <- old %>%
  select(
    all_of(
      required_plot_columns
    )
  ) %>%
  mutate(

    plot_order =
      as.integer(plot_order) + 2L,

    source_group =
      "GlobDB genomes"
  )


## ================================================================== ##
## Read metagenome neighborhoods
## ================================================================== ##

message("Reading metagenome neighborhood table...")

mg <- read_tsv(
  mg_neighborhood_file,
  show_col_types = FALSE
)


required_mg_columns <- c(
  "genome",
  "focal_id",
  "focal_protein_id",
  "focal_contig",
  "focal_start",
  "focal_end",
  "focal_strand",
  "neighbor_protein_id",
  "neighbor_start",
  "neighbor_end",
  "neighbor_strand",
  "oriented_gene_offset",
  "oriented_midpoint_offset_bp",
  "local_family_id"
)


require_columns(
  mg,
  required_mg_columns,
  "Metagenome neighborhood table"
)


target_focals <- c(
  "contig_25_317",
  "contig_26_118"
)


mg <- mg %>%
  filter(
    focal_protein_id %in%
      target_focals
  )


if (
  n_distinct(
    mg$focal_protein_id
  ) != 2
) {

  stop(
    "Did not recover both Cluster_00050 MAG focal proteins."
  )
}


## ================================================================== ##
## Old GlobDB-cluster assignments for metagenome heme proteins
##
## Prefer an old_global_mmseq_cluster column if Stage 60 already
## transferred it into the neighborhood table. Otherwise use the
## standalone Methylococcales mapping table.
## ================================================================== ##

if (
  "old_global_mmseq_cluster" %in%
    names(mg)
) {

  mg$mapped_old_cluster <-
    as.character(
      mg$old_global_mmseq_cluster
    )

} else {

  mg$mapped_old_cluster <- ""


  if (
    file.exists(
      cluster_assignment_file
    )
  ) {

    assignments <- read_tsv(
      cluster_assignment_file,
      show_col_types = FALSE
    )


    protein_col <- first_existing_column(
      assignments,
      c(
        "protein_id",
        "query_protein_id",
        "query",
        "query_id"
      ),
      regex = "protein|query",
      required = TRUE
    )


    cluster_col <- first_existing_column(
      assignments,
      c(
        "top_cluster",
        "top_cluster_id",
        "assigned_cluster",
        "best_cluster",
        "old_global_mmseq_cluster"
      ),
      regex = "(top|best|assigned).*cluster|cluster.*(top|best)",
      required = TRUE
    )


    cluster_map <- assignments %>%
      transmute(

        neighbor_protein_id =
          as.character(
            .data[[protein_col]]
          ),

        mapped_old_cluster =
          as.character(
            .data[[cluster_col]]
          )
      ) %>%
      filter(
        !is.na(neighbor_protein_id),
        neighbor_protein_id != ""
      ) %>%
      distinct(
        neighbor_protein_id,
        .keep_all = TRUE
      )


    mg <- mg %>%
      select(
        -mapped_old_cluster
      ) %>%
      left_join(
        cluster_map,
        by = "neighbor_protein_id"
      )
  }
}


mg$mapped_old_cluster[
  is.na(
    mg$mapped_old_cluster
  )
] <- ""


## Force the two focal proteins to their validated old family.

mg$mapped_old_cluster[
  mg$neighbor_protein_id %in%
    target_focals
] <- "Cluster_00050"


## ================================================================== ##
## Extract annotation fields
## ================================================================== ##

mg$fegenie_positive <- numeric_column(
  mg,
  c(
    "fegenie_positive",
    "neighbor_fegenie_positive"
  )
)


mg$findmehemes_positive <- numeric_column(
  mg,
  c(
    "findmehemes_positive",
    "neighbor_findmehemes_positive"
  )
)


mg$heme_count <- numeric_column(
  mg,
  c(
    "number_of_hemes",
    "neighbor_number_of_hemes"
  )
)


mg$fegenie_hmms <- character_column(
  mg,
  c(
    "fegenie_HMMs",
    "neighbor_fegenie_HMMs"
  )
)


mg$cog_for_plot <- first_nonempty(
  mg,
  c(
    "comparison_cog",
    "globdb_cog",
    "direct_COG20_cog",
    "neighbor_globdb_cog"
  )
)


mg$product_for_plot <- first_nonempty(
  mg,
  c(
    "comparison_product",
    "globdb_product",
    "direct_COG20_function",
    "neighbor_globdb_product"
  )
)


## ================================================================== ##
## Reorient coordinates around the focal Cyc2 midpoint
##
## All focal Cyc2 genes point to the right.
## ================================================================== ##

mg_plot <- mg %>%
  mutate(

    focal_start =
      as.numeric(focal_start),

    focal_end =
      as.numeric(focal_end),

    neighbor_start =
      as.numeric(neighbor_start),

    neighbor_end =
      as.numeric(neighbor_end),


    focal_midpoint =
      (
        focal_start +
        focal_end
      ) / 2,


    plot_start_bp = case_when(

      focal_strand == "+" ~
        neighbor_start -
        focal_midpoint,

      focal_strand == "-" ~
        focal_midpoint -
        neighbor_end,

      TRUE ~
        NA_real_
    ),


    plot_end_bp = case_when(

      focal_strand == "+" ~
        neighbor_end -
        focal_midpoint,

      focal_strand == "-" ~
        focal_midpoint -
        neighbor_start,

      TRUE ~
        NA_real_
    ),


    plot_strand = case_when(

      "oriented_neighbor_strand" %in%
        names(mg) ~

        as.character(
          oriented_neighbor_strand
        ),

      focal_strand == "+" ~
        as.character(
          neighbor_strand
        ),

      focal_strand == "-" &
        neighbor_strand == "+" ~
        "-",

      focal_strand == "-" &
        neighbor_strand == "-" ~
        "+",

      TRUE ~
        ""
    )
  )


if (
  anyNA(mg_plot$plot_start_bp) ||
  anyNA(mg_plot$plot_end_bp)
) {

  stop(
    "Failed to calculate MAG plotting coordinates."
  )
}


## ================================================================== ##
## MAG display metadata
## ================================================================== ##

mg_plot <- mg_plot %>%
  mutate(

    plot_region_id = case_when(

      focal_protein_id ==
        "contig_25_317" ~

        "MAG_C00050_barcode10_SemiBin_1",

      focal_protein_id ==
        "contig_26_118" ~

        "MAG_C00050_barcode10_SemiBin_8"
    ),


    plot_order = case_when(

      focal_protein_id ==
        "contig_25_317" ~
        1L,

      focal_protein_id ==
        "contig_26_118" ~
        2L
    ),


    taxonomy_display =
      "Methylumidiphilus",


    track_label = case_when(

      focal_protein_id ==
        "contig_25_317" ~

        "Methylumidiphilus [barcode10 SemiBin_1]",

      focal_protein_id ==
        "contig_26_118" ~

        "Methylumidiphilus [barcode10 SemiBin_8]"
    ),


    display_class = case_when(

      neighbor_protein_id ==
        focal_protein_id ~

        "focal",

      fegenie_positive == 1 ~

        "fegenie_other",

      findmehemes_positive == 1 ~

        "findmehemes_other",

      TRUE ~

        "background_gene"
    ),


    source_group =
      "Metagenome MAGs"
  )


## ================================================================== ##
## Compact MAG table matching the old plotting schema
## ================================================================== ##

mg_plot <- mg_plot %>%
  transmute(

    genome =
      as.character(genome),

    protein_id =
      as.character(
        neighbor_protein_id
      ),

    plot_region_id,

    plot_order,

    track_label,

    taxonomy_display,

    plot_start_bp,

    plot_end_bp,

    plot_strand,

    display_class,

    local_family_id =
      as.character(
        local_family_id
      ),

    neighbor_mmseq_cluster =
      mapped_old_cluster,

    neighbor_fegenie_HMMs =
      fegenie_hmms,

    neighbor_number_of_hemes =
      heme_count,

    neighbor_globdb_cog =
      cog_for_plot,

    neighbor_globdb_product =
      product_for_plot,

    source_group
  )


## ================================================================== ##
## Combine
## ================================================================== ##

combined <- bind_rows(
  mg_plot,
  old_plot
) %>%
  arrange(
    plot_order,
    plot_start_bp
  )


if (
  n_distinct(
    combined$plot_region_id
  ) != 14
) {

  stop(
    "Expected 14 total Module-10 regions."
  )
}


if (
  n_distinct(
    combined$track_label
  ) != 14
) {

  stop(
    "Expected 14 unique track labels."
  )
}


write_tsv(
  combined,
  output_table
)


## ================================================================== ##
## Contig availability
##
## Keep the original validated GlobDB contig-QC information.
##
## For the MAGs:
##   - if the extracted neighborhood reaches ±20 genes OR ±20 kb,
##     that side is treated as extending beyond the displayed window;
##   - otherwise the observed edge is retained as the available
##     boundary.
##
## This keeps the visual convention of the old Module-10 figure
## without inventing additional sequence outside the observed data.
## ================================================================== ##

old_qc <- read_tsv(
  old_contig_qc,
  show_col_types = FALSE
) %>%
  filter(
    figure == "Module_10"
  ) %>%
  transmute(

    figure =
      "Module_10_with_MAGs",

    track_label =
      as.character(
        track_label
      ),

    plot_left_available_bp =
      as.numeric(
        plot_left_available_bp
      ),

    plot_right_available_bp =
      as.numeric(
        plot_right_available_bp
      )
  )


mag_qc <- mg %>%
  mutate(

    focal_midpoint =
      (
        as.numeric(focal_start) +
        as.numeric(focal_end)
      ) / 2,


    oriented_gene_offset =
      as.numeric(
        oriented_gene_offset
      ),


    oriented_midpoint_offset_bp =
      as.numeric(
        oriented_midpoint_offset_bp
      ),


    track_label = case_when(

      focal_protein_id ==
        "contig_25_317" ~

        "Methylumidiphilus [barcode10 SemiBin_1]",

      focal_protein_id ==
        "contig_26_118" ~

        "Methylumidiphilus [barcode10 SemiBin_8]"
    )
  ) %>%
  group_by(
    focal_protein_id,
    track_label
  ) %>%
  summarise(

    left_gene_reach =
      abs(
        min(
          oriented_gene_offset,
          na.rm = TRUE
        )
      ),

    right_gene_reach =
      max(
        oriented_gene_offset,
        na.rm = TRUE
      ),

    left_bp_reach =
      abs(
        min(
          oriented_midpoint_offset_bp,
          na.rm = TRUE
        )
      ),

    right_bp_reach =
      max(
        oriented_midpoint_offset_bp,
        na.rm = TRUE
      ),

    observed_left =
      abs(
        min(
          oriented_midpoint_offset_bp,
          na.rm = TRUE
        )
      ),

    observed_right =
      max(
        oriented_midpoint_offset_bp,
        na.rm = TRUE
      ),

    .groups = "drop"
  ) %>%
  mutate(

    left_complete =
      left_gene_reach >= 20 |
      left_bp_reach >= 20000,

    right_complete =
      right_gene_reach >= 20 |
      right_bp_reach >= 20000,


    ## Large value means the contig line will simply be clipped to the
    ## common plotting window, with no boundary mark visible.

    plot_left_available_bp =
      if_else(
        left_complete,
        1e9,
        observed_left
      ),

    plot_right_available_bp =
      if_else(
        right_complete,
        1e9,
        observed_right
      ),


    figure =
      "Module_10_with_MAGs"
  ) %>%
  select(
    figure,
    track_label,
    plot_left_available_bp,
    plot_right_available_bp
  )


combined_qc <- bind_rows(
  mag_qc,
  old_qc
)


if (
  nrow(
    combined_qc
  ) != 14
) {

  stop(
    "Expected 14 contig-QC tracks; found ",
    nrow(combined_qc),
    "."
  )
}


write_tsv(
  combined_qc,
  output_contig_qc
)


## ================================================================== ##
## Report
## ================================================================== ##

message("")
message("Stage 80A complete.")
message("")

message(
  "MAG regions:    ",
  n_distinct(
    mg_plot$plot_region_id
  )
)

message(
  "GlobDB regions: ",
  n_distinct(
    old_plot$plot_region_id
  )
)

message(
  "Total regions:  ",
  n_distinct(
    combined$plot_region_id
  )
)

message(
  "Genes drawn:    ",
  nrow(combined)
)

message("")
message("Combined plotting table:")
message("  ", output_table)

message("")
message("Combined contig QC:")
message("  ", output_contig_qc)

message("")
message("MAG old-cluster assignments present:")

print(
  mg_plot %>%
    filter(
      neighbor_mmseq_cluster != ""
    ) %>%
    distinct(
      protein_id,
      neighbor_mmseq_cluster
    ) %>%
    arrange(
      neighbor_mmseq_cluster,
      protein_id
    )
)
