#!/usr/bin/env Rscript

## ================================================================== ##
## STAGE 80B - MODULE 10 / Cluster_00050 Cyc2 + METAGENOME MAGs
##
## Based directly on the validated Stage-16B Module-10 plotting script.
##
## Changes:
##   - 12 -> 14 focal regions
##   - adds two Methylumidiphilus MAG neighborhoods
##   - separates Metagenome MAGs / GlobDB genomes into blocks
##
## All original gene-map aesthetics are otherwise retained.
## ================================================================== ##

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(gggenes)
  library(ggrepel)
})


## ================================================================== ##
## Paths
## ================================================================== ##

project <- path.expand(
  "~/methanotrophs/methanotroph_project/jeppe/metagenome"
)

input_dir <- file.path(
  project,
  "comparative_analysis",
  "Module_10_with_MAGs"
)

input_file <- file.path(
  input_dir,
  "Module_10_selected_gene_map_data_with_MAGs.tsv"
)

contig_qc_file <- file.path(
  input_dir,
  "Module_10_with_MAGs_contig_boundary_qc.tsv"
)

output_pdf <- file.path(
  input_dir,
  "Module_10_selected_gene_map_with_MAGs.pdf"
)

output_png <- file.path(
  input_dir,
  "Module_10_selected_gene_map_with_MAGs.png"
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


kb_labels <- function(x) {

  ifelse(
    x == 0,
    "0",
    paste0(
      x / 1000,
      " kb"
    )
  )
}


zero_centred_breaks <- function(
  x_limit,
  step = 5000
) {

  break_max <-
    floor(
      x_limit /
        step
    ) *
    step

  seq(
    -break_max,
    break_max,
    by = step
  )
}


build_contig_spans <- function(
  contig_qc,
  figure_name,
  track_table,
  x_limit
) {

  qc <- contig_qc %>%
    filter(
      figure ==
        figure_name
    ) %>%
    transmute(

      track_label =
        as.character(
          track_label
        ),

      left_available_bp =
        as.numeric(
          plot_left_available_bp
        ),

      right_available_bp =
        as.numeric(
          plot_right_available_bp
        )
    )


  spans <- track_table %>%
    inner_join(
      qc,
      by = "track_label"
    ) %>%
    mutate(

      sequence_start_bp =
        pmax(
          -left_available_bp,
          -x_limit
        ),

      sequence_end_bp =
        pmin(
          right_available_bp,
          x_limit
        ),

      left_boundary_visible =
        left_available_bp <
        x_limit,

      right_boundary_visible =
        right_available_bp <
        x_limit
    )


  if (
    nrow(spans) !=
    nrow(track_table)
  ) {

    stop(
      "Contig-QC matching failed for ",
      figure_name,
      ": expected ",
      nrow(track_table),
      " tracks, matched ",
      nrow(spans),
      "."
    )
  }


  spans
}


## ================================================================== ##
## Read and validate
## ================================================================== ##

message(
  "Reading combined Module-10 + MAG gene-map data..."
)


genes <- read_tsv(
  input_file,
  show_col_types = FALSE
)


contig_qc <- read_tsv(
  contig_qc_file,
  show_col_types = FALSE
)


require_columns(
  genes,
  c(
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
    "neighbor_globdb_product",
    "source_group"
  ),
  "Module-10 + MAG gene-map table"
)


require_columns(
  contig_qc,
  c(
    "figure",
    "track_label",
    "plot_left_available_bp",
    "plot_right_available_bp"
  ),
  "Contig-boundary QC table"
)


if (
  n_distinct(
    genes$plot_region_id
  ) != 14
) {

  stop(
    "Expected 14 focal regions, found ",
    n_distinct(
      genes$plot_region_id
    ),
    "."
  )
}


if (
  n_distinct(
    genes$genome
  ) != 14
) {

  stop(
    "Expected 14 genomes, found ",
    n_distinct(
      genes$genome
    ),
    "."
  )
}


if (
  any(
    duplicated(
      genes[
        c(
          "plot_region_id",
          "protein_id"
        )
      ]
    )
  )
) {

  stop(
    "Duplicate plot_region_id + protein_id rows detected."
  )
}


genes <- genes %>%
  mutate(

    plot_start_bp =
      as.numeric(
        plot_start_bp
      ),

    plot_end_bp =
      as.numeric(
        plot_end_bp
      ),

    plot_order =
      as.integer(
        plot_order
      ),

    gene_midpoint =
      (
        plot_start_bp +
        plot_end_bp
      ) / 2,

    forward =
      plot_strand ==
      "+",

    source_group =
      factor(
        source_group,
        levels = c(
          "Metagenome MAGs",
          "GlobDB genomes"
        )
      )
  )


if (
  anyNA(
    genes$plot_start_bp
  ) ||
  anyNA(
    genes$plot_end_bp
  )
) {

  stop(
    "Missing/non-numeric plotting coordinates detected."
  )
}


## ================================================================== ##
## Biological classes
##
## UNCHANGED from original Stage 16B.
## ================================================================== ##

genes <- genes %>%
  mutate(

    architecture_class = case_when(

      neighbor_mmseq_cluster ==
        "Cluster_00050" ~
        "cyc2",

      neighbor_mmseq_cluster ==
        "Cluster_00091" ~
        "c00091",

      ## ---------------------------------------------------------- ##
      ## Cluster_00064
      ##
      ## Original GlobDB classification:
      ##   M10_local_002 = C00064-A
      ##   M10_local_007 = C00064-B
      ##
      ## Both new MAG Cluster_00064 proteins have their strongest
      ## old-family hit to MOTU40_079870_000000000028_8, which is
      ## M10_local_002. Therefore they are displayed as C00064-A.
      ## ---------------------------------------------------------- ##

      neighbor_mmseq_cluster ==
        "Cluster_00064" &
        (
          local_family_id ==
            "M10_local_002" |
          source_group ==
            "Metagenome MAGs"
        ) ~
        "c00064_A",

      neighbor_mmseq_cluster ==
        "Cluster_00064" &
        local_family_id ==
        "M10_local_007" ~
        "c00064_B",

      neighbor_mmseq_cluster ==
        "Cluster_00069" ~
        "c00069",

      neighbor_mmseq_cluster ==
        "Cluster_00222" &
        local_family_id ==
        "M10_local_031" ~
        "c00222_A",

      neighbor_mmseq_cluster ==
        "Cluster_00222" &
        local_family_id ==
        "M10_local_077" ~
        "c00222_B",

      ## ---------------------------------------------------------- ##
      ## Accessory CytC553-like
      ##
      ## The original Module-10 figure defines M10_local_008/013 as
      ## accessory CytC553-like proteins.
      ##
      ## Cluster_00141 is represented in the selected GlobDB figure
      ## by GCF_006175985_000000000019_54, assigned to M10_local_008.
      ## The two MAG proteins mapping to Cluster_00141 therefore use
      ## the same display class.
      ## ---------------------------------------------------------- ##

      local_family_id %in%
        c(
          "M10_local_008",
          "M10_local_013"
        ) |
        neighbor_mmseq_cluster ==
          "Cluster_00141" ~
        "accessory_cytc553",

      local_family_id ==
        "M10_local_005" ~
        "cytb_like",

      display_class ==
        "fegenie_other" ~
        "other_fegenie",

      display_class ==
        "findmehemes_other" ~
        "other_heme",

      TRUE ~
        "background"
    ),


    display_label_plot = case_when(

      architecture_class ==
        "cyc2" ~
        "Cyc2",

      architecture_class ==
        "c00091" ~
        "C00091",

      architecture_class ==
        "c00064_A" ~
        "C00064-A",

      architecture_class ==
        "c00064_B" ~
        "C00064-B",

      architecture_class ==
        "c00069" ~
        "C00069",

      architecture_class ==
        "c00222_A" ~
        "C00222-A",

      architecture_class ==
        "c00222_B" ~
        "C00222-B",

      architecture_class ==
        "accessory_cytc553" ~
        "Accessory CytC553-like",

      architecture_class ==
        "cytb_like" ~
        "CytB-like",

      architecture_class ==
        "other_fegenie" ~
        "Other FeGenie hit",

      architecture_class ==
        "other_heme" ~
        "Other heme protein",

      TRUE ~
        ""
    ),


    label_side = case_when(

      architecture_class %in%
        c(
          "cyc2",
          "c00091",
          "c00064_A",
          "c00064_B",
          "c00069",
          "c00222_A",
          "c00222_B"
        ) ~
        "above",

      architecture_class %in%
        c(
          "accessory_cytc553",
          "cytb_like",
          "other_fegenie",
          "other_heme"
        ) ~
        "below",

      TRUE ~
        "none"
    )
  )


## ================================================================== ##
## Track order
##
## MAGs first, then preserve the original 12 GlobDB orders.
## ================================================================== ##

track_table <- genes %>%
  distinct(
    source_group,
    plot_order,
    track_label
  ) %>%
  arrange(
    source_group,
    plot_order
  )


if (
  nrow(
    track_table
  ) != 14
) {

  stop(
    "Expected exactly 14 unique plotting tracks."
  )
}


track_levels <-
  track_table$track_label


genes <- genes %>%
  mutate(

    track_factor =
      factor(
        track_label,
        levels =
          rev(
            track_levels
          )
      )
  )


## ================================================================== ##
## Sparse labels
## ================================================================== ##

labels <- genes %>%
  filter(
    display_label_plot != "",
    label_side != "none"
  ) %>%
  group_by(
    plot_region_id,
    display_label_plot
  ) %>%
  arrange(
    abs(
      gene_midpoint
    ),
    gene_midpoint,
    .by_group = TRUE
  ) %>%
  slice(1) %>%
  ungroup()


labels_above <- labels %>%
  filter(
    label_side ==
      "above"
  )


labels_below <- labels %>%
  filter(
    label_side ==
      "below"
  )


## ================================================================== ##
## X scale and contig spans
## ================================================================== ##

x_limit <- ceiling(
  (
    max(
      abs(
        c(
          genes$plot_start_bp,
          genes$plot_end_bp
        )
      ),
      na.rm = TRUE
    ) +
      500
  ) /
    1000
) *
  1000


x_breaks <-
  zero_centred_breaks(
    x_limit
  )


contig_tracks <- genes %>%
  distinct(
    source_group,
    track_label,
    track_factor
  )


contig_spans <- build_contig_spans(
  contig_qc =
    contig_qc,
  figure_name =
    "Module_10_with_MAGs",
  track_table =
    contig_tracks,
  x_limit =
    x_limit
)


## ================================================================== ##
## Visual scales
##
## UNCHANGED.
## ================================================================== ##

architecture_colors <- c(

  "cyc2" =
    "#0072B2",

  "c00091" =
    "#E69F00",

  "c00064_A" =
    "#009E73",

  "c00064_B" =
    "#56B4E9",

  "c00069" =
    "#D8A600",

  "c00222_A" =
    "#8C6BB1",

  "c00222_B" =
    "#CC79A7",

  "accessory_cytc553" =
    "#B8B0A5",

  "cytb_like" =
    "#D55E00",

  "other_fegenie" =
    "#C77CFF",

  "other_heme" =
    "#73A857",

  "background" =
    "#D9D9D9"
)


legend_breaks <- c(
  "cyc2",
  "c00091",
  "c00064_A",
  "c00064_B",
  "c00069",
  "c00222_A",
  "c00222_B",
  "accessory_cytc553",
  "cytb_like",
  "other_fegenie",
  "other_heme"
)


legend_labels <- c(

  "cyc2" =
    "Cyc2 (Cluster_00050)",

  "c00091" =
    "C00091 (CytC553-like)",

  "c00064_A" =
    "C00064-A (CytC553-like)",

  "c00064_B" =
    "C00064-B (CytC553-like)",

  "c00069" =
    "C00069 (1-heme)",

  "c00222_A" =
    "C00222-A (CytC553-like)",

  "c00222_B" =
    "C00222-B (CytC553-like)",

  "accessory_cytc553" =
    "Accessory CytC553-like",

  "cytb_like" =
    "CytB-like",

  "other_fegenie" =
    "Other FeGenie hit",

  "other_heme" =
    "Other heme protein"
)


## ================================================================== ##
## Plot
## ================================================================== ##

p <- ggplot() +

  geom_segment(
    data =
      contig_spans,
    aes(
      x =
        sequence_start_bp,
      xend =
        sequence_end_bp,
      y =
        track_factor,
      yend =
        track_factor
    ),
    inherit.aes = FALSE,
    colour = "grey78",
    linewidth = 0.45
  ) +

  geom_point(
    data =
      contig_spans %>%
      filter(
        left_boundary_visible
      ),
    aes(
      x =
        sequence_start_bp,
      y =
        track_factor
    ),
    inherit.aes = FALSE,
    shape = 124,
    size = 2.6,
    colour = "grey55"
  ) +

  geom_point(
    data =
      contig_spans %>%
      filter(
        right_boundary_visible
      ),
    aes(
      x =
        sequence_end_bp,
      y =
        track_factor
    ),
    inherit.aes = FALSE,
    shape = 124,
    size = 2.6,
    colour = "grey55"
  ) +

  geom_gene_arrow(
    data =
      genes,
    aes(
      xmin =
        plot_start_bp,
      xmax =
        plot_end_bp,
      y =
        track_factor,
      fill =
        architecture_class,
      forward =
        forward
    ),
    arrowhead_height =
      grid::unit(
        3,
        "mm"
      ),
    arrowhead_width =
      grid::unit(
        1.5,
        "mm"
      ),
    arrow_body_height =
      grid::unit(
        4,
        "mm"
      ),
    colour =
      "grey25",
    size =
      0.25
  ) +

  geom_text_repel(
    data =
      labels_above,
    aes(
      x =
        gene_midpoint,
      y =
        track_factor,
      label =
        display_label_plot
    ),
    inherit.aes =
      FALSE,
    direction =
      "x",
    nudge_y =
      0.32,
    box.padding =
      0.10,
    point.padding =
      0.05,
    min.segment.length =
      0,
    segment.colour =
      "grey50",
    segment.size =
      0.25,
    segment.alpha =
      0.75,
    max.overlaps =
      Inf,
    seed =
      42,
    size =
      3,
    fontface =
      "italic",
    colour =
      "grey20"
  ) +

  geom_text_repel(
    data =
      labels_below,
    aes(
      x =
        gene_midpoint,
      y =
        track_factor,
      label =
        display_label_plot
    ),
    inherit.aes =
      FALSE,
    direction =
      "x",
    nudge_y =
      -0.32,
    box.padding =
      0.10,
    point.padding =
      0.05,
    min.segment.length =
      0,
    segment.colour =
      "grey50",
    segment.size =
      0.25,
    segment.alpha =
      0.75,
    max.overlaps =
      Inf,
    seed =
      42,
    size =
      3,
    fontface =
      "italic",
    colour =
      "grey20"
  ) +

  ## -------------------------------------------------------------- ##
  ## Only substantive presentation change:
  ## MAG block above GlobDB block.
  ## -------------------------------------------------------------- ##

  facet_grid(
    source_group ~ .,
    scales =
      "free_y",
    space =
      "free_y",
    switch =
      "y"
  ) +

  scale_fill_manual(
    values =
      architecture_colors,
    breaks =
      legend_breaks,
    labels =
      legend_labels,
    drop =
      TRUE,
    name =
      NULL
  ) +

  scale_x_continuous(
    limits =
      c(
        -x_limit,
        x_limit
      ),
    breaks =
      x_breaks,
    labels =
      kb_labels,
    expand =
      expansion(
        mult =
          c(
            0.01,
            0.01
          )
      )
  ) +

  labs(
    x =
      "Position relative to focal Cyc2 (Cluster_00050)",
    y =
      NULL
  ) +

  coord_cartesian(
    clip =
      "off"
  ) +

  guides(
    fill =
      guide_legend(
        nrow =
          2,
        byrow =
          TRUE,
        override.aes =
          list(
            colour =
              "grey25"
          )
      )
  ) +

  theme_classic(
    base_size =
      11
  ) +

  theme(

    axis.line.y =
      element_blank(),

    axis.ticks.y =
      element_blank(),

    axis.text.y =
      element_text(
        face =
          "italic",
        size =
          10,
        colour =
          "grey15",
        margin =
          margin(
            r = 8
          )
      ),

    axis.text.x =
      element_text(
        size =
          9,
        colour =
          "grey20"
      ),

    axis.title.x =
      element_text(
        margin =
          margin(
            t = 8
          )
      ),

    ## Cluster_00121-style group block

    strip.background =
      element_rect(
        fill =
          "grey95",
        colour =
          "grey45",
        linewidth =
          0.4
      ),

    strip.text.y.left =
      element_text(
        angle =
          0,
        face =
          "bold",
        size =
          10,
        colour =
          "grey15"
      ),

    legend.position =
      "bottom",

    legend.box =
      "horizontal",

    legend.text =
      element_text(
        size =
          9
      ),

    legend.key.width =
      grid::unit(
        1.1,
        "cm"
      ),

    legend.spacing.x =
      grid::unit(
        0.20,
        "cm"
      ),

    panel.grid =
      element_blank(),

    plot.margin =
      margin(
        t = 12,
        r = 20,
        b = 10,
        l = 10
      )
  )


## ================================================================== ##
## Save
## ================================================================== ##

ggsave(
  output_pdf,
  p,
  width =
    14,
  height =
    12,
  units =
    "in",
  device =
    cairo_pdf
)


ggsave(
  output_png,
  p,
  width =
    14,
  height =
    12,
  units =
    "in",
  dpi =
    300
)


## ================================================================== ##
## Report
## ================================================================== ##

message("")
message(
  "Module 10 + MAG gene-map rendering complete."
)

message(
  "MAG regions: ",
  n_distinct(
    genes$plot_region_id[
      genes$source_group ==
        "Metagenome MAGs"
    ]
  )
)

message(
  "GlobDB regions: ",
  n_distinct(
    genes$plot_region_id[
      genes$source_group ==
        "GlobDB genomes"
    ]
  )
)

message(
  "Total regions: ",
  n_distinct(
    genes$plot_region_id
  )
)

message(
  "Genes drawn: ",
  nrow(
    genes
  )
)

message(
  "Labels drawn: ",
  nrow(
    labels
  )
)

message(
  "PDF: ",
  output_pdf
)

message(
  "PNG: ",
  output_png
)
