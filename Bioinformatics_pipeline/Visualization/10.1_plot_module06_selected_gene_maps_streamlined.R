#!/usr/bin/env Rscript

## ================================================================== ##
## STAGE 16B - MODULE 06 / CLUSTER_00063 STACKED GENE MAP
##
## Focal family:
##   Cluster_00063 = MCA0421 / c553O-like multiheme cytochrome
##
## Figure:
##   12 representative genomes
##   5 Methylococcus
##   5 Methylobacter lineage
##   2 additional lineages
##
## Highlights:
##   - Cluster_00063
##   - other Module_06 clusters
##   - FeGenie-positive proteins
##   - other heme proteins
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

workflow_dir <- path.expand(
  "~/methanotrophs/methanotroph_project/jeppe/genome_analysis_workflow"
)

input_file <- file.path(
  workflow_dir,
  "16_visualization",
  "16A3_module06_representative_selection",
  "Module_06_selected_gene_map_data.tsv"
)

observability_file <- file.path(
  workflow_dir,
  "12_neighborhood_extraction",
  "focal_context_observability.tsv"
)

output_dir <- file.path(
  workflow_dir,
  "16_visualization",
  "16B_module06_gene_maps"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

output_pdf <- file.path(
  output_dir,
  "Module_06_selected_gene_map.pdf"
)

output_png <- file.path(
  output_dir,
  "Module_06_selected_gene_map.png"
)

output_labels <- file.path(
  output_dir,
  "Module_06_selected_gene_map_labels.tsv"
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
      paste(
        missing,
        collapse = ", "
      )
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

  break_max <- floor(
    x_limit / step
  ) * step

  seq(
    -break_max,
    break_max,
    by = step
  )
}


clean_text <- function(x) {

  x <- as.character(x)

  x[
    is.na(x) |
      x == "" |
      x == "nan"
  ] <- NA_character_

  x
}


## ================================================================== ##
## Read data
## ================================================================== ##

message(
  "Reading Module-06 selected gene-map data..."
)

genes <- read_tsv(
  input_file,
  show_col_types = FALSE
)

obs <- read_tsv(
  observability_file,
  show_col_types = FALSE
)


require_columns(
  genes,
  c(
    "plot_region_id",
    "genome",
    "focal_protein_id",
    "protein_id",
    "taxonomy_display",
    "taxonomy_family",
    "taxonomy_genus",
    "taxonomy_species",
    "plot_start_bp",
    "plot_end_bp",
    "plot_strand",
    "display_class",
    "display_label",
    "neighbor_cluster",
    "neighbor_fegenie_HMMs",
    "neighbor_number_of_hemes",
    "selection_order"
  ),
  "Module-06 selected gene-map data"
)


if (
  n_distinct(
    genes$plot_region_id
  ) != 12
) {

  stop(
    "Expected 12 Module-06 focal regions, found ",
    n_distinct(
      genes$plot_region_id
    ),
    "."
  )
}


## ================================================================== ##
## Numeric fields
## ================================================================== ##

genes <- genes %>%
  mutate(
    plot_start_bp = as.numeric(
      plot_start_bp
    ),
    plot_end_bp = as.numeric(
      plot_end_bp
    ),
    selection_order = as.numeric(
      selection_order
    ),
    neighbor_number_of_hemes = as.numeric(
      neighbor_number_of_hemes
    ),
    gene_midpoint = (
      plot_start_bp +
        plot_end_bp
    ) / 2,
    forward = (
      plot_strand == "+"
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
## ================================================================== ##

genes <- genes %>%
  mutate(

    architecture_class = case_when(

      neighbor_cluster ==
        "Cluster_00063" ~
        "core_mca0421",

      neighbor_cluster ==
        "Cluster_00154" ~
        "module06_00154",

      neighbor_cluster ==
        "Cluster_00163" ~
        "module06_00163",

      display_class ==
        "module06_other" ~
        "module06_other",

      display_class ==
        "fegenie_other" ~
        "fegenie_other",

      display_class ==
        "findmehemes_other" ~
        "findmehemes_other",

      TRUE ~
        "background_gene"
    )
  )


## ================================================================== ##
## Labels
## ================================================================== ##

genes <- genes %>%
  mutate(

    neighbor_fegenie_HMMs = clean_text(
      neighbor_fegenie_HMMs
    ),

    display_label_plot = case_when(

      architecture_class ==
        "core_mca0421" ~
        "MCA0421-like",

      architecture_class ==
        "module06_00154" ~
        "Cluster_00154",

      architecture_class ==
        "module06_00163" ~
        "Cluster_00163",

      architecture_class ==
        "module06_other" &
        !is.na(
          neighbor_cluster
        ) ~
        neighbor_cluster,

      architecture_class ==
        "fegenie_other" &
        !is.na(
          neighbor_fegenie_HMMs
        ) ~
        neighbor_fegenie_HMMs,

      architecture_class ==
        "fegenie_other" ~
        "FeGenie",

      architecture_class ==
        "findmehemes_other" ~
        "Other heme protein",

      TRUE ~
        NA_character_
    )
  )


architecture_levels <- c(
  "core_mca0421",
  "module06_00154",
  "module06_00163",
  "module06_other",
  "fegenie_other",
  "findmehemes_other",
  "background_gene"
)


genes <- genes %>%
  mutate(
    architecture_class = factor(
      architecture_class,
      levels = architecture_levels
    )
  )


## ================================================================== ##
## Track order
## ================================================================== ##

genome_order <- genes %>%
  distinct(
    plot_region_id,
    genome,
    focal_protein_id,
    taxonomy_family,
    taxonomy_genus,
    taxonomy_species,
    taxonomy_display,
    selection_order
  ) %>%
  arrange(
    selection_order
  ) %>%
  mutate(
    y_track = rev(
      seq_len(
        n()
      )
    )
  )


genes <- genes %>%
  left_join(
    genome_order %>%
      select(
        plot_region_id,
        y_track
      ),
    by = "plot_region_id"
  )


## ================================================================== ##
## Labels above / below
## ================================================================== ##

label_genes <- genes %>%
  filter(
    !is.na(
      display_label_plot
    ),
    display_label_plot != ""
  ) %>%
  group_by(
    plot_region_id
  ) %>%
  arrange(
    gene_midpoint,
    .by_group = TRUE
  ) %>%
  mutate(
    label_side = if_else(
      row_number() %% 2 == 1,
      "above",
      "below"
    )
  ) %>%
  ungroup()


labels_above <- label_genes %>%
  filter(
    label_side == "above"
  )

labels_below <- label_genes %>%
  filter(
    label_side == "below"
  )


write_tsv(
  label_genes,
  output_labels
)


## ================================================================== ##
## Contig availability
##
## Since every focal locus was oriented so Cluster_00063 points right:
##
## plot left  = oriented upstream
## plot right = oriented downstream
## ================================================================== ##

obs06 <- obs %>%
  filter(
    focal_cluster ==
      "Cluster_00063"
  ) %>%
  select(
    genome,
    focal_protein_id,
    oriented_upstream_bp_available,
    oriented_downstream_bp_available
  ) %>%
  mutate(
    oriented_upstream_bp_available =
      as.numeric(
        oriented_upstream_bp_available
      ),
    oriented_downstream_bp_available =
      as.numeric(
        oriented_downstream_bp_available
      )
  )


contig_spans <- genome_order %>%
  left_join(
    obs06,
    by = c(
      "genome",
      "focal_protein_id"
    )
  )


if (
  nrow(
    contig_spans
  ) != 12
) {

  stop(
    "Expected 12 contig spans, found ",
    nrow(
      contig_spans
    ),
    "."
  )
}


## ================================================================== ##
## X-axis
## ================================================================== ##

x_max <- 20000

x_breaks <- zero_centred_breaks(
  x_max
)


contig_spans <- contig_spans %>%
  mutate(

    sequence_start_bp = pmax(
      -oriented_upstream_bp_available,
      -x_max
    ),

    sequence_end_bp = pmin(
      oriented_downstream_bp_available,
      x_max
    ),

    left_boundary_visible =
      oriented_upstream_bp_available <
      x_max,

    right_boundary_visible =
      oriented_downstream_bp_available <
      x_max
  )


## ================================================================== ##
## Colours
##
## Uses the same broad palette conventions as the previous maps.
## ================================================================== ##

architecture_colors <- c(

  "core_mca0421" =
    "#0072B2",

  "module06_00154" =
    "#E69F00",

  "module06_00163" =
    "#8C6BB1",

  "module06_other" =
    "#56B4E9",

  "fegenie_other" =
    "#CC79A7",

  "findmehemes_other" =
    "#73A857",

  "background_gene" =
    "#D9D9D9"
)


legend_breaks <- c(
  "core_mca0421",
  "module06_00154",
  "module06_00163",
  "module06_other",
  "fegenie_other",
  "findmehemes_other"
)


legend_labels <- c(

  "core_mca0421" =
    "MCA0421-like (Cluster_00063)",

  "module06_00154" =
    "Module 06: Cluster_00154",

  "module06_00163" =
    "Module 06: Cluster_00163",

  "module06_other" =
    "Other Module 06 protein",

  "fegenie_other" =
    "Other FeGenie protein",

  "findmehemes_other" =
    "Other heme protein"
)


## ================================================================== ##
## Plot
## ================================================================== ##

p <- ggplot() +

  ## Available assembled sequence
  geom_segment(
    data = contig_spans,
    aes(
      x = sequence_start_bp,
      xend = sequence_end_bp,
      y = y_track,
      yend = y_track
    ),
    inherit.aes = FALSE,
    colour = "grey78",
    linewidth = 0.45
  ) +

  ## Left contig boundary
  geom_point(
    data = contig_spans %>%
      filter(
        left_boundary_visible
      ),
    aes(
      x = sequence_start_bp,
      y = y_track
    ),
    inherit.aes = FALSE,
    shape = 124,
    size = 2.6,
    colour = "grey55"
  ) +

  ## Right contig boundary
  geom_point(
    data = contig_spans %>%
      filter(
        right_boundary_visible
      ),
    aes(
      x = sequence_end_bp,
      y = y_track
    ),
    inherit.aes = FALSE,
    shape = 124,
    size = 2.6,
    colour = "grey55"
  ) +

  ## Gene arrows
  geom_gene_arrow(
    data = genes,
    aes(
      xmin = plot_start_bp,
      xmax = plot_end_bp,
      y = y_track,
      fill = architecture_class,
      forward = forward
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
    colour = "grey25",
    size = 0.25
  ) +

  ## Labels above
  geom_text_repel(
    data = labels_above,
    aes(
      x = gene_midpoint,
      y = y_track,
      label = display_label_plot
    ),
    inherit.aes = FALSE,
    direction = "x",
    nudge_y = 0.32,
    box.padding = 0.10,
    point.padding = 0.05,
    min.segment.length = 0,
    segment.colour = "grey50",
    segment.size = 0.25,
    segment.alpha = 0.75,
    max.overlaps = Inf,
    seed = 42,
    size = 2.7,
    fontface = "italic"
  ) +

  ## Labels below
  geom_text_repel(
    data = labels_below,
    aes(
      x = gene_midpoint,
      y = y_track,
      label = display_label_plot
    ),
    inherit.aes = FALSE,
    direction = "x",
    nudge_y = -0.32,
    box.padding = 0.10,
    point.padding = 0.05,
    min.segment.length = 0,
    segment.colour = "grey50",
    segment.size = 0.25,
    segment.alpha = 0.75,
    max.overlaps = Inf,
    seed = 42,
    size = 2.7,
    fontface = "italic"
  ) +

  scale_fill_manual(
    values = architecture_colors,
    breaks = legend_breaks,
    labels = legend_labels,
    drop = TRUE,
    name = NULL
  ) +

  scale_x_continuous(
    limits = c(
      -x_max,
      x_max
    ),
    breaks = x_breaks,
    labels = kb_labels,
    expand = expansion(
      mult = c(
        0.01,
        0.01
      )
    )
  ) +

  scale_y_continuous(
    breaks =
      genome_order$y_track,
    labels =
      genome_order$taxonomy_display,
    expand = expansion(
      add = c(
        0.7,
        0.7
      )
    )
  ) +

  labs(
    x = paste0(
      "Position relative to focal ",
      "MCA0421-like protein (Cluster_00063)"
    ),
    y = NULL
  ) +

  coord_cartesian(
    clip = "off"
  ) +

  guides(
    fill = guide_legend(
      nrow = 2,
      byrow = TRUE,
      override.aes = list(
        colour = "grey25"
      )
    )
  ) +

  theme_classic(
    base_size = 11
  ) +

  theme(

    axis.line.y =
      element_blank(),

    axis.ticks.y =
      element_blank(),

    axis.text.y =
      element_text(
        size = 9.5,
        face = "italic",
        colour = "black",
        margin = margin(
          r = 8
        )
      ),

    axis.text.x =
      element_text(
        size = 9
      ),

    axis.title.x =
      element_text(
        margin = margin(
          t = 8
        )
      ),

    legend.position =
      "bottom",

    legend.box =
      "horizontal",

    legend.text =
      element_text(
        size = 9
      ),

    panel.grid =
      element_blank(),

    plot.margin =
      margin(
        t = 15,
        r = 15,
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
  width = 12,
  height = 9,
  units = "in"
)

ggsave(
  output_png,
  p,
  width = 12,
  height = 9,
  units = "in",
  dpi = 400
)


message("")
message(
  "Module 06 gene-map rendering complete."
)

message(
  "Focal regions: ",
  n_distinct(
    genes$plot_region_id
  )
)

message(
  "Genomes: ",
  n_distinct(
    genes$genome
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
    label_genes
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
