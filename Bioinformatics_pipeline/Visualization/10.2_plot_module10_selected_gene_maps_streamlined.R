#!/usr/bin/env Rscript

## ================================================================== ##
## STAGE 16B - MODULE 10 / Cluster_00050 Cyc2 GENE MAP
##
## Streamlined plotting script.
## - Preserves all actual genes.
## - Uses compact on-gene labels and descriptive legend labels.
## - Uses zero-centred 5-kb x-axis ticks.
## - Adds subtle contig-availability lines and true contig-boundary marks.
## - Uses the same arrow/label geometry as the Module-20/35 figures.
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
  "16A3_module10_representative_selection",
  "Module_10_selected_gene_map_data.tsv"
)

contig_qc_file <- file.path(
  workflow_dir,
  "16_visualization",
  "16C_cross_figure_qc",
  "selected_gene_map_contig_boundary_qc.tsv"
)

output_dir <- file.path(
  workflow_dir,
  "16_visualization",
  "16B_module10_gene_maps"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

output_pdf <- file.path(output_dir, "Module_10_selected_gene_map.pdf")
output_png <- file.path(output_dir, "Module_10_selected_gene_map.png")

## ================================================================== ##
## Helpers
## ================================================================== ##

require_columns <- function(df, required, table_name) {
  missing <- setdiff(required, names(df))
  if (length(missing) > 0) {
    stop(
      table_name,
      " is missing required columns: ",
      paste(missing, collapse = ", ")
    )
  }
}

kb_labels <- function(x) {
  ifelse(x == 0, "0", paste0(x / 1000, " kb"))
}

zero_centred_breaks <- function(x_limit, step = 5000) {
  break_max <- floor(x_limit / step) * step
  seq(-break_max, break_max, by = step)
}

build_contig_spans <- function(contig_qc, figure_name, track_table, x_limit) {
  qc <- contig_qc %>%
    filter(figure == figure_name) %>%
    transmute(
      track_label = as.character(track_label),
      left_available_bp = as.numeric(plot_left_available_bp),
      right_available_bp = as.numeric(plot_right_available_bp)
    )

  spans <- track_table %>%
    inner_join(qc, by = "track_label") %>%
    mutate(
      sequence_start_bp = pmax(-left_available_bp, -x_limit),
      sequence_end_bp = pmin(right_available_bp, x_limit),
      left_boundary_visible = left_available_bp < x_limit,
      right_boundary_visible = right_available_bp < x_limit
    )

  if (nrow(spans) != nrow(track_table)) {
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
## Read and validate data
## ================================================================== ##

message("Reading selected Module-10 gene-map data...")

genes <- read_tsv(input_file, show_col_types = FALSE)
contig_qc <- read_tsv(contig_qc_file, show_col_types = FALSE)

require_columns(
  genes,
  c(
    "genome", "protein_id",
    "plot_region_id", "plot_order", "track_label",
    "taxonomy_display",
    "plot_start_bp", "plot_end_bp", "plot_strand",
    "display_class",
    "local_family_id", "neighbor_mmseq_cluster",
    "neighbor_fegenie_HMMs", "neighbor_number_of_hemes",
    "neighbor_globdb_cog", "neighbor_globdb_product"
  ),
  "Module-10 gene-map table"
)

require_columns(
  contig_qc,
  c(
    "figure", "track_label",
    "plot_left_available_bp", "plot_right_available_bp"
  ),
  "Contig-boundary QC table"
)

if (n_distinct(genes$plot_region_id) != 12) {
  stop("Expected 12 focal regions, found ", n_distinct(genes$plot_region_id), ".")
}

if (n_distinct(genes$genome) != 12) {
  stop("Expected 12 genomes, found ", n_distinct(genes$genome), ".")
}

if (any(duplicated(genes[c("plot_region_id", "protein_id")]))) {
  stop("Duplicate plot_region_id + protein_id rows detected.")
}

genes <- genes %>%
  mutate(
    plot_start_bp = as.numeric(plot_start_bp),
    plot_end_bp = as.numeric(plot_end_bp),
    plot_order = as.integer(plot_order),
    gene_midpoint = (plot_start_bp + plot_end_bp) / 2,
    forward = plot_strand == "+"
  )

if (anyNA(genes$plot_start_bp) || anyNA(genes$plot_end_bp)) {
  stop("Missing/non-numeric plotting coordinates detected.")
}

## ================================================================== ##
## Biological classes
##
## C00064-A/B and C00222-A/B remain figure-level shorthand for local
## sequence variants; global MMseq cluster identities remain unchanged.
## ================================================================== ##

genes <- genes %>%
  mutate(
    architecture_class = case_when(
      neighbor_mmseq_cluster == "Cluster_00050" ~ "cyc2",
      neighbor_mmseq_cluster == "Cluster_00091" ~ "c00091",

      neighbor_mmseq_cluster == "Cluster_00064" &
        local_family_id == "M10_local_002" ~ "c00064_A",
      neighbor_mmseq_cluster == "Cluster_00064" &
        local_family_id == "M10_local_007" ~ "c00064_B",

      neighbor_mmseq_cluster == "Cluster_00069" ~ "c00069",

      neighbor_mmseq_cluster == "Cluster_00222" &
        local_family_id == "M10_local_031" ~ "c00222_A",
      neighbor_mmseq_cluster == "Cluster_00222" &
        local_family_id == "M10_local_077" ~ "c00222_B",

      local_family_id %in% c("M10_local_008", "M10_local_013") ~
        "accessory_cytc553",

      local_family_id == "M10_local_005" ~ "cytb_like",

      display_class == "fegenie_other" ~ "other_fegenie",
      display_class == "findmehemes_other" ~ "other_heme",

      TRUE ~ "background"
    ),

    ## Compact on-gene labels; biological detail belongs in the legend.
    display_label_plot = case_when(
      architecture_class == "cyc2" ~ "Cyc2",
      architecture_class == "c00091" ~ "C00091",
      architecture_class == "c00064_A" ~ "C00064-A",
      architecture_class == "c00064_B" ~ "C00064-B",
      architecture_class == "c00069" ~ "C00069",
      architecture_class == "c00222_A" ~ "C00222-A",
      architecture_class == "c00222_B" ~ "C00222-B",
      architecture_class == "accessory_cytc553" ~ "Accessory CytC553-like",
      architecture_class == "cytb_like" ~ "CytB-like",
      architecture_class == "other_fegenie" ~ "Other FeGenie hit",
      architecture_class == "other_heme" ~ "Other heme protein",
      TRUE ~ ""
    ),

    label_side = case_when(
      architecture_class %in% c(
        "cyc2",
        "c00091",
        "c00064_A",
        "c00064_B",
        "c00069",
        "c00222_A",
        "c00222_B"
      ) ~ "above",

      architecture_class %in% c(
        "accessory_cytc553",
        "cytb_like",
        "other_fegenie",
        "other_heme"
      ) ~ "below",

      TRUE ~ "none"
    )
  )

## ================================================================== ##
## Track and sparse-label preparation
## ================================================================== ##

track_table <- genes %>%
  distinct(plot_order, track_label) %>%
  arrange(plot_order)

if (nrow(track_table) != 12) {
  stop("Expected exactly 12 unique plotting tracks.")
}

track_levels <- track_table$track_label

genes <- genes %>%
  mutate(
    track_factor = factor(
      track_label,
      levels = rev(track_levels)
    )
  )

## If the same displayed label occurs more than once in a locus,
## label only the copy closest to the focal Cyc2.
labels <- genes %>%
  filter(display_label_plot != "", label_side != "none") %>%
  group_by(plot_region_id, display_label_plot) %>%
  arrange(abs(gene_midpoint), gene_midpoint, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup()

labels_above <- labels %>% filter(label_side == "above")
labels_below <- labels %>% filter(label_side == "below")

## ================================================================== ##
## X scale and contig spans
## ================================================================== ##

x_limit <- ceiling(
  (
    max(abs(c(genes$plot_start_bp, genes$plot_end_bp)), na.rm = TRUE) +
      500
  ) / 1000
) * 1000

x_breaks <- zero_centred_breaks(x_limit)

contig_tracks <- genes %>%
  distinct(track_label, track_factor)

contig_spans <- build_contig_spans(
  contig_qc = contig_qc,
  figure_name = "Module_10",
  track_table = contig_tracks,
  x_limit = x_limit
)

## ================================================================== ##
## Visual scales
## ================================================================== ##

architecture_colors <- c(
  "cyc2" = "#0072B2",
  "c00091" = "#E69F00",
  "c00064_A" = "#009E73",
  "c00064_B" = "#56B4E9",
  "c00069" = "#D8A600",
  "c00222_A" = "#8C6BB1",
  "c00222_B" = "#CC79A7",
  "accessory_cytc553" = "#B8B0A5",
  "cytb_like" = "#D55E00",
  "other_fegenie" = "#C77CFF",
  "other_heme" = "#73A857",
  "background" = "#D9D9D9"
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
  "cyc2" = "Cyc2 (Cluster_00050)",
  "c00091" = "C00091 (CytC553-like)",
  "c00064_A" = "C00064-A (CytC553-like)",
  "c00064_B" = "C00064-B (CytC553-like)",
  "c00069" = "C00069 (1-heme)",
  "c00222_A" = "C00222-A (CytC553-like)",
  "c00222_B" = "C00222-B (CytC553-like)",
  "accessory_cytc553" = "Accessory CytC553-like",
  "cytb_like" = "CytB-like",
  "other_fegenie" = "Other FeGenie hit",
  "other_heme" = "Other heme protein"
)

## ================================================================== ##
## Plot
## ================================================================== ##

p <- ggplot() +

  geom_segment(
    data = contig_spans,
    aes(
      x = sequence_start_bp,
      xend = sequence_end_bp,
      y = track_factor,
      yend = track_factor
    ),
    inherit.aes = FALSE,
    colour = "grey78",
    linewidth = 0.45
  ) +

  geom_point(
    data = contig_spans %>% filter(left_boundary_visible),
    aes(x = sequence_start_bp, y = track_factor),
    inherit.aes = FALSE,
    shape = 124,
    size = 2.6,
    colour = "grey55"
  ) +

  geom_point(
    data = contig_spans %>% filter(right_boundary_visible),
    aes(x = sequence_end_bp, y = track_factor),
    inherit.aes = FALSE,
    shape = 124,
    size = 2.6,
    colour = "grey55"
  ) +

  geom_gene_arrow(
    data = genes,
    aes(
      xmin = plot_start_bp,
      xmax = plot_end_bp,
      y = track_factor,
      fill = architecture_class,
      forward = forward
    ),
    arrowhead_height = grid::unit(3, "mm"),
    arrowhead_width = grid::unit(1.5, "mm"),
    arrow_body_height = grid::unit(4, "mm"),
    colour = "grey25",
    size = 0.25
  ) +

  geom_text_repel(
    data = labels_above,
    aes(x = gene_midpoint, y = track_factor, label = display_label_plot),
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
    size = 3,
    fontface = "italic",
    colour = "grey20"
  ) +

  geom_text_repel(
    data = labels_below,
    aes(x = gene_midpoint, y = track_factor, label = display_label_plot),
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
    size = 3,
    fontface = "italic",
    colour = "grey20"
  ) +

  scale_fill_manual(
    values = architecture_colors,
    breaks = legend_breaks,
    labels = legend_labels,
    drop = TRUE,
    name = NULL
  ) +

  scale_x_continuous(
    limits = c(-x_limit, x_limit),
    breaks = x_breaks,
    labels = kb_labels,
    expand = expansion(mult = c(0.01, 0.01))
  ) +

  labs(
    x = "Position relative to focal Cyc2 (Cluster_00050)",
    y = NULL
  ) +

  coord_cartesian(clip = "off") +

  guides(
    fill = guide_legend(
      nrow = 2,
      byrow = TRUE,
      override.aes = list(colour = "grey25")
    )
  ) +

  theme_classic(base_size = 11) +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(
      face = "italic",
      size = 10,
      colour = "grey15",
      margin = margin(r = 8)
    ),
    axis.text.x = element_text(size = 9, colour = "grey20"),
    axis.title.x = element_text(margin = margin(t = 8)),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.text = element_text(size = 9),
    legend.key.width = grid::unit(1.1, "cm"),
    legend.spacing.x = grid::unit(0.20, "cm"),
    panel.grid = element_blank(),
    plot.margin = margin(t = 12, r = 20, b = 10, l = 10)
  )

## ================================================================== ##
## Save
## ================================================================== ##

ggsave(
  output_pdf,
  p,
  width = 14,
  height = 10.5,
  units = "in",
  device = cairo_pdf
)

ggsave(
  output_png,
  p,
  width = 14,
  height = 10.5,
  units = "in",
  dpi = 300
)

message("")
message("Module 10 gene-map rendering complete.")
message("Regions: ", n_distinct(genes$plot_region_id))
message("Genomes: ", n_distinct(genes$genome))
message("Genes drawn: ", nrow(genes))
message("Labels drawn: ", nrow(labels))
message("PDF: ", output_pdf)
message("PNG: ", output_png)
