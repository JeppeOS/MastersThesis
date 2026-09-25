#!/usr/bin/env Rscript

## ================================================================== ##
## STAGE 16B - MODULE 20 STACKED GENE MAP
##
## Streamlined plotting script.
## - Preserves the established Module-20 biological interpretation.
## - Uses zero-centred 5-kb x-axis ticks.
## - Adds subtle contig-availability lines and true contig-boundary marks.
## - Keeps every actual gene visible.
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
  "16A3_module20_representative_selection",
  "Module_20_selected_gene_map_data.tsv"
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
  "16B_module20_gene_maps"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

output_pdf <- file.path(output_dir, "Module_20_selected_gene_map.pdf")
output_png <- file.path(output_dir, "Module_20_selected_gene_map.png")

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

message("Reading selected Module-20 architecture data...")

genes <- read_tsv(input_file, show_col_types = FALSE)
contig_qc <- read_tsv(contig_qc_file, show_col_types = FALSE)

require_columns(
  genes,
  c(
    "plot_region_id",
    "genome", "protein_id",
    "taxonomy_display", "taxonomy_family", "taxonomy_genus", "taxonomy_species",
    "plot_start_bp", "plot_end_bp", "plot_strand",
    "display_class", "display_label",
    "local_family_id", "neighbor_mmseq_cluster",
    "focal_call", "is_mtrb", "selection_order"
  ),
  "Module-20 gene-map table"
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
  stop(
    "Expected 12 Module-20 focal regions, found ",
    n_distinct(genes$plot_region_id),
    "."
  )
}

genes <- genes %>%
  mutate(
    plot_start_bp = as.numeric(plot_start_bp),
    plot_end_bp = as.numeric(plot_end_bp),
    selection_order = as.numeric(selection_order),
    is_mtrb = as.numeric(is_mtrb),
    gene_midpoint = (plot_start_bp + plot_end_bp) / 2,
    forward = plot_strand == "+"
  )

if (anyNA(genes$plot_start_bp) || anyNA(genes$plot_end_bp)) {
  stop("Missing/non-numeric plotting coordinates detected.")
}

## ================================================================== ##
## Biological classes and labels
## ================================================================== ##

genes <- genes %>%
  mutate(
    architecture_class = case_when(
      local_family_id == "M20_local_001" & focal_call == "MtoA" ~
        "core_mtoa",
      local_family_id == "M20_local_001" & focal_call == "MtrA" ~
        "core_mtra",
      local_family_id == "M20_local_001" &
        focal_call == "FeGenie_negative" ~
        "core_fegenie_negative",

      is_mtrb == 1 ~ "core_mtrb",
      local_family_id == "M20_local_002" ~ "core_cytc552",
      local_family_id == "M20_local_003" ~ "recurrent_unknown",

      local_family_id %in% c(
        "M20_local_007",
        "M20_local_008",
        "M20_local_013",
        "M20_local_014"
      ) ~ "accessory_redox_block",

      local_family_id == "M20_local_015" ~
        "one_heme_paralogue_block",

      display_class == "fegenie_other" ~ "fegenie_other",
      display_class == "findmehemes_other" ~ "findmehemes_other",
      display_class == "other_mmseq_family" ~ "other_mmseq_family",

      TRUE ~ "background_gene"
    ),

    display_label_plot = case_when(
      local_family_id == "M20_local_001" & focal_call == "MtoA" ~ "MtoA",
      local_family_id == "M20_local_001" & focal_call == "MtrA" ~ "MtrA",
      local_family_id == "M20_local_001" &
        focal_call == "FeGenie_negative" ~ "Cluster_00035 homolog",

      is_mtrb == 1 ~ "MtrB",
      local_family_id == "M20_local_002" ~ "CytC551/552-like",
      local_family_id == "M20_local_003" ~ "M20-U1",

      local_family_id == "M20_local_007" ~ "4-heme NapC-like",
      local_family_id == "M20_local_008" ~ "2-heme cyt.c",
      local_family_id == "M20_local_013" ~ "CytC553-like",
      local_family_id == "M20_local_014" ~ "Multicopper oxidase",

      local_family_id == "M20_local_015" ~ NA_character_,

      architecture_class == "findmehemes_other" ~ "Other heme protein",

      architecture_class %in% c(
        "fegenie_other",
        "other_mmseq_family"
      ) &
        !is.na(display_label) &
        display_label != "" ~ display_label,

      TRUE ~ NA_character_
    )
  )

## One shared text label for each local-015 paralogue block.
genes <- genes %>%
  group_by(plot_region_id) %>%
  arrange(gene_midpoint, .by_group = TRUE) %>%
  mutate(
    display_label_plot = case_when(
      local_family_id == "M20_local_015" &
        cumsum(local_family_id == "M20_local_015") == 1 ~
        "1-heme cyt.c paralogue block",
      TRUE ~ display_label_plot
    )
  ) %>%
  ungroup()

architecture_levels <- c(
  "core_mtoa",
  "core_mtra",
  "core_fegenie_negative",
  "core_mtrb",
  "core_cytc552",
  "recurrent_unknown",
  "accessory_redox_block",
  "one_heme_paralogue_block",
  "fegenie_other",
  "findmehemes_other",
  "other_mmseq_family",
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
## Track and label preparation
## ================================================================== ##

genome_order <- genes %>%
  distinct(
    plot_region_id,
    genome,
    taxonomy_family,
    taxonomy_genus,
    taxonomy_species,
    taxonomy_display,
    selection_order
  ) %>%
  arrange(selection_order) %>%
  mutate(y_track = rev(seq_len(n())))

genes <- genes %>%
  left_join(
    genome_order %>% select(plot_region_id, y_track),
    by = "plot_region_id"
  )

label_genes <- genes %>%
  filter(!is.na(display_label_plot), display_label_plot != "") %>%
  group_by(plot_region_id) %>%
  arrange(gene_midpoint, .by_group = TRUE) %>%
  mutate(
    label_side = if_else(row_number() %% 2 == 1, "above", "below")
  ) %>%
  ungroup()

labels_above <- label_genes %>% filter(label_side == "above")
labels_below <- label_genes %>% filter(label_side == "below")

## ================================================================== ##
## X scale and contig spans
## ================================================================== ##

x_max <- ceiling(
  max(abs(c(genes$plot_start_bp, genes$plot_end_bp)), na.rm = TRUE) / 1000
) * 1000

x_breaks <- zero_centred_breaks(x_max)

contig_tracks <- genome_order %>%
  transmute(
    track_label = taxonomy_display,
    y_track
  )

contig_spans <- build_contig_spans(
  contig_qc = contig_qc,
  figure_name = "Module_20",
  track_table = contig_tracks,
  x_limit = x_max
)

## ================================================================== ##
## Visual scales
## ================================================================== ##

architecture_colors <- c(
  "core_mtoa" = "#0072B2",
  "core_mtra" = "#8C6BB1",
  "core_fegenie_negative" = "#999999",
  "core_mtrb" = "#E69F00",
  "core_cytc552" = "#009E73",
  "recurrent_unknown" = "#D8A600",
  "accessory_redox_block" = "#CC79A7",
  "one_heme_paralogue_block" = "#B8B0A5",
  "fegenie_other" = "#D55E00",
  "findmehemes_other" = "#73A857",
  "other_mmseq_family" = "#9ECAE1",
  "background_gene" = "#D9D9D9"
)

legend_breaks <- c(
  "core_mtoa",
  "core_mtra",
  "core_fegenie_negative",
  "core_mtrb",
  "core_cytc552",
  "recurrent_unknown",
  "accessory_redox_block",
  "one_heme_paralogue_block",
  "findmehemes_other"
)

legend_labels <- c(
  "core_mtoa" = "MtoA",
  "core_mtra" = "MtrA",
  "core_fegenie_negative" = "Cluster_00035 homolog",
  "core_mtrb" = "MtrB",
  "core_cytc552" = "CytC551/552-like",
  "recurrent_unknown" = "M20-U1",
  "accessory_redox_block" = "Accessory cytochrome/redox block",
  "one_heme_paralogue_block" = "1-heme cyt.c paralogue block",
  "findmehemes_other" = "Other heme protein"
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
      y = y_track,
      yend = y_track
    ),
    inherit.aes = FALSE,
    colour = "grey78",
    linewidth = 0.45
  ) +

  geom_point(
    data = contig_spans %>% filter(left_boundary_visible),
    aes(x = sequence_start_bp, y = y_track),
    inherit.aes = FALSE,
    shape = 124,
    size = 2.6,
    colour = "grey55"
  ) +

  geom_point(
    data = contig_spans %>% filter(right_boundary_visible),
    aes(x = sequence_end_bp, y = y_track),
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
      y = y_track,
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
    aes(x = gene_midpoint, y = y_track, label = display_label_plot),
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
    fontface = "italic"
  ) +

  geom_text_repel(
    data = labels_below,
    aes(x = gene_midpoint, y = y_track, label = display_label_plot),
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
    limits = c(-x_max, x_max),
    breaks = x_breaks,
    labels = kb_labels,
    expand = expansion(mult = c(0.01, 0.01))
  ) +

  scale_y_continuous(
    breaks = genome_order$y_track,
    labels = genome_order$taxonomy_display,
    expand = expansion(add = c(0.7, 0.7))
  ) +

  labs(
    x = "Position relative to focal Cluster_00035 protein",
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
      size = 10,
      face = "italic",
      colour = "black",
      margin = margin(r = 8)
    ),
    axis.text.x = element_text(size = 9),
    axis.title.x = element_text(margin = margin(t = 8)),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.text = element_text(size = 9),
    panel.grid = element_blank(),
    plot.margin = margin(t = 15, r = 15, b = 10, l = 10)
  )

## ================================================================== ##
## Save
## ================================================================== ##

ggsave(output_pdf, p, width = 12, height = 9, units = "in")
ggsave(output_png, p, width = 12, height = 9, units = "in", dpi = 400)

message("")
message("Module 20 gene-map rendering complete.")
message("Focal regions: ", n_distinct(genes$plot_region_id))
message("Genomes: ", n_distinct(genes$genome))
message("Genes drawn: ", nrow(genes))
message("Labels drawn: ", nrow(label_genes))
message("PDF: ", output_pdf)
message("PNG: ", output_png)
