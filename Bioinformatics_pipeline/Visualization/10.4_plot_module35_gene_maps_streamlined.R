#!/usr/bin/env Rscript

## ================================================================== ##
## STAGE 16B - MODULE 35 STACKED GENE MAP
##
## Streamlined plotting script.
## - Preserves the established Module-35 biological classes and style.
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
  "16A2_local_architecture_pipeline",
  "runs",
  "Module_35",
  "Module_35_gene_map_data_with_local_families.tsv"
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
  "16B_module35_gene_maps"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

output_pdf <- file.path(output_dir, "Module_35_all_genomes_gene_map.pdf")
output_png <- file.path(output_dir, "Module_35_all_genomes_gene_map.png")

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

message("Reading Module-35 gene-map data...")

genes <- read_tsv(input_file, show_col_types = FALSE)
contig_qc <- read_tsv(contig_qc_file, show_col_types = FALSE)

require_columns(
  genes,
  c(
    "genome", "protein_id",
    "taxonomy_display", "taxonomy_family", "taxonomy_genus", "taxonomy_species",
    "plot_start_bp", "plot_end_bp", "plot_strand",
    "display_class", "display_label",
    "local_family_id", "neighbor_mmseq_cluster"
  ),
  "Module-35 gene-map table"
)

require_columns(
  contig_qc,
  c(
    "figure", "track_label",
    "plot_left_available_bp", "plot_right_available_bp"
  ),
  "Contig-boundary QC table"
)

if (n_distinct(genes$genome) != 8) {
  stop("Expected 8 Module-35 genomes, found ", n_distinct(genes$genome), ".")
}

if (nrow(genes) != 174) {
  stop("Expected 174 Module-35 genes, found ", nrow(genes), ".")
}

genes <- genes %>%
  mutate(
    plot_start_bp = as.numeric(plot_start_bp),
    plot_end_bp = as.numeric(plot_end_bp),
    gene_midpoint = (plot_start_bp + plot_end_bp) / 2,
    forward = plot_strand == "+"
  )

if (anyNA(genes$plot_start_bp) || anyNA(genes$plot_end_bp)) {
  stop("Missing/non-numeric plotting coordinates detected.")
}

## ================================================================== ##
## Biological classes and labels
## ================================================================== ##

lineage_specific_families <- sprintf("M35_local_%03d", 14:22)

genes <- genes %>%
  mutate(
    architecture_class = case_when(
      local_family_id == "M35_local_001" ~ "core_cyc2",
      local_family_id == "M35_local_002" ~ "core_cytc5",
      local_family_id == "M35_local_008" ~ "core_unknown",

      local_family_id %in% c("M35_local_010", "M35_local_012") ~
        "regulatory_cassette",
      local_family_id == "M35_local_013" ~ "accessory_cytochrome",

      local_family_id %in% c("M35_local_007", "M35_local_009") ~
        "efflux_pair",

      local_family_id %in% lineage_specific_families ~
        "lineage_specific_block",

      display_class == "fegenie_other" ~ "fegenie_other",
      display_class == "findmehemes_other" ~ "findmehemes_other",
      display_class == "other_mmseq_family" ~ "other_mmseq_family",

      TRUE ~ "background_gene"
    ),

    display_label_plot = case_when(
      local_family_id == "M35_local_001" ~ "Cyc2",
      local_family_id == "M35_local_002" ~ "CytC5-like",
      local_family_id == "M35_local_008" ~ "M35-U1",

      local_family_id == "M35_local_010" ~ "BaeS-like",
      local_family_id == "M35_local_012" ~ "OmpR-like",
      local_family_id == "M35_local_013" ~ "CccA/CytC553-like",

      local_family_id == "M35_local_007" ~ "TolC-like",
      local_family_id == "M35_local_009" ~ "AcrA-like",

      local_family_id == "M35_local_014" ~ "PcoB-like",
      local_family_id == "M35_local_017" ~ "DUF411-like",
      local_family_id == "M35_local_019" ~ "TolC-like",
      local_family_id == "M35_local_022" ~ "Multicopper oxidase",

      neighbor_mmseq_cluster == "Cluster_00033" ~ "TrxA-like",

      architecture_class %in% c(
        "fegenie_other",
        "findmehemes_other",
        "other_mmseq_family"
      ) &
        !is.na(display_label) &
        display_label != "" ~ display_label,

      TRUE ~ NA_character_
    )
  )

architecture_levels <- c(
  "core_cyc2",
  "core_cytc5",
  "core_unknown",
  "regulatory_cassette",
  "accessory_cytochrome",
  "efflux_pair",
  "lineage_specific_block",
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
    genome,
    taxonomy_family,
    taxonomy_genus,
    taxonomy_species,
    taxonomy_display
  ) %>%
  arrange(
    taxonomy_family,
    taxonomy_genus,
    taxonomy_species,
    taxonomy_display
  ) %>%
  mutate(y_track = rev(seq_len(n())))

genes <- genes %>%
  left_join(
    genome_order %>% select(genome, y_track),
    by = "genome"
  )

label_genes <- genes %>%
  filter(!is.na(display_label_plot), display_label_plot != "") %>%
  group_by(genome) %>%
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
  figure_name = "Module_35",
  track_table = contig_tracks,
  x_limit = x_max
)

## ================================================================== ##
## Visual scales
## ================================================================== ##

architecture_colors <- c(
  "core_cyc2" = "#0072B2",
  "core_cytc5" = "#E69F00",
  "core_unknown" = "#D8A600",
  "regulatory_cassette" = "#8C6BB1",
  "accessory_cytochrome" = "#009E73",
  "efflux_pair" = "#56B4E9",
  "lineage_specific_block" = "#B8B0A5",
  "fegenie_other" = "#CC79A7",
  "findmehemes_other" = "#73A857",
  "other_mmseq_family" = "#9ECAE1",
  "background_gene" = "#D9D9D9"
)

legend_breaks <- c(
  "core_cyc2",
  "core_cytc5",
  "core_unknown",
  "regulatory_cassette",
  "accessory_cytochrome",
  "efflux_pair",
  "lineage_specific_block"
)

legend_labels <- c(
  "core_cyc2" = "Cyc2 (Cluster_00206)",
  "core_cytc5" = "CytC5-like (Cluster_00228)",
  "core_unknown" = "M35-U1",
  "regulatory_cassette" = "BaeS/OmpR cassette",
  "accessory_cytochrome" = "CccA/CytC553",
  "efflux_pair" = "TolC/AcrA pair",
  "lineage_specific_block" = "3-genome block"
)

## ================================================================== ##
## Plot
## ================================================================== ##

p <- ggplot() +

  ## Available assembled sequence behind the genes.
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

  ## True contig boundaries that fall inside the displayed window.
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
    x = "Position relative to focal Cyc2 (Cluster_00206)",
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

ggsave(output_pdf, p, width = 12, height = 7, units = "in")
ggsave(output_png, p, width = 12, height = 7, units = "in", dpi = 400)

message("")
message("Module 35 gene-map rendering complete.")
message("Genomes: ", n_distinct(genes$genome))
message("Genes drawn: ", nrow(genes))
message("Labels drawn: ", nrow(label_genes))
message("PDF: ", output_pdf)
message("PNG: ", output_png)
