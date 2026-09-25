#!/usr/bin/env Rscript

## ================================================================== ##
## Comparative genomics workflow figure
##
## Purpose:
##   Clean, generic, publication-style vertical flowchart for a
##   protein-family comparative genomics pipeline.
##
## Notes:
##   - No dataset-specific counts are shown.
##   - The four predictor/annotation tools are displayed as parallel
##     evidence streams.
##   - The final figure is produced entirely in R using ggplot2/grid.
##   - Output is vector PDF plus high-resolution PNG.
## ================================================================== ##

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(grid)
})

## ================================================================== ##
## Output directory
## ================================================================== ##

workflow_dir <- path.expand(
  "~/methanotrophs/methanotroph_project/jeppe/genome_analysis_workflow"
)

out_dir <- file.path(
  workflow_dir,
  "16_visualization",
  "16D_pipeline_overview"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

pdf_out <- file.path(
  out_dir,
  "comparative_genomics_pipeline_vertical.pdf"
)

png_out <- file.path(
  out_dir,
  "comparative_genomics_pipeline_vertical.png"
)

## ================================================================== ##
## Palette
## ================================================================== ##

cols <- c(
  input_fill       = "#F0F0F0",
  input_border     = "#7A7A7A",

  evidence_fill    = "#FCE9D5",
  evidence_border  = "#E88C2A",

  sequence_fill    = "#EAF2FC",
  sequence_border  = "#4C82D4",

  network_fill     = "#F0EAFB",
  network_border   = "#8C63D9",

  context_fill     = "#EAF5EE",
  context_border   = "#4D9A72",

  arrow            = "#666666",
  text             = "#1F1F1F",
  title            = "#111111"
)

## ================================================================== ##
## Layout constants
## ================================================================== ##

main_x <- 0

main_w <- 5.6
main_h <- 0.84

tool_w <- 2.35
tool_h <- 0.82

tool_x <- c(
  -4.05,
  -1.35,
   1.35,
   4.05
)

## Vertical positions
y_title        <- 16.2
y_genome       <- 15.0
y_catalogue    <- 13.6
y_tools        <- 11.9
y_integrated   <- 10.4
y_exported     <- 9.1
y_mmseqs       <- 7.9
y_families     <- 6.7
y_prevalence   <- 5.5
y_retained     <- 4.3
y_jaccard      <- 3.0
y_network      <- 1.8
y_mcl          <- 0.6
y_communities  <- -0.6
y_neighborhood <- -1.9
y_architecture <- -3.2

## ================================================================== ##
## Helpers
## ================================================================== ##

make_box <- function(
    id,
    x,
    y,
    width,
    height,
    label,
    fill,
    border,
    text_size = 4.3,
    fontface = "plain"
) {

  tibble(
    id = id,
    x = x,
    y = y,
    xmin = x - width / 2,
    xmax = x + width / 2,
    ymin = y - height / 2,
    ymax = y + height / 2,
    label = label,
    fill = fill,
    border = border,
    text_size = text_size,
    fontface = fontface
  )
}

## ================================================================== ##
## Boxes
## ================================================================== ##

boxes <- bind_rows(

  make_box(
    "genome",
    main_x,
    y_genome,
    main_w,
    main_h,
    "Genome dataset",
    cols["input_fill"],
    cols["input_border"],
    fontface = "bold"
  ),

  make_box(
    "catalogue",
    main_x,
    y_catalogue,
    main_w,
    main_h,
    "Protein catalogue",
    cols["input_fill"],
    cols["input_border"],
    fontface = "bold"
  ),

  make_box(
    "fegenie",
    tool_x[1],
    y_tools,
    tool_w,
    tool_h,
    "FeGenie",
    cols["evidence_fill"],
    cols["evidence_border"],
    fontface = "bold"
  ),

  make_box(
    "findmehemes",
    tool_x[2],
    y_tools,
    tool_w,
    tool_h,
    "FindMeHemes",
    cols["evidence_fill"],
    cols["evidence_border"],
    fontface = "bold"
  ),

  make_box(
    "signalp",
    tool_x[3],
    y_tools,
    tool_w,
    tool_h,
    "SignalP",
    cols["evidence_fill"],
    cols["evidence_border"],
    fontface = "bold"
  ),

  make_box(
    "deeptmhmm",
    tool_x[4],
    y_tools,
    tool_w,
    tool_h,
    "DeepTMHMM",
    cols["evidence_fill"],
    cols["evidence_border"],
    fontface = "bold"
  ),

  make_box(
    "integrated",
    main_x,
    y_integrated,
    main_w,
    main_h,
    "Integrated evidence",
    cols["evidence_fill"],
    cols["evidence_border"],
    fontface = "bold"
  ),

  make_box(
    "exported",
    main_x,
    y_exported,
    main_w,
    main_h,
    "Exported heme proteins",
    cols["sequence_fill"],
    cols["sequence_border"],
    fontface = "bold"
  ),

  make_box(
    "mmseqs",
    main_x,
    y_mmseqs,
    main_w,
    main_h,
    "MMseqs2 clustering",
    cols["sequence_fill"],
    cols["sequence_border"],
    fontface = "bold"
  ),

  make_box(
    "families",
    main_x,
    y_families,
    main_w,
    main_h,
    "Protein families",
    cols["sequence_fill"],
    cols["sequence_border"],
    fontface = "bold"
  ),

  make_box(
    "prevalence",
    main_x,
    y_prevalence,
    main_w,
    main_h,
    "Prevalence filtering",
    cols["sequence_fill"],
    cols["sequence_border"],
    fontface = "bold"
  ),

  make_box(
    "retained",
    main_x,
    y_retained,
    main_w,
    main_h,
    "Retained families",
    cols["sequence_fill"],
    cols["sequence_border"],
    fontface = "bold"
  ),

  make_box(
    "jaccard",
    main_x,
    y_jaccard,
    main_w,
    main_h,
    "Jaccard co-occurrence",
    cols["network_fill"],
    cols["network_border"],
    fontface = "bold"
  ),

  make_box(
    "network",
    main_x,
    y_network,
    main_w,
    main_h,
    "Co-occurrence network",
    cols["network_fill"],
    cols["network_border"],
    fontface = "bold"
  ),

  make_box(
    "mcl",
    main_x,
    y_mcl,
    main_w,
    main_h,
    "MCL clustering",
    cols["network_fill"],
    cols["network_border"],
    fontface = "bold"
  ),

  make_box(
    "communities",
    main_x,
    y_communities,
    main_w,
    main_h,
    "Communities",
    cols["network_fill"],
    cols["network_border"],
    fontface = "bold"
  ),

  make_box(
    "neighborhood",
    main_x,
    y_neighborhood,
    main_w,
    main_h,
    "Neighborhood analysis",
    cols["context_fill"],
    cols["context_border"],
    fontface = "bold"
  ),

  make_box(
    "architecture",
    main_x,
    y_architecture,
    main_w,
    main_h,
    "Representative architectures",
    cols["context_fill"],
    cols["context_border"],
    fontface = "bold"
  )
)

## ================================================================== ##
## Lookup helpers
## ================================================================== ##

lookup <- boxes %>%
  select(id, x, y, xmin, xmax, ymin, ymax)

bx <- function(id) lookup$x[lookup$id == id]
by <- function(id) lookup$y[lookup$id == id]
btop <- function(id) lookup$ymax[lookup$id == id]
bbottom <- function(id) lookup$ymin[lookup$id == id]

## ================================================================== ##
## Main vertical arrows
## ================================================================== ##

main_pairs <- tribble(
  ~from,           ~to,
  "genome",        "catalogue",
  "integrated",    "exported",
  "exported",      "mmseqs",
  "mmseqs",        "families",
  "families",      "prevalence",
  "prevalence",    "retained",
  "retained",      "jaccard",
  "jaccard",       "network",
  "network",       "mcl",
  "mcl",           "communities",
  "communities",   "neighborhood",
  "neighborhood",  "architecture"
)

main_edges <- main_pairs %>%
  rowwise() %>%
  mutate(
    x = bx(from),
    xend = bx(to),
    y = bbottom(from),
    yend = btop(to)
  ) %>%
  ungroup()

## ================================================================== ##
## Branch geometry: catalogue -> four tools -> integrated evidence
##
## This is intentionally drawn as clean orthogonal connector geometry.
## ================================================================== ##

branch_top_y <- 12.75
branch_bottom_y <- 11.00

catalogue_drop <- tibble(
  x = main_x,
  y = bbottom("catalogue"),
  xend = main_x,
  yend = branch_top_y
)

top_bus <- tibble(
  x = min(tool_x),
  y = branch_top_y,
  xend = max(tool_x),
  yend = branch_top_y
)

top_drops <- tibble(
  x = tool_x,
  y = branch_top_y,
  xend = tool_x,
  yend = y_tools + tool_h / 2
)

bottom_risers <- tibble(
  x = tool_x,
  y = y_tools - tool_h / 2,
  xend = tool_x,
  yend = branch_bottom_y
)

bottom_bus <- tibble(
  x = min(tool_x),
  y = branch_bottom_y,
  xend = max(tool_x),
  yend = branch_bottom_y
)

integrated_drop <- tibble(
  x = main_x,
  y = branch_bottom_y,
  xend = main_x,
  yend = btop("integrated")
)

## ================================================================== ##
## Plot
## ================================================================== ##

p <- ggplot() +

  ## Title
  annotate(
    "text",
    x = 0,
    y = y_title,
    label = "Protein-family comparative genomics pipeline",
    size = 6.6,
    fontface = "bold",
    colour = cols["title"]
  ) +

  ## Main vertical arrows
  geom_segment(
    data = main_edges,
    aes(
      x = x,
      y = y,
      xend = xend,
      yend = yend
    ),
    linewidth = 0.55,
    colour = cols["arrow"],
    lineend = "round",
    arrow = arrow(
      length = unit(2.2, "mm"),
      type = "closed"
    )
  ) +

  ## Branch: catalogue down to top bus
  geom_segment(
    data = catalogue_drop,
    aes(
      x = x,
      y = y,
      xend = xend,
      yend = yend
    ),
    linewidth = 0.55,
    colour = cols["arrow"],
    lineend = "round"
  ) +

  ## Branch: horizontal top bus
  geom_segment(
    data = top_bus,
    aes(
      x = x,
      y = y,
      xend = xend,
      yend = yend
    ),
    linewidth = 0.55,
    colour = cols["arrow"],
    lineend = "round"
  ) +

  ## Branch: top bus down into tools
  geom_segment(
    data = top_drops,
    aes(
      x = x,
      y = y,
      xend = xend,
      yend = yend
    ),
    linewidth = 0.55,
    colour = cols["arrow"],
    lineend = "round",
    arrow = arrow(
      length = unit(2.0, "mm"),
      type = "closed"
    )
  ) +

  ## Branch: tools down to bottom bus
  geom_segment(
    data = bottom_risers,
    aes(
      x = x,
      y = y,
      xend = xend,
      yend = yend
    ),
    linewidth = 0.55,
    colour = cols["arrow"],
    lineend = "round"
  ) +

  ## Branch: bottom bus
  geom_segment(
    data = bottom_bus,
    aes(
      x = x,
      y = y,
      xend = xend,
      yend = yend
    ),
    linewidth = 0.55,
    colour = cols["arrow"],
    lineend = "round"
  ) +

  ## Branch: bottom bus to integrated evidence
  geom_segment(
    data = integrated_drop,
    aes(
      x = x,
      y = y,
      xend = xend,
      yend = yend
    ),
    linewidth = 0.55,
    colour = cols["arrow"],
    lineend = "round",
    arrow = arrow(
      length = unit(2.2, "mm"),
      type = "closed"
    )
  ) +

  ## Boxes
  geom_rect(
    data = boxes,
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax,
      fill = fill,
      colour = border
    ),
    linewidth = 0.65
  ) +

  ## Labels
  geom_text(
    data = boxes,
    aes(
      x = x,
      y = y,
      label = label,
      size = text_size,
      fontface = fontface
    ),
    colour = cols["text"],
    lineheight = 0.98
  ) +

  scale_fill_identity() +
  scale_colour_identity() +
  scale_size_identity() +

  coord_cartesian(
    xlim = c(-5.7, 5.7),
    ylim = c(-4.0, 16.8),
    clip = "off"
  ) +

  theme_void(base_size = 11) +

  theme(
    plot.background = element_rect(
      fill = "white",
      colour = NA
    ),
    panel.background = element_rect(
      fill = "white",
      colour = NA
    ),
    plot.margin = margin(
      t = 16,
      r = 16,
      b = 16,
      l = 16
    )
  )

## ================================================================== ##
## Save
## ================================================================== ##

ggsave(
  filename = pdf_out,
  plot = p,
  width = 8.6,
  height = 12.2,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = png_out,
  plot = p,
  width = 8.6,
  height = 12.2,
  units = "in",
  dpi = 400
)

message("")
message("Vertical comparative-genomics workflow written:")
message("  PDF: ", pdf_out)
message("  PNG: ", png_out)
