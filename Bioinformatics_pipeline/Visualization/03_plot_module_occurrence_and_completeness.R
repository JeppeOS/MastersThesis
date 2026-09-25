#!/usr/bin/env Rscript

## ================================================================== ##
## Occurrence and completeness of MCL communities
##
## x-axis:
##   Number of genomes containing at least one member family
##
## y-axis:
##   Median module completeness among genomes containing the module
##
## Point size:
##   Number of member protein families in the MCL community
##
## Modules 10, 20 and 35 are highlighted.
##
## No error bars are shown because this figure is intended as a
## global descriptive overview. Within-module variation can instead
## be reported using IQRs in a supplementary table or targeted
## distribution plot.
## ================================================================== ##

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
})


## ================================================================== ##
## Paths
## ================================================================== ##

project_dir <- path.expand(
  "~/methanotrophs/methanotroph_project/jeppe/genome_analysis_workflow"
)

input_file <- file.path(
  project_dir,
  "09_module_occurrence",
  "module_completeness_summary.tsv"
)

output_dir <- file.path(
  project_dir,
  "16_visualization",
  "16F_module_occurrence"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

output_png <- file.path(
  output_dir,
  "module_occurrence_vs_median_completeness.png"
)

output_pdf <- file.path(
  output_dir,
  "module_occurrence_vs_median_completeness.pdf"
)


## ================================================================== ##
## Check input
## ================================================================== ##

if (!file.exists(input_file)) {

  stop(
    "Input file does not exist:\n",
    input_file
  )
}


## ================================================================== ##
## Read authoritative Stage-09 summary
## ================================================================== ##

x <- read.delim(
  input_file,
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

message(
  "Rows read: ",
  nrow(x)
)

message(
  "Columns detected: ",
  paste(
    names(x),
    collapse = ", "
  )
)


## ================================================================== ##
## Validate required columns
## ================================================================== ##

required_columns <- c(
  "module",
  "module_size_clusters",
  "n_genomes_present",
  "median_completeness_present"
)

missing_columns <- setdiff(
  required_columns,
  names(x)
)

if (length(missing_columns) > 0) {

  stop(
    "Missing required Stage-09 columns: ",
    paste(
      missing_columns,
      collapse = ", "
    )
  )
}


## ================================================================== ##
## Prepare plotting table
## ================================================================== ##

plot_df <- data.frame(

  module =
    as.character(
      x$module
    ),

  n_member_families =
    as.numeric(
      x$module_size_clusters
    ),

  n_genomes =
    as.numeric(
      x$n_genomes_present
    ),

  median_completeness_pct =
    100 *
    as.numeric(
      x$median_completeness_present
    ),

  stringsAsFactors = FALSE
)


## ================================================================== ##
## Focal modules
## ================================================================== ##

focal_modules <- c(
  "Module_10",
  "Module_20",
  "Module_35"
)

plot_df$highlight_group <- ifelse(
  plot_df$module %in% focal_modules,
  plot_df$module,
  "Other"
)

plot_df$highlight_group <- factor(
  plot_df$highlight_group,
  levels = c(
    "Other",
    "Module_10",
    "Module_20",
    "Module_35"
  )
)

plot_df$label <- ifelse(
  plot_df$module %in% focal_modules,
  plot_df$module,
  NA_character_
)


## ================================================================== ##
## QC
## ================================================================== ##

if (nrow(plot_df) != 35) {

  warning(
    "Expected 35 MCL communities but found ",
    nrow(plot_df),
    "."
  )
}

if (
  any(
    plot_df$median_completeness_pct < 0 |
      plot_df$median_completeness_pct > 100,
    na.rm = TRUE
  )
) {

  stop(
    "Median completeness contains values outside 0-100%."
  )
}


message("")
message("Focal modules:")

print(
  plot_df[
    plot_df$module %in% focal_modules,
    c(
      "module",
      "n_member_families",
      "n_genomes",
      "median_completeness_pct"
    )
  ]
)


## ================================================================== ##
## Colours
## ================================================================== ##

module_colours <- c(

  "Other" =
    "#D0D0D0",

  "Module_10" =
    "#009E73",

  "Module_20" =
    "#E69F00",

  "Module_35" =
    "#7570B3"
)


## ================================================================== ##
## Plot
## ================================================================== ##

p <- ggplot(
  plot_df,
  aes(
    x = n_genomes,
    y = median_completeness_pct
  )
) +

  ## -------------------------------------------------------------- ##
  ## Communities
  ## -------------------------------------------------------------- ##

  geom_point(
    aes(
      size = n_member_families,
      fill = highlight_group
    ),
    shape = 21,
    colour = "grey35",
    stroke = 0.45,
    alpha = 0.95
  ) +

  ## -------------------------------------------------------------- ##
  ## Labels only for focal modules
  ## -------------------------------------------------------------- ##

  geom_text_repel(
    data = plot_df[
      !is.na(plot_df$label),
    ],
    aes(
      label = label
    ),
    size = 4.2,
    fontface = "bold",
    colour = "grey15",
    box.padding = 0.45,
    point.padding = 0.4,
    segment.colour = "grey55",
    segment.size = 0.3,
    min.segment.length = 0,
    seed = 20260821,
    show.legend = FALSE
  ) +

  ## -------------------------------------------------------------- ##
  ## Scales
  ## -------------------------------------------------------------- ##

  scale_fill_manual(
    values = module_colours,
    name = "MCL community"
  ) +

  scale_size_continuous(
    name = "Member families",
    range = c(
      2.5,
      9
    ),
    breaks = c(
      2,
      5,
      10,
      15,
      20
    )
  ) +

  scale_x_continuous(
    name = "Number of genomes",
    breaks = scales::pretty_breaks(
      n = 6
    ),
    expand = expansion(
      mult = c(
        0.02,
        0.05
      )
    )
  ) +

  scale_y_continuous(
    name = "Median module completeness (%)",
    limits = c(
      0,
      100
    ),
    breaks = seq(
      0,
      100,
      by = 20
    ),
    expand = expansion(
      mult = c(
        0.01,
        0.02
      )
    )
  ) +

  ## -------------------------------------------------------------- ##
  ## Titles
  ## -------------------------------------------------------------- ##

  labs(
    title =
      "Occurrence and completeness of MCL communities",

    subtitle =
      "Point size represents the number of member protein families"
  ) +

  ## -------------------------------------------------------------- ##
  ## Theme
  ## -------------------------------------------------------------- ##

  theme_classic(
    base_size = 12
  ) +

  theme(

    plot.title =
      element_text(
        size = 15,
        face = "bold",
        margin = margin(
          b = 4
        )
      ),

    plot.subtitle =
      element_text(
        size = 10.5,
        colour = "grey30",
        margin = margin(
          b = 10
        )
      ),

    axis.title =
      element_text(
        size = 11
      ),

    axis.text =
      element_text(
        colour = "grey25"
      ),

    legend.position =
      "right",

    legend.title =
      element_text(
        size = 10,
        face = "bold"
      ),

    legend.text =
      element_text(
        size = 9
      ),

    plot.margin =
      margin(
        12,
        18,
        12,
        12
      )
  ) +

  guides(

    size =
      guide_legend(
        order = 1
      ),

    fill =
      guide_legend(
        order = 2,
        override.aes = list(
          size = 5
        )
      )
  )


## ================================================================== ##
## Save
## ================================================================== ##

ggsave(
  filename = output_png,
  plot = p,
  width = 9,
  height = 6.5,
  units = "in",
  dpi = 400,
  bg = "white"
)

ggsave(
  filename = output_pdf,
  plot = p,
  width = 9,
  height = 6.5,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)


## ================================================================== ##
## Finish
## ================================================================== ##

message("")
message("Figure written:")
message("  PNG: ", output_png)
message("  PDF: ", output_pdf)
message("")
message("Done.")