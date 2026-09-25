#!/usr/bin/env Rscript

## ================================================================== ##
## Jaccard-threshold sensitivity analysis
##
## Shows how network topology changes across tested Jaccard cutoffs.
##
## Metrics:
##   - retained network nodes
##   - retained network edges
##   - connected components
##   - size of largest connected component
##
## Jaccard = 0.25 is the selected threshold.
## ================================================================== ##

suppressPackageStartupMessages({
  library(ggplot2)
})


## ================================================================== ##
## Paths
## ================================================================== ##

project_dir <- path.expand(
  "~/methanotrophs/methanotroph_project/jeppe/genome_analysis_workflow"
)

output_dir <- file.path(
  project_dir,
  "16_visualization",
  "16G_jaccard_sensitivity"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

output_png <- file.path(
  output_dir,
  "jaccard_threshold_sensitivity.png"
)

output_pdf <- file.path(
  output_dir,
  "jaccard_threshold_sensitivity.pdf"
)


## ================================================================== ##
## Authoritative sensitivity results
## ================================================================== ##

sensitivity <- data.frame(

  jaccard = c(
    0.15,
    0.20,
    0.25,
    0.30,
    0.35,
    0.40
  ),

  nodes = c(
    198,
    177,
    156,
    126,
    105,
    86
  ),

  edges = c(
    841,
    506,
    332,
    226,
    157,
    118
  ),

  subnetworks = c(
    6,
    18,
    26,
    28,
    27,
    24
  ),

  largest_subnetwork = c(
    188,
    127,
    36,
    23,
    19,
    14
  )
)


## ================================================================== ##
## Convert to plotting format
## ================================================================== ##

plot_df <- rbind(

  data.frame(
    jaccard = sensitivity$jaccard,
    metric = "Nodes retained",
    value = sensitivity$nodes
  ),

  data.frame(
    jaccard = sensitivity$jaccard,
    metric = "Edges retained",
    value = sensitivity$edges
  ),

  data.frame(
    jaccard = sensitivity$jaccard,
    metric = "Disconnected modules",
    value = sensitivity$subnetworks
  ),

  data.frame(
    jaccard = sensitivity$jaccard,
    metric = "Nodes in largest module",
    value = sensitivity$largest_subnetwork
  )
)


## ================================================================== ##
## Set panel order
## ================================================================== ##

plot_df$metric <- factor(
  plot_df$metric,
  levels = c(
    "Nodes retained",
    "Edges retained",
    "Disconnected modules",
    "Nodes in largest module"
  )
)


## ================================================================== ##
## Mark selected threshold
## ================================================================== ##

selected_threshold <- 0.25

plot_df$selected <- (
  plot_df$jaccard == selected_threshold
)

selected_df <- plot_df[
  plot_df$selected,
]


## ================================================================== ##
## QC
## ================================================================== ##

message("")
message("Jaccard sensitivity table:")
print(sensitivity)

message("")
message("Selected threshold:")
print(selected_df)


## ================================================================== ##
## Plot
## ================================================================== ##

p <- ggplot(
  plot_df,
  aes(
    x = jaccard,
    y = value
  )
) +

  ## -------------------------------------------------------------- ##
  ## Selected threshold
  ## -------------------------------------------------------------- ##

  geom_vline(
    xintercept = selected_threshold,
    linetype = "dashed",
    linewidth = 0.45,
    colour = "grey55"
  ) +

  ## -------------------------------------------------------------- ##
  ## Sensitivity trajectory
  ## -------------------------------------------------------------- ##

  geom_line(
    linewidth = 0.8,
    colour = "grey40"
  ) +

  geom_point(
    size = 2.8,
    shape = 21,
    fill = "grey78",
    colour = "grey35",
    stroke = 0.45
  ) +

  ## -------------------------------------------------------------- ##
  ## Selected Jaccard = 0.25 points
  ## -------------------------------------------------------------- ##

  geom_point(
    data = selected_df,
    size = 4.2,
    shape = 21,
    fill = "#0072B2",
    colour = "grey20",
    stroke = 0.6
  ) +

  ## -------------------------------------------------------------- ##
  ## Label value at selected threshold
  ## -------------------------------------------------------------- ##

  geom_text(
    data = selected_df,
    aes(
      label = value
    ),
    nudge_x = 0.009,
    hjust = 0,
    vjust = 0.5,
    size = 3.5,
    fontface = "bold",
    colour = "grey20"
  ) +

  ## -------------------------------------------------------------- ##
  ## Separate scale for each network property
  ## -------------------------------------------------------------- ##

  facet_wrap(
    ~ metric,
    ncol = 2,
    scales = "free_y"
  ) +

  ## -------------------------------------------------------------- ##
  ## Axes
  ## -------------------------------------------------------------- ##

  scale_x_continuous(
    name = "Jaccard threshold",
    breaks = c(
      0.15,
      0.20,
      0.25,
      0.30,
      0.35,
      0.40
    ),
    limits = c(
      0.14,
      0.41
    ),
    expand = expansion(
      mult = c(
        0,
        0
      )
    )
  ) +

  scale_y_continuous(
    name = "Count",
    breaks = scales::pretty_breaks(
      n = 5
    ),
    expand = expansion(
      mult = c(
        0.05,
        0.12
      )
    )
  ) +

  expand_limits(
    y = 0
  ) +

  ## -------------------------------------------------------------- ##
  ## Titles
  ## -------------------------------------------------------------- ##

  labs(
    title =
      "Sensitivity of network structure to Jaccard threshold",

    subtitle =
      "Dashed line and blue points indicate the selected threshold (Jaccard = 0.25)"
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
          b = 12
        )
      ),

    strip.background =
      element_blank(),

    strip.text =
      element_text(
        size = 11,
        face = "bold",
        margin = margin(
          b = 5
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

    panel.spacing =
      grid::unit(
        1.2,
        "lines"
      ),

    plot.margin =
      margin(
        12,
        15,
        12,
        12
      )
  )


## ================================================================== ##
## Save
## ================================================================== ##

ggsave(
  filename = output_png,
  plot = p,
  width = 9,
  height = 7,
  units = "in",
  dpi = 400,
  bg = "white"
)

ggsave(
  filename = output_pdf,
  plot = p,
  width = 9,
  height = 7,
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
