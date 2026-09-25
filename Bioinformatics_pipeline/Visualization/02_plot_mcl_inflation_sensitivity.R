#!/usr/bin/env Rscript

## ================================================================== ##
## MCL inflation sensitivity analysis
##
## Tests stability of MCL community assignments across inflation
## parameters.
##
## Panel A:
##   Number of MCL communities
##
## Panel B:
##   Adjusted Rand Index (ARI) relative to the selected
##   inflation = 2.0 solution
##
## Inflation = 2.0 is the selected value.
## ================================================================== ##

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
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
  "16H_mcl_sensitivity"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

output_png <- file.path(
  output_dir,
  "mcl_inflation_sensitivity.png"
)

output_pdf <- file.path(
  output_dir,
  "mcl_inflation_sensitivity.pdf"
)


## ================================================================== ##
## Authoritative sensitivity results
## ================================================================== ##

sensitivity <- data.frame(

  inflation = c(
    1.5,
    2.0,
    2.5,
    3.0
  ),

  n_communities = c(
    35,
    35,
    38,
    39
  ),

  ari_vs_selected = c(
    1.000,
    1.000,
    0.928,
    0.918
  ),

  n_assigned_nodes = c(
    156,
    156,
    156,
    155
  )
)


## ================================================================== ##
## Selected inflation
## ================================================================== ##

selected_inflation <- 2.0

selected_df <- sensitivity[
  sensitivity$inflation == selected_inflation,
]


## ================================================================== ##
## QC
## ================================================================== ##

if (nrow(selected_df) != 1) {

  stop(
    "Could not identify exactly one selected inflation value."
  )
}

if (
  any(
    sensitivity$ari_vs_selected < 0 |
      sensitivity$ari_vs_selected > 1
  )
) {

  stop(
    "ARI values must fall between 0 and 1."
  )
}


message("")
message("MCL inflation sensitivity:")
print(sensitivity)

message("")
message(
  "Selected inflation: ",
  selected_inflation
)

message(
  "Communities at selected inflation: ",
  selected_df$n_communities
)

message(
  "Nodes assigned at selected inflation: ",
  selected_df$n_assigned_nodes
)


## ================================================================== ##
## Shared theme
## ================================================================== ##

shared_theme <- theme_classic(
  base_size = 12
) +

  theme(

    plot.title =
      element_text(
        size = 11.5,
        face = "bold",
        hjust = 0.5,
        margin = margin(
          b = 7
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

    plot.margin =
      margin(
        8,
        12,
        8,
        8
      )
  )


## ================================================================== ##
## Panel A
## Number of communities
## ================================================================== ##

p_communities <- ggplot(
  sensitivity,
  aes(
    x = inflation,
    y = n_communities
  )
) +

  ## selected inflation
  geom_vline(
    xintercept = selected_inflation,
    linetype = "dashed",
    linewidth = 0.45,
    colour = "grey55"
  ) +

  ## trajectory
  geom_line(
    linewidth = 0.85,
    colour = "grey40"
  ) +

  geom_point(
    size = 3.0,
    shape = 21,
    fill = "grey78",
    colour = "grey35",
    stroke = 0.5
  ) +

  ## selected point
  geom_point(
    data = selected_df,
    size = 4.4,
    shape = 21,
    fill = "#0072B2",
    colour = "grey20",
    stroke = 0.65
  ) +

  ## values
  geom_text(
    aes(
      label = n_communities
    ),
    nudge_y = 0.55,
    size = 3.5,
    colour = "grey25"
  ) +

  scale_x_continuous(
    name = "MCL inflation",
    breaks = sensitivity$inflation,
    limits = c(
      1.4,
      3.1
    ),
    expand = expansion(
      mult = c(
        0,
        0
      )
    )
  ) +

  scale_y_continuous(
    name = "Number of communities",
    breaks = seq(
      34,
      40,
      by = 1
    ),
    limits = c(
      34,
      40
    ),
    expand = expansion(
      mult = c(
        0,
        0
      )
    )
  ) +

  labs(
    title = "A   Community number"
  ) +

  shared_theme


## ================================================================== ##
## Panel B
## ARI relative to selected solution
## ================================================================== ##

p_ari <- ggplot(
  sensitivity,
  aes(
    x = inflation,
    y = ari_vs_selected
  )
) +

  ## selected inflation
  geom_vline(
    xintercept = selected_inflation,
    linetype = "dashed",
    linewidth = 0.45,
    colour = "grey55"
  ) +

  ## trajectory
  geom_line(
    linewidth = 0.85,
    colour = "grey40"
  ) +

  geom_point(
    size = 3.0,
    shape = 21,
    fill = "grey78",
    colour = "grey35",
    stroke = 0.5
  ) +

  ## selected point
  geom_point(
    data = selected_df,
    size = 4.4,
    shape = 21,
    fill = "#0072B2",
    colour = "grey20",
    stroke = 0.65
  ) +

  ## values
  geom_text(
    aes(
      label = sprintf(
        "%.3f",
        ari_vs_selected
      )
    ),
    nudge_y = 0.006,
    size = 3.5,
    colour = "grey25"
  ) +

  scale_x_continuous(
    name = "MCL inflation",
    breaks = sensitivity$inflation,
    limits = c(
      1.4,
      3.1
    ),
    expand = expansion(
      mult = c(
        0,
        0
      )
    )
  ) +

  scale_y_continuous(
    name = "ARI relative to inflation = 2.0",
    breaks = seq(
      0.90,
      1.00,
      by = 0.02
    ),
    limits = c(
      0.90,
      1.01
    ),
    expand = expansion(
      mult = c(
        0,
        0
      )
    )
  ) +

  labs(
    title = "B   Assignment stability"
  ) +

  shared_theme


## ================================================================== ##
## Combine panels
## ================================================================== ##

combined_plot <-

  p_communities +
  p_ari +

  plot_layout(
    ncol = 2,
    widths = c(
      1,
      1
    )
  ) +

  plot_annotation(

    title =
      "Sensitivity of MCL clustering to inflation parameter",

    subtitle =
      paste0(
        "Dashed line and blue points indicate the selected inflation ",
        "(I = 2.0)"
      ),

    caption =
      paste0(
        "ARI = adjusted Rand index relative to the selected I = 2.0 ",
        "partition. At I = 3.0, 155 of 156 network nodes were assigned."
      ),

    theme =
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

        plot.caption =
          element_text(
            size = 8.5,
            colour = "grey40",
            hjust = 0,
            margin = margin(
              t = 8
            )
          )
      )
  )


## ================================================================== ##
## Save
## ================================================================== ##

ggsave(
  filename = output_png,
  plot = combined_plot,
  width = 9,
  height = 5.5,
  units = "in",
  dpi = 400,
  bg = "white"
)

ggsave(
  filename = output_pdf,
  plot = combined_plot,
  width = 9,
  height = 5.5,
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
