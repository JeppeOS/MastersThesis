#!/usr/bin/env Rscript


## ------------------------------------------------------------
## Compare:
##
## A) SemiBin2 MAG-derived mean coverage
## B) Sylph read-based sequence abundance
##
## Important:
##
## - MAG coverage is NOT normalized to 100%.
## - Sylph abundance IS a percentage.
## - Same-genus MAGs are retained as distinct populations
##   upstream, then summed at genus level for this figure.
## - No manual fragmentation correction is applied.
## - barcode09 is excluded.
## ------------------------------------------------------------


library(tidyverse)
library(patchwork)


## ------------------------------------------------------------
## Paths
## ------------------------------------------------------------

project <- path.expand(
  "~/methanotrophs/methanotroph_project/jeppe/metagenome"
)

semibin_dir <- file.path(
  project,
  "semibin2_analysis"
)

sylph_dir <- file.path(
  project,
  "abundance_comparison"
)

figure_dir <- file.path(
  semibin_dir,
  "figures"
)

dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


mag_file <- file.path(
  semibin_dir,
  "SemiBin2_MAG_coverage_by_category.tsv"
)

sylph_file <- file.path(
  sylph_dir,
  "Sylph_sequence_abundance_by_category.tsv"
)


## ------------------------------------------------------------
## Read data
## ------------------------------------------------------------

mag <- read_tsv(
  mag_file,
  show_col_types = FALSE
)

sylph <- read_tsv(
  sylph_file,
  show_col_types = FALSE
)


## ------------------------------------------------------------
## Sample order
## ------------------------------------------------------------

sample_order <- paste0(
  "barcode",
  10:16,
  "_seqs"
)

sample_labels <- paste0(
  "barcode",
  10:16
)


## ------------------------------------------------------------
## Preferred taxon order
##
## Anything unexpected is retained automatically.
## ------------------------------------------------------------

preferred_order <- c(
  "Methylobacter_A",
  "Methylobacter_C",
  "Methylovulum",
  "Methylomonas",
  "Methylumidiphilus",
  "Methylicorpusculum",
  "UBA4132",
  "UBA6136",
  "Methylococcales_unclassified"
)


all_categories <- union(
  mag$Category,
  sylph$Category
)


extra_categories <- setdiff(
  all_categories,
  c(
    preferred_order,
    "Other"
  )
)


category_order <- c(
  preferred_order[
    preferred_order %in% all_categories
  ],
  sort(
    extra_categories
  ),
  "Other"
)


## ------------------------------------------------------------
## Common colour palette
##
## Same taxon = same colour in both panels.
## ------------------------------------------------------------

taxon_colours <- c(

  "Methylobacter_A" =
    "#1F78B4",

  "Methylobacter_C" =
    "#66C2FF",

  "Methylovulum" =
    "#E6AB02",

  "Methylomonas" =
    "#1B9E77",

  "Methylumidiphilus" =
    "#CC79A7",

  "Methylicorpusculum" =
    "#D95F02",

  "UBA4132" =
    "#B15928",

  "UBA6136" =
    "#8E63CE",

  "Methylococcales_unclassified" =
    "#F0E442",

  "Other" =
    "#BDBDBD"
)


## ------------------------------------------------------------
## Check that every category has a defined colour
## ------------------------------------------------------------

missing_colours <- setdiff(
  category_order,
  names(
    taxon_colours
  )
)


if (
  length(
    missing_colours
  ) > 0
) {

  stop(
    "Missing colours for: ",
    paste(
      missing_colours,
      collapse = ", "
    )
  )
}


## ------------------------------------------------------------
## Complete MAG table
##
## Missing sample × category combinations become zero.
## ------------------------------------------------------------

mag_plot <- mag %>%

  filter(
    Sample %in% sample_order
  ) %>%

  complete(
    Sample = sample_order,
    Category = category_order,
    fill = list(
      Summed_MAG_mean_coverage_x = 0
    )
  ) %>%

  mutate(

    Sample = factor(
      Sample,
      levels = sample_order,
      labels = sample_labels
    ),

    Category = factor(
      Category,
      levels = category_order
    )
  )


## ------------------------------------------------------------
## Complete Sylph table
## ------------------------------------------------------------

sylph_plot <- sylph %>%

  filter(
    Sample %in% sample_order
  ) %>%

  complete(
    Sample = sample_order,
    Category = category_order,
    fill = list(
      Relative_sequence_abundance_pct = 0
    )
  ) %>%

  mutate(

    Sample = factor(
      Sample,
      levels = sample_order,
      labels = sample_labels
    ),

    Category = factor(
      Category,
      levels = category_order
    )
  )


## ------------------------------------------------------------
## QC: Sylph totals
##
## Small deviations such as 100.0002 are rounding.
## ------------------------------------------------------------

cat(
  "\nSylph totals used for plotting:\n"
)

print(

  sylph_plot %>%

    group_by(
      Sample
    ) %>%

    summarise(
      Total_pct = sum(
        Relative_sequence_abundance_pct,
        na.rm = TRUE
      ),
      .groups = "drop"
    )
)


## ------------------------------------------------------------
## QC: print Methylococcales MAG coverage
## ------------------------------------------------------------

cat(
  "\nSemiBin2 Methylococcales genus-level coverage:\n"
)

print(

  mag %>%

    filter(
      Sample %in% sample_order,
      Category != "Other"
    ) %>%

    arrange(
      Sample,
      Category
    )
)


## ------------------------------------------------------------
## Common fill scale
## ------------------------------------------------------------

common_fill_scale <- scale_fill_manual(

  name = "Taxon",

  values = taxon_colours,

  breaks = category_order,

  limits = category_order,

  drop = FALSE,

  na.value = "#000000"
)


## ------------------------------------------------------------
## Explicit legend colours
##
## Prevents taxa absent from one panel from getting blank keys.
## ------------------------------------------------------------

legend_guide <- guide_legend(

  title = "Taxon",

  override.aes = list(

    fill = unname(
      taxon_colours[
        category_order
      ]
    ),

    colour = NA,

    alpha = 1
  )
)


## ------------------------------------------------------------
## Panel A:
## SemiBin2 MAG coverage
##
## This panel provides the single legend.
## ------------------------------------------------------------

p_mag <- ggplot(
  mag_plot,
  aes(
    x = Sample,
    y = Summed_MAG_mean_coverage_x,
    fill = Category
  )
) +

  geom_col(
    width = 0.8
  ) +

  common_fill_scale +

  guides(
    fill = legend_guide
  ) +

  scale_y_continuous(

    expand = expansion(
      mult = c(
        0,
        0.05
      )
    )
  ) +

  labs(
    x = NULL,
    y = "Summed MAG mean coverage (×)",
    title = "A  SemiBin2 MAG-derived population coverage"
  ) +

  theme_classic(
    base_size = 12
  ) +

  theme(

    plot.title = element_text(
      size = 14,
      face = "plain"
    ),

    axis.text.x = element_blank(),

    axis.ticks.x = element_blank(),

    legend.position = "right",

    legend.title = element_text(
      face = "bold"
    )
  )


## ------------------------------------------------------------
## Panel B:
## Sylph read-based abundance
##
## coord_cartesian() is deliberately used rather than
## scale limits so values such as 100.0002 are not censored.
## ------------------------------------------------------------

p_sylph <- ggplot(
  sylph_plot,
  aes(
    x = Sample,
    y = Relative_sequence_abundance_pct,
    fill = Category
  )
) +

  geom_col(
    width = 0.8
  ) +

  common_fill_scale +

  scale_y_continuous(

    breaks = seq(
      0,
      100,
      20
    ),

    labels = function(x) {

      paste0(
        x,
        "%"
      )
    },

    expand = expansion(
      mult = c(
        0,
        0
      )
    )
  ) +

  coord_cartesian(
    ylim = c(
      0,
      100
    )
  ) +

  labs(
    x = NULL,
    y = "Sylph sequence abundance (%)",
    title = "B  Sylph read-based abundance"
  ) +

  theme_classic(
    base_size = 12
  ) +

  theme(

    plot.title = element_text(
      size = 14,
      face = "plain"
    ),

    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),

    axis.ticks.x = element_blank(),

    legend.position = "none"
  )


## ------------------------------------------------------------
## Combine
##
## Panel A provides the common legend.
## ------------------------------------------------------------

combined <- (
  p_mag
  /
  p_sylph
) +

  plot_layout(
    heights = c(
      1,
      1
    )
  )


## ------------------------------------------------------------
## Save
## ------------------------------------------------------------

pdf_file <- file.path(
  figure_dir,
  "SemiBin2_MAG_coverage_vs_Sylph_abundance.pdf"
)

png_file <- file.path(
  figure_dir,
  "SemiBin2_MAG_coverage_vs_Sylph_abundance.png"
)


ggsave(
  pdf_file,
  combined,
  width = 11,
  height = 7.5
)


ggsave(
  png_file,
  combined,
  width = 11,
  height = 7.5,
  dpi = 400
)


print(
  combined
)


cat(
  "\nSaved:\n",
  pdf_file,
  "\n",
  png_file,
  "\n"
)
