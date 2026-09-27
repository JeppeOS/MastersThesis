#!/usr/bin/env Rscript


## ================================================================== ##
## STAGE 79B
##
## FINAL Cluster_00121 NEIGHBORHOOD FIGURE
##
## Displays:
##
##   - all 4 metagenome MAG loci
##   - upper 6 GlobDB loci from Stage-79 ordering
##   - recurrent Cluster_00121 local families
##   - all genes
##   - subtle grey neighborhood / contig lines
##
## Additional heme-aware annotation:
##
##   MOTU40_034781 / UBA6136 contains immediately downstream:
##
##      +2  Cluster_00069
##      +3  Cluster_00091
##      +4  Cluster_00222
##      +5  Cluster_00050 / Cyc2
##
##   C00069/C00091/C00222 share one color as:
##
##      "Module 10-associated 1-heme cytochromes"
##
##   Cyc2 / Cluster_00050 receives its own color.
##
##   Other non-focal FindMeHemes-positive proteins are shown with a
##   separate generic color but are not assigned specific functional
##   interpretations.
##
## Only the four UBA6136 Module-10-associated genes receive gene-level
## labels with leader lines.
## ================================================================== ##


suppressPackageStartupMessages({

  library(tidyverse)
  library(gggenes)
  library(ggrepel)
  library(scales)

})


## ================================================================== ##
## Paths
## ================================================================== ##

project <- file.path(
  Sys.getenv("HOME"),
  "methanotrophs",
  "methanotroph_project",
  "jeppe",
  "metagenome"
)


outdir <- file.path(
  project,
  "comparative_analysis",
  "Cluster_00121",
  "final_neighborhood_figure"
)


plot_data_file <- file.path(
  outdir,
  "Cluster_00121_neighborhood_gene_map_data.tsv"
)


highlight_file <- file.path(
  outdir,
  "Cluster_00121_highlighted_local_families.tsv"
)


region_file <- file.path(
  outdir,
  "Cluster_00121_region_plot_order.tsv"
)


figure_png <- file.path(
  outdir,
  "Cluster_00121_neighborhood_gene_map_FINAL.png"
)


figure_pdf <- file.path(
  outdir,
  "Cluster_00121_neighborhood_gene_map_FINAL.pdf"
)


selected_regions_file <- file.path(
  outdir,
  "Cluster_00121_final_displayed_regions.tsv"
)


## ================================================================== ##
## Settings
## ================================================================== ##

N_GLOBDB <- 6


FOCAL_FAMILY <- "C121_local_001"


UBA_FOCAL_ID <- "GlobDB_C121_17"


UBA_MODULE10_CYTOCHROMES <- c(
  "Cluster_00069",
  "Cluster_00091",
  "Cluster_00222"
)


UBA_CYC2_CLUSTER <- "Cluster_00050"


## ================================================================== ##
## Read Stage-79 outputs
##
## Read plot data as character first to avoid vroom type-guessing
## warnings caused by mixed annotation columns.
## ================================================================== ##

plot_genes <- read_tsv(
  plot_data_file,
  col_types =
    cols(
      .default =
        col_character()
    ),
  show_col_types =
    FALSE
)


highlighted <- read_tsv(
  highlight_file,
  show_col_types =
    FALSE
)


region_order <- read_tsv(
  region_file,
  show_col_types =
    FALSE
)


## ================================================================== ##
## Convert plotting columns
## ================================================================== ##

plot_genes <- plot_genes %>%
  mutate(

    xmin_kb =
      as.numeric(
        xmin_kb
      ),

    xmax_kb =
      as.numeric(
        xmax_kb
      ),

    forward =
      as.logical(
        forward
      ),

    is_focal =
      as.numeric(
        is_focal
      ),

    number_of_hemes =
      as.numeric(
        number_of_hemes
      ),

    findmehemes_positive =
      as.numeric(
        findmehemes_positive
      ),

    old_global_mmseq_cluster =
      replace_na(
        old_global_mmseq_cluster,
        ""
      )
  )


## ================================================================== ##
## QC
## ================================================================== ##

required_columns <- c(

  "focal_id",

  "source_group",

  "row_label",

  "neighbor_protein_id",

  "xmin_kb",

  "xmax_kb",

  "forward",

  "local_family_id",

  "is_focal",

  "number_of_hemes",

  "findmehemes_positive",

  "old_global_mmseq_cluster"
)


missing_columns <- setdiff(
  required_columns,
  names(
    plot_genes
  )
)


if (
  length(
    missing_columns
  ) >
    0
) {

  stop(
    "Stage-79 plot table is missing required columns: ",
    paste(
      missing_columns,
      collapse = ", "
    )
  )
}


## ================================================================== ##
## Select final displayed regions
##
##   4 MAGs
##   upper 6 GlobDB regions
## ================================================================== ##

mag_regions <- region_order %>%
  filter(
    source_group ==
      "Metagenome MAGs"
  )


globdb_regions <- region_order %>%
  filter(
    source_group ==
      "GlobDB genomes"
  ) %>%
  slice_head(
    n =
      N_GLOBDB
  )


selected_regions <- bind_rows(
  mag_regions,
  globdb_regions
)


plot_genes <- plot_genes %>%
  filter(
    focal_id %in%
      selected_regions$focal_id
  )


write_tsv(
  selected_regions,
  selected_regions_file
)


## ================================================================== ##
## Row ordering
## ================================================================== ##

row_levels <- rev(
  selected_regions$row_label
)


plot_genes <- plot_genes %>%
  mutate(

    row_label =
      factor(
        row_label,
        levels =
          row_levels
      ),

    source_group =
      factor(
        source_group,
        levels =
          c(
            "Metagenome MAGs",
            "GlobDB genomes"
          )
      )
  )


## ================================================================== ##
## Subtle neighborhood lines
## ================================================================== ##

contig_lines <- plot_genes %>%
  group_by(
    focal_id,
    source_group,
    row_label
  ) %>%
  summarise(

    xmin =
      min(
        xmin_kb,
        na.rm =
          TRUE
      ),

    xmax =
      max(
        xmax_kb,
        na.rm =
          TRUE
      ),

    .groups =
      "drop"
  )


## ================================================================== ##
## Recurrent Cluster_00121 families
## ================================================================== ##

highlight_ids <-
  highlighted$local_family_id


## ================================================================== ##
## Final display classification
##
## Priority:
##
##   1. Cluster_00121 focal protein
##   2. special UBA6136 Module-10-associated 1-heme cytochromes
##   3. special UBA6136 Cyc2 / Cluster_00050
##   4. recurrent C121 local families
##   5. other non-focal FindMeHemes-positive proteins
##   6. background genes
##
## The focal family has highest priority so Cluster_00121 itself is
## never accidentally recolored as a generic heme protein.
## ================================================================== ##

plot_genes <- plot_genes %>%
  mutate(

    final_display_class = case_when(

      is_focal ==
        1 |

      local_family_id ==
        FOCAL_FAMILY ~

        FOCAL_FAMILY,


      focal_id ==
        UBA_FOCAL_ID &

      old_global_mmseq_cluster %in%
        UBA_MODULE10_CYTOCHROMES ~

        "UBA_Module10_1heme",


      focal_id ==
        UBA_FOCAL_ID &

      old_global_mmseq_cluster ==
        UBA_CYC2_CLUSTER ~

        "UBA_Cyc2",


      local_family_id %in%
        highlight_ids ~

        local_family_id,


      is_focal ==
        0 &

      (
        coalesce(
          findmehemes_positive,
          0
        ) ==
          1 |

        coalesce(
          number_of_hemes,
          0
        ) >
          0
      ) ~

        "Other_heme_positive",


      TRUE ~

        "Other genes"
    )
  )


## ================================================================== ##
## Special UBA6136 labels
## ================================================================== ##

special_cytochrome_labels <- plot_genes %>%
  filter(

    focal_id ==
      UBA_FOCAL_ID,

    old_global_mmseq_cluster %in%
      c(
        UBA_MODULE10_CYTOCHROMES,
        UBA_CYC2_CLUSTER
      )
  ) %>%
  mutate(

    special_label =
      recode(

        old_global_mmseq_cluster,

        "Cluster_00069" =
          "C00069",

        "Cluster_00091" =
          "C00091",

        "Cluster_00222" =
          "C00222",

        "Cluster_00050" =
          "Cyc2"
      ),

    gene_midpoint_kb =
      (
        xmin_kb +
        xmax_kb
      ) /
        2
  ) %>%
  arrange(
    gene_midpoint_kb
  )


if (
  nrow(
    special_cytochrome_labels
  ) !=
    4
) {

  stop(
    "Expected four Module-10-associated genes in UBA6136; found ",
    nrow(
      special_cytochrome_labels
    ),
    "."
  )
}


## ================================================================== ##
## Colors
##
## Existing Cluster_00121 colors are unchanged.
##
## New colors:
##
##   dark magenta = Module-10-associated 1-heme cytochromes
##   indigo       = Cyc2 / Cluster_00050
##   olive        = other FindMeHemes-positive proteins
## ================================================================== ##

family_colors <- c(

  "C121_local_001" =
    "#D55E00",

  "C121_local_002" =
    "#0072B2",

  "C121_local_003" =
    "#009E73",

  "C121_local_004" =
    "#CC79A7",

  "C121_local_005" =
    "#E69F00",

  "C121_local_006" =
    "#56B4E9",


  "UBA_Module10_1heme" =
    "#882255",

  "UBA_Cyc2" =
    "#332288",

  "Other_heme_positive" =
    "#999933",


  "Other genes" =
    "#D9D9D9"
)


## ================================================================== ##
## Legend labels
## ================================================================== ##

legend_labels <- c()


for (
  fam
  in
  highlight_ids
) {

  hit <- highlighted %>%
    filter(
      local_family_id ==
        fam
    )


  if (
    nrow(
      hit
    ) ==
      1
  ) {

    legend_labels[
      fam
    ] <-
      paste0(
        hit$figure_label,
        "  (",
        hit$n_regions,
        "/23)"
      )

  } else {

    legend_labels[
      fam
    ] <-
      fam
  }
}


legend_labels[
  "UBA_Module10_1heme"
] <-
  "Module 10-associated 1-heme cytochromes\n(C00069 / C00091 / C00222)"


legend_labels[
  "UBA_Cyc2"
] <-
  "Cyc2 (Cluster_00050)"


legend_labels[
  "Other_heme_positive"
] <-
  "Other FindMeHemes-positive protein"


legend_labels[
  "Other genes"
] <-
  "Other local families"


legend_breaks <- c(

  highlight_ids,

  "UBA_Module10_1heme",

  "UBA_Cyc2",

  "Other_heme_positive",

  "Other genes"
)


## ================================================================== ##
## X-axis
## ================================================================== ##

max_abs_kb <- max(

  abs(
    c(
      plot_genes$xmin_kb,
      plot_genes$xmax_kb
    )
  ),

  na.rm =
    TRUE
)


x_limit <- ceiling(
  max_abs_kb /
    5
) *
  5


## ================================================================== ##
## Plot
## ================================================================== ##

p <- ggplot() +

  ## Focal midpoint

  geom_vline(
    xintercept =
      0,
    linewidth =
      0.35,
    linetype =
      "dashed",
    color =
      "grey55"
  ) +


  ## Subtle contig / neighborhood line

  geom_segment(
    data =
      contig_lines,
    aes(
      x =
        xmin,
      xend =
        xmax,
      y =
        row_label,
      yend =
        row_label
    ),
    inherit.aes =
      FALSE,
    linewidth =
      0.45,
    color =
      "grey72",
    lineend =
      "round"
  ) +


  ## Genes

  geom_gene_arrow(
    data =
      plot_genes,
    aes(
      xmin =
        xmin_kb,
      xmax =
        xmax_kb,
      y =
        row_label,
      fill =
        final_display_class,
      forward =
        forward
    ),
    arrowhead_height =
      grid::unit(
        3.0,
        "mm"
      ),
    arrowhead_width =
      grid::unit(
        1.4,
        "mm"
      ),
    arrow_body_height =
      grid::unit(
        3.0,
        "mm"
      ),
    colour =
      "grey25",
    linewidth =
      0.18
  ) +


  ## UBA6136 Module-10-style labels

  geom_text_repel(
    data =
      special_cytochrome_labels,
    aes(
      x =
        gene_midpoint_kb,
      y =
        row_label,
      label =
        special_label
    ),
    inherit.aes =
      FALSE,
    direction =
      "x",
    nudge_y =
      0.34,
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


  ## MAG / GlobDB blocks

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
      family_colors,
    breaks =
      legend_breaks,
    labels =
      legend_labels[
        legend_breaks
      ],
    drop =
      FALSE
  ) +


  scale_x_continuous(
    limits =
      c(
        -x_limit,
        x_limit
      ),
    breaks =
      pretty_breaks(
        n =
          9
      ),
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
      "Position relative to Cluster_00121 focal midpoint (kb)",
    y =
      NULL,
    fill =
      "Protein family / annotation"
  ) +


  guides(
    fill =
      guide_legend(
        ncol =
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


  theme_bw(
    base_size =
      10
  ) +

  theme(

    panel.grid.major.y =
      element_blank(),

    panel.grid.minor =
      element_blank(),

    panel.grid.major.x =
      element_line(
        linewidth =
          0.25,
        color =
          "grey90"
      ),

    axis.text.y =
      element_text(
        size =
          7.4,
        color =
          "black"
      ),

    axis.text.x =
      element_text(
        size =
          8.5,
        color =
          "black"
      ),

    axis.title.x =
      element_text(
        size =
          10,
        margin =
          margin(
            t =
              8
          )
      ),

    strip.background =
      element_rect(
        fill =
          "grey95",
        color =
          "grey40"
      ),

    strip.text.y.left =
      element_text(
        angle =
          0,
        face =
          "bold",
        size =
          9
      ),

    legend.position =
      "bottom",

    legend.title =
      element_text(
        face =
          "bold"
      ),

    legend.text =
      element_text(
        size =
          8
      ),

    legend.key.width =
      grid::unit(
        1.4,
        "cm"
      ),

    plot.margin =
      margin(
        10,
        15,
        10,
        10
      )
  )


## ================================================================== ##
## Save
## ================================================================== ##

ggsave(
  filename =
    figure_png,
  plot =
    p,
  width =
    15,
  height =
    8.8,
  units =
    "in",
  dpi =
    300,
  bg =
    "white"
)


ggsave(
  filename =
    figure_pdf,
  plot =
    p,
  width =
    15,
  height =
    8.8,
  units =
    "in",
  device =
    cairo_pdf,
  bg =
    "white"
)


## ================================================================== ##
## Report
## ================================================================== ##

cat(
  "\n",
  paste(
    rep(
      "=",
      80
    ),
    collapse =
      ""
  ),
  "\n",
  sep =
    ""
)


cat(
  "STAGE 79B COMPLETE - FINAL Cluster_00121 DISPLAY\n"
)


cat(
  paste(
    rep(
      "=",
      80
    ),
    collapse =
      ""
  ),
  "\n\n",
  sep =
    ""
)


cat(
  "Metagenome MAG regions: ",
  nrow(
    mag_regions
  ),
  "\n",
  sep =
    ""
)


cat(
  "GlobDB regions:         ",
  nrow(
    globdb_regions
  ),
  "\n",
  sep =
    ""
)


cat(
  "Total displayed:        ",
  nrow(
    selected_regions
  ),
  "\n\n",
  sep =
    ""
)


cat(
  "Special UBA6136 architecture:\n"
)


print(
  special_cytochrome_labels %>%
    select(
      neighbor_protein_id,
      old_global_mmseq_cluster,
      special_label,
      gene_midpoint_kb
    ),
  n =
    Inf
)


cat(
  "\nOther displayed non-focal heme-positive proteins:\n"
)


print(
  plot_genes %>%
    filter(
      final_display_class ==
        "Other_heme_positive"
    ) %>%
    select(
      focal_id,
      neighbor_protein_id,
      number_of_hemes,
      old_global_mmseq_cluster,
      comparison_cog,
      comparison_product
    ),
  n =
    Inf,
  width =
    Inf
)


cat(
  "\nPNG:\n  ",
  figure_png,
  "\n",
  sep =
    ""
)


cat(
  "\nPDF:\n  ",
  figure_pdf,
  "\n",
  sep =
    ""
)
