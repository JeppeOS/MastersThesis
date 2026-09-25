#!/usr/bin/env Rscript


## ================================================================== ##
## STAGE 79
##
## FINAL CLUSTER_00121 NEIGHBORHOOD GENE MAP
##
## Input:
##
##   Cluster_00121_combined_neighborhood_genes_local_families.tsv
##   C121_local_family_summary.tsv
##   C121_shared_recurrent_architecture_families.tsv
##
## Design:
##
##   - observational unit = focal region
##   - all 23 Cluster_00121 regions retained
##   - all genes retained
##   - focal gene centered at 0 bp
##   - all neighborhoods reoriented so focal gene points right
##   - C121_local_001 highlighted as the focal 9-heme family
##   - recurrent neighboring families occurring in >=5 regions colored
##   - all other genes shown in grey
##   - MAG and GlobDB regions separated
##   - likely contig-censored ends marked with "x"
##
## IMPORTANT:
##
## The recurrent neighboring families are NOT interpreted as a
## universal operon/cassette. Their prevalence and positional
## variability are retained explicitly.
## ================================================================== ##


suppressPackageStartupMessages({

  library(tidyverse)

  library(gggenes)

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


old_project <- file.path(
  dirname(project),
  "genome_analysis_workflow"
)


indir <- file.path(
  project,
  "comparative_analysis",
  "Cluster_00121"
)


local_dir <- file.path(
  indir,
  "local_clustering"
)


outdir <- file.path(
  indir,
  "final_neighborhood_figure"
)


dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)


neighborhood_file <- file.path(
  indir,
  "Cluster_00121_combined_neighborhood_genes_local_families.tsv"
)


family_summary_file <- file.path(
  local_dir,
  "C121_local_family_summary.tsv"
)


shared_file <- file.path(
  local_dir,
  "C121_shared_recurrent_architecture_families.tsv"
)


taxonomy_file <- file.path(
  old_project,
  "16_visualization",
  "globdb_r226_taxonomy.tsv"
)


## ================================================================== ##
## Output
## ================================================================== ##

figure_png <- file.path(
  outdir,
  "Cluster_00121_neighborhood_gene_map.png"
)


figure_pdf <- file.path(
  outdir,
  "Cluster_00121_neighborhood_gene_map.pdf"
)


plot_data_file <- file.path(
  outdir,
  "Cluster_00121_neighborhood_gene_map_data.tsv"
)


highlight_file <- file.path(
  outdir,
  "Cluster_00121_highlighted_local_families.tsv"
)


censor_file <- file.path(
  outdir,
  "Cluster_00121_neighborhood_censoring_qc.tsv"
)


region_file <- file.path(
  outdir,
  "Cluster_00121_region_plot_order.tsv"
)


## ================================================================== ##
## Constants
## ================================================================== ##

FOCAL_FAMILY <- "C121_local_001"


## Neighborhood extraction originally used:
##
##   ±20 genes OR CDS intersecting focal ±20 kb.
##
## These thresholds are also used below to infer whether a displayed
## side terminated before the expected observation window.

GENE_RADIUS <- 20

BP_RADIUS <- 20000


## ================================================================== ##
## Input checks
## ================================================================== ##

required_files <- c(
  neighborhood_file,
  family_summary_file,
  shared_file
)


missing_files <- required_files[
  !file.exists(required_files)
]


if (length(missing_files) > 0) {

  stop(
    paste(
      "Missing input file(s):\n",
      paste(
        missing_files,
        collapse = "\n"
      )
    )
  )
}


## ================================================================== ##
## Read data
## ================================================================== ##

genes <- read_tsv(
  neighborhood_file,
  show_col_types = FALSE
)


families <- read_tsv(
  family_summary_file,
  show_col_types = FALSE
)


shared <- read_tsv(
  shared_file,
  show_col_types = FALSE
)


## ================================================================== ##
## Basic QC
## ================================================================== ##

required_gene_columns <- c(
  "source_dataset",
  "focal_id",
  "genome",
  "focal_protein_id",
  "focal_start",
  "focal_end",
  "focal_strand",
  "neighbor_protein_id",
  "neighbor_start",
  "neighbor_end",
  "neighbor_strand",
  "oriented_gene_offset",
  "oriented_midpoint_offset_bp",
  "oriented_neighbor_strand",
  "is_focal",
  "local_family_id"
)


missing_columns <- setdiff(
  required_gene_columns,
  names(genes)
)


if (length(missing_columns) > 0) {

  stop(
    paste(
      "Neighborhood table missing columns:",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  )
}


n_regions <- n_distinct(
  genes$focal_id
)


if (n_regions != 23) {

  stop(
    paste0(
      "Expected 23 Cluster_00121 focal regions; found ",
      n_regions,
      "."
    )
  )
}


focal_qc <- genes %>%
  filter(
    is_focal == 1
  ) %>%
  count(
    focal_id
  )


if (
  nrow(focal_qc) != 23 ||
  any(focal_qc$n != 1)
) {

  stop(
    "Each of the 23 regions must contain exactly one focal gene."
  )
}


focal_family_qc <- genes %>%
  filter(
    is_focal == 1
  ) %>%
  distinct(
    focal_id,
    local_family_id
  )


if (
  any(
    focal_family_qc$local_family_id !=
      FOCAL_FAMILY
  )
) {

  stop(
    "Not all focal genes map to C121_local_001."
  )
}


## ================================================================== ##
## Determine highlighted families
##
## Data-driven rule:
##
##   - focal family always highlighted
##   - non-focal shared recurrent families occurring in >=5/23 regions
##
## In the present dataset this corresponds to:
##
##   C121_local_001
##   C121_local_002
##   C121_local_003
##   C121_local_004
##   C121_local_005
##   C121_local_006
##
## Everything else remains visible but grey.
## ================================================================== ##

highlighted <- shared %>%
  filter(
    local_family_id == FOCAL_FAMILY |
      (
        is_focal_family == 0 &
          n_regions >= 5
      )
  ) %>%
  arrange(
    desc(
      is_focal_family
    ),
    local_family_id
  )


highlight_ids <- highlighted$local_family_id


if (
  !FOCAL_FAMILY %in%
    highlight_ids
) {

  stop(
    "C121_local_001 was not recovered as a highlighted family."
  )
}


## ================================================================== ##
## Concise biological labels
##
## Do NOT use the transferred COG1271/AppC annotation for
## Cluster_00121 itself because that transfer is incongruent with
## the conserved 9-heme architecture.
## ================================================================== ##

family_labels <- c(

  "C121_local_001" =
    "Cluster_00121\n9-heme cytochrome",

  "C121_local_002" =
    "AraJ-like\nMFS permease",

  "C121_local_003" =
    "Recurrent\nunknown",

  "C121_local_004" =
    "GspE/PilB-like\nATPase",

  "C121_local_005" =
    "NtrC-family\nregulator",

  "C121_local_006" =
    "Recurrent\nunknown"
)


highlighted <- highlighted %>%
  mutate(

    figure_label = case_when(

      local_family_id %in%
        names(family_labels) ~

        family_labels[
          local_family_id
        ],

      TRUE ~
        local_family_id
    )
  )


write_tsv(
  highlighted,
  highlight_file
)


## ================================================================== ##
## GlobDB taxonomy
## ================================================================== ##

normalize_genome_id <- function(x) {

  x %>%
    str_remove(
      "^(GB_|RS_)"
    ) %>%
    str_remove(
      "\\.\\d+$"
    )
}


best_taxon <- function(x) {

  if (
    is.na(x) ||
    x == ""
  ) {

    return(
      NA_character_
    )
  }


  pieces <- str_split(
    x,
    ";",
    simplify = FALSE
  )[[1]]


  pieces <- str_trim(
    pieces
  )


  get_rank <- function(prefix) {

    hit <- pieces[
      str_starts(
        pieces,
        prefix
      )
    ]


    if (
      length(hit) == 0
    ) {

      return(
        NA_character_
      )
    }


    out <- str_remove(
      hit[1],
      fixed(prefix)
    )


    if (
      is.na(out) ||
      out == ""
    ) {

      return(
        NA_character_
      )
    }


    out
  }


  species <- get_rank(
    "s__"
  )


  genus <- get_rank(
    "g__"
  )


  family <- get_rank(
    "f__"
  )


  order <- get_rank(
    "o__"
  )


  phylum <- get_rank(
    "p__"
  )


  for (
    candidate
    in
    c(
      species,
      genus,
      family,
      order,
      phylum
    )
  ) {

    if (
      !is.na(candidate) &&
      candidate != ""
    ) {

      return(
        candidate
      )
    }
  }


  NA_character_
}


taxonomy <- tibble(
  genome_norm = character(),
  taxon = character()
)


if (
  file.exists(
    taxonomy_file
  )
) {

  first_line <- readLines(
    taxonomy_file,
    n = 1
  )


  first_fields <- str_split(
    first_line,
    "\t"
  )[[1]]


  header_like <- any(
    str_to_lower(
      first_fields
    ) %in%
      c(
        "genome",
        "genome_id",
        "accession",
        "taxonomy",
        "gtdb_taxonomy",
        "classification"
      )
  )


  if (
    header_like
  ) {

    tx <- read_tsv(
      taxonomy_file,
      show_col_types = FALSE
    )


    lower_names <- str_to_lower(
      names(tx)
    )


    id_i <- which(
      lower_names %in%
        c(
          "genome",
          "genome_id",
          "accession",
          "id"
        )
    )[1]


    tax_i <- which(
      lower_names %in%
        c(
          "taxonomy",
          "gtdb_taxonomy",
          "classification"
        )
    )[1]


    if (
      !is.na(id_i) &&
      !is.na(tax_i)
    ) {

      taxonomy <- tibble(

        genome_norm =
          normalize_genome_id(
            tx[[id_i]]
          ),

        taxonomy_string =
          tx[[tax_i]]
      )
    }

  } else {

    tx <- read_tsv(
      taxonomy_file,
      col_names = FALSE,
      show_col_types = FALSE
    )


    taxonomy <- tibble(

      genome_norm =
        normalize_genome_id(
          tx[[1]]
        ),

      taxonomy_string =
        tx[[2]]
    )
  }


  if (
    "taxonomy_string" %in%
      names(taxonomy)
  ) {

    taxonomy <- taxonomy %>%
      mutate(

        taxon =
          map_chr(
            taxonomy_string,
            best_taxon
          )
      ) %>%
      select(
        genome_norm,
        taxon
      ) %>%
      distinct()
  }
}


## ================================================================== ##
## Region metadata
## ================================================================== ##

regions <- genes %>%
  filter(
    is_focal == 1
  ) %>%
  transmute(

    focal_id,

    source_dataset,

    genome,

    focal_protein_id,

    focal_strand
  ) %>%
  distinct()


regions <- regions %>%
  mutate(

    genome_norm =
      normalize_genome_id(
        genome
      )
  ) %>%
  left_join(
    taxonomy,
    by = "genome_norm"
  )


## ================================================================== ##
## Explicit MAG labels
##
## These four are the Cluster_00121 MAG carriers established in the
## metagenome analysis.
## ================================================================== ##

mag_taxon <- c(

  "contig_25_1899" =
    "Methylumidiphilus",

  "contig_1431_201" =
    "Methylumidiphilus",

  "contig_122_3517" =
    "Methylovulum",

  "contig_871_20" =
    "Methylomonas"
)


regions <- regions %>%
  mutate(

    source_group = case_when(

      str_detect(
        str_to_lower(
          source_dataset
        ),
        "meta|mag"
      ) ~
        "Metagenome MAGs",

      TRUE ~
        "GlobDB genomes"
    ),


    display_taxon = case_when(

      focal_protein_id %in%
        names(
          mag_taxon
        ) ~

        unname(
          mag_taxon[
            focal_protein_id
          ]
        ),

      !is.na(
        taxon
      ) &
        taxon != "" ~

        taxon,

      TRUE ~
        genome
    ),


    row_label = case_when(

      source_group ==
        "Metagenome MAGs" ~

        paste0(
          display_taxon,
          " [",
          focal_id,
          "]"
        ),

      TRUE ~

        paste0(
          display_taxon,
          " [",
          genome,
          "]"
        )
    )
  )


## ================================================================== ##
## Convert genomic coordinates into focal-oriented coordinates
##
## Focal midpoint = 0.
##
## If focal strand is "-", genomic coordinates are reversed so every
## Cluster_00121 focal gene points to the right in the figure.
## ================================================================== ##

plot_genes <- genes %>%
  mutate(

    focal_midpoint =
      (
        as.numeric(
          focal_start
        ) +
          as.numeric(
            focal_end
          )
      ) /
        2,


    neighbor_start =
      as.numeric(
        neighbor_start
      ),


    neighbor_end =
      as.numeric(
        neighbor_end
      ),


    x1_bp = case_when(

      focal_strand == "+" ~

        neighbor_start -
        focal_midpoint,

      focal_strand == "-" ~

        focal_midpoint -
        neighbor_end
    ),


    x2_bp = case_when(

      focal_strand == "+" ~

        neighbor_end -
        focal_midpoint,

      focal_strand == "-" ~

        focal_midpoint -
        neighbor_start
    ),


    xmin_bp =
      pmin(
        x1_bp,
        x2_bp
      ),


    xmax_bp =
      pmax(
        x1_bp,
        x2_bp
      ),


    xmin_kb =
      xmin_bp /
      1000,


    xmax_kb =
      xmax_bp /
      1000,


    forward =
      oriented_neighbor_strand ==
      "+",


    plot_family = case_when(

      local_family_id %in%
        highlight_ids ~

        local_family_id,

      TRUE ~
        "Other genes"
    )
  ) %>%
  left_join(

    regions %>%
      select(
        focal_id,
        source_group,
        row_label
      ),

    by = "focal_id"
  )


## ================================================================== ##
## Count highlighted neighboring families per region
##
## Used only to order rows so regions carrying more of the recurrent
## architecture are visually adjacent.
## ================================================================== ##

region_order <- plot_genes %>%
  filter(
    local_family_id %in%
      highlight_ids
  ) %>%
  distinct(
    focal_id,
    local_family_id
  ) %>%
  count(
    focal_id,
    name = "n_highlighted_families"
  ) %>%
  right_join(
    regions,
    by = "focal_id"
  ) %>%
  mutate(

    n_highlighted_families =
      replace_na(
        n_highlighted_families,
        0L
      )
  ) %>%
  arrange(

    factor(
      source_group,
      levels = c(
        "Metagenome MAGs",
        "GlobDB genomes"
      )
    ),

    desc(
      n_highlighted_families
    ),

    display_taxon,

    genome
  )


## ================================================================== ##
## Factor ordering
##
## ggplot draws the first factor level at the bottom, hence rev().
## ================================================================== ##

row_levels <- rev(
  region_order$row_label
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
        levels = c(
          "Metagenome MAGs",
          "GlobDB genomes"
        )
      )
  )


region_order <- region_order %>%
  mutate(

    row_label =
      factor(
        row_label,
        levels =
          row_levels
      )
  )


write_tsv(
  region_order,
  region_file
)


## ================================================================== ##
## Infer contig censoring
##
## Extraction was ±20 genes OR ±20 kb.
##
## Therefore, if one side stops before BOTH:
##
##   20 genes
##   AND
##   20 kb
##
## the most plausible reason is that the assembled contig ended.
##
## This is explicitly treated as an inferred censoring flag rather
## than a directly measured contig coordinate.
## ================================================================== ##

censoring <- plot_genes %>%
  group_by(
    focal_id,
    source_group,
    row_label
  ) %>%
  summarise(

    left_gene_reach =
      abs(
        min(
          oriented_gene_offset,
          na.rm = TRUE
        )
      ),

    right_gene_reach =
      max(
        oriented_gene_offset,
        na.rm = TRUE
      ),


    left_bp_reach =
      abs(
        min(
          oriented_midpoint_offset_bp,
          na.rm = TRUE
        )
      ),

    right_bp_reach =
      max(
        oriented_midpoint_offset_bp,
        na.rm = TRUE
      ),


    observed_xmin =
      min(
        xmin_kb,
        na.rm = TRUE
      ),

    observed_xmax =
      max(
        xmax_kb,
        na.rm = TRUE
      ),

    .groups = "drop"
  ) %>%
  mutate(

    left_censored =
      left_gene_reach <
        GENE_RADIUS &
      left_bp_reach <
        BP_RADIUS,


    right_censored =
      right_gene_reach <
        GENE_RADIUS &
      right_bp_reach <
        BP_RADIUS
  )


write_tsv(
  censoring,
  censor_file
)


censor_points <- bind_rows(

  censoring %>%
    filter(
      left_censored
    ) %>%
    transmute(

      focal_id,

      source_group,

      row_label,

      side =
        "left",

      x =
        observed_xmin -
        0.35
    ),


  censoring %>%
    filter(
      right_censored
    ) %>%
    transmute(

      focal_id,

      source_group,

      row_label,

      side =
        "right",

      x =
        observed_xmax +
        0.35
    )
)


## ================================================================== ##
## Save complete plotting table
## ================================================================== ##

write_tsv(
  plot_genes,
  plot_data_file
)


## ================================================================== ##
## Colors
##
## Focal family is deliberately visually dominant.
## Other recurrent families use a color-blind-friendly palette.
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

  "Other genes" =
    "#D9D9D9"
)


## Any unexpected highlighted family gets a neutral dark tone rather
## than silently disappearing.

missing_colors <- setdiff(
  highlight_ids,
  names(
    family_colors
  )
)


if (
  length(
    missing_colors
  ) >
    0
) {

  extra_colors <- rep(
    "#7F7F7F",
    length(
      missing_colors
    )
  )


  names(
    extra_colors
  ) <- missing_colors


  family_colors <- c(
    family_colors,
    extra_colors
  )
}


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
    nrow(hit) ==
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
  "Other genes"
] <-
  "Other local families"


## ================================================================== ##
## Symmetric x-axis
## ================================================================== ##

max_abs_kb <- max(
  abs(
    c(
      plot_genes$xmin_kb,
      plot_genes$xmax_kb
    )
  ),
  na.rm = TRUE
)


x_limit <- ceiling(
  max_abs_kb /
    5
) *
  5


## ================================================================== ##
## Plot
## ================================================================== ##

p <- ggplot(
  plot_genes,
  aes(
    xmin = xmin_kb,
    xmax = xmax_kb,
    y = row_label,
    fill = plot_family,
    forward = forward
  )
) +

  ## Focal center
  geom_vline(
    xintercept = 0,
    linewidth = 0.35,
    linetype = "dashed",
    color = "grey55"
  ) +

  ## Genes
  geom_gene_arrow(
    colour = "grey25",
    linewidth = 0.18,
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
      )
  ) +

  ## Inferred contig-censored ends
  geom_point(
    data = censor_points,
    aes(
      x = x,
      y = row_label
    ),
    inherit.aes = FALSE,
    shape = 4,
    stroke = 0.8,
    size = 2.2,
    color = "black"
  ) +

  facet_grid(
    source_group ~ .,
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +

  scale_fill_manual(
    values =
      family_colors,
    breaks =
      c(
        highlight_ids,
        "Other genes"
      ),
    labels =
      legend_labels[
        c(
          highlight_ids,
          "Other genes"
        )
      ],
    drop = FALSE
  ) +

  scale_x_continuous(
    limits =
      c(
        -x_limit,
        x_limit
      ),
    breaks =
      pretty_breaks(
        n = 9
      ),
    expand =
      expansion(
        mult = c(
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
      "Local protein family"
  ) +

  guides(
    fill =
      guide_legend(
        ncol = 2,
        byrow = TRUE,
        override.aes =
          list(
            colour =
              "grey25"
          )
      )
  ) +

  theme_bw(
    base_size = 10
  ) +

  theme(

    panel.grid.major.y =
      element_blank(),

    panel.grid.minor =
      element_blank(),

    panel.grid.major.x =
      element_line(
        linewidth = 0.25,
        color = "grey90"
      ),

    axis.text.y =
      element_text(
        size = 7.4,
        color = "black"
      ),

    axis.text.x =
      element_text(
        size = 8.5,
        color = "black"
      ),

    axis.title.x =
      element_text(
        size = 10,
        margin =
          margin(
            t = 8
          )
      ),

    strip.background =
      element_rect(
        fill = "grey95",
        color = "grey40"
      ),

    strip.text.y.left =
      element_text(
        angle = 0,
        face = "bold",
        size = 9
      ),

    legend.position =
      "bottom",

    legend.title =
      element_text(
        face = "bold"
      ),

    legend.text =
      element_text(
        size = 8
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
    13.5,
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
    13.5,
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
    collapse = ""
  ),
  "\n",
  sep = ""
)


cat(
  "STAGE 79 COMPLETE - CLUSTER_00121 NEIGHBORHOOD FIGURE\n"
)


cat(
  paste(
    rep(
      "=",
      80
    ),
    collapse = ""
  ),
  "\n\n",
  sep = ""
)


cat(
  "Focal regions:            ",
  n_distinct(
    plot_genes$focal_id
  ),
  "\n",
  sep = ""
)


cat(
  "Neighborhood genes:       ",
  nrow(
    plot_genes
  ),
  "\n",
  sep = ""
)


cat(
  "Local families:           ",
  n_distinct(
    plot_genes$local_family_id
  ),
  "\n",
  sep = ""
)


cat(
  "Highlighted families:     ",
  length(
    highlight_ids
  ),
  "\n",
  sep = ""
)


cat(
  "Inferred left-censored:   ",
  sum(
    censoring$left_censored
  ),
  "\n",
  sep = ""
)


cat(
  "Inferred right-censored:  ",
  sum(
    censoring$right_censored
  ),
  "\n",
  sep = ""
)


cat(
  "Displayed x-range:        ±",
  x_limit,
  " kb\n\n",
  sep = ""
)


cat(
  "Highlighted families:\n"
)


print(
  highlighted %>%
    select(
      local_family_id,
      n_regions,
      prevalence_all_23_position,
      median_oriented_gene_offset,
      figure_label
    )
)


cat(
  "\nPNG:\n  ",
  figure_png,
  "\n",
  sep = ""
)


cat(
  "\nPDF:\n  ",
  figure_pdf,
  "\n",
  sep = ""
)


cat(
  "\nPlot data:\n  ",
  plot_data_file,
  "\n",
  sep = ""
)


cat(
  "\nCensoring QC:\n  ",
  censor_file,
  "\n",
  sep = ""
)
