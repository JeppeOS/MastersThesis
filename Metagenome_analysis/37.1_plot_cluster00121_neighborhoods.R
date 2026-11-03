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

##   - BCRBG_22465 explicitly included because the genome contains

##     both Cluster_00121 and an MtoA/Module-20 locus

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

## Mto/Mtr highlighting:

##

##   - FeGenie MtoA is highlighted wherever it occurs in a displayed

##     Cluster_00121 neighborhood.

##

##   - literal MtoB, if present, is highlighted.

##

##   - FeGenie MtrB_TIGR03509 / MtrB is also highlighted as the

##     MtoB/MtrB-like outer-membrane partner used in this workflow.

##

## IMPORTANT:

##

## BCRBG_22465 contains Cluster_00121 and MtoA at the GENOME level,

## but the established proteins occur on different contigs. MtoA is

## therefore NOT artificially inserted into the Cluster_00121 locus.

## It is colored only if an MtoA/MtoB/MtrB annotation genuinely occurs

## inside the plotted Cluster_00121 neighborhood.

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



root <- file.path(

  Sys.getenv("HOME"),

  "methanotrophs",

  "methanotroph_project",

  "jeppe"

)





project <- file.path(

  root,

  "metagenome"

)





old_project <- file.path(

  root,

  "genome_analysis_workflow"

)





outdir <- file.path(

  project,

  "comparative_analysis",

  "Cluster_00121",

  "final_neighborhood_figure"

)





dir.create(

  outdir,

  recursive = TRUE,

  showWarnings = FALSE

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





## Authoritative GlobDB integrated annotation table.

## This is used to restore FeGenie HMM calls such as MtrB that may

## have been lost from intermediate neighborhood annotation tables.



globdb_master_file <- file.path(

  old_project,

  "05_integrated_annotations",

  "MASTER_PROTEIN_ANNOTATIONS.tsv"

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





mto_mtr_file <- file.path(

  outdir,

  "Cluster_00121_displayed_MtoA_MtoB_MtrB_hits.tsv"

)





## ================================================================== ##

## Settings

## ================================================================== ##



N_GLOBDB <- 6





FOCAL_FAMILY <- "C121_local_001"





## UBA6136 region used for the explicitly labelled downstream

## Module-10-associated cytochromes.



UBA_FOCAL_ID <- "GlobDB_C121_17"





UBA_MODULE10_CYTOCHROMES <- c(

  "Cluster_00069",

  "Cluster_00091",

  "Cluster_00222"

)





UBA_CYC2_CLUSTER <- "Cluster_00050"





## Explicitly requested additional GlobDB genome.



EXTRA_GLOBDB_GENOME <- "BCRBG_22465"





## ================================================================== ##

## Input checks

## ================================================================== ##



required_files <- c(

  plot_data_file,

  highlight_file,

  region_file,

  globdb_master_file

)





missing_files <- required_files[

  !file.exists(

    required_files

  )

]





if (

  length(

    missing_files

  ) >

    0

) {



  stop(

    paste0(

      "Missing required input file(s):\n",

      paste(

        missing_files,

        collapse = "\n"

      )

    )

  )

}





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

      ),



    fegenie_HMMs =

      replace_na(

        fegenie_HMMs,

        ""

      )



  )





## ================================================================== ##

## QC

## ================================================================== ##



required_columns <- c(



  "focal_id",



  "source_group",



  "genome",



  "row_label",



  "neighbor_protein_id",



  "xmin_kb",



  "xmax_kb",



  "forward",



  "local_family_id",



  "is_focal",



  "number_of_hemes",



  "findmehemes_positive",



  "fegenie_HMMs",



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

## Supplement GlobDB FeGenie annotations

##

## Earlier analysis established that some MtrB_TIGR03509 calls can be

## absent from intermediate neighborhood annotation tables.

##

## The authoritative GlobDB MASTER_PROTEIN_ANNOTATIONS.tsv is therefore

## used as an annotation supplement.

##

## Existing Stage-79 annotations are preserved for MAG proteins.

## ================================================================== ##



globdb_master <- read_tsv(

  globdb_master_file,

  show_col_types = FALSE

)





globdb_required <- c(

  "genome",

  "protein_id",

  "fegenie_HMMs"

)





missing_globdb_columns <- setdiff(

  globdb_required,

  names(

    globdb_master

  )

)





if (

  length(

    missing_globdb_columns

  ) >

    0

) {



  stop(

    "GlobDB master annotation table is missing columns: ",

    paste(

      missing_globdb_columns,

      collapse = ", "

    )

  )

}





globdb_lookup <- globdb_master %>%

  transmute(



    genome =

      as.character(

        genome

      ),



    neighbor_protein_id =

      as.character(

        protein_id

      ),



    master_fegenie_HMMs =

      coalesce(

        as.character(

          fegenie_HMMs

        ),

        ""

      )



  ) %>%

  distinct(

    genome,

    neighbor_protein_id,

    .keep_all = TRUE

  )





plot_genes <- plot_genes %>%

  left_join(

    globdb_lookup,

    by =

      c(

        "genome",

        "neighbor_protein_id"

      )

  ) %>%

  mutate(



    master_fegenie_HMMs =

      replace_na(

        master_fegenie_HMMs,

        ""

      ),



    ## Prefer the authoritative master annotation when available.

    ## Otherwise retain the Stage-79 annotation, which is especially

    ## important for metagenome MAG proteins.



    effective_fegenie_HMMs =

      case_when(



        master_fegenie_HMMs != "" ~

          master_fegenie_HMMs,



        TRUE ~

          fegenie_HMMs

      )



  )





## ================================================================== ##

## Helper functions for FeGenie HMM recognition

## ================================================================== ##



has_mtoa <- function(x) {



  str_detect(

    coalesce(

      x,

      ""

    ),

    regex(

      "(^|;)\\s*MtoA\\s*(;|$)",

      ignore_case = TRUE

    )

  )

}





has_mtob_mtrb <- function(x) {



  str_detect(

    coalesce(

      x,

      ""

    ),

    regex(

      "(^|;)\\s*(MtoB|MtrB|MtrB_TIGR03509)\\s*(;|$)",

      ignore_case = TRUE

    )

  )

}





## ================================================================== ##

## Select final displayed regions

##

##   4 MAGs

##   upper 6 GlobDB regions

##   + BCRBG_22465

## ================================================================== ##



mag_regions <- region_order %>%

  filter(

    source_group ==

      "Metagenome MAGs"

  )





## -------------------------------------------------------------- ##

## Original upper six GlobDB loci

## -------------------------------------------------------------- ##



globdb_top <- region_order %>%

  filter(

    source_group ==

      "GlobDB genomes"

  ) %>%

  slice_head(

    n =

      N_GLOBDB

  )





## -------------------------------------------------------------- ##

## Explicitly include BCRBG_22465

## -------------------------------------------------------------- ##



globdb_extra <- region_order %>%

  filter(

    source_group ==

      "GlobDB genomes",

    genome ==

      EXTRA_GLOBDB_GENOME

  )





if (

  nrow(

    globdb_extra

  ) !=

    1

) {



  stop(

    paste0(

      "Expected exactly one Cluster_00121 region for ",

      EXTRA_GLOBDB_GENOME,

      ", but found ",

      nrow(

        globdb_extra

      ),

      "."

    )

  )

}





message(

  "\nExplicitly added GlobDB region:"

)





print(

  globdb_extra %>%

    select(

      focal_id,

      genome,

      row_label

    )

)





## -------------------------------------------------------------- ##

## Combine.

##

## distinct() prevents duplication if BCRBG_22465 ever moves into

## the upper-six Stage-79 ordering.

## -------------------------------------------------------------- ##



globdb_regions <- bind_rows(

  globdb_top,

  globdb_extra

) %>%

  distinct(

    focal_id,

    .keep_all = TRUE

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

## Subtle neighborhood / contig lines

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

## Identify MtoA and MtoB/MtrB-like proteins in displayed regions

##

## These are based on FeGenie annotation, with GlobDB calls

## supplemented from the authoritative master annotation table.

## ================================================================== ##



plot_genes <- plot_genes %>%

  mutate(



    is_MtoA =

      has_mtoa(

        effective_fegenie_HMMs

      ),



    is_MtoB_MtrB =

      has_mtob_mtrb(

        effective_fegenie_HMMs

      )



  )





mto_mtr_hits <- plot_genes %>%

  filter(

    is_MtoA |

      is_MtoB_MtrB

  ) %>%

  mutate(



    redox_annotation =

      case_when(



        is_MtoA ~

          "MtoA",



        is_MtoB_MtrB ~

          "MtoB/MtrB-like",



        TRUE ~

          NA_character_

      )



  ) %>%

  select(

    focal_id,

    source_group,

    genome,

    row_label,

    neighbor_protein_id,

    effective_fegenie_HMMs,

    redox_annotation,

    xmin_kb,

    xmax_kb

  )





write_tsv(

  mto_mtr_hits,

  mto_mtr_file

)





message(

  "\nMtoA / MtoB-MtrB-like proteins inside DISPLAYED neighborhoods:"

)





if (

  nrow(

    mto_mtr_hits

  ) ==

    0

) {



  message(

    "  None detected."

  )



} else {



  print(

    mto_mtr_hits,

    n =

      Inf,

    width =

      Inf

  )

}





## ================================================================== ##

## Final display classification

##

## Priority:

##

##   1. Cluster_00121 focal protein

##   2. MtoA

##   3. MtoB/MtrB-like

##   4. special UBA6136 Module-10-associated 1-heme cytochromes

##   5. special UBA6136 Cyc2 / Cluster_00050

##   6. recurrent C121 local families

##   7. other non-focal FindMeHemes-positive proteins

##   8. background genes

##

## The focal family has highest priority so Cluster_00121 itself is

## never accidentally recolored by another annotation.

## ================================================================== ##



plot_genes <- plot_genes %>%

  mutate(



    final_display_class =

      case_when(



        is_focal ==

          1 |



        local_family_id ==

          FOCAL_FAMILY ~



          FOCAL_FAMILY,





        is_MtoA ~



          "MtoA",





        is_MtoB_MtrB ~



          "MtoB_MtrB",





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
## Module-10-style sparse gene labels
##
## The highlighted Cluster_00121 families receive compact labels.
## As in the Module-10 figure:
##   - labels are split above/below the genes;
##   - repeated copies receive only one label per locus / displayed label;
##   - ordinary background genes remain unlabeled.
##
## C121-L3 and C121-L6 are figure shorthand for the two recurrent
## unannotated local families C121_local_003 and C121_local_006.
## ================================================================== ##


plot_genes <- plot_genes %>%
  mutate(

    gene_midpoint_kb =
      (
        xmin_kb +
          xmax_kb
      ) /
        2,

    display_label_plot =
      case_when(

        final_display_class ==
          "C121_local_001" ~
          "C00121",

        final_display_class ==
          "C121_local_002" ~
          "AraJ-like",

        final_display_class ==
          "C121_local_003" ~
          "C121-L3",

        final_display_class ==
          "C121_local_004" ~
          "GspE/PilB-like",

        final_display_class ==
          "C121_local_005" ~
          "NtrC-like",

        final_display_class ==
          "C121_local_006" ~
          "C121-L6",

        final_display_class ==
          "MtoA" ~
          "MtoA",

        final_display_class ==
          "MtoB_MtrB" ~
          "MtoB/MtrB-like",

        final_display_class ==
          "UBA_Module10_1heme" &
          old_global_mmseq_cluster ==
            "Cluster_00069" ~
          "C00069",

        final_display_class ==
          "UBA_Module10_1heme" &
          old_global_mmseq_cluster ==
            "Cluster_00091" ~
          "C00091",

        final_display_class ==
          "UBA_Module10_1heme" &
          old_global_mmseq_cluster ==
            "Cluster_00222" ~
          "C00222",

        final_display_class ==
          "UBA_Cyc2" ~
          "Cyc2",

        final_display_class ==
          "Other_heme_positive" ~
          "Other heme protein",

        TRUE ~
          ""
      ),

    label_side =
      case_when(

        final_display_class %in%
          c(
            "C121_local_001",
            "C121_local_003",
            "C121_local_005",
            "UBA_Module10_1heme",
            "UBA_Cyc2"
          ) ~
          "above",

        final_display_class %in%
          c(
            "C121_local_002",
            "C121_local_004",
            "C121_local_006",
            "MtoA",
            "MtoB_MtrB",
            "Other_heme_positive"
          ) ~
          "below",

        TRUE ~
          "none"
      )
  )


## If the same displayed label occurs more than once in a locus,
## label only the copy closest to the Cluster_00121 focal protein.

labels <- plot_genes %>%
  filter(
    display_label_plot != "",
    label_side != "none"
  ) %>%
  group_by(
    focal_id,
    display_label_plot
  ) %>%
  arrange(
    abs(
      gene_midpoint_kb
    ),
    gene_midpoint_kb,
    .by_group = TRUE
  ) %>%
  slice(
    1
  ) %>%
  ungroup()


labels_above <- labels %>%
  filter(
    label_side ==
      "above"
  )


labels_below <- labels %>%
  filter(
    label_side ==
      "below"
  )


## Retain the explicit UBA6136 four-gene QC used by the original script.

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

      )

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

## Additional colors:

##

##   brown         = MtoA

##   purple        = MtoB/MtrB-like outer-membrane partner

##   dark magenta  = Module-10-associated 1-heme cytochromes

##   indigo        = Cyc2 / Cluster_00050

##   olive         = other FindMeHemes-positive proteins

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





  "MtoA" =

    "#A6761D",



  "MtoB_MtrB" =

    "#7570B3",





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

  "MtoA"

] <-

  "MtoA"





legend_labels[

  "MtoB_MtrB"

] <-

  "MtoB/MtrB-like outer-membrane partner"





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





## Only show MtoA / MtoB-MtrB legend entries if at least one such

## protein is actually present in the displayed neighborhoods.



mto_legend_breaks <- character()





if (

  any(

    plot_genes$final_display_class ==

      "MtoA",

    na.rm = TRUE

  )

) {



  mto_legend_breaks <- c(

    mto_legend_breaks,

    "MtoA"

  )

}





if (

  any(

    plot_genes$final_display_class ==

      "MtoB_MtrB",

    na.rm = TRUE

  )

) {



  mto_legend_breaks <- c(

    mto_legend_breaks,

    "MtoB_MtrB"

  )

}





legend_breaks <- c(

  highlight_ids,

  mto_legend_breaks,

  "UBA_Module10_1heme",

  "UBA_Cyc2",

  "Other_heme_positive"

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





  ## -------------------------------------------------------------- ##

  ## Subtle contig / neighborhood line

  ## -------------------------------------------------------------- ##



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





  ## -------------------------------------------------------------- ##

  ## Genes

  ## -------------------------------------------------------------- ##



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





  ## -------------------------------------------------------------- ##
  ## Module-10-style sparse labels above genes
  ## -------------------------------------------------------------- ##

  geom_text_repel(

    data =
      labels_above,

    aes(
      x =
        gene_midpoint_kb,

      y =
        row_label,

      label =
        display_label_plot
    ),

    inherit.aes =
      FALSE,

    direction =
      "x",

    nudge_y =
      0.32,

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


  ## -------------------------------------------------------------- ##
  ## Module-10-style sparse labels below genes
  ## -------------------------------------------------------------- ##

  geom_text_repel(

    data =
      labels_below,

    aes(
      x =
        gene_midpoint_kb,

      y =
        row_label,

      label =
        display_label_plot
    ),

    inherit.aes =
      FALSE,

    direction =
      "x",

    nudge_y =
      -0.32,

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


  ## -------------------------------------------------------------- ##

  ## MAG / GlobDB blocks

  ## -------------------------------------------------------------- ##



  facet_grid(



    source_group ~ .,



    scales =

      "free_y",



    space =

      "free_y",



    switch =

      "y"



  ) +





  ## -------------------------------------------------------------- ##

  ## Fill scale

  ## -------------------------------------------------------------- ##



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
      TRUE,

    name =
      NULL



  ) +





  ## -------------------------------------------------------------- ##

  ## X scale

  ## -------------------------------------------------------------- ##



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

      NULL



  ) +





  guides(



    fill =

      guide_legend(



        nrow =
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
  base_size = 10
) +

theme(

  panel.border =
    element_blank(),

  panel.background =
    element_blank(),

  axis.line.x.bottom =
    element_line(
      color = "black",
      linewidth = 0.4
    ),

  axis.line.x.top =
    element_blank(),

  axis.line.y.left =
    element_line(
      color = "black",
      linewidth = 0.4
    ),

  axis.line.y.right =
    element_blank(),

  panel.grid.major.y =
    element_blank(),

  panel.grid.minor =
    element_blank(),

  panel.grid.major.x =
    element_blank(),

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
      margin = margin(t = 8)
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

    9.6,



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

    9.6,



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

  "Explicit extra GlobDB genome:\n"

)





print(

  globdb_extra %>%

    select(

      focal_id,

      genome,

      row_label

    ),

  n =

    Inf

)





cat(

  "\nMtoA / MtoB-MtrB-like proteins in displayed neighborhoods:\n"

)





if (

  nrow(

    mto_mtr_hits

  ) ==

    0

) {



  cat(

    "None detected.\n"

  )



} else {



  print(

    mto_mtr_hits,

    n =

      Inf,

    width =

      Inf

  )

}





cat(

  "\nSpecial UBA6136 architecture:\n"

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

      any_of(

        c(

          "comparison_cog",

          "comparison_product"

        )

      )

    ),

  n =

    Inf,

  width =

    Inf

)





cat(

  "\nLabels drawn: ",

  nrow(
    labels
  ),

  "\n",

  sep =
    ""

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





cat(

  "\nDisplayed-region table:\n  ",

  selected_regions_file,

  "\n",

  sep =

    ""

)





cat(

  "\nMtoA / MtoB-MtrB neighborhood QC:\n  ",

  mto_mtr_file,

  "\n",

  sep =

    ""

)


