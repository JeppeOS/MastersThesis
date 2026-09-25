#!/usr/bin/env Rscript


suppressPackageStartupMessages({
  library(ape)
})


project <- file.path(
  Sys.getenv("HOME"),
  "methanotrophs",
  "methanotroph_project",
  "jeppe",
  "metagenome"
)


base <- file.path(
  project,
  "comparative_analysis",
  "final_species_phylogeny"
)


input_tree <- file.path(
  base,
  "04_fasttree",
  "final_methylococcales_GlobDB_MAGs_Umezawa_bac120.tree"
)


outdir <- file.path(
  base,
  "05_rooted"
)


dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)


output_tree <- file.path(
  outdir,
  "final_methylococcales_GlobDB_MAGs_Umezawa_bac120_ROOTED.tree"
)


outgroup <- "OUT_Methylophaga_nitratireducenticrescens"


tree <- read.tree(
  input_tree
)


if (
  !outgroup %in%
    tree$tip.label
) {

  stop(
    "Outgroup not present in FastTree tree: ",
    outgroup
  )
}


rooted <- root(
  tree,
  outgroup = outgroup,
  resolve.root = TRUE
)


write.tree(
  rooted,
  file = output_tree
)


## ================================================================== ##
## Root verification
## ================================================================== ##

tip_number <- which(
  rooted$tip.label ==
    outgroup
)


parent <- rooted$edge[
  rooted$edge[, 2] ==
    tip_number,
  1
]


root_node <- Ntip(
  rooted
) + 1


cat(
  "\n============================================================\n"
)

cat(
  "FINAL ROOTED SPECIES TREE\n"
)

cat(
  "============================================================\n\n"
)


cat(
  "Tips:      ",
  Ntip(rooted),
  "\n",
  sep = ""
)


cat(
  "Outgroup:  ",
  outgroup,
  "\n",
  sep = ""
)


cat(
  "Root node: ",
  root_node,
  "\n",
  sep = ""
)


cat(
  "Outgroup parent: ",
  parent,
  "\n",
  sep = ""
)


cat(
  "Outgroup parent == root: ",
  identical(
    as.integer(parent),
    as.integer(root_node)
  ),
  "\n",
  sep = ""
)


cat(
  "\nTree:\n  ",
  output_tree,
  "\n",
  sep = ""
)
