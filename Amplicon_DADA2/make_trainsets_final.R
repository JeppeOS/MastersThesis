###############################################################
## Build GTDB r226 representative-SSU training sets for DADA2
##
## Input:
##   bac120_ssu_reps_r226.fna.gz
##   ar53_ssu_reps_r226.fna.gz
##
## Output:
##   GTDB_r226_bac120_reps_DADA2_genus_trainset.fa.gz
##   GTDB_r226_ar53_reps_DADA2_genus_trainset.fa.gz
##
## Taxonomy retained:
##   Domain -> Phylum -> Class -> Order -> Family -> Genus
###############################################################

library(Biostrings)


###############################################################
## Function to convert GTDB SSU FASTA -> DADA2 training FASTA
###############################################################

build_gtdb_dada2_trainset <- function(
    input_fasta,
    output_fasta,
    expected_domain,
    max_rank = "Genus"
) {
  
  ## Standard GTDB ranks ##
  gtdb_ranks <- c(
    "Domain",
    "Phylum",
    "Class",
    "Order",
    "Family",
    "Genus",
    "Species"
  )
  
  ## Determine how many ranks to retain ##
  n_keep <- match(max_rank, gtdb_ranks)
  
  if (is.na(n_keep)) {
    stop(
      "max_rank must be one of: ",
      paste(gtdb_ranks, collapse = ", ")
    )
  }
  
  ## Check input ##
  if (!file.exists(input_fasta)) {
    stop("Input file does not exist: ", input_fasta)
  }
  
  cat("\n==================================================\n")
  cat("Input:  ", input_fasta, "\n")
  cat("Output: ", output_fasta, "\n")
  cat("Keeping taxonomy through: ", max_rank, "\n")
  cat("==================================================\n\n")
  
  
  #############################################################
  ## Read GTDB sequences
  #############################################################
  
  gtdb <- readDNAStringSet(
    input_fasta,
    format = "fasta"
  )
  
  cat("Reference SSU sequences:", length(gtdb), "\n")
  
  raw_headers <- names(gtdb)
  
  
  #############################################################
  ## Extract GTDB taxonomy
  #############################################################
  
  ## Raw GTDB header example:
  ##
  ## RS_GCF_010645065.1
  ## d__Bacteria;p__...;g__Flavobacterium;s__...
  ## [locus_tag=...] [location=...] ...
  
  ## Remove genome/accession token at beginning ##
  taxonomy <- sub(
    "^[^[:space:]]+[[:space:]]+",
    "",
    raw_headers
  )
  
  ## Remove GTDB sequence metadata beginning with [ ... ##
  taxonomy <- sub(
    "[[:space:]]+\\[.*$",
    "",
    taxonomy
  )
  
  ## Split taxonomy into ranks ##
  tax_parts <- strsplit(
    taxonomy,
    ";",
    fixed = TRUE
  )
  
  
  #############################################################
  ## Validate original GTDB taxonomy
  #############################################################
  
  n_ranks <- lengths(tax_parts)
  
  cat("\nOriginal number of taxonomy ranks:\n")
  print(table(n_ranks))
  
  if (any(n_ranks < n_keep)) {
    stop(
      "At least one reference has fewer than ",
      n_keep,
      " taxonomy ranks."
    )
  }
  
  
  #############################################################
  ## Remove GTDB rank prefixes
  ##
  ## d__Bacteria       -> Bacteria
  ## p__Pseudomonadota -> Pseudomonadota
  ## g__Methylobacter  -> Methylobacter
  #############################################################
  
  tax_parts <- lapply(
    tax_parts,
    function(x) sub("^[a-z]__", "", x)
  )
  
  
  #############################################################
  ## Retain only requested ranks
  #############################################################
  
  taxonomy_clean <- vapply(
    tax_parts,
    function(x) {
      
      paste(
        x[seq_len(n_keep)],
        collapse = ";"
      )
      
    },
    character(1)
  )
  
  
  #############################################################
  ## Check domain
  #############################################################
  
  detected_domain <- vapply(
    tax_parts,
    function(x) x[1],
    character(1)
  )
  
  cat("\nDomains detected:\n")
  print(table(detected_domain))
  
  if (!all(detected_domain == expected_domain)) {
    
    stop(
      "Unexpected domain found. Expected only: ",
      expected_domain
    )
  }
  
  
  #############################################################
  ## Make DADA2-compatible headers
  #############################################################
  
  ## DADA2 format:
  ##
  ## >Domain;Phylum;Class;Order;Family;Genus;
  ## ACTG...
  
  names(gtdb) <- paste0(
    taxonomy_clean,
    ";"
  )
  
  
  #############################################################
  ## Validate final training set
  #############################################################
  
  final_parts <- strsplit(
    sub(";$", "", names(gtdb)),
    ";",
    fixed = TRUE
  )
  
  cat("\nFinal number of taxonomy ranks:\n")
  print(table(lengths(final_parts)))
  
  cat(
    "\nReference sequences:",
    length(gtdb),
    "\n"
  )
  
  cat(
    "Unique terminal taxonomy strings:",
    length(unique(names(gtdb))),
    "\n"
  )
  
  cat("\nExample formatted headers:\n")
  print(head(names(gtdb), 5))
  
  
  #############################################################
  ## Write compressed DADA2 training FASTA
  #############################################################
  
  writeXStringSet(
    gtdb,
    filepath = output_fasta,
    format = "fasta",
    compress = TRUE
  )
  
  cat(
    "\nTraining set written successfully:\n",
    output_fasta,
    "\n"
  )
  
  ## Return summary invisibly ##
  invisible(
    list(
      sequences = length(gtdb),
      unique_terminal_taxa = length(unique(names(gtdb))),
      ranks = gtdb_ranks[seq_len(n_keep)],
      output = output_fasta
    )
  )
}


###############################################################
## BACTERIA
###############################################################

bac_summary <- build_gtdb_dada2_trainset(
  input_fasta =
    "bac120_ssu_reps_r226.fna.gz",
  
  output_fasta =
    "GTDB_r226_bac120_reps_DADA2_genus_trainset.fa.gz",
  
  expected_domain =
    "Bacteria",
  
  max_rank =
    "Genus"
)

bac_summary


###############################################################
## Clear memory before Archaea
###############################################################

gc()


###############################################################
## ARCHAEA
###############################################################

arch_summary <- build_gtdb_dada2_trainset(
  input_fasta =
    "ar53_ssu_reps_r226.fna.gz",
  
  output_fasta =
    "GTDB_r226_ar53_reps_DADA2_genus_trainset.fa.gz",
  
  expected_domain =
    "Archaea",
  
  max_rank =
    "Genus"
)

arch_summary

