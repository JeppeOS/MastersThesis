library(phyloseq)
library(dplyr)

ps <- readRDS("results_FH/phyloseq_gtdb.rds")

ps.rel <- transform_sample_counts(ps, function(x) x / sum(x))

psmelt(ps.rel) %>%
  filter(
    Sample == "JO-Bac2_S70",   # FH2
    Genus == "UBA4132"
  ) %>%
  summarise(
    Percent = sum(Abundance) * 100
  )
