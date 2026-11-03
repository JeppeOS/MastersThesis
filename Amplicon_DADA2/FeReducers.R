## ============================================================
## SCREEN UL_SA_3 FOR DESCRIBED Fe(III)-REDUCING BACTERIA
## ============================================================

library(phyloseq)
library(dplyr)

## Sample corresponding to UL_SA_3 ##
target_sample <- "Bac21_S159"

## Genera with well-described Fe(III)-reducing representatives ##
iron_reducing_genera <- c(
  "Geobacter",
  "Geomonas",
  "Geopsychrobacter",
  "Geothrix",
  "Shewanella",
  "Anaeromyxobacter",
  "Rhodoferax",
  "Ferribacterium",
  "Ferrimonas",
  "Desulfuromonas",
  "Pelobacter",
  "Deferribacter"
)

## Allow GTDB suffixes such as Geobacter_A ##
iron_regex <- paste0(
  "^(",
  paste(iron_reducing_genera, collapse = "|"),
  ")(_[A-Z0-9]+)?$"
)

## Convert to relative abundance ##
ps.rel <- transform_sample_counts(
  ps,
  function(x) x / sum(x)
)

## Extract UL_SA_3 ##
df.iron <- psmelt(ps.rel) %>%
  filter(Sample == target_sample) %>%
  mutate(Genus = as.character(Genus)) %>%
  filter(
    !is.na(Genus),
    grepl(iron_regex, Genus)
  ) %>%
  group_by(Genus) %>%
  summarise(
    Abundance = sum(Abundance),
    .groups = "drop"
  ) %>%
  mutate(
    Percent = 100 * Abundance
  ) %>%
  arrange(desc(Percent))

## Show detected iron-reducing genera ##
df.iron
