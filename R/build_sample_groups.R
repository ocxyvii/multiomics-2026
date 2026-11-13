# Build sample groups for notebook 03-ordination.
# Defines "community" from Lactobacillus fraction > 0.5.
source(here::here("R/count_tables.R"))

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
})

d <- read_counts("instructional")

# taxonomy table (already has a 'taxon' column)
tax <- as_tibble(d$taxonomy)

# identify Lactobacillus taxa
lacto_taxa <- tax |>
  filter(Genus == "Lactobacillus") |>
  pull(taxon)

# compute per-sample total and Lactobacillus sum
spp <- d$species
total_counts <- colSums(spp)
lacto_counts <- colSums(spp[lacto_taxa, , drop = FALSE])

sample_groups <- tibble(
  sample_id = names(total_counts),
  total_reads = as.integer(total_counts),
  lacto_reads = as.integer(lacto_counts),
  lacto_fraction = lacto_reads / total_reads,
  community = if_else(lacto_fraction > 0.5, "high_lacto", "low_lacto")
)

write_csv(sample_groups, here("data/my_sample_groups.csv"))
message("Saved ", nrow(sample_groups), " rows to data/my_sample_groups.csv")
