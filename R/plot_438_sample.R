# Label the 438-read sample and redraw PCoA without it
source(here::here("R/count_tables.R"))

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(phyloseq)
  library(microViz)
  library(ggrepel)
})

# Build phyloseq and attach community/lacto_fraction
ps_species <- as_counts("instructional", table = "microbiome")

sample_groups <- read_csv(here("data/my_sample_groups.csv"), show_col_types = FALSE)

meta <- read_csv(here("data/2-workshop2/instructional_dataset/01_metadata.csv"), show_col_types = FALSE) |>
  left_join(sample_groups |> select(sample_id, community, lacto_fraction), by = "sample_id") |>
  mutate(library_size = sample_groups$total_reads[match(sample_id, sample_groups$sample_id)]) |>
  column_to_rownames("sample_id")

sample_data(ps_species) <- sample_data(meta[sample_names(ps_species), , drop = FALSE])

# Identify the 438-read sample
lib_size <- sample_data(ps_species)$library_size
sample_438 <- sample_names(ps_species)[which.min(lib_size)]
cat("438-read sample:", sample_438, "\n")
cat("Library size:", min(lib_size), "\n")

# Build the ordination
ord <- ps_species |>
  tax_transform("identity") |>
  dist_calc("bray") |>
  ord_calc("PCoA")

# Extract coordinates
axis_scores <- as.data.frame(vegan::scores(ord_get(ord), display = "sites"))
axis_scores$sample_id <- rownames(axis_scores)
axis_scores <- axis_scores |>
  mutate(
    is_438 = sample_id == sample_438,
    label = if_else(is_438, sample_id, NA_character_)
  )

# Base PCoA plot
p_base <- ord |>
  ord_plot(colour = "community", size = 1.5) +
  theme_bw()

# Add label for the 438 sample using the underlying data
p_original <- p_base +
  geom_point(data = axis_scores |> filter(is_438),
             aes(x = MDS1, y = MDS2), colour = "black", size = 3, shape = 21, stroke = 1.5) +
  geom_text_repel(data = axis_scores |> filter(is_438),
                  aes(x = MDS1, y = MDS2, label = sample_id),
                  colour = "black", nudge_x = 0.2, nudge_y = 0.2, size = 3) +
  ggtitle("Original (438-read sample highlighted)")

# Remove the 438-read sample
ps_dropped <- prune_samples(setdiff(sample_names(ps_species), sample_438), ps_species)
ps_dropped <- prune_taxa(taxa_sums(otu_table(ps_dropped)) > 0, ps_dropped)

ord_dropped <- ps_dropped |>
  tax_transform("identity") |>
  dist_calc("bray") |>
  ord_calc("PCoA")

p_dropped <- ord_dropped |>
  ord_plot(colour = "community", size = 1.5) +
  theme_bw() +
  ggtitle("After dropping 438-read sample")

# Side by side
combined <- cowplot::plot_grid(p_original, p_dropped, ncol = 2)

ggsave(here("figures/ordination/sample_438_comparison.png"), combined,
       width = 14, height = 6, dpi = 150)

cat("Plot saved to: figures/ordination/sample_438_comparison.png\n")
