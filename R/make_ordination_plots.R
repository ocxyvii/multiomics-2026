# Generate ordination plots for notebook 03-ordination
source(here::here("R/count_tables.R"))

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(microViz)
  library(phyloseq)
  library(patchwork)
  library(vegan)
})

# --- Build phyloseq objects and attach community / metadata -------------------
ps_species <- as_counts("instructional", table = "microbiome")
ps_genes   <- as_counts("instructional", table = "rnaseq",
                        min_count = 10, min_samples = 100)

sample_groups <- read_csv(here("data/my_sample_groups.csv"), show_col_types = FALSE)

meta <- read_csv(here("data/2-workshop2/instructional_dataset/01_metadata.csv"),
                 show_col_types = FALSE) |>
  left_join(sample_groups |> select(sample_id, community, lacto_fraction),
            by = "sample_id") |>
  mutate(library_size = sample_groups$total_reads[match(sample_id, sample_groups$sample_id)]) |>
  column_to_rownames("sample_id")

# Reorder to match phyloseq sample order
sample_data(ps_species) <- sample_data(meta[sample_names(ps_species), , drop = FALSE])
sample_data(ps_genes)   <- sample_data(meta[sample_names(ps_genes), , drop = FALSE])

# --- Helper: build ordination and return object + eigenvalues ---------------
make_ord <- function(ps, trans, dist) {
  ps |>
    tax_transform(trans) |>
    dist_calc(dist) |>
    ord_calc("PCoA")
}

# --- Plot 1: three distances side by side with scree -------------------------
make_pcoa_scree <- function(ord, title) {
  p_pcoa <- ord |>
    ord_plot(colour = "community", size = 1.2) +
    ggtitle(title) +
    theme_bw()

  eig <- vegan::eigenvals(ord_get(ord))
  eig <- sort(eig, decreasing = TRUE)
  n <- min(8, length(eig))
  scree_df <- tibble(
    axis = seq_len(n),
    variance_pct = (eig[seq_len(n)] / sum(eig[eig > 0])) * 100
  )

  p_scree <- ggplot(scree_df, aes(x = axis, y = variance_pct)) +
    geom_col() +
    scale_x_continuous(breaks = 1:8) +
    labs(x = "Axis", y = "% variation") +
    theme_bw()

  p_pcoa / p_scree
}

ord_bray  <- make_ord(ps_species, "identity", "bray")
ord_jacc  <- make_ord(ps_species, "binary",   "jaccard")
ord_eucl  <- make_ord(ps_species, "clr",      "euclidean")

p_bray  <- make_pcoa_scree(ord_bray,  "Bray-Curtis")
p_jacc  <- make_pcoa_scree(ord_jacc,  "Jaccard")
p_eucl  <- make_pcoa_scree(ord_eucl,  "Euclidean (clr)")

three_panel <- wrap_plots(p_bray, p_jacc, p_eucl, ncol = 3)
ggsave(here("figures/ordination/three_distances.png"), three_panel,
       width = 15, height = 6, dpi = 150)

# --- Plot 2: multi-coloured ordinations ---------------------------------------
make_coloured_ord <- function(ord, colour_var) {
  ord |>
    ord_plot(colour = colour_var, size = 1.2) +
    theme_bw() +
    labs(title = colour_var)
}

colour_vars <- c("community", "outcome", "time_point", "nugent", "lacto_fraction", "library_size")

colour_plots <- map(colour_vars, ~ make_coloured_ord(ord_bray, .x))
names(colour_plots) <- colour_vars

# Arrange in a 2x3 grid
colour_grid <- wrap_plots(colour_plots, ncol = 3)
ggsave(here("figures/ordination/multi_colour_grid.png"), colour_grid,
       width = 15, height = 10, dpi = 150)

# --- Plot 3: single plot with community (the main one) -----------------------
p_community <- ord_bray |>
  ord_plot(colour = "community", size = 1.5) +
  theme_bw() +
  ggtitle("Bray-Curtis PCoA coloured by community")
ggsave(here("figures/ordination/community_pcoa.png"), p_community,
       width = 7, height = 5, dpi = 150)

# --- Print some diagnostics --------------------------------------------------
cat("n samples in ps_species:", nsamples(ps_species), "\n")
cat("n features in ps_species:", ntaxa(ps_species), "\n")
cat("Missing community values:", sum(is.na(sample_data(ps_species)$community)), "\n")
cat("Saved plots to figures/ordination/\n")
