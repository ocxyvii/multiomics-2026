# Interpret the species ordination:
# - correlation of axis 1 with lacto_fraction
# - five genera that track axis 1 most closely
source(here::here("R/count_tables.R"))

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(phyloseq)
  library(microViz)
})

# Build phyloseq and attach community/lacto_fraction
ps_species <- as_counts("instructional", table = "microbiome")

sample_groups <- read_csv(here("data/my_sample_groups.csv"), show_col_types = FALSE)

meta <- read_csv(here("data/2-workshop2/instructional_dataset/01_metadata.csv"), show_col_types = FALSE) |>
  left_join(sample_groups |> select(sample_id, community, lacto_fraction), by = "sample_id") |>
  mutate(library_size = sample_groups$total_reads[match(sample_id, sample_groups$sample_id)]) |>
  column_to_rownames("sample_id")

sample_data(ps_species) <- sample_data(meta[sample_names(ps_species), , drop = FALSE])

# PCoA with Bray-Curtis
ord <- ps_species |>
  tax_transform("identity") |>
  dist_calc("bray") |>
  ord_calc("PCoA")

# Axis 1 scores (MDS1)
axis_scores <- as.data.frame(vegan::scores(ord_get(ord), display = "sites"))
axis1 <- axis_scores[, 1]
names(axis1) <- rownames(axis_scores)
axis_df <- tibble(
  sample_id = names(axis1),
  axis1 = as.numeric(axis1)
)

# Attach lacto_fraction
axis_df <- axis_df |>
  left_join(
    tibble(
      sample_id = rownames(sample_data(ps_species)),
      lacto_fraction = as.numeric(sample_data(ps_species)$lacto_fraction)
    ),
    by = "sample_id"
  )

# Correlation of axis 1 with lacto_fraction
cor_pearson <- cor(axis_df$axis1, axis_df$lacto_fraction, use = "complete.obs")
cor_spearman <- cor(axis_df$axis1, axis_df$lacto_fraction, method = "spearman", use = "complete.obs")

cat("Pearson correlation (axis1 ~ lacto_fraction):", round(cor_pearson, 3), "\n")
cat("Spearman correlation (axis1 ~ lacto_fraction):", round(cor_spearman, 3), "\n")

# Aggregate to genus-level relative abundance per sample
ps_genus <- tax_glom(ps_species, taxrank = "Genus")
# Convert to relative abundance (sample sums to 1)
genus_relab <- as(otu_table(ps_genus), "matrix")
genus_relab <- sweep(genus_relab, 2, colSums(genus_relab), "/")

# Get genus names
genus_names <- tax_table(ps_genus)[, "Genus"]
rownames(genus_relab) <- genus_names

# Compute correlation of each genus with axis 1
sample_order <- names(axis1)
genus_relab <- genus_relab[, sample_order, drop = FALSE]

genus_cor <- apply(genus_relab, 1, function(x) {
  cor(x, axis1, use = "complete.obs")
})

genus_df <- tibble(
  genus = names(genus_cor),
  axis1_cor = as.numeric(genus_cor)
) |>
  arrange(desc(abs(axis1_cor)))

cat("\nTop 10 genera by absolute correlation with axis 1:\n")
print(genus_df |> head(10), n = 10)

cat("\nTop 5 positively correlated with axis 1:\n")
print(genus_df |> arrange(desc(axis1_cor)) |> head(5), n = 5)

cat("\nTop 5 negatively correlated with axis 1:\n")
print(genus_df |> arrange(axis1_cor) |> head(5), n = 5)
