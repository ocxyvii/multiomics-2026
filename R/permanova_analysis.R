# PERMANOVA for species ordination: real variables + shuffled outcome
source(here::here("R/count_tables.R"))

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(phyloseq)
  library(microViz)
  library(vegan)
})

# Build phyloseq and attach community/lacto_fraction
ps_species <- as_counts("instructional", table = "microbiome")

sample_groups <- read_csv(here("data/my_sample_groups.csv"), show_col_types = FALSE)

meta <- read_csv(here("data/2-workshop2/instructional_dataset/01_metadata.csv"), show_col_types = FALSE) |>
  left_join(sample_groups |> select(sample_id, community, lacto_fraction), by = "sample_id") |>
  mutate(library_size = sample_groups$total_reads[match(sample_id, sample_groups$sample_id)]) |>
  column_to_rownames("sample_id")

sample_data(ps_species) <- sample_data(meta[sample_names(ps_species), , drop = FALSE])

# Compute Bray-Curtis distance
ord <- ps_species |>
  tax_transform("identity") |>
  dist_calc("bray")

d <- dist_get(ord)
sd_df <- as(sample_data(ord), "data.frame") |>
  rownames_to_column("sample_id")

# --- Real PERMANOVAs ---------------------------------------------------------
run_permanova <- function(formula_str, perm = 999) {
  mod <- adonis2(as.formula(formula_str), data = sd_df, permutations = perm, by = "terms")
  tibble(
    term = rownames(mod),
    R2 = mod$R2,
    p = mod$`Pr(>F)`
  ) |>
    filter(term != "Residual" & term != "Total")
}

cat("=== Real PERMANOVA results (999 permutations) ===\n")
real_vars <- c("community", "outcome", "time_point", "rnaseq_batch")
real_results <- map_dfr(real_vars, ~ run_permanova(paste0("d ~ ", .x))) |>
  mutate(model = real_vars, .before = everything())
print(real_results)

# --- Shuffled outcome PERMANOVA (20 times, 99 permutations each) -------------
set.seed(2026)
shuffled_R2 <- map_dfr(1:20, function(i) {
  df <- sd_df |> mutate(outcome_shuffled = sample(outcome))
  mod <- adonis2(d ~ outcome_shuffled, data = df, permutations = 99, by = "terms")
  tibble(
    shuffle = i,
    R2 = mod$R2[1],
    p = mod$`Pr(>F)`[1]
  )
})

real_outcome_R2 <- real_results |> filter(model == "outcome") |> pull(R2)

cat("\n=== Shuffled outcome PERMANOVA (20 times, 99 permutations) ===\n")
print(shuffled_R2, n = 20)

# --- Plot --------------------------------------------------------------------
shuffled_R2 <- shuffled_R2 |>
  mutate(
    above_real = R2 >= real_outcome_R2,
    label_y = if_else(R2 >= real_outcome_R2, 0.05, -0.05)
  )

p <- ggplot(shuffled_R2, aes(x = R2)) +
  geom_histogram(bins = 10, fill = "grey70", colour = "white") +
  geom_vline(xintercept = real_outcome_R2, colour = "red", linewidth = 1) +
  annotate("text", x = real_outcome_R2, y = Inf, label = "Real outcome R²",
           vjust = 1.5, hjust = -0.1, colour = "red") +
  labs(
    title = "Outcome PERMANOVA R²: 20 shuffles (99 permutations each)",
    subtitle = paste("Real R² =", round(real_outcome_R2, 4), "(red line)"),
    x = expression("PERMANOVA R"^2),
    y = "Count"
  ) +
  theme_bw()

ggsave(here("figures/ordination/shuffled_outcome_permanova.png"), p,
       width = 8, height = 5, dpi = 150)

cat("\nPlot saved to: figures/ordination/shuffled_outcome_permanova.png\n")

# Count how many shuffled R2 are >= real
cat("\nNumber of shuffled R2 >= real R2:", sum(shuffled_R2$above_real), "out of 20\n")
