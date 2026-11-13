# Plots for enrichment analysis
library(tidyverse)
library(here)
library(fgsea)
library(msigdbr)

source(here("R/count_tables.R"))

read_result <- function(file, where) {
  path <- here("data", file)
  if (!file.exists(path)) {
    stop("Could not find ", path, "\nGo back to notebook 04, ", where,
         ", and save the results table there.", call. = FALSE)
  }
  read_csv(path, show_col_types = FALSE)
}

genes_by_outcome <- read_result("genes_by_outcome.csv", "Exercise 2, step 5")

# Rank vector
gene_stat <- genes_by_outcome %>%
  filter(!is.na(stat)) %>%
  pull(stat, name = feature)

# Hallmark sets
hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>%
  select(gs_name, ensembl_gene)
hallmark_sets <- split(hallmark$ensembl_gene, hallmark$gs_name)

set.seed(42)
fgsea_hall <- fgsea(
  pathways = hallmark_sets,
  stats = gene_stat,
  minSize = 15,
  maxSize = 500,
  nproc = 1
) %>% arrange(padj, pval)

# 1. Bar chart of top 12 NES
plot_data <- fgsea_hall %>%
  mutate(
    direction = ifelse(NES > 0, "Higher in preterm", "Lower in preterm"),
    pathway = str_remove(pathway, "HALLMARK_"),
    pathway = str_replace_all(pathway, "_", " ")
  ) %>%
  slice_max(order_by = abs(NES), n = 12)

p_bar <- ggplot(plot_data, aes(x = reorder(pathway, NES), y = NES, fill = direction)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("Higher in preterm" = "#D55E00", "Lower in preterm" = "#0072B2")) +
  labs(
    title = "Top 12 Hallmark gene sets by NES",
    subtitle = "Preterm vs term birth, instructional dataset",
    x = NULL,
    y = "Normalized enrichment score (NES)",
    fill = "Direction"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")

ggsave(here("data", "hallmark_nes_bar.png"), p_bar, width = 10, height = 7, dpi = 150)

# 2. Running score plot for top set
top_pathway <- fgsea_hall$pathway[1]
plotEnrichment(hallmark_sets[[top_pathway]], gene_stat) +
  labs(
    title = str_remove(top_pathway, "HALLMARK_"),
    subtitle = paste("NES =", round(fgsea_hall$NES[1], 2), "| padj =", format(fgsea_hall$padj[1], digits = 2)),
    x = "Rank",
    y = "Enrichment score"
  ) +
  theme_minimal(base_size = 14)

ggsave(here("data", "hallmark_running_score.png"), width = 8, height = 5, dpi = 150)

# 3. Dot plot of all 50 sets
p_dot <- fgsea_hall %>%
  mutate(
    direction = ifelse(NES > 0, "Higher in preterm", "Lower in preterm"),
    pathway = str_remove(pathway, "HALLMARK_"),
    pathway = str_replace_all(pathway, "_", " ")
  ) %>%
  ggplot(aes(x = NES, y = reorder(pathway, NES), color = -log10(padj), size = size)) +
  geom_point() +
  scale_color_gradient(low = "grey80", high = "#D55E00") +
  labs(
    title = "All 50 Hallmark sets",
    x = "NES",
    y = NULL,
    color = "-log10 padj",
    size = "Set size"
  ) +
  theme_minimal(base_size = 12)

ggsave(here("data", "hallmark_dotplot.png"), p_dot, width = 10, height = 10, dpi = 150)

cat("Plots saved to data/:\n")
cat("- hallmark_nes_bar.png\n")
cat("- hallmark_running_score.png\n")
cat("- hallmark_dotplot.png\n")
