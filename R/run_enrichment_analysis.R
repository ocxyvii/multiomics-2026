# Module 5 enrichment analysis
# Run from project root
library(tidyverse)
library(here)
library(fgsea)
library(msigdbr)

source(here("R/count_tables.R"))

# ---- Load result tables from notebook 04 ----
read_result <- function(file, where) {
  path <- here("data", file)
  if (!file.exists(path)) {
    stop("Could not find ", path, "\nGo back to notebook 04, ", where,
         ", and save the results table there.", call. = FALSE)
  }
  read_csv(path, show_col_types = FALSE)
}

genes_by_outcome   <- read_result("genes_by_outcome.csv",   "Exercise 2, step 5")
genes_by_community <- read_result("genes_by_community.csv", "Exercise 1, step 1")
species_by_outcome <- read_result("species_by_outcome.csv", "Exercise 2, step 3")

cat("Dimensions:\n")
cat("genes_by_outcome:", nrow(genes_by_outcome), "\n")
cat("genes_by_community:", nrow(genes_by_community), "\n")
cat("species_by_outcome:", nrow(species_by_outcome), "\n")

# ---- 1. Rank genes by outcome ----
gene_stat <- genes_by_outcome %>%
  filter(!is.na(stat)) %>%
  pull(stat, name = feature)

cat("\nGenes in ranking:", length(gene_stat), "\n")
cat("Range of stat:", round(min(gene_stat), 3), "to", round(max(gene_stat), 3), "\n")
cat("Positive stat (up in preterm):", sum(gene_stat > 0), 
    "(", round(100*sum(gene_stat > 0)/length(gene_stat), 1), "%)\n")

# ---- 2. Hallmark gene sets ----
hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>%
  select(gs_name, ensembl_gene)

hallmark_sets <- split(hallmark$ensembl_gene, hallmark$gs_name)

cat("\nHallmark sets:", length(hallmark_sets), "\n")
cat("Unique genes in Hallmark:", length(unique(hallmark$ensembl_gene)), "\n")

# Identifier matching check
in_ranking <- names(gene_stat)
in_sets <- unique(hallmark$ensembl_gene)
cat("\n--- Identifier matching ---\n")
cat("Genes in ranking:", length(in_ranking), "\n")
cat("Genes in Hallmark sets:", length(in_sets), "\n")
cat("Overlap:", length(intersect(in_ranking, in_sets)), "\n")
cat("Ranking genes not in any Hallmark set:", length(setdiff(in_ranking, in_sets)), "\n")
cat("Hallmark genes not in ranking:", length(setdiff(in_sets, in_ranking)), "\n")

# ---- 3. Run fgsea on Hallmark ----
set.seed(42)
fgsea_hall <- fgsea(
  pathways = hallmark_sets,
  stats = gene_stat,
  minSize = 15,
  maxSize = 500,
  nproc = 1
)

fgsea_hall <- fgsea_hall %>% arrange(padj, pval)

cat("\n--- Hallmark fgsea results ---\n")
cat("Significant sets (padj < 0.05):", sum(fgsea_hall$padj < 0.05, na.rm = TRUE), "\n")
cat("Positive NES (higher in preterm):", sum(fgsea_hall$padj < 0.05 & fgsea_hall$NES > 0, na.rm = TRUE), "\n")
cat("Negative NES (lower in preterm):", sum(fgsea_hall$padj < 0.05 & fgsea_hall$NES < 0, na.rm = TRUE), "\n")

print(fgsea_hall %>% select(pathway, pval, padj, NES, size, leadingEdge) %>% head(12))

# ---- 4. Shared genes between top sets ----
cat("\n--- Shared genes between top 5 Hallmark sets ---\n")
top5 <- fgsea_hall$pathway[1:5]
shared_matrix <- matrix(0, nrow = 5, ncol = 5, dimnames = list(top5, top5))
for (i in 1:5) {
  for (j in 1:5) {
    shared_matrix[i, j] <- length(intersect(hallmark_sets[[top5[i]]], hallmark_sets[[top5[j]]]))
  }
}
print(shared_matrix)

# Jaccard overlap
cat("\nJaccard overlap (shared / union) between top 5 sets:\n")
jaccard_matrix <- matrix(0, nrow = 5, ncol = 5, dimnames = list(top5, top5))
for (i in 1:5) {
  for (j in 1:5) {
    jaccard_matrix[i, j] <- length(intersect(hallmark_sets[[top5[i]]], hallmark_sets[[top5[j]]])) /
      length(union(hallmark_sets[[top5[i]]], hallmark_sets[[top5[j]]]))
  }
}
print(round(jaccard_matrix, 3))

# ---- 5. Check ranking imbalance ----
cat("\n--- Ranking imbalance check ---\n")
cat("All genes: positive stat =", sum(gene_stat > 0), 
    "(", round(100*sum(gene_stat > 0)/length(gene_stat), 1), "%)\n")

sig_genes <- genes_by_outcome %>% filter(padj < 0.05 & !is.na(stat))
cat("Significant genes (padj < 0.05):", nrow(sig_genes),
    "; positive stat =", sum(sig_genes$stat > 0),
    "(", round(100*sum(sig_genes$stat > 0)/nrow(sig_genes), 1), "%)\n")

strong500 <- sort(abs(gene_stat), decreasing = TRUE)[1:min(500, length(gene_stat))]
strong500_names <- names(strong500)
strong500_signs <- sign(gene_stat[strong500_names])
cat("Top 500 strongest |stat|: positive =", sum(strong500_signs > 0),
    "(", round(100*sum(strong500_signs > 0)/length(strong500_signs), 1), "%)\n")

# Compare to community ranking
community_stat <- genes_by_community %>%
  filter(!is.na(stat)) %>%
  pull(stat, name = feature)
cat("\n--- Community ranking imbalance ---\n")
cat("All genes: positive stat =", sum(community_stat > 0),
    "(", round(100*sum(community_stat > 0)/length(community_stat), 1), "%)\n")
community_strong500 <- sort(abs(community_stat), decreasing = TRUE)[1:min(500, length(community_stat))]
cat("Top 500 strongest |stat|: positive =", sum(community_stat[names(community_strong500)] > 0),
    "(", round(100*sum(community_stat[names(community_strong500)] > 0)/length(community_strong500), 1), "%)\n")

# ---- 6. Cutoff-based test (hypergeometric / over-representation) ----
cat("\n--- Cutoff-based hypergeometric test ---\n")

# Universe: genes in ranking and in at least one Hallmark set
universe <- intersect(names(gene_stat), in_sets)
gene_stat_in_universe <- gene_stat[universe]

# Up-regulated significant genes
up_genes <- genes_by_outcome %>%
  filter(padj < 0.05, stat > 0, feature %in% universe) %>%
  pull(feature)

down_genes <- genes_by_outcome %>%
  filter(padj < 0.05, stat < 0, feature %in% universe) %>%
  pull(feature)

cat("Up-regulated significant genes in universe:", length(up_genes), "\n")
cat("Down-regulated significant genes in universe:", length(down_genes), "\n")

# Run fgsea-simple (ORA equivalent) for up and down separately
ora_up <- fora(
  pathways = hallmark_sets,
  genes = up_genes,
  universe = universe,
  minSize = 15,
  maxSize = 500
) %>% arrange(padj)

ora_down <- fora(
  pathways = hallmark_sets,
  genes = down_genes,
  universe = universe,
  minSize = 15,
  maxSize = 500
) %>% arrange(padj)

cat("\nORA up-regulated: significant sets (padj < 0.05):", sum(ora_up$padj < 0.05, na.rm = TRUE), "\n")
print(ora_up %>% select(pathway, overlap, size, pval, padj) %>% head(8))

cat("\nORA down-regulated: significant sets (padj < 0.05):", sum(ora_down$padj < 0.05, na.rm = TRUE), "\n")
print(ora_down %>% select(pathway, overlap, size, pval, padj) %>% head(8))

# ---- 7. Second set collection: Reactome (C2 CP:REACTOME) ----
cat("\n--- Reactome fgsea ---\n")
reactome <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") %>%
  select(gs_name, ensembl_gene)
reactome_sets <- split(reactome$ensembl_gene, reactome$gs_name)

set.seed(42)
fgsea_react <- fgsea(
  pathways = reactome_sets,
  stats = gene_stat,
  minSize = 15,
  maxSize = 500,
  nproc = 1
) %>% arrange(padj, pval)

cat("Reactome sets:", length(reactome_sets), "\n")
cat("Significant Reactome sets (padj < 0.05):", sum(fgsea_react$padj < 0.05, na.rm = TRUE), "\n")
cat("Positive NES:", sum(fgsea_react$padj < 0.05 & fgsea_react$NES > 0, na.rm = TRUE), "\n")
cat("Negative NES:", sum(fgsea_react$padj < 0.05 & fgsea_react$NES < 0, na.rm = TRUE), "\n")
print(fgsea_react %>% select(pathway, pval, padj, NES, size) %>% head(12))

# Save results
write_csv(fgsea_hall %>% select(-leadingEdge), here("data", "fgsea_hallmark_outcome.csv"))
write_csv(fgsea_react %>% select(-leadingEdge), here("data", "fgsea_reactome_outcome.csv"))
write_csv(ora_up, here("data", "ora_hallmark_up.csv"))
write_csv(ora_down, here("data", "ora_hallmark_down.csv"))

# ---- 8. Taxonomic sets: genus and family from species table ----
cat("\n\n=== TAXONOMIC SETS ===\n")

species_stat <- species_by_outcome %>%
  filter(!is.na(stat)) %>%
  pull(stat, name = feature)

cat("Species in ranking:", length(species_stat), "\n")

# Genus sets
genus_sets <- species_by_outcome %>%
  filter(!is.na(Genus)) %>%
  group_by(Genus) %>%
  filter(n() >= 2) %>%
  summarise(species = list(feature), .groups = "drop") %>%
  deframe()

cat("\nGenus sets (>=2 species):", length(genus_sets), "\n")
cat("Set sizes:", paste(sapply(genus_sets, length), collapse = ", "), "\n")
cat("Largest genus size:", max(sapply(genus_sets, length)), "\n")
print(tibble(Genus = names(genus_sets), n_species = sapply(genus_sets, length)))

# Check if any species appears in more than one genus set
all_species_in_sets <- unlist(genus_sets)
cat("Species in multiple genus sets:", sum(duplicated(all_species_in_sets)), "\n")

set.seed(42)
fgsea_genus <- fgsea(
  pathways = genus_sets,
  stats = species_stat,
  minSize = 2,
  maxSize = 50,
  nproc = 1
) %>% arrange(padj, pval)

cat("\n--- Genus fgsea results ---\n")
cat("Significant genus sets (padj < 0.05):", sum(fgsea_genus$padj < 0.05, na.rm = TRUE), "\n")
print(fgsea_genus %>% select(pathway, pval, padj, NES, size, leadingEdge))

# Family sets
family_sets <- species_by_outcome %>%
  filter(!is.na(Family)) %>%
  group_by(Family) %>%
  filter(n() >= 2) %>%
  summarise(species = list(feature), .groups = "drop") %>%
  deframe()

cat("\nFamily sets (>=2 species):", length(family_sets), "\n")
cat("Set sizes:", paste(sapply(family_sets, length), collapse = ", "), "\n")

set.seed(42)
fgsea_family <- fgsea(
  pathways = family_sets,
  stats = species_stat,
  minSize = 2,
  maxSize = 50,
  nproc = 1
) %>% arrange(padj, pval)

cat("\n--- Family fgsea results ---\n")
cat("Significant family sets (padj < 0.05):", sum(fgsea_family$padj < 0.05, na.rm = TRUE), "\n")
print(fgsea_family %>% select(pathway, pval, padj, NES, size, leadingEdge))

# Save taxonomic results
write_csv(fgsea_genus %>% select(-leadingEdge), here("data", "fgsea_genus_outcome.csv"))
write_csv(fgsea_family %>% select(-leadingEdge), here("data", "fgsea_family_outcome.csv"))

# Simple table: species sorted by stat
cat("\n--- Species sorted by stat ---\n")
print(species_by_outcome %>%
        select(feature, stat, padj, Genus, Family) %>%
        arrange(desc(stat)) %>%
        head(15))

cat("\nAnalysis complete.\n")
