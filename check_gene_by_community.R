library(tidyverse)
library(here)
library(DESeq2)

source(here("R/count_tables.R"))

# Load the data -----------------------------------------------------------
d <- read_counts("instructional")

# The grouping you built and saved in notebook 02
sample_groups <- read_csv(here("data/my_sample_groups.csv"), show_col_types = FALSE)

metadata <- d$metadata |>
  left_join(sample_groups, by = "sample_id")

# Keep only week 18 samples ------------------------------------------------
meta18 <- metadata |>
  filter(time_point == "week_18") |>
  mutate(
    community = case_when(
      community == "high_lacto" ~ "lactobacillus",
      community == "low_lacto"  ~ "anaerobe",
      TRUE                       ~ community
    ),
    community = factor(community, levels = c("lactobacillus", "anaerobe"))
  )

message("Samples at visit 18: ", nrow(meta18))
table(meta18$community)

# Get the gene count matrix and keep only week 18 columns in metadata order
counts_all <- d$genes
counts18 <- counts_all[, meta18$sample_id, drop = FALSE]

# Check 1: sample order -----------------------------------------------------
cat("\n--- Sample order check ---\n")
cat("First 5 sample_id in metadata:\n")
print(head(meta18$sample_id, 5))

cat("First 5 column names of count matrix:\n")
print(head(colnames(counts18), 5))

cat("All match, in order:", all(colnames(counts18) == meta18$sample_id), "\n")

# Filter genes: count >= 10 in at least 25% of samples ----------------------
min_samples <- ceiling(0.25 * ncol(counts18))
keep <- rowSums(counts18 >= 10) >= min_samples
counts18_filt <- counts18[keep, , drop = FALSE]

message("\nGenes kept after filtering: ", nrow(counts18_filt),
        " (out of ", nrow(counts18), ")")
message("Filter: >= 10 counts in >= ", min_samples, " samples")

# Run DESeq2 with design ~ community ----------------------------------------
dds <- DESeqDataSetFromMatrix(
  countData = counts18_filt,
  colData   = meta18,
  design    = ~ community
)

dds <- DESeq(dds)

# Results: contrast anaerobe vs lactobacillus -------------------------------
# In DESeq2, contrast = c("community", "anaerobe", "lactobacillus") means
# the log2 fold change is log2(anaerobe / lactobacillus)
res <- results(dds, contrast = c("community", "anaerobe", "lactobacillus"))

# Add gene info
res_df <- as.data.frame(res) |>
  rownames_to_column("feature") |>
  left_join(d$gene_info, by = c("feature" = "ensembl_gene_id")) |>
  select(feature, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj,
         symbol, gene_name)

# Save
write_csv(res_df, here("data/genes_by_community.csv"))
message("\nSaved results to data/genes_by_community.csv")

# Summary -------------------------------------------------------------------
cat("\n--- Result summary ---\n")
cat("Raw p < 0.05:", sum(res_df$pvalue < 0.05, na.rm = TRUE), "\n")
cat("Adjusted p (padj) < 0.05:", sum(res_df$padj < 0.05, na.rm = TRUE), "\n")
cat("padj < 0.05 & |log2FC| > 1:", sum(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, na.rm = TRUE), "\n")
cat("padj < 0.05 & log2FC > 1 (higher in anaerobe):",
    sum(res_df$padj < 0.05 & res_df$log2FoldChange > 1, na.rm = TRUE), "\n")
cat("padj < 0.05 & log2FC < -1 (higher in lactobacillus):",
    sum(res_df$padj < 0.05 & res_df$log2FoldChange < -1, na.rm = TRUE), "\n")

# Top 20 by padj
cat("\n--- Top 20 genes by padj ---\n")
print(
  res_df |>
    dplyr::arrange(padj) |>
    head(20) |>
    dplyr::select(feature, symbol, gene_name, log2FoldChange, padj)
)

# Direction check: show the sign convention
# The first level of the factor is lactobacillus (reference).
# A positive log2FoldChange means the gene is higher in anaerobe.
cat("\n--- Contrast direction ---\n")
cat("Factor levels: ", paste(levels(meta18$community), collapse = ", "), "\n")
cat("Contrast used: c('community', 'anaerobe', 'lactobacillus')\n")
cat("Positive log2FoldChange = higher in anaerobe-dominant community.\n")
cat("Negative log2FoldChange = higher in lactobacillus-dominant community.\n")
