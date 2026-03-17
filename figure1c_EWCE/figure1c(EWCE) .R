# ==============================================================================
# Single-Cell Expression Weighted Celltype Enrichment (EWCE) Analysis
#
# Workflow:
#   1. Compute per-cell-type average expression from an integrated Seurat object
#   2. CPM-normalise and calculate per-gene specificity scores
#   3. Run a permutation-based EWCE test for up- and down-regulated gene sets
#   4. Visualise enrichment results as a diverging bar chart
# ==============================================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(Matrix)


# ==============================================================================
# Section 1: Compute per-cell-type average expression
#
# Input:  scRNA_combined — integrated Seurat object with a metadata column
#         that identifies the cells cell-type label for each cell.
# Output: exp_lvl5 — long-format tibble with columns:
#           Gene | cells | Expr_sum_mean | specificity
# ==============================================================================

# FIX: Replace "cells" with the actual metadata column that holds the cells
#      cell-type labels (e.g. "cells", "cell_type").
#      "cells" is a Seurat reserved keyword for cell barcodes, not cell types.
cells_META_COL <- "cells"   # <-- change to match your Seurat metadata column name

avg_exp <- AverageExpression(
  scRNA_combined,
  assays   = "RNA",
  group.by = cells_META_COL,   # Group cells by cells cell-type annotation
  return.seurat = TRUE
)

# Extract the log-normalised expression matrix (genes × cell-types)
exp_matrix <- Seurat::GetAssayData(avg_exp, assay = "RNA", layer = "data")

# Reshape to long format so each row is one gene–cell-type observation
exp_lvl5 <- exp_matrix %>%
  as.data.frame() %>%
  rownames_to_column("Gene") %>%
  pivot_longer(
    cols      = -Gene,
    names_to  = "cells",
    values_to = "Expr_sum_mean"
  )

# Remove genes that are completely unexpressed across all cell types,
# as they carry no specificity information
unexpressed_genes <- exp_lvl5 %>%
  group_by(Gene) %>%
  summarise(total = sum(Expr_sum_mean), .groups = "drop") %>%
  filter(total == 0) %>%
  pull(Gene)

exp_lvl5 <- exp_lvl5 %>% filter(!Gene %in% unexpressed_genes)

# CPM normalisation within each cell type:
# rescale expression values so they sum to 1e6 per cell type,
# making cell types with different total counts comparable
exp_lvl5 <- exp_lvl5 %>%
  group_by(cells) %>%
  mutate(Expr_sum_mean = Expr_sum_mean * 1e6 / sum(Expr_sum_mean)) %>%
  ungroup()

# Compute per-gene specificity: the fraction of a gene's total expression
# that is attributable to each cell type (sums to 1 across cell types per gene)
exp_lvl5 <- exp_lvl5 %>%
  group_by(Gene) %>%
  mutate(specificity = Expr_sum_mean / sum(Expr_sum_mean)) %>%
  ungroup()

# Column order expected by ewce_analysis(): Gene | CellType | Expr_sum_mean | specificity
exp_lvl5 <- exp_lvl5 %>%
  rename(CellType = cells) %>%               # Standardise column name up front
  dplyr::select(Gene, CellType, Expr_sum_mean, specificity)


# ==============================================================================
# Section 2: EWCE permutation test
#
# For each cell type, tests whether the mean specificity of the target gene set
# is significantly higher than expected by chance (drawn from the background).
#
# Arguments:
#   target_genes     — character vector of genes of interest
#   background_genes — universe of expressed genes (target set must be a subset)
#   specificity_df   — long-format data frame with columns:
#                        Gene | CellType | Expr_sum_mean | specificity
#   n_perm           — number of random permutations (default 1000)
#   seed             — random seed for reproducibility (default 42)
#   p_adjust_method  — multiple-testing correction method passed to p.adjust()
#
# Returns a tibble: CellType | obs_score | p_value | p_adj
# ==============================================================================
ewce_analysis <- function(target_genes,
                          background_genes,
                          specificity_df,
                          n_perm          = 1000,
                          seed            = 42,
                          p_adjust_method = "BH") {
  
  set.seed(seed)   # Ensure permutation results are reproducible
  
  # Restrict both gene sets to genes present in the specificity matrix
  target_genes     <- intersect(target_genes,     specificity_df$Gene)
  background_genes <- intersect(background_genes, specificity_df$Gene)
  
  # FIX: exclude target genes from the background pool so that random draws
  #      cannot contain target genes, keeping the null distribution unbiased
  background_only <- setdiff(background_genes, target_genes)
  
  if (length(target_genes) == 0)
    stop("No target genes found in specificity_df after intersection.")
  if (length(background_only) < length(target_genes))
    stop("Background gene pool (after excluding targets) is smaller than the target set.")
  
  # Observed mean specificity per cell type for the actual target genes
  obs_scores <- specificity_df %>%
    filter(Gene %in% target_genes) %>%
    group_by(CellType) %>%
    summarise(obs_score = mean(specificity), .groups = "drop")
  
  cell_types <- obs_scores$CellType
  
  # Build a permutation null distribution:
  # sample the same number of genes as the target set from the background,
  # and record the mean specificity per cell type for each random draw
  perm_matrix <- matrix(
    nrow     = n_perm,
    ncol     = length(cell_types),
    dimnames = list(NULL, cell_types)
  )
  
  # FIX: use seq_len() instead of 1:n_perm to avoid the c(1,0) edge case
  #      when n_perm = 0
  for (i in seq_len(n_perm)) {
    random_genes <- sample(background_only, length(target_genes))
    
    rand_scores <- specificity_df %>%
      filter(Gene %in% random_genes) %>%
      group_by(CellType) %>%
      summarise(score = mean(specificity), .groups = "drop")
    
    # Align random scores to the fixed cell-type column order
    perm_matrix[i, ] <- rand_scores$score[match(cell_types, rand_scores$CellType)]
  }
  
  # Empirical one-tailed p-value: proportion of random draws that are at least
  # as extreme as the observed score (right-tail test for positive enrichment)
  results <- obs_scores %>%
    mutate(
      p_value = sapply(seq_len(nrow(obs_scores)), function(j) {
        ct <- CellType[j]
        mean(perm_matrix[, ct] >= obs_score[j])
      }),
      p_adj = p.adjust(p_value, method = p_adjust_method)
    )
  
  return(results)
}


# ==============================================================================
# Section 3: Identify differentially expressed genes (PD vs control)
# ==============================================================================

Idents(scRNA_combined) <- "group"

# Wilcoxon rank-sum test comparing PD cells ("pd") to control cells ("con").
# Only genes detected in ≥ 25 % of cells and with |log2FC| > 0.25 are tested.
markers_after <- FindMarkers(
  object          = scRNA_combined,
  ident.1         = "pd",
  ident.2         = "con",
  only.pos        = FALSE,
  min.pct         = 0.25,
  logfc.threshold = 0.25,
  test.use        = "wilcox",
  verbose         = TRUE
)
markers_after$gene <- rownames(markers_after)

# Select the top 600 up- and down-regulated genes by log2 fold-change.
# NA values in avg_log2FC are excluded first to prevent silent mis-sorting.
target_genes_up <- markers_after %>%
  filter(!is.na(avg_log2FC)) %>%
  slice_max(order_by = avg_log2FC, n = 600) %>%
  pull(gene)

target_genes_down <- markers_after %>%
  filter(!is.na(avg_log2FC)) %>%
  slice_min(order_by = avg_log2FC, n = 600) %>%
  pull(gene)

# Use all genes present in the specificity matrix as the background universe
background_genes <- unique(exp_lvl5$Gene)


# ==============================================================================
# Section 4: Run EWCE for up- and down-regulated gene sets
# ==============================================================================
results_up   <- ewce_analysis(target_genes_up,   background_genes, exp_lvl5)
results_down <- ewce_analysis(target_genes_down, background_genes, exp_lvl5)


# ==============================================================================
# Section 5: Visualise enrichment as a diverging bar chart
#
# Convention: upregulated enrichment scores are shown as positive bars (red),
#             downregulated enrichment scores are mirrored as negative bars
#             (blue) for visual symmetry ONLY — the raw obs_score values are
#             stored separately and the negation is applied only to the plot
#             column so that the original values remain accessible.
# ==============================================================================

# FIX: Keep the original obs_score intact; use a separate plot_score column
#      for the sign-flipped downregulated values so that the source data is
#      not overwritten.
results_combined <- rbind(
  results_up   %>% mutate(Direction = "Upregulated",   plot_score =  obs_score),
  results_down %>% mutate(Direction = "Downregulated", plot_score = -obs_score)
)

# FIX: Derive cell-type display order from upregulated scores only,
#      so the order reflects true positive enrichment magnitude rather than
#      a mixture of positive and sign-flipped values.
cell_order <- results_up %>%
  arrange(desc(obs_score)) %>%
  pull(CellType)

# Compute a dynamic y-axis range so breaks stay meaningful regardless of
# the actual range of specificity scores in this dataset
max_val <- max(abs(results_combined$plot_score), na.rm = TRUE) * 1.1
break_step <- round(max_val / 4, 2)   # ~5 evenly spaced breaks
y_breaks <- seq(-max_val, max_val, by = break_step)

plot_data <- results_combined %>%
  mutate(
    CellType   = factor(CellType, levels = cell_order),
    Direction  = factor(Direction, levels = c("Upregulated", "Downregulated")),
    # Significance asterisk shown only when FDR-adjusted p < 0.05
    sig_label  = ifelse(p_adj < 0.05, "*", ""),
    color_group = case_when(
      p_adj < 0.05 & Direction == "Upregulated"   ~ "up_sig",
      p_adj < 0.05 & Direction == "Downregulated" ~ "down_sig",
      TRUE ~ "ns"
    )
  )

ggplot(plot_data, aes(x = CellType, y = plot_score, fill = color_group)) +
  geom_col(width = 0.7) +
  # Dashed reference line at y = 0 separating up- and down-regulated bars
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.5, linetype = "dashed") +
  # Place significance asterisks just beyond the tip of each bar
  geom_text(
    aes(
      label = sig_label,
      y     = plot_score + sign(plot_score) * max_val * 0.02,
      vjust = ifelse(Direction == "Upregulated", 0, 1)
    ),
    colour = "black",
    size   = 5
  ) +
  scale_y_continuous(
    name   = "Mean Specificity Score",
    breaks = y_breaks,
    # Display absolute values on the axis so both directions read positively
    labels = abs(y_breaks),
    limits = c(-max_val, max_val)
  ) +
  scale_fill_manual(
    values = c("up_sig" = "#E41A1C", "down_sig" = "#377EB8", "ns" = "grey80"),
    breaks = c("up_sig", "down_sig", "ns"),
    labels = c("Upregulated (FDR < 0.05)",
               "Downregulated (FDR < 0.05)",
               "Not significant")
  ) +
  labs(
    title    = "Cell-Type Enrichment Analysis (EWCE)",
    subtitle = "Top 600 upregulated vs downregulated PD genes",
    fill     = "Significance"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x        = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.title.x       = element_blank(),
    plot.title         = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle      = element_text(hjust = 0.5),
    legend.position    = "top"
  )