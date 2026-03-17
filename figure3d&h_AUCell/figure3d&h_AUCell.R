# ============================================================
# AUCell Gene Set Activity Scoring on Integrated snRNA-seq Data
# Input: scRNA_combined.RDS — fully annotated, Harmony-integrated object
# Gene sets: oxidative stress signatures from GSE7621 and GSE49036
# ============================================================

# ── Dependencies ─────────────────────────────────────────────
library(Seurat)           
library(AUCell)           
library(GSEABase)        
library(dplyr)           
library(ggplot2)          
library(tidydr)           
library(gplots)           
library(tidyverse)        
library(cowplot)         
library(ggthemes)        
library(clusterProfiler)  

# ============================================================
# STEP 1 — Load Data and Subset to PD and Control Cells
# ============================================================

scRNA_combined <- readRDS("/data/h01005/hc/RE/scRNA_combined.RDS")

# Retain only PD and healthy control cells; exclude LBD, aged, or other groups
scRNA_pdcon <- subset(scRNA_combined, subset = group %in% c("pd", "con"))
rm(scRNA_combined)
gc()

# Cache metadata for fast barcode-level lookups in the loop below
metadata <- scRNA_pdcon@meta.data
message("Cells retained (PD + CON): ", ncol(scRNA_pdcon))

# Confirm metadata column names before referencing them in the loop
message("Metadata columns: ", paste(colnames(metadata), collapse = ", "))

# ============================================================
# STEP 2 — Extract Expression Matrix for AUCell
# ============================================================

# AUCell requires a genes × cells matrix.
# Use the SCT log-normalized data layer (Seurat v5 compatible accessor).
# NOTE: Seurat v4 used @assays$SCT@data; Seurat v5 requires GetAssayData().
expr_matrix <- GetAssayData(scRNA_pdcon, assay = "SCT", layer = "data")
message("Expression matrix dimensions: ",
        nrow(expr_matrix), " genes × ", ncol(expr_matrix), " cells")
gc()

# ============================================================
# STEP 3 — Build Gene Rankings Across All Cells
# ============================================================

# AUCell_buildRankings ranks all genes by expression in each cell.
# These rankings are reused for every gene set scored in STEP 5.
# nCores = 4: parallel processing; adjust to available CPU cores.
# plotStats = FALSE: suppresses the interactive ranking distribution plot.
cells_rankings <- AUCell_buildRankings(expr_matrix, nCores = 4, plotStats = FALSE)
gc()

# ============================================================
# STEP 4 — Define Oxidative Stress Gene Sets
# ============================================================

# Gene sets derived from the two bulk RNA-seq analyses:
#   oxidative_genes_vector1: DEGs from GSE7621 (stress response, apoptosis,
#                             transcription factors, NF-κB pathway)
#   oxidative_genes_vector2: DEGs from GSE49036 (Nrf2 pathway, oxidative
#                             stress markers, immune signaling)
#
# IMPORTANT: Both gene sets are scored in the same loop below.
# Inspect per-set results separately — do not combine scores across sets.

oxidative_genes_vector1 <- c(
  "ATF4",  "CAPN2", "BANF1", "PNPT1", "PTPRK",
  "FOXO3", "SIRPA", "STAU1", "ETV5",  "RELA",
  "KEAP1", "ATM",   "DHFR",  "ABL1"
)

oxidative_genes_vector2 <- c(
  "NFE2L2", "CYP1B1", "HMOX1", "PTGS2", "TXNIP",
  "ATM",    "CRYAB",  "PTPRK", "KAT2B", "DDR2",
  "CAT",    "GPR37"
)

# Combine into a named list; loop variable `i` will be each list name
oxidative_genes_vector <- list(
  oxidative_genes_GSE7621  = oxidative_genes_vector1,
  oxidative_genes_GSE49036 = oxidative_genes_vector2
)

# ============================================================
# STEP 5 — AUCell Scoring Loop: One Iteration Per Gene Set
# ============================================================

# Set output directory once — update here to change all output paths
OUTPUT_DIR <- "/data/h01005/hc/AUCell/20251221/"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
message("Output directory: ", OUTPUT_DIR)

# Known cell types present in the annotated object
# (used for consistent color mapping across all plots)
cell_type_colors <- c(
  "Oligodendrocyte"               = "red",
  "Astrocyte"                     = "#F78000",
  "Oligodendrocyte precursor cell" = "yellow",
  "Neuron"                        = "pink",
  "Microglial cell"               = "#3176B7",
  "Fibroblast"                    = "brown"
)

for (i in names(oxidative_genes_vector)) {
  
  message("\n── Processing gene set: ", i, " ──")
  
  current_genes <- oxidative_genes_vector[[i]]
  
  # ── Gene presence check ──────────────────────────────────
  # AUCell silently scores on whatever subset of genes it finds.
  # Report coverage so the analyst knows if the gene set is well-represented.
  genes_found   <- intersect(current_genes, rownames(expr_matrix))
  genes_missing <- setdiff(current_genes, rownames(expr_matrix))
  message("  Genes in set: ",    length(current_genes),
          " | Found in matrix: ", length(genes_found),
          " | Missing: ",         length(genes_missing))
  if (length(genes_missing) > 0) {
    warning("Missing genes for ", i, ": ",
            paste(genes_missing, collapse = ", "))
  }
  
  # ── AUC Calculation ───────────────────────────────────────
  # aucMaxRank = top 5% of all ranked genes: defines the "expressed" threshold.
  # The gene set name is passed dynamically using setNames() so the AUC matrix
  # row is named after the gene set (not hardcoded as "geneset").
  cells_AUC <- AUCell_calcAUC(
    geneSets   = setNames(list(genes_found), i),
    rankings   = cells_rankings,
    aucMaxRank = nrow(cells_rankings) * 0.05
  )
  
  # Extract AUC scores into a data.frame (1 row per cell)
  auc_matrix <- getAUC(cells_AUC)
  auc_df <- data.frame(
    Cell = colnames(auc_matrix),
    AUC  = auc_matrix[1, ]
  )
  
  # ── Attach cell type and group labels from metadata ───────
  # metadata is indexed by barcode; any unmatched barcodes return NA
  auc_df$Cluster <- metadata[auc_df$Cell, "cells"]  # Cell type annotation
  auc_df$Group   <- metadata[auc_df$Cell, "group"]  # Disease group (pd / con)
  
  # Validate that no barcodes were silently lost in the metadata join
  n_na_cluster <- sum(is.na(auc_df$Cluster))
  n_na_group   <- sum(is.na(auc_df$Group))
  if (n_na_cluster > 0 || n_na_group > 0) {
    warning("Metadata join produced NAs for gene set ", i,
            ": Cluster NAs = ", n_na_cluster,
            ", Group NAs = ",   n_na_group,
            ". Check barcode consistency between AUC results and metadata.")
  }
  
  # ── Sort clusters by mean AUC for ordered boxplot display ─
  # .groups = "drop" prevents retained grouping from affecting arrange()
  auc_summary <- auc_df %>%
    group_by(Cluster) %>%
    summarise(mean_AUC = mean(AUC, na.rm = TRUE), .groups = "drop")
  
  cluster_order  <- auc_summary %>% arrange(mean_AUC) %>% pull(Cluster)
  auc_df$Cluster <- factor(auc_df$Cluster, levels = cluster_order)
  
  # ── Save per-cell AUC scores to CSV ───────────────────────
  out_csv <- paste0(OUTPUT_DIR, i, ".csv")
  write.csv(auc_df, out_csv, row.names = FALSE)
  message("  AUC scores saved to: ", out_csv)
  
  # ── Boxplot: AUC distribution per cell type ───────────────
  p <- ggplot(auc_df, aes(x = Cluster, y = AUC, fill = Cluster)) +
    geom_boxplot(outlier.size = 0.5) +
    labs(
      title = paste0("AUCell Oxidative Stress Score — ", i),
      x     = "Cell Type",
      y     = "AUC Score"
    ) +
    scale_fill_manual(values = cell_type_colors) +
    theme_minimal() +
    theme(
      axis.text.x  = element_text(angle = 45, hjust = 1),
      plot.title   = element_text(hjust = 0.5, face = "bold"),
      legend.position = "none"   # Legend is redundant when x-axis shows cell types
    )
  
  out_pdf <- paste0(OUTPUT_DIR, i, ".pdf")
  ggsave(out_pdf, plot = p, width = 10, height = 6)
  message("  Plot saved to: ", out_pdf)
}

message("\nAll gene sets processed. Outputs written to: ", OUTPUT_DIR)