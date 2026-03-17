# ============================================================
# Single-Cell RNA-seq Analysis Pipeline
# Dataset: GSE243639 (Parkinson's Disease snRNA-seq)
# Reference: Martirosyan et al., Mol Neurodegener 2024
# ============================================================

library(Seurat)
library(harmony)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(patchwork)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(data.table)

# Increase the maximum allowed size for global variables passed to parallel workers.
# Required when using future-based parallelism (e.g., via BiocParallel or plan(multisession)).
options(future.globals.maxSize = 50 * 1024^3)  # 50 GB limit

# ============================================================
# 1. Load and Clean Raw Count Matrix
# ============================================================

# Read the filtered count table (genes x cells); fread is faster than read.csv for large files
data <- fread("D://R//GEO//GSE243639//Filtered_count_table.csv")

# Convert to data.frame and set gene names as row names
data2 <- as.data.frame(data)
rownames(data2) <- data2$V1
data2 <- data2[, -1]  # Remove the gene-name column now stored as rownames

# Remove genes with zero total counts across all cells (uninformative features)
data2 <- data2[rowSums(data2) > 0, ]

# Remove duplicated gene names to prevent downstream errors in Seurat
data2 <- data2[!duplicated(rownames(data2)), ]

# ============================================================
# 2. Split Count Matrix by Sample
# ============================================================

# List available sample directories (used as ground-truth sample IDs)
dir_name <- list.files('D://R//GEO//GSE243639//GSE243639/')

# Extract sample prefixes from cell barcodes (everything before the first underscore)
# Assumes barcodes follow the pattern: <SampleID>_<Barcode>
level <- unique(gsub("_.*", "", colnames(data2)))

# Validate that sample IDs derived from barcodes exactly match directory names.
# NOTE: This check is order-sensitive. Sort both vectors before comparing
# to avoid false failures caused by ordering differences.
stopifnot(all(sort(level) == sort(dir_name)))

# Subset the count matrix for each sample and store in a named list
filtered_dfs <- list()
for (i in seq_along(level)) {
  filtered_dfs[[dir_name[i]]] <- data2 %>%
    dplyr::select(starts_with(level[i]))  # Select columns (cells) belonging to sample i
}

# ============================================================
# 3. Create Seurat Objects Per Sample
# ============================================================

se_list <- list()
for (i in seq_along(dir_name)) {
  # Create a Seurat object; project name is set to the sample ID for traceability
  se_list[[i]] <- CreateSeuratObject(filtered_dfs[[dir_name[i]]],
                                     project = dir_name[i])
}
names(se_list) <- dir_name

# ============================================================
# 4. Normalize and Reduce Each Sample (SCTransform + PCA)
# ============================================================

for (i in seq_along(dir_name)) {
  sc <- se_list[[i]]
  
  # SCTransform: regularized negative-binomial regression to normalize and stabilize variance.
  # nfeatures = 6000 retains more variable features for better integration across samples.
  sc <- SCTransform(sc, nfeatures = 6000)
  
  # PCA on the top variable features from SCTransform
  sc <- RunPCA(sc, npcs = 30)
  
  se_list[[i]] <- sc
  rm(sc)
  gc()  # Free memory after processing each sample
}

# ============================================================
# 5. SCT-based Integration Across Samples
# ============================================================

# Select the most informative features shared across all samples for integration
integration_features <- SelectIntegrationFeatures(object.list = se_list,
                                                  nfeatures = 3000)

# Prepare each Seurat object for SCT-based integration (scales residuals across samples)
se_list <- PrepSCTIntegration(object.list = se_list,
                              anchor.features = integration_features)

# Verify sample identities before selecting the reference sample
all_idents <- unlist(lapply(se_list, function(x) unique(x$orig.ident)))
print(all_idents)  # Diagnostic: confirm GSM7792123 is present in orig.ident

# Identify the list index containing the reference sample (GSM7792123).
# Stop early with an informative error if the reference sample is not found.
reference_index <- which(sapply(se_list, function(x) "GSM7792123" %in% unique(x$orig.ident)))
if (length(reference_index) == 0) {
  stop("Reference sample 'GSM7792123' not found in any Seurat object. Check sample IDs.")
}

# Find integration anchors using the reference-based approach
combined_anchors <- FindIntegrationAnchors(
  object.list = se_list,
  normalization.method = "SCT",
  anchor.features = integration_features,
  reference = reference_index
)
gc()

# Integrate data across all samples using the identified anchors
combin_integrated <- IntegrateData(
  anchorset = combined_anchors,
  normalization.method = "SCT",
  dims = 1:10
)
gc()

# ============================================================
# 6. Dimensionality Reduction and Clustering
# ============================================================

# PCA on the integrated assay
combin_integrated <- RunPCA(combin_integrated, npcs = 30)

# Use top 10 PCs as per the original publication (Martirosyan et al., 2024)
combin_integrated <- FindNeighbors(combin_integrated, dims = 1:10)
combin_integrated <- FindClusters(combin_integrated, resolution = 0.6)
combin_integrated <- RunUMAP(combin_integrated, dims = 1:10)

# Visualize clusters on UMAP embedding
DimPlot(combin_integrated, reduction = "umap")

# ============================================================
# 7. Marker Gene Detection
# ============================================================

# Switch to the SCT assay for marker gene detection (required after integration)
DefaultAssay(combin_integrated) <- "SCT"

# Re-prepare SCT model parameters across all samples before calling FindAllMarkers
combin_integrated <- PrepSCTFindMarkers(combin_integrated)

# Find positive marker genes for each cluster
# min.pct = 0.25: gene must be detected in ≥25% of cells in the cluster
# logfc.threshold = 0.25: minimum average log2 fold-change filter
markers_all <- FindAllMarkers(
  combin_integrated,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  verbose = TRUE
)

# Extract the top 100 marker genes per cluster ranked by average log2 fold-change.
# NOTE: top_n() is deprecated as of dplyr 1.0.0; use slice_max() instead.
top100 <- markers_all %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 100)

# Print top markers per cluster in a compact format for manual cell type annotation
tmp <- split(top100$gene, top100$cluster)
cat(unlist(lapply(names(tmp), function(x) {
  paste0("cluster", x, ": ", paste(tmp[[x]], collapse = ", "))
})), sep = "\n")

# ============================================================
# 8. Cell Type Annotation
# ============================================================

# Manually assigned cell type labels for each cluster (based on marker gene inspection)
new.cluster.ids <- c(
  "Oligodendrocyte",             # Cluster 0
  "Oligodendrocyte",             # Cluster 1
  "Oligodendrocyte",             # Cluster 2
  "Microglial cell",             # Cluster 3
  "Astrocyte",                   # Cluster 4
  "Oligodendrocyte",             # Cluster 5
  "Oligodendrocyte precursor cell", # Cluster 6
  "Astrocyte",                   # Cluster 7
  "Neurons",                     # Cluster 8
  "Oligodendrocyte",             # Cluster 9
  "Microglial cell",             # Cluster 10
  "Astrocyte",                   # Cluster 11
  "Microglial cell",             # Cluster 12
  "Astrocyte",                   # Cluster 13
  "Oligodendrocyte",             # Cluster 14
  "Neurons",                     # Cluster 15
  "Astrocyte",                   # Cluster 16
  "Vascular cell",               # Cluster 17
  "Neurons",                     # Cluster 18
  "Oligodendrocyte precursor cell", # Cluster 19
  "Astrocyte",                   # Cluster 20
  "Microglial cell",             # Cluster 21
  "Microglial cell",             # Cluster 22
  "Oligodendrocyte precursor cell"  # Cluster 23
)

# Sanity check: number of labels must match number of clusters
message("Number of clusters detected: ", length(levels(combin_integrated)))
stopifnot(length(new.cluster.ids) == length(levels(combin_integrated)))

# Rename cluster identities to cell type labels
names(new.cluster.ids) <- levels(combin_integrated)
combin_integrated <- RenameIdents(combin_integrated, new.cluster.ids)

# Store the cell type annotation in the metadata for downstream use
combin_integrated@meta.data$cells <- combin_integrated@active.ident