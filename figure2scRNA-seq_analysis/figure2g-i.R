# ============================================================
# Mouse MPTP Parkinson's Disease Model snRNA-seq Pipeline
# Dataset: CNP0000892 MPTP mouse model
# Note: Mouse mitochondrial genes use lowercase prefix (^mt-)
# ============================================================

# ── Dependencies ─────────────────────────────────────────────
library(Seurat)           
library(harmony)          
library(tidyverse)        
library(cowplot)          
library(patchwork)        
library(pheatmap)         
library(ggthemes)         
library(tidydr)           
library(SingleR)          
library(celldex)          
library(clusterProfiler)  
library(org.Mm.eg.db)     
library(org.Hs.eg.db)     
library(R.utils)          
library(data.table)      
library(readr)           
library(Matrix)         

# Increase global variable size limit for parallel workers
options(future.globals.maxSize = 150 * 1024^3)  # 150 GB

# ============================================================
# STEP 1 — Load 10X Data and Create Seurat Objects
# ============================================================

# List all sample subdirectories in the raw data folder
dir_name <- list.files('/data/h01005/hc/MPTP/rowdata/')
print(dir_name)  # Diagnostic: confirm expected sample directories are present

# Build one Seurat object per sample from 10X CellRanger output
scRNAlist <- list()
for (i in seq_along(dir_name)) {
  counts <- Read10X(
    data.dir = paste0('/data/h01005/hc/MPTP/rowdata/', dir_name[i])
  )
  scRNAlist[[i]] <- CreateSeuratObject(
    counts       = counts,
    project      = dir_name[i],
    min.cells    = 3,    # Retain genes detected in at least 3 cells
    min.features = 200   # Retain cells with at least 200 detected genes
  )
}

# ============================================================
# STEP 2 — Quality Control Filtering
# ============================================================

for (i in seq_along(scRNAlist)) {
  sc <- scRNAlist[[i]]
  
  # Compute mitochondrial transcript percentage per cell.
  # Mouse mitochondrial genes use lowercase prefix "mt-" (unlike human "MT-")
  sc[["mt_percent"]] <- PercentageFeatureSet(sc, pattern = "^mt-")
  
  # Filter cells by total UMI count using Seurat's pre-computed nCount_RNA.
  # 500–6000 UMIs: removes empty droplets (too low) and likely doublets (too high)
  sc <- subset(sc, nCount_RNA >= 500 & nCount_RNA <= 6000)
  
  scRNAlist[[i]] <- sc
  rm(sc)
  gc()
}

# Remove cells with high mitochondrial content (>10%) — likely damaged or dying
scRNAlist <- lapply(scRNAlist, function(x) subset(x, mt_percent < 10))

# ============================================================
# STEP 3 — Per-Sample Normalization and PCA
# ============================================================

for (i in seq_along(scRNAlist)) {
  scRNA <- scRNAlist[[i]]
  
  # SCTransform: regularized negative-binomial regression that handles
  # normalization, variable gene selection, and scaling in one step
  scRNA <- SCTransform(scRNA)
  
  # PCA on top variable genes identified by SCTransform
  scRNA <- RunPCA(scRNA, npcs = 30)
  
  scRNAlist[[i]] <- scRNA
  rm(scRNA)
  gc()
}

# ============================================================
# STEP 4 — SCT-based Cross-Sample Integration
# ============================================================

# Select the most informative shared features across all samples
integration_features <- SelectIntegrationFeatures(
  object.list = scRNAlist,
  nfeatures   = 3000
)

# Prepare SCT residuals for anchor-based integration
scRNAlist <- PrepSCTIntegration(
  object.list     = scRNAlist,
  anchor.features = integration_features
)

# Find integration anchors using SCT-normalized data
combined_anchors <- FindIntegrationAnchors(
  object.list          = scRNAlist,
  normalization.method = "SCT",
  anchor.features      = integration_features
)

# Integrate all samples into a single batch-corrected object
combin_integrated <- IntegrateData(
  anchorset            = combined_anchors,
  normalization.method = "SCT",
  dims                 = 1:30
)

# ============================================================
# STEP 5 — Dimensionality Reduction and Clustering
# ============================================================

# PCA on integrated assay; 50 PCs computed to allow flexible selection
combin_integrated <- RunPCA(combin_integrated, npcs = 50)

# Use top 30 PCs for neighborhood graph and clustering
# (selected based on elbow plot inspection)
combin_integrated <- FindNeighbors(combin_integrated, dims = 1:30)
combin_integrated <- FindClusters(combin_integrated, resolution = 0.4)
combin_integrated <- RunUMAP(combin_integrated, dims = 1:30)

# Preliminary UMAP view before annotation
p_unannotated <- DimPlot(
  combin_integrated,
  reduction = "umap",
  label     = TRUE,
  repel     = TRUE,
  pt.size   = 0.8
) +
  theme_dr(
    xlength = 0.22,
    ylength = 0.22,
    arrow   = grid::arrow(length = unit(0.15, "inches"), type = "closed")
  ) +
  theme(panel.grid = element_blank())
print(p_unannotated)

# ============================================================
# STEP 6 — Initial Marker Detection for Cluster Annotation
# ============================================================

# Switch to SCT assay and prepare models before differential testing
DefaultAssay(combin_integrated) <- "SCT"
combin_integrated <- PrepSCTFindMarkers(combin_integrated)

# Find positive marker genes for each cluster
# min.pct = 0.25: gene expressed in ≥ 25% of cluster cells
# logfc.threshold = 0.25: minimum average log2 fold-change filter
markers_tmp <- FindAllMarkers(
  combin_integrated,
  only.pos        = TRUE,
  min.pct         = 0.25,
  logfc.threshold = 0.25,
  verbose         = TRUE
)

# Extract top 60 markers per cluster for manual cell type annotation.
# NOTE: top_n() is deprecated since dplyr 1.0.0; use slice_max() instead.
top60 <- markers_tmp %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 60)

# Print marker genes per cluster for review on the ACT annotation portal
tmp <- split(top60$gene, top60$cluster)
cat(unlist(lapply(names(tmp), function(x) {
  paste0("cluster", x, ": ", paste(tmp[[x]], collapse = ", "))
})), sep = "\n")

# ============================================================
# STEP 7 — Manual Cell Type Annotation (27 clusters)
# ============================================================

# Cell type labels assigned in cluster order (verified against levels())
# MSCs = Mesenchymal Stem Cells
new.cluster.ids <- c(
  "Oligodendrocyte",               # Cluster 0
  "Neuron",                        # Cluster 1
  "Astrocyte",                     # Cluster 2
  "Oligodendrocyte",               # Cluster 3
  "Neuron",                        # Cluster 4
  "Neuron",                        # Cluster 5
  "Microglial cell",               # Cluster 6
  "Oligodendrocyte precursor cell",# Cluster 7
  "Neuron",                        # Cluster 8
  "Oligodendrocyte",               # Cluster 9
  "Neuron",                        # Cluster 10
  "Neuron",                        # Cluster 11
  "Neuron",                        # Cluster 12
  "Neuron",                        # Cluster 13
  "Neuron",                        # Cluster 14
  "Astrocyte",                     # Cluster 15
  "Neuron",                        # Cluster 16
  "MSCs",                          # Cluster 17
  "Neuron",                        # Cluster 18
  "Neuron",                        # Cluster 19
  "Neuron",                        # Cluster 20
  "Astrocyte",                     # Cluster 21
  "MSCs",                          # Cluster 22
  "Microglial cell",               # Cluster 23
  "Oligodendrocyte",               # Cluster 24
  "Neuron",                        # Cluster 25
  "Neuron"                         # Cluster 26
)

# Sanity check: label count must match cluster count before renaming
message("Clusters detected: ", length(levels(combin_integrated)),
        " | Labels provided: ", length(new.cluster.ids))
stopifnot(length(new.cluster.ids) == length(levels(combin_integrated)))

names(new.cluster.ids) <- levels(combin_integrated)
combin_integrated      <- RenameIdents(combin_integrated, new.cluster.ids)

# Store cell type annotation in metadata for downstream use
combin_integrated@meta.data$cells <- combin_integrated@active.ident

# UMAP after annotation
p_annotated <- DimPlot(
  combin_integrated,
  reduction = "umap",
  label     = TRUE,
  repel     = TRUE,
  pt.size   = 0.8
) +
  theme_dr(
    xlength = 0.22,
    ylength = 0.22,
    arrow   = grid::arrow(length = unit(0.15, "inches"), type = "closed")
  ) +
  theme(panel.grid = element_blank())
print(p_annotated)

# ============================================================
# STEP 8 — Re-run Marker Detection Per Annotated Cell Type
# ============================================================

# NOTE: This second FindAllMarkers call uses the renamed cell type identities
# (not original cluster numbers), so it returns markers that distinguish
# broad cell types — useful for validating annotation and downstream analysis.
markers <- FindAllMarkers(
  combin_integrated,
  only.pos        = TRUE,
  min.pct         = 0.25,
  logfc.threshold = 0.25,
  verbose         = TRUE
)

# Remove any rows with NA values (can occur in edge-case clusters)
markers <- na.omit(markers)

# ============================================================
# STEP 9 — Assign Disease Group Labels
# ============================================================

combin_integrated@meta.data$group <- NA
con_samples <- c("blank_1", "blank_2", "blank_3")   # Saline control samples
pd_samples  <- c("model_1", "model_2", "model_3")   # MPTP-treated PD model samples

combin_integrated@meta.data$group[
  combin_integrated@meta.data$orig.ident %in% pd_samples]  <- "pd"
combin_integrated@meta.data$group[
  combin_integrated@meta.data$orig.ident %in% con_samples] <- "con"

# Validate all cells received a group label
n_unassigned <- sum(is.na(combin_integrated@meta.data$group))
if (n_unassigned > 0) {
  warning(n_unassigned, " cell(s) have no group label. Check orig.ident values.")
} else {
  message("All cells successfully assigned to a group.")
}

# ============================================================
# STEP 10 — Cell Type Proportion Stacked Bar Chart (PD vs Control)
# ============================================================

meta_data <- combin_integrated@meta.data

# Define all cell types to include in proportion analysis
cell_types_to_compare <- c(
  "Oligodendrocyte precursor cell",
  "Astrocyte", "Microglial cell", "Neuron",
  "Oligodendrocyte", "MSCs"
)

# Filter metadata to retain only the target cell types
filtered_meta       <- meta_data[meta_data$cells %in% cell_types_to_compare, ]
filtered_meta$group <- as.factor(filtered_meta$group)
filtered_meta$cells <- as.factor(filtered_meta$cells)

# Count cells per group–cell type combination
# .groups = "drop" prevents retained grouping from affecting downstream mutate()
cell_counts <- filtered_meta %>%
  group_by(group, cells) %>%
  summarise(count = n(), .groups = "drop")

# Compute within-group proportions.
# NOTE: use an explicit named column (total) — avoids the fragile backtick
# column name created by the bare mutate(sum(count)) pattern.
bar_per <- cell_counts %>%
  group_by(group) %>%
  mutate(
    total   = sum(count),
    percent = count / total
  ) %>%
  ungroup()

# Define a consistent color palette for all cell types
cell_type_colors <- c(
  "Oligodendrocyte"               = "red",
  "Astrocyte"                     = "#F78000",
  "Oligodendrocyte precursor cell" = "yellow",
  "Neuron"                        = "pink",
  "Microglial cell"               = "#3176B7",
  "MSCs"                          = "brown"
)

# Stacked bar chart: cell type proportions per disease group.
# NOTE: theme_few() is applied first as the base theme; all custom theme()
# overrides are consolidated into a single block afterward — this ensures
# no settings are silently discarded by a later theme call.
p_bar <- ggplot(bar_per, aes(x = group, y = percent, fill = cells)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  # Optionally display raw counts inside bar segments (currently disabled):
  # geom_text(aes(label = count), position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = cell_type_colors) +
  labs(
    title = "Cell Type Proportions by Group",
    x     = "Group",
    y     = "Proportion",
    fill  = "Cell Type"
  ) +
  theme_few() +                                    # Apply base theme first
  theme(                                           # Consolidated custom overrides
    plot.title       = element_text(size = 12, hjust = 0.5),
    axis.ticks       = element_line(linetype = "blank"),
    legend.position  = "top",
    panel.grid.minor = element_line(colour = NA, linetype = "blank"),
    panel.background = element_rect(fill = NA),
    plot.background  = element_rect(colour = NA)
  )
print(p_bar)

# ============================================================
# STEP 11 — Marker Gene Dot Plot and Heatmaps
# ============================================================

# Curated mouse marker genes for major cell types (mouse gene symbols, mixed case):
# Opalin, Mal       → Oligodendrocytes
# Gm10754, Gm30382  → Oligodendrocyte-associated (predicted/novel transcripts)
# Aqp4, Slc39a12    → Astrocytes
# Dock8, P2ry12     → Microglial cells
# Vcan, Vxn         → Oligodendrocyte precursor cells
# Lama1, Cped1      → MSCs / Vascular-associated cells
marker_genes <- c(
  "Opalin", "Mal",
  "Gm10754", "Gm30382",
  "Aqp4", "Slc39a12",
  "Dock8", "P2ry12",
  "Vcan", "Vxn",
  "Lama1", "Cped1"
)

# Dot plot: dot size = fraction of cells expressing gene; color = scaled mean expression
top1_dotplot <- DotPlot(combin_integrated, features = marker_genes) +
  theme(axis.text.x = element_text(size = 6, angle = 45, vjust = 1, hjust = 1))
print(top1_dotplot)

# ── Heatmap Preparation ──────────────────────────────────────

# Extract scaled expression values from dot plot and reshape to wide format
# Result: rows = cell type clusters, columns = marker genes
dotplot_data <- top1_dotplot$data

heatmap_data <- dotplot_data %>%
  dplyr::select(features.plot, id, avg.exp.scaled) %>%
  pivot_wider(names_from = features.plot, values_from = avg.exp.scaled)

# Move cluster ID labels from a column into row names for pheatmap compatibility
heatmap_data <- column_to_rownames(heatmap_data, var = "id")

# ── Heatmap 1: Cell types (rows) × Genes (columns) ──────────
heatmap1 <- pheatmap(
  heatmap_data,
  cellwidth     = 10,
  cellheight    = 15,
  cluster_rows  = FALSE,
  cluster_cols  = FALSE,
  legend_breaks = c(-2, -1, 0, 1, 2),
  color         = colorRampPalette(c("#d3eafa", "white", "#e85b3f"))(100)
)

# ── Heatmap 2: Genes (rows) × Cell types (columns) ──────────
# Transpose so genes become rows and cell types become columns.
# The id column was already moved to row names above — no further
# column removal is needed before transposing.
heatmap_data_transposed <- t(heatmap_data)

heatmap2 <- pheatmap(
  heatmap_data_transposed,
  cellwidth     = 15,
  cellheight    = 10,
  cluster_rows  = FALSE,
  cluster_cols  = FALSE,
  angle_col     = 45,
  legend_breaks = c(-2, -1, 0, 1, 2),
  color         = colorRampPalette(c("#d3eafa", "white", "#e85b3f"))(100)
)