# ============================================================
# Mouse Midbrain snRNA-seq Analysis Pipeline
# Dataset: CNP0000892 (CNS midbrain samples)
# ============================================================

# ── Dependencies ─────────────────────────────────────────────
library(Seurat)           
library(harmony)          
library(tidyverse)        
library(cowplot)         
library(patchwork)        
library(R.utils)          
library(celldex)          
library(SingleR)          
library(pheatmap)        
library(clusterProfiler)  
library(org.Mm.eg.db)     
library(org.Hs.eg.db)    

# Organism variable for enrichment analysis (set to human here)
# NOTE: defined here for use in downstream enrichment steps
organism <- "org.Hs.eg.db"

# Increase global variable size limit for parallel workers
options(future.globals.maxSize = 5000 * 1024^2)  # 5 GB

# ============================================================
# STEP 1 — Load 10X Data and Create Seurat Objects
# ============================================================

# List all sample subdirectories in the midbrain data folder
dir_name <- list.files('D:/R/GEO/CNP0000892/CNS/mid')
print(dir_name)  # Diagnostic: confirm expected samples are found

# Build one Seurat object per sample from 10X output directories
scRNAlist <- list()
for (i in seq_along(dir_name)) {
  counts <- Read10X(
    data.dir = paste0("D:/R/GEO/CNP0000892/CNS/mid/", dir_name[i],
                      "/filtered_feature_bc_matrix")
  )
  scRNAlist[[i]] <- CreateSeuratObject(
    counts       = counts,
    project      = dir_name[i],
    min.cells    = 3,    # Retain genes detected in at least 3 cells
    min.features = 200   # Retain cells with at least 200 detected genes
  )
}
print(scRNAlist)  # Diagnostic: review sample sizes and feature counts

# ============================================================
# STEP 2 — Quality Control Filtering
# ============================================================

for (i in seq_along(scRNAlist)) {
  sc <- scRNAlist[[i]]
  
  # Compute mitochondrial transcript percentage (high % indicates damaged cells)
  sc[["mt_percent"]] <- PercentageFeatureSet(sc, pattern = "^MT-")
  
  # Filter cells by total UMI count using Seurat's pre-computed nCount_RNA
  # 500–6000 UMIs: removes empty droplets (low) and likely doublets (high)
  sc <- subset(sc, nCount_RNA >= 500 & nCount_RNA <= 6000)
  
  scRNAlist[[i]] <- sc
  rm(sc)
}

# Remove cells with high mitochondrial content (>2%) — likely dying or damaged
scRNAlist <- lapply(scRNAlist, function(x) subset(x, mt_percent < 2))

# ============================================================
# STEP 3 — Per-Sample Normalization and PCA
# ============================================================

for (i in seq_along(dir_name)) {
  sc <- scRNAlist[[i]]
  
  # SCTransform: regularized negative-binomial regression for normalization,
  # variable gene selection, and scaling in a single step
  sc <- SCTransform(sc)
  
  # PCA on top variable genes; 50 PCs captured for flexible downstream selection
  sc <- RunPCA(sc, npcs = 50)
  
  scRNAlist[[i]] <- sc
  rm(sc)
  gc()
}

# Backup the processed list before integration (preserved for potential reuse)
scrnalist2 <- scRNAlist

# ============================================================
# STEP 4 — SCT-based Cross-Sample Integration
# ============================================================

# Select the most informative shared features across all samples for integration
features <- SelectIntegrationFeatures(object.list = scRNAlist, nfeatures = 3000)

# Prepare SCT residuals for integration
scRNAlist <- PrepSCTIntegration(
  object.list    = scRNAlist,
  assay          = "SCT",
  anchor.features = features
)

# Find integration anchors using SCT-normalized data
combined_anchors <- FindIntegrationAnchors(
  object.list        = scRNAlist,
  normalization.method = "SCT",
  anchor.features    = features
)

# Integrate all samples into a single corrected object
combin_integrated <- IntegrateData(
  anchorset          = combined_anchors,
  dims               = 1:20,
  normalization.method = "SCT"
)

# ============================================================
# STEP 5 — Dimensionality Reduction and Clustering
# ============================================================

combin_integrated <- RunPCA(combin_integrated, npcs = 20)
combin_integrated <- FindNeighbors(combin_integrated, dims = 1:20)
combin_integrated <- FindClusters(combin_integrated, resolution = 1)

# UMAP with custom neighbor count for midbrain dataset
combin_integrated <- RunUMAP(
  combin_integrated,
  dims       = 1:20,
  reduction  = "pca",
  n.neighbors = 60,
  min.dist   = 0.3
)

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
# STEP 6 — Assign Disease Group Labels
# ============================================================

# Inspect all unique sample IDs before group assignment
print(unique(combin_integrated@meta.data$orig.ident))

combin_integrated@meta.data$group <- NA
CON_samples <- "CNS0195744"  # Healthy control sample
PD_samples  <- "CNS0195747"  # Parkinson's Disease sample

combin_integrated@meta.data$group[
  combin_integrated@meta.data$orig.ident %in% PD_samples]  <- "pd"
combin_integrated@meta.data$group[
  combin_integrated@meta.data$orig.ident %in% CON_samples] <- "con"

# Validate all cells received a group label
n_unassigned <- sum(is.na(combin_integrated@meta.data$group))
if (n_unassigned > 0) {
  warning(n_unassigned, " cell(s) have no group label. Check orig.ident values.")
} else {
  message("All cells successfully assigned to a group.")
}

# ============================================================
# STEP 7 — Marker Gene Detection
# ============================================================

# PrepSCTFindMarkers is REQUIRED before FindAllMarkers when SCT normalization
# was used — it re-corrects counts across samples for valid statistical testing
combin_integrated <- PrepSCTFindMarkers(combin_integrated)

# Find positive marker genes per cluster
# min.pct = 0.25: gene expressed in ≥ 25% of cluster cells
# logfc.threshold = 0.25: minimum average log2 fold-change filter
markers <- FindAllMarkers(
  combin_integrated,
  only.pos       = TRUE,
  min.pct        = 0.25,
  logfc.threshold = 0.25,
  verbose        = TRUE
)

# Extract top 100 markers per cluster by avg log2FC for annotation review.
# NOTE: top_n() is deprecated since dplyr 1.0.0; use slice_max() instead.
top100 <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 100)

# Print top markers per cluster for manual review on the ACT annotation portal
tmp <- split(top100$gene, top100$cluster)
cat(unlist(lapply(names(tmp), function(x) {
  paste0("cluster", x, ": ", paste(tmp[[x]], collapse = ", "))
})), sep = "\n")

# ============================================================
# STEP 8 — Manual Cell Type Annotation
# ============================================================

# Cell type labels assigned in cluster order (matching levels of combin_integrated)
new.cluster.ids <- c(
  "Glutamatergic Neuron",            # Cluster 0
  "Oligodendrocyte",                 # Cluster 1
  "GABAergic Neuron",                # Cluster 2
  "GABAergic Neuron",                # Cluster 3
  "Astrocyte",                       # Cluster 4
  "Glutamatergic Neuron",            # Cluster 5
  "Oligodendrocyte",                 # Cluster 6
  "Microglial cell",                 # Cluster 7
  "Oligodendrocyte",                 # Cluster 8
  "GABAergic Neuron",                # Cluster 9
  "Glutamatergic Neuron",            # Cluster 10
  "Oligodendrocyte",                 # Cluster 11
  "Glutamatergic Neuron",            # Cluster 12
  "GABAergic Neuron",                # Cluster 13
  "Glutamatergic Neuron",            # Cluster 14
  "Oligodendrocyte precursor cell",  # Cluster 15
  "Dopaminergic Neuron",             # Cluster 16
  "Oligodendrocyte precursor cell",  # Cluster 17
  "Glutamatergic Neuron",            # Cluster 18
  "Glutamatergic Neuron",            # Cluster 19
  "Vascular cell",                   # Cluster 20
  "GABAergic Neuron",                # Cluster 21
  "AP",                              # Cluster 22 (Area Postrema / Astrocyte Precursor)
  "GABAergic Neuron",                # Cluster 23
  "Oligodendrocyte",                 # Cluster 24
  "Pericyte",                        # Cluster 25
  "Glutamatergic Neuron",            # Cluster 26
  "Glutamatergic Neuron",            # Cluster 27
  "Astrocyte",                       # Cluster 28
  "Serotonergic Neuron"              # Cluster 29
)

# Sanity check: label count must match cluster count before renaming
message("Clusters detected: ", length(levels(combin_integrated)),
        " | Labels provided: ", length(new.cluster.ids))
stopifnot(length(new.cluster.ids) == length(levels(combin_integrated)))

# Assign cluster-number-indexed names before renaming (used for marker mapping below)
names(new.cluster.ids) <- levels(combin_integrated)
combin_integrated      <- RenameIdents(combin_integrated, new.cluster.ids)

# Store cell type annotation in metadata for downstream subsetting and plotting
combin_integrated@meta.data$cells <- combin_integrated@active.ident

# Map cell type labels onto the markers table using original cluster IDs as keys.
# names(new.cluster.ids) were set to original numeric cluster levels above,
# so indexing by markers$cluster correctly retrieves the corresponding cell type.
markers_all             <- markers
markers_all$new_cluster <- new.cluster.ids[as.character(markers_all$cluster)]
markers                 <- markers_all

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
# STEP 9 — Cell Type Proportion Stacked Bar Chart (PD vs Control)
# ============================================================

meta_data <- combin_integrated@meta.data

# Define all cell types to include in proportion analysis
cell_types_to_compare <- c(
  "Oligodendrocyte", "Oligodendrocyte precursor cell",
  "Glutamatergic Neuron", "Astrocyte", "Pericyte",
  "GABAergic Neuron", "Microglial cell", "Serotonergic Neuron",
  "Vascular cell", "AP", "Dopaminergic Neuron"
)

# Filter metadata to retain only the target cell types
filtered_meta        <- meta_data[meta_data$cells %in% cell_types_to_compare, ]
filtered_meta$group  <- as.factor(filtered_meta$group)
filtered_meta$cells  <- as.factor(filtered_meta$cells)

# Count cells per group–cell type combination
# .groups = "drop" prevents retained grouping from affecting downstream mutate()
cell_counts <- filtered_meta %>%
  group_by(group, cells) %>%
  summarise(count = n(), .groups = "drop")

# Compute within-group proportions
# NOTE: use an explicit named column (total) to avoid fragile backtick column names
bar_per <- cell_counts %>%
  group_by(group) %>%
  mutate(
    total   = sum(count),
    percent = count / total
  ) %>%
  ungroup()

# Define a consistent color palette for all cell types
cell_type_colors <- c(
  "Oligodendrocyte"              = "red",
  "Astrocyte"                    = "#F78000",
  "Oligodendrocyte precursor cell" = "yellow",
  "Glutamatergic Neuron"         = "green",
  "Pericyte"                     = "#41BED1",
  "GABAergic Neuron"             = "pink",
  "Microglial cell"              = "#3176B7",
  "Serotonergic Neuron"          = "#e85b3f",
  "Vascular cell"                = "#d3eafa",
  "AP"                           = "brown",
  "Dopaminergic Neuron"          = "purple"
)

# Stacked bar chart: cell type proportions per disease group
# NOTE: theme_few() is applied first as the base theme; custom theme() overrides
# follow afterward — reversing this order would discard all custom settings.
p2 <- ggplot(bar_per, aes(x = group, y = percent, fill = cells)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  # Optionally display raw counts inside bar segments (currently disabled):
  # geom_text(aes(label = count), position = position_stack(vjust = 0.5)) +
  labs(
    title = "Cell Type Proportions by Group",
    x     = "Group",
    y     = "Proportion",
    fill  = "Cell Type"
  ) +
  scale_fill_manual(values = cell_type_colors) +
  theme_few() +                                    # Apply base theme first
  theme(                                           # Then apply custom overrides
    plot.title       = element_text(size = 12, hjust = 0.5),
    axis.ticks       = element_line(linetype = "blank"),
    legend.position  = "top",
    panel.grid.minor = element_line(colour = NA, linetype = "blank"),
    panel.background = element_rect(fill = NA),
    plot.background  = element_rect(colour = NA)
  )
print(p2)

# ============================================================
# STEP 10 — Top Marker Gene Dot Plot and Heatmap
# ============================================================

# Extract top 10 marker genes per cluster for dot plot visualization
top10 <- markers_all %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10)

# Export top 10 gene table (cluster + gene columns) for record-keeping
write.table(
  top10[, c("cluster", "gene")],
  "D:/R/GEO/CNP_RE/new/top10gene.txt",
  row.names = FALSE
)

# Dot plot: dot size = fraction of cells expressing gene; color = scaled mean expression
top10_dotplot <- DotPlot(combin_integrated, features = unique(top10$gene)) +
  theme(axis.text.x = element_text(size = 6, angle = 45, vjust = 1, hjust = 1))
print(top10_dotplot)

# ── Heatmap: Genes (columns) × Cell Types (rows) ─────────────

# Extract scaled expression values from dot plot and reshape to wide format
dotplot_data <- top10_dotplot$data

heatmap_data <- dotplot_data %>%
  dplyr::select(features.plot, id, avg.exp.scaled) %>%
  pivot_wider(names_from = features.plot, values_from = avg.exp.scaled)

# Move cluster ID labels from a column into row names for pheatmap compatibility
heatmap_data <- column_to_rownames(heatmap_data, var = "id")

# Plot heatmap: rows = cell type clusters, columns = top marker genes
pheatmap(
  heatmap_data,
  cellwidth     = 10,
  cellheight    = 15,
  cluster_rows  = FALSE,
  cluster_cols  = FALSE,
  legend_breaks = c(-2, -1, 0, 1, 2),
  color         = colorRampPalette(c("#d3eafa", "white", "#e85b3f"))(100)
)