# ============================================================
# Oligodendrocyte Subcluster Analysis
# Input: scRNA_combined.RDS — fully annotated, integrated Seurat object
# Goal: Re-cluster oligodendrocytes, compare subcluster composition
#       between PD and CON, and identify subcluster marker genes
# ============================================================

# ── Dependencies ─────────────────────────────────────────────
library(Seurat)      
library(ggplot2)     
library(tidyverse)   
library(cowplot)     
library(ggthemes)    
library(dplyr)      

# ============================================================
# STEP 1 — Load Integrated Object and Extract Oligodendrocytes
# ============================================================

scRNA_combined <- readRDS("/data/h01005/hc/RE/scRNA_combined.RDS")

# Confirm the default assay — should be "SCT" if SCTransform was run
message("Default assay: ", DefaultAssay(scRNA_combined))

# Subset to oligodendrocytes only (annotation stored in 'cells' metadata column)
oligo <- subset(scRNA_combined, subset = cells %in% "Oligodendrocyte")

# Retain only PD and healthy control cells for the comparison
oligo_conpd <- subset(oligo, subset = group %in% c("con", "pd"))
message("Oligodendrocytes retained (CON + PD): ", ncol(oligo_conpd))

# Free the full integrated object from memory
rm(scRNA_combined, oligo)
gc()

# ============================================================
# STEP 2 — Re-normalize the Oligodendrocyte Subset
# ============================================================

# DELIBERATE CHOICE: SCTransform is re-run on the oligodendrocyte subset
# rather than reusing the whole-object SCT assay. This ensures that the
# variance stabilization and variable feature selection reflect only the
# biology of oligodendrocytes, not the full heterogeneous cell population.
#
# CAUTION: If scRNA_combined was already SCT-normalized, this step applies
# SCTransform a second time. This is intentional for subcluster analysis
# but means the new SCT assay replaces the inherited one for this object.
oligo_obj <- SCTransform(oligo_conpd, verbose = TRUE)

# Confirm SCT is now the default assay
DefaultAssay(oligo_obj) <- "SCT"
message("Default assay after SCTransform: ", DefaultAssay(oligo_obj))

# ============================================================
# STEP 3 — Dimensionality Reduction and Clustering
# ============================================================

# Use top 30 PCs in PCA; downstream steps use the top N_DIMS of these.
# N_DIMS = 20 was selected based on elbow plot inspection.
# Adjust N_DIMS if the elbow shifts for different subsets.
N_PCS  <- 30
N_DIMS <- 20  # Number of significant PCs used for neighbors and UMAP

oligo_obj <- RunPCA(oligo_obj, npcs = N_PCS, verbose = TRUE)

# Optionally inspect the elbow plot to validate N_DIMS before proceeding:
# ElbowPlot(oligo_obj, ndims = N_PCS)

oligo_obj <- FindNeighbors(oligo_obj, dims = 1:N_DIMS, verbose = TRUE)

# resolution = 0.3: relatively low to produce broad subclusters.
# Increase to 0.5–0.8 for finer subpopulation resolution.
oligo_obj <- FindClusters(oligo_obj, resolution = 0.3, verbose = TRUE)
oligo_obj <- RunUMAP(oligo_obj, dims = 1:N_DIMS, verbose = TRUE)

# ============================================================
# STEP 4 — Validate Cluster Count and Define Color Palette
# ============================================================

n_clusters <- length(levels(oligo_obj$seurat_clusters))
message("Number of subclusters detected: ", n_clusters)

# Color palette for clusters 0–7 (ordered consistently)
# Update this vector if the number of clusters changes
cluster_colors <- c(
  "0" = "#f57c6e",
  "1" = "#f2b56f",
  "2" = "#fae69e",
  "3" = "#84c3b7",
  "4" = "#71b7ed",
  "5" = "#b8aeeb",
  "6" = "#f2a7da",
  "7" = "#88d8db"
)

# Warn if the detected cluster IDs do not match the color palette keys
detected_clusters <- as.character(levels(oligo_obj$seurat_clusters))
missing_colors    <- setdiff(detected_clusters, names(cluster_colors))
if (length(missing_colors) > 0) {
  warning("No color defined for cluster(s): ",
          paste(missing_colors, collapse = ", "),
          ". Add entries to cluster_colors.")
}

# ============================================================
# STEP 5 — UMAP Visualization of Subclusters
# ============================================================

p_umap <- DimPlot(
  oligo_obj,
  reduction = "umap",
  group.by  = "seurat_clusters",
  label     = TRUE,
  repel     = TRUE,
  pt.size   = 0.5,
  cols      = cluster_colors
) +
  ggtitle("Oligodendrocyte Subclusters (CON + PD)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
print(p_umap)

# ============================================================
# STEP 6 — Compute Subcluster Proportions Per Group
# ============================================================

# NOTE: The object was already filtered to CON and PD in STEP 1.
# The group filter below is retained as an explicit safety check
# to guard against accidental reuse of this code with a wider object.
meta_for_plot <- oligo_obj@meta.data

# Count cells per group–cluster combination
# .groups = "drop" prevents retained grouping from affecting downstream mutate()
cell_counts <- meta_for_plot %>%
  group_by(group, seurat_clusters) %>%
  summarise(count = n(), .groups = "drop")

# Compute within-group proportions for each subcluster
bar_per <- cell_counts %>%
  group_by(group) %>%
  mutate(
    total   = sum(count),
    percent = count / total
  ) %>%
  ungroup() %>%
  filter(group %in% c("con", "pd"))  # Safety filter (see note above)

# ============================================================
# STEP 7 — Stacked Bar Chart: Subcluster Proportions
# ============================================================

# NOTE: theme_few() is applied first as the base theme; all custom
# theme() settings follow in a single consolidated block.
# Reversing this order causes theme_few() to silently overwrite all
# custom settings (legend position, background, axis ticks).
p_bar <- ggplot(bar_per, aes(x = group, y = percent, fill = seurat_clusters)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  # Optionally display raw counts inside bar segments (currently disabled):
  # geom_text(aes(label = count), position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = cluster_colors) +
  labs(
    title = "Oligodendrocyte Subcluster Proportions by Group",
    x     = "Group",
    y     = "Proportion",
    fill  = "Subcluster"
  ) +
  theme_few() +                                    # Base theme first
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
# STEP 8 — Identify Subcluster Marker Genes
# ============================================================

# PrepSCTFindMarkers is REQUIRED before FindAllMarkers when SCTransform
# was used — it re-corrects counts across samples for valid differential testing
oligo_obj <- PrepSCTFindMarkers(oligo_obj, assay = "SCT")

# Find marker genes for each subcluster.
# only.pos = FALSE: returns both up- and down-regulated markers per cluster.
# This differs from pipeline-wide convention (only.pos = TRUE) because
# subcluster analysis benefits from bidirectional contrast to characterize
# each population relative to all others.
markers <- FindAllMarkers(
  oligo_obj,
  assay           = "SCT",
  only.pos        = FALSE,   # Return both up- and down-regulated markers
  min.pct         = 0.25,    # Gene expressed in ≥ 25% of cells in the cluster
  logfc.threshold = 0.25,    # Minimum |log2FC| threshold
  verbose         = TRUE
)

# Report marker gene counts per cluster
marker_summary <- markers %>%
  group_by(cluster) %>%
  summarise(
    n_up   = sum(avg_log2FC > 0),
    n_down = sum(avg_log2FC < 0),
    .groups = "drop"
  )
message("Marker gene summary per subcluster:")
print(marker_summary)