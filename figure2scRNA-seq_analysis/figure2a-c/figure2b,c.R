# ============================================================
# Single-Cell RNA-seq Downstream Visualization Pipeline
# Input: scRNA_combined — fully annotated, integrated Seurat object
# ============================================================

# ── Dependencies ─────────────────────────────────────────────
library(Seurat)           
library(harmony)          
library(tidyverse)       
library(dplyr)            
library(R.utils)          
library(tidydr)          
library(glmGamPoi)        
library(patchwork)        
library(clustree)         
library(SingleR)          
library(tidyr)           
library(networkD3)        
library(clusterProfiler)  
library(future)           
library(ggthemes)         
library(pheatmap)         

# ============================================================
# SECTION 1 — Annotated UMAP Plot
# ============================================================

# scRNA_combined: fully integrated and annotated human snRNA-seq object
# Visualize cell type clusters on the UMAP embedding
p1 <- DimPlot(
  scRNA_combined,
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
print(p1)

# ============================================================
# SECTION 2 — Cell Type Proportion Stacked Bar Chart (PD vs Control)
# ============================================================

# Extract cell-level metadata for proportion calculations
meta_data <- scRNA_combined@meta.data

# Define the cell types to include in the comparison
cell_types_to_compare <- c(
  "Oligodendrocyte", "Oligodendrocyte precursor cell",
  "Astrocyte", "Microglial cell", "Neuron", "Fibroblast"
)

# Filter metadata to retain only the target cell types
filtered_meta <- meta_data[meta_data$cells %in% cell_types_to_compare, ]
filtered_meta$group <- as.factor(filtered_meta$group)
filtered_meta$cells <- as.factor(filtered_meta$cells)

# Count cells per group–cell type combination
# .groups = "drop" prevents retained grouping from affecting downstream mutate()
cell_counts <- filtered_meta %>%
  group_by(group, cells) %>%
  summarise(count = n(), .groups = "drop")

# Compute within-group proportions for each cell type
bar_per <- cell_counts %>%
  group_by(group) %>%
  mutate(
    total   = sum(count),
    percent = count / total
  ) %>%
  ungroup()

# Retain only PD and healthy control samples for the final comparison
bar_per <- bar_per %>%
  filter(group %in% c("con", "pd"))

# Define a consistent color palette for cell types
cell_type_colors <- c(
  "Oligodendrocyte"              = "red",
  "Astrocyte"                    = "#F78000",
  "Oligodendrocyte precursor cell" = "yellow",
  "Neuron"                       = "pink",
  "Microglial cell"              = "#3176B7",
  "Fibroblast"                   = "brown"
)

# Stacked bar chart showing cell type proportions per group
# NOTE: theme_few() is applied first; custom theme() overrides follow afterward
# to prevent theme_few() from overwriting the custom settings.
p2 <- ggplot(bar_per, aes(x = group, y = percent, fill = cells)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  # Optionally display raw counts in bar segments (currently disabled):
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
    plot.title      = element_text(size = 12, hjust = 0.5),
    axis.ticks      = element_line(linetype = "blank"),
    legend.position = "top",
    panel.grid.minor = element_line(colour = NA, linetype = "blank"),
    panel.background = element_rect(fill = NA),
    plot.background  = element_rect(colour = NA)
  )
print(p2)

# ============================================================
# SECTION 3 — Marker Gene Dot Plot and Heatmaps
# ============================================================

# Curated marker genes for major cell types present in the dataset:
# OPALIN, MOG        → Oligodendrocytes
# DOCK8, P2RY12      → Microglial cells
# MAP2, DCX          → Neurons
# AQP4, SLC39A12     → Astrocytes
# VCAN, OLIG2        → Oligodendrocyte precursor cells
# PDGFRB, FOXC1      → Vascular/Fibroblast cells
marker_genes <- c(
  "OPALIN", "MOG",
  "DOCK8",  "P2RY12",
  "MAP2",   "DCX",
  "AQP4",   "SLC39A12",
  "VCAN",   "OLIG2",
  "PDGFRB", "FOXC1"
)

# Dot plot: dot size = fraction of cells expressing gene; color = scaled mean expression
top1_dotplot <- DotPlot(scRNA_combined, features = marker_genes) +
  theme(axis.text.x = element_text(size = 6, angle = 45, vjust = 1, hjust = 1))
print(top1_dotplot)

# ── Heatmap Preparation ───────────────────────────────────────

# Extract scaled expression values from dot plot data and reshape to wide format
# Rows = cell types (id), Columns = genes
dotplot_data <- top1_dotplot$data

heatmap_data <- dotplot_data %>%
  dplyr::select(features.plot, id, avg.exp.scaled) %>%
  pivot_wider(names_from = features.plot, values_from = avg.exp.scaled)

# Move cell type labels from column to row names for pheatmap compatibility
heatmap_data <- tibble::column_to_rownames(heatmap_data, var = "id")

# ── Heatmap 1: Cell types (rows) × Genes (columns) ───────────
heatmap1 <- pheatmap(
  heatmap_data,
  cellwidth     = 10,
  cellheight    = 15,
  cluster_rows  = FALSE,
  cluster_cols  = FALSE,
  legend_breaks = c(-2, -1, 0, 1, 2),
  color         = colorRampPalette(c("#d3eafa", "white", "#e85b3f"))(100)
)

# ── Heatmap 2: Genes (rows) × Cell types (columns) ───────────
# Transpose so genes become rows and cell types become columns.
# The id column was already converted to row names above,
# so no additional column removal is needed before transposing.
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