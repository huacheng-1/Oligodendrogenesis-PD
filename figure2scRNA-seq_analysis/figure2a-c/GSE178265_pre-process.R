# ==============================================================================
# Single-Sample scRNA-seq Processing Pipeline — GSE178265 (Human Brain)
#
# Workflow:
#   1. Load raw 10x count matrix
#   2. Quality-control filtering (UMI count, mitochondrial content)
#   3. Normalisation, HVG selection, PCA, clustering, UMAP
#   4. Marker identification and cluster annotation
#   5. Map per-sample disease group labels (Ctrl / PD ) to cell metadata
# ==============================================================================

# ── Core single-cell analysis framework ──────────────────────────────────────
library(Seurat)
library(harmony)        
library(glmGamPoi)      
library(clustree)       
library(celldex)        
library(SingleR)        

# ── Data wrangling & visualisation ───────────────────────────────────────────
library(tidyverse)      
library(ggplot2)        
library(cowplot)        
library(dplyr)         
library(patchwork)      
library(reshape2)       
library(ggthemes)       
library(tibble)         
library(pheatmap)       
library(networkD3)      

# ── Pathway & gene-set analysis ───────────────────────────────────────────────
library(clusterProfiler) 
library(org.Mm.eg.db)  
library(org.Hs.eg.db)   

# ── Utility ───────────────────────────────────────────────────────────────────
library(R.utils)       
library(tidydr)         
library(future)       

# Allow up to 150 GB of global data to be exported to parallel workers
options(future.globals.maxSize = 150 * 1024^3)


# ==============================================================================
# Section 1: Load raw 10x Genomics count matrix
# ==============================================================================
counts <- Read10X(data.dir = "D:/R/GSE/GSE178265-HOMO/GSE178265")

# Create Seurat object; retain only:
#   - genes detected in ≥ 3 cells         (min.cells   = 3)
#   - cells expressing ≥ 200 distinct genes (min.features = 200)
scRNAlist <- CreateSeuratObject(
  counts       = counts,
  project      = "GSE178265",
  min.cells    = 3,
  min.features = 200
)


# ==============================================================================
# Section 2: Quality-control filtering
# ==============================================================================

# Calculate the percentage of counts from mitochondrial genes (prefix "MT-").
# High mt_percent indicates damaged / dying cells that have lost cytoplasmic RNA.
scRNAlist[["mt_percent"]] <- PercentageFeatureSet(scRNAlist, pattern = "^MT-")

# Filter cells by total UMI count to remove:
#   - Low-count cells  (< 500 UMIs) — likely empty droplets or low-quality cells
#   - High-count cells (> 6000 UMIs) — likely doublets (two cells in one droplet)
umi_counts    <- Matrix::colSums(scRNAlist@assays$RNA$counts)
low_threshold  <- 500
high_threshold <- 6000
keep_cells     <- umi_counts >= low_threshold & umi_counts <= high_threshold
scRNAlist      <- subset(scRNAlist, cells = colnames(scRNAlist)[keep_cells])

scRNAlist <- subset(scRNAlist, mt_percent < 2)


# ==============================================================================
# Section 3: Normalisation, variable-feature selection, and scaling
# ==============================================================================

scRNAlist <- NormalizeData(scRNAlist)

# Identify the top 6000 highly variable genes using the VST method.
# NOTE: The Seurat default is 2000. 6000 is intentionally increased here to
# maintain feature consistency with the companion dataset GSE243639.
scRNAlist <- FindVariableFeatures(
  scRNAlist,
  selection.method = "vst",
  nfeatures        = 6000
)

scRNAlist <- ScaleData(scRNAlist, vars.to.regress = "mt_percent")


# ==============================================================================
# Section 4: Dimensionality reduction and clustering
# ==============================================================================

# Run PCA; restrict to 10 PCs to match the GSE243639 integration settings
# and reduce computation time.
scRNAlist <- RunPCA(scRNAlist, npcs = 10)

# Build a shared-nearest-neighbour graph on the first 10 PCs, then partition
# it into clusters at resolution 0.8 (higher = more, smaller clusters).
scRNAlist <- FindNeighbors(scRNAlist, dims = 1:10)
scRNAlist <- FindClusters(scRNAlist,  resolution = 0.8)

# Embed cells in 2-D UMAP space for visualisation
scRNAlist <- RunUMAP(scRNAlist, dims = 1:10)
DimPlot(scRNAlist, reduction = "umap")


# ==============================================================================
# Section 5: Marker gene identification
# ==============================================================================

DefaultAssay(scRNAlist) <- "RNA"

# Find positive marker genes for every cluster simultaneously.
# Genes must be detected in ≥ 25 % of cells and show log2FC > 0.25.
markers_all <- FindAllMarkers(
  scRNAlist,
  only.pos        = TRUE,
  min.pct         = 0.25,
  logfc.threshold = 0.25,
  verbose         = TRUE
)

# Extract the top 200 markers per cluster (by log2FC) and print a
# compact cluster → gene list for manual inspection / annotation
top200 <- markers_all %>%
  group_by(cluster) %>%
  top_n(n = 200, wt = avg_log2FC)

tmp <- split(top200$gene, top200$cluster)
cat(
  unlist(lapply(names(tmp), function(x) {
    paste0("cluster", x, ": ", paste(tmp[[x]], collapse = ", "))
  })),
  sep = "\n"
)


# ==============================================================================
# Section 6: Assign cell-type labels to clusters
# ==============================================================================

# FIX: Guard that the label vector length matches the number of clusters
#      BEFORE calling RenameIdents() to get an informative error early.
n_clusters <- length(levels(scRNAlist))
message("Total clusters to rename: ", n_clusters)

new.cluster.ids <- c(
  "Oligodendrocyte", "Oligodendrocyte", "Astrocyte",       "Microglial cell", "Oligodendrocyte",
  "Oligodendrocyte", "Neuron",          "Oligodendrocyte precursor cell", "Neuron", "Oligodendrocyte",
  "Oligodendrocyte", "Neuron",          "Neuron",           "Epithelial cell",
  "Oligodendrocyte", "Neuron",          "Neuron",           "Oligodendrocyte precursor cell", "Astrocyte",
  "Fibroblast",      "Neuron",          "Macrophage",       "Epithelial cell", "Neuron", "Microglial cell",
  "Neuron",          "Oligodendrocyte precursor cell", "Neuron", "Epithelial cell", "Epithelial cell", "Neuron",
  "Neuron",          "Neuron",          "Astrocyte"
)

stopifnot(
  "new.cluster.ids length must equal the number of Seurat clusters" =
    length(new.cluster.ids) == n_clusters
)

names(new.cluster.ids) <- levels(scRNAlist)
scRNAlist <- RenameIdents(scRNAlist, new.cluster.ids)

# FIX: Renamed "cells" → "cell_type" to avoid confusion with the Seurat
#      reserved parameter name "cells" used in subset(), group.by =, etc.
scRNAlist@meta.data$cell_type <- scRNAlist@active.ident


# ==============================================================================
# Section 7: Map disease-group labels (Ctrl / PD / LBD) to each sample
# ==============================================================================

# One label per unique orig.ident value; order must exactly match
# levels(scRNAlist$orig.ident).
group <- c(
  "Ctrl","Ctrl","Ctrl","Ctrl","Ctrl","Ctrl","PD",
  "Ctrl","PD",  "Ctrl","PD",  "Ctrl","Ctrl","Ctrl",
  "Ctrl","Ctrl","Ctrl","Ctrl","Ctrl","Ctrl","Ctrl",
  "Ctrl","LBD", "LBD", "PD",  "PD",  "PD",  "Ctrl","Ctrl",
  "PD",  "Ctrl","Ctrl","PD",  "PD",  "PD",  "PD",  "PD",
  "PD",  "LBD", "PD",  "LBD","Ctrl","LBD", "PD",  "Ctrl",
  "PD",  "LBD", "LBD","Ctrl","PD",  "Ctrl","Ctrl","PD",
  "PD",  "Ctrl","Ctrl","PD",  "PD",  "PD",  "LBD","PD",
  "PD",  "PD",  "LBD","PD",  "LBD","PD",  "Ctrl","Ctrl",
  "PD",  "Ctrl","LBD","LBD","Ctrl","LBD","Ctrl","Ctrl",
  "LBD","LBD","Ctrl","LBD","Ctrl","Ctrl","Ctrl","PD",
  "Ctrl","PD",  "LBD","PD",  "LBD","PD",
  "Ctrl","PD",  "Ctrl","Ctrl","PD",  "Ctrl"
)

stopifnot(
  "group vector length must equal the number of unique orig.ident levels" =
    length(group) == length(levels(scRNAlist$orig.ident))
)

# Build a sample → group mapping table
mapping_df <- data.frame(
  orig.ident = levels(scRNAlist$orig.ident),
  group      = group,
  stringsAsFactors = FALSE
)

# Store cell barcodes before merging so we can restore row order afterwards
scRNAlist@meta.data$cell_barcode <- rownames(scRNAlist@meta.data)

# Left-join the group label onto the per-cell metadata.
# all.x = TRUE keeps every cell even if its orig.ident has no match in mapping_df.
updated_meta <- merge(
  x    = scRNAlist@meta.data,
  y    = mapping_df,
  by   = "orig.ident",
  all.x = TRUE
)

# FIX: merge() sorts rows by the join key and does NOT preserve the original
#      cell order. Reindex to the original barcode order before writing back,
#      otherwise every metadata column is misaligned with the wrong cell.
rownames(updated_meta) <- updated_meta$cell_barcode
updated_meta           <- updated_meta[rownames(scRNAlist@meta.data), ]

scRNAlist@meta.data <- updated_meta