# ==============================================================================
# Multi-Sample scRNA-seq Integration Pipeline — GSE193688 (Human Brain)
#
# Workflow:
#   1. Load per-sample 10x H5 files and apply per-sample QC filtering
#   2. Merge all samples, normalise with SCTransform, and reduce with PCA
#   3. Correct batch effects with Harmony, cluster, and visualise with UMAP
#   4. Identify marker genes and annotate clusters with cell-type labels
#   5. Attach disease-group labels (aged / pd / young) to cell metadata
# ==============================================================================

# ── Core single-cell analysis ─────────────────────────────────────────────────
library(Seurat)
library(harmony)        
library(glmGamPoi)      
library(clustree)       
library(celldex)        
library(SingleR)        

# ── Data wrangling & visualisation ────────────────────────────────────────────
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
library(hdf5r)          
library(future)

# Allow up to 50 GB of data to be exported to parallel workers
options(future.globals.maxSize = 50000 * 1024^2)

# FIX: Scope the parallel plan only to the per-sample QC loop.
#      Leaving multisession active globally can cause unexpected behaviour
#      and memory exhaustion in SCTransform / FindAllMarkers.
plan(multisession, workers = 4)


# ==============================================================================
# Section 1: Load per-sample 10x H5 files
# ==============================================================================
data_dir <- "D:/R/GSE/GSE193688-HOMO/GSE193688_RAW/"

# FIX: use seq_along() instead of 1:length() to avoid producing c(1,0)
#      when the directory is empty, which would iterate in the wrong direction.
samples <- list.files(data_dir)
message("Samples found: ", length(samples))

scRNAlist <- vector("list", length(samples))

for (i in seq_along(samples)) {
  h5_path <- file.path(data_dir, samples[i])
  counts  <- Read10X_h5(filename = h5_path)
  
  # FIX: Read10X_h5() returns a list only for multi-modal data (e.g. CITE-seq).
  #      Detect whether the "Gene Expression" modality exists before indexing;
  #      fall back to using the matrix directly for single-modality files.
  if (is.list(counts) && "Gene Expression" %in% names(counts)) {
    rna_counts <- counts[["Gene Expression"]]
  } else {
    rna_counts <- counts
  }
  
  scRNAlist[[i]] <- CreateSeuratObject(
    counts       = rna_counts,
    project      = samples[i],
    min.cells    = 3,     # Retain genes detected in ≥ 3 cells
    min.features = 200    # Retain cells expressing ≥ 200 distinct genes
  )
}


# ==============================================================================
# Section 2: Per-sample quality-control filtering
# ==============================================================================

for (i in seq_along(scRNAlist)) {
  sc <- scRNAlist[[i]]
  
  # Flag mitochondrial reads — high mt_percent indicates damaged / dying cells
  sc[["mt_percent"]] <- PercentageFeatureSet(sc, pattern = "^MT-")
  
  # Remove likely empty droplets (< 500 UMIs) and doublets (> 6000 UMIs)
  umi_counts     <- Matrix::colSums(sc@assays$RNA$counts)
  keep_cells     <- umi_counts >= 500 & umi_counts <= 6000
  sc             <- subset(sc, cells = colnames(sc)[keep_cells])
  
  # NOTE: A 2 % mt_percent cutoff is stricter than the typical 5–20 % range
  # for brain tissue. Validate against VlnPlot() QC distributions before use.
  sc <- subset(sc, mt_percent < 2)
  
  # NOTE: nfeatures = 6000 intentionally exceeds the Seurat default (2000)
  # to maintain HVG consistency with the companion dataset GSE243639.
  sc <- NormalizeData(sc)
  sc <- FindVariableFeatures(sc, nfeatures = 6000)
  
  # Regress out mitochondrial content to remove QC-driven variation.
  # Add "S.Score" / "G2M.Score" here if cell-cycle effects are present.
  sc <- ScaleData(sc, vars.to.regress = "mt_percent")
  
  # FIX: RunPCA() here is intentionally removed. After merge() + SCTransform()
  #      a new SCT assay is created and RunPCA() is re-run on the merged object.
  #      Computing per-sample PCA at this stage wastes time without contributing
  #      to the final integrated embeddings.
  
  scRNAlist[[i]] <- sc
  rm(sc)
}

# Reset to sequential execution after the parallelised QC loop
plan(sequential)


# ==============================================================================
# Section 3: Merge samples and integrate with SCTransform + Harmony
# ==============================================================================

# Concatenate all per-sample Seurat objects; prefix cell barcodes with the
# sample name to prevent barcode collisions across samples.
scRNA_combined <- merge(
  x           = scRNAlist[[1]],
  y           = scRNAlist[-1],
  add.cell.ids = samples
)

# SCTransform: regularised negative binomial regression that simultaneously
# normalises, scales, and selects variable features across the merged object.
# This replaces the NormalizeData → ScaleData → FindVariableFeatures sequence
# for the integrated analysis.
scRNA_combined <- SCTransform(scRNA_combined, assay = "RNA", verbose = FALSE)

# FIX: RunPCA() MUST be called here on the SCT assay before RunHarmony().
#      RunHarmony() corrects existing PCA embeddings; it does NOT create them.
#      Omitting this step causes RunHarmony() to crash or use stale embeddings.
scRNA_combined <- RunPCA(scRNA_combined, assay = "SCT", npcs = 30)

# Harmony removes batch effects by iteratively adjusting PCA embeddings so
# that cells cluster by biology rather than by sample of origin.
scRNA_combined <- RunHarmony(scRNA_combined, group.by.vars = "orig.ident")

# Build an SNN graph on the first 10 Harmony dimensions, then cluster.
# NOTE: 10 dims chosen to maintain consistency with GSE243639.
scRNA_combined <- FindNeighbors(scRNA_combined, reduction = "harmony", dims = 1:10)
scRNA_combined <- FindClusters(scRNA_combined,  resolution = 0.8)

# Project cells into 2-D UMAP space using the batch-corrected Harmony embedding
scRNA_combined <- RunUMAP(scRNA_combined, reduction = "harmony", dims = 1:10)
DimPlot(scRNA_combined, reduction = "umap")


# ==============================================================================
# Section 4: Marker gene identification
# ==============================================================================

# PrepSCTFindMarkers() re-corrects the SCT model for each cluster comparison,
# which is required before calling FindAllMarkers() on an SCT-normalised object.
scRNA_combined <- PrepSCTFindMarkers(scRNA_combined)

# Find positive marker genes for every cluster (genes expressed in ≥ 25 % of
# cells in the cluster, log2FC > 0.25 vs all other clusters).
markers_all <- FindAllMarkers(
  scRNA_combined,
  only.pos        = TRUE,
  min.pct         = 0.25,
  logfc.threshold = 0.25,
  verbose         = TRUE
)

# Print a compact cluster → top-400-marker-gene summary for manual annotation
top400 <- markers_all %>%
  group_by(cluster) %>%
  top_n(n = 400, wt = avg_log2FC)

tmp <- split(top400$gene, top400$cluster)
cat(
  unlist(lapply(names(tmp), function(x) {
    paste0("cluster", x, ": ", paste(tmp[[x]], collapse = ", "))
  })),
  sep = "\n"
)


# ==============================================================================
# Section 5: Assign cell-type labels to clusters
# ==============================================================================

n_clusters <- length(levels(scRNA_combined))
message("Total clusters to rename: ", n_clusters)

new.cluster.ids <- c(
  "Oligodendrocyte", "Oligodendrocyte", "Oligodendrocyte",
  "Oligodendrocyte", "Oligodendrocyte", "Oligodendrocyte", "Oligodendrocyte precursor cell",
  "Oligodendrocyte", "Oligodendrocyte", "Astrocyte",       "Oligodendrocyte",
  "Oligodendrocyte", "Oligodendrocyte", "Oligodendrocyte",
  "Microglial cell", "Microglial cell", "Astrocyte",       "Oligodendrocyte",
  "Microglial cell", "Oligodendrocyte", "Oligodendrocyte precursor cell",
  "Neuron",          "Oligodendrocyte precursor cell", "Oligodendrocyte",
  "Oligodendrocyte", "Oligodendrocyte", "Oligodendrocyte", "Oligodendrocyte",
  "Neuron",          "Oligodendrocyte precursor cell",
  "Oligodendrocyte", "Pericyte"
)

# FIX: validate label count BEFORE RenameIdents() to get an informative error
stopifnot(
  "new.cluster.ids length must equal the number of Seurat clusters" =
    length(new.cluster.ids) == n_clusters
)

names(new.cluster.ids) <- levels(scRNA_combined)
scRNA_combined <- RenameIdents(scRNA_combined, new.cluster.ids)

# FIX: renamed "cells" → "cell_type" to avoid collision with the Seurat
#      reserved keyword "cells" used internally by subset(), group.by =, etc.
scRNA_combined@meta.data$cell_type <- scRNA_combined@active.ident


# ==============================================================================
# Section 6: Map disease-group labels to each sample
# ==============================================================================

# One label per unique orig.ident level; order must exactly match
# levels(scRNA_combined$orig.ident).
group <- c(
  "aged","aged","aged","aged","aged","aged","aged","aged","aged",
  "pd",  "pd",  "pd",  "pd",  "pd",  "pd",  "pd",  "pd",  "pd",  "pd",
  "pd",  "pd",  "pd",  "pd",  "pd",
  "young","young","young","young","young","young","young","young","young","young"
)

stopifnot(
  "group vector length must equal the number of unique orig.ident levels" =
    length(group) == length(levels(scRNA_combined$orig.ident))
)

# Build a sample-level orig.ident → group mapping table
mapping_df <- data.frame(
  orig.ident       = levels(scRNA_combined@meta.data$orig.ident),
  group            = group,
  stringsAsFactors = FALSE
)

# Store cell barcodes before merging so we can restore row order afterwards
scRNA_combined@meta.data$cell_barcode <- rownames(scRNA_combined@meta.data)

# Left-join the group label onto the per-cell metadata.
# all.x = TRUE retains every cell even if its orig.ident has no mapping entry.
updated_meta <- merge(
  x     = scRNA_combined@meta.data,
  y     = mapping_df,
  by    = "orig.ident",
  all.x = TRUE
)

# FIX: merge() sorts rows by the join key and does NOT preserve the original
#      cell order. Without re-sorting, every metadata column is misaligned
#      with the wrong cell, silently corrupting the entire Seurat object.
rownames(updated_meta) <- updated_meta$cell_barcode
updated_meta           <- updated_meta[rownames(scRNA_combined@meta.data), ]

scRNA_combined@meta.data <- updated_meta