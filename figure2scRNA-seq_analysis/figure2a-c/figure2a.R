# ============================================================
# Multi-Dataset Single-Cell RNA-seq Integration Pipeline
# Datasets: GSE243639, GSE178265, GSE193688
# Disease context: Parkinson's Disease (PD), LBD, Control
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

# Increase global variable size limit for parallel workers (required for large Seurat objects)
options(future.globals.maxSize = 300000 * 1024^2)  # 300 GB

# ============================================================
# STEP 1 — Reconstruct Per-Sample Count Matrices
#          from Pre-integrated Object (com3965)
# ============================================================
# com3965 is a previously integrated Seurat object combining GSE243639 and GSE178265.
# We extract raw counts and re-split by sample to allow fresh QC and re-normalization.

filtered_dfs <- list()
se_list      <- list()

load("~/hc/RE/GSE243639-GSE178265/com3965.RData")
stopifnot(exists("com3965"))  # Abort if the file did not load correctly

# Retrieve the ordered sample IDs from metadata
level    <- levels(com3965@meta.data$orig.ident)
dir_name <- level

# Extract the raw count matrix (genes × cells) as a plain data.frame
data2 <- as.data.frame(com3965@assays$RNA$counts)

# Split the count matrix into per-sample subsets using barcode prefixes
for (i in seq_along(level)) {
  filtered_dfs[[dir_name[i]]] <- data2 %>%
    dplyr::select(starts_with(level[i]))  # Select cells belonging to sample i
  gc()
}
gc()

# ============================================================
# STEP 2 — Create Per-Sample Seurat Objects and Apply QC Filters
# ============================================================

# Build one Seurat object per sample
for (i in seq_along(dir_name)) {
  se_list[[i]] <- CreateSeuratObject(
    counts  = filtered_dfs[[dir_name[i]]],
    project = dir_name[i]   # Store sample ID in orig.ident
  )
}

# Apply quality control filters to each sample independently
for (i in seq_along(se_list)) {
  sc <- se_list[[i]]
  
  # Compute mitochondrial transcript percentage per cell
  sc[["mt_percent"]] <- PercentageFeatureSet(sc, pattern = "^MT-")
  
  # Filter cells by total UMI count (nCount_RNA already computed by Seurat)
  # Retain cells with 200–6000 UMIs to exclude empty droplets and doublets
  sc <- subset(sc, nCount_RNA >= 200 & nCount_RNA <= 6000)
  
  # Exclude cells with high mitochondrial content (likely damaged/dying)
  sc <- subset(sc, mt_percent < 10)
  
  # NOTE: Normalization (NormalizeData / ScaleData) is intentionally omitted here.
  # SCTransform will be applied after merging, making per-sample normalization redundant.
  
  se_list[[i]] <- sc
  rm(sc)
  gc()
}
names(se_list) <- dir_name

# ============================================================
# STEP 3 — Merge, Normalize, Cluster (GSE243639 + GSE178265)
# ============================================================

# Merge all per-sample objects into one; add sample ID prefix to each barcode
scRNA_combined <- merge(
  x            = se_list[[1]],
  y            = se_list[-1],
  add.cell.ids = names(se_list)
)

# Normalize using SCTransform with mitochondrial regression
# variable.features.n = 4000 captures more biology before cross-dataset integration
scRNA_combined <- SCTransform(
  scRNA_combined,
  vars.to.regress     = "mt_percent",
  variable.features.n = 4000,
  verbose             = FALSE
)
gc()

# Dimensionality reduction and clustering (preliminary; will be redone after Harmony)
scRNA_combined <- RunPCA(scRNA_combined, npcs = 30, verbose = FALSE)
gc()
scRNA_combined <- FindNeighbors(scRNA_combined, dims = 1:30)
gc()
scRNA_combined <- FindClusters(scRNA_combined, resolution = 0.8)
gc()
scRNA_combined <- RunUMAP(scRNA_combined, dims = 1:30, min.dist = 0.5, n.neighbors = 30)
gc()

# ============================================================
# STEP 4 — Load GSE193688 and Harmonize Gene Space
# ============================================================

load("~/hc/RE/GSE193688-HOMO/Before_annotation.RData")
stopifnot(exists("gse193688"))  # Abort if object is missing

# Restrict both datasets to their shared gene set to enable integration
common_genes           <- intersect(rownames(gse193688), rownames(scRNA_combined))
gse193688_filtered     <- gse193688[common_genes, ]
scRNA_combined_filtered <- scRNA_combined[common_genes, ]
gc()

# Keep only the two filtered objects to free memory before re-normalization
rm(list = ls()[!ls() %in% c("gse193688_filtered", "scRNA_combined_filtered")])
gc()

# ============================================================
# STEP 5 — Re-run SCTransform on Each Dataset Before Merging
# ============================================================
# NOTE: These two SCTransform calls prepare each dataset independently.
# A final SCTransform will be run on the merged object (Step 6) — 
# that final run is what matters for downstream analysis.
# Unifying variable.features.n and vars.to.regress ensures consistency.

scRNA_filtered_repro <- SCTransform(
  scRNA_combined_filtered,
  vars.to.regress     = "mt_percent",
  variable.features.n = 3000,
  verbose             = FALSE
)
gc()

gse193688_filtered_repro <- SCTransform(
  gse193688_filtered,
  vars.to.regress     = "mt_percent",
  variable.features.n = 3000,
  verbose             = FALSE
)
gc()

# ============================================================
# STEP 6 — Merge Datasets and Run Final SCTransform
# ============================================================

combined <- merge(
  x            = gse193688_filtered_repro,
  y            = scRNA_filtered_repro,
  add.cell.ids = c("GSE193688", "scRNA_Combined"),
  project      = "Integrated_Project"
)

# Final normalization on the merged object — this is the SCT assay used downstream
combined <- SCTransform(
  combined,
  vars.to.regress     = "mt_percent",
  variable.features.n = 3000,
  verbose             = FALSE
)
gc()

# Confirm SCTransform set the default assay correctly (should already be "SCT")
DefaultAssay(combined) <- "SCT"
message("Variable features count: ", length(VariableFeatures(combined)))  # Expected: 3000
head(VariableFeatures(combined))  # Diagnostic: confirm gene names are present

# ============================================================
# STEP 7 — PCA and Harmony Batch Correction
# ============================================================

combined <- RunPCA(combined, assay = "SCT", npcs = 30, verbose = FALSE)

# Harmony corrects for sample-level batch effects using orig.ident as the batch variable.
# theta controls diversity penalty (higher = more mixing); lambda controls correction strength.
combined <- combined %>%
  RunHarmony(
    group.by.vars    = "orig.ident",
    dims.use         = 1:30,
    theta            = 2,
    lambda           = 0.5,
    plot_convergence = TRUE   # Display convergence diagnostic plot
  )

# ============================================================
# STEP 8 — Clustering and UMAP on Harmony-Corrected Embedding
# ============================================================

combined <- FindNeighbors(combined, reduction = "harmony", dims = 1:30)
combined <- FindClusters(combined, resolution = 0.8)
combined <- RunUMAP(
  combined,
  reduction  = "harmony",
  dims       = 1:30,
  n.neighbors = 30,
  min.dist   = 0.3
)

# Rename to reflect this is now the final integrated object
harmony_integrated <- combined
rm(combined)
gc()

# Visualize clusters before annotation
p_unannotated <- DimPlot(
  harmony_integrated,
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
# STEP 9 — Identify Cluster Marker Genes
# ============================================================

# Prepare SCT model parameters across all samples before differential testing
harmony_integrated <- PrepSCTFindMarkers(harmony_integrated)

# Find positive marker genes per cluster
# min.pct = 0.25: gene expressed in ≥ 25% of cluster cells
# logfc.threshold = 0.25: minimum log2 fold-change filter
markers <- FindAllMarkers(
  harmony_integrated,
  only.pos       = TRUE,
  min.pct        = 0.25,
  logfc.threshold = 0.25,
  verbose        = FALSE
)

# Extract top 300 markers per cluster by average log2 fold-change for manual annotation.
# NOTE: top_n() is deprecated since dplyr 1.0.0; use slice_max() instead.
top300 <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 300)

# Print marker genes per cluster for manual review on the ACT annotation portal
tmp <- split(top300$gene, top300$cluster)
cat(unlist(lapply(names(tmp), function(x) {
  paste0("cluster", x, ": ", paste(tmp[[x]], collapse = ", "))
})), sep = "\n")

# ============================================================
# STEP 10 — Manual Cell Type Annotation
# ============================================================

# Sanity check: confirm expected cluster count before assigning labels
message("Total clusters detected: ", length(levels(harmony_integrated)))

# Cell type labels assigned in cluster order (0-indexed, left to right)
new.cluster.ids <- c(
  "Oligodendrocyte",               # Cluster 0
  "Oligodendrocyte",               # Cluster 1
  "Microglial cell",               # Cluster 2
  "Oligodendrocyte",               # Cluster 3
  "Neuron",                        # Cluster 4
  "Neuron",                        # Cluster 5
  "Astrocyte",                     # Cluster 6
  "Oligodendrocyte",               # Cluster 7
  "Oligodendrocyte",               # Cluster 8
  "Astrocyte",                     # Cluster 9
  "Oligodendrocyte precursor cell",# Cluster 10
  "Oligodendrocyte",               # Cluster 11
  "Fibroblast",                    # Cluster 12
  "Oligodendrocyte",               # Cluster 13
  "Microglial cell",               # Cluster 14
  "Neuron",                        # Cluster 15
  "Neuron",                        # Cluster 16
  "Microglial cell",               # Cluster 17
  "Oligodendrocyte",               # Cluster 18
  "Oligodendrocyte",               # Cluster 19
  "Oligodendrocyte",               # Cluster 20
  "Neuron",                        # Cluster 21
  "Oligodendrocyte",               # Cluster 22
  "Astrocyte",                     # Cluster 23
  "Neuron",                        # Cluster 24
  "Neuron",                        # Cluster 25
  "Neuron",                        # Cluster 26
  "Astrocyte",                     # Cluster 27
  "Oligodendrocyte",               # Cluster 28
  "Neuron",                        # Cluster 29
  "Oligodendrocyte",               # Cluster 30
  "Neuron",                        # Cluster 31
  "Neuron",                        # Cluster 32
  "Neuron",                        # Cluster 33
  "Neuron",                        # Cluster 34
  "Neuron",                        # Cluster 35
  "Oligodendrocyte",               # Cluster 36
  "Neuron",                        # Cluster 37
  "Fibroblast",                    # Cluster 38
  "Neuron",                        # Cluster 39
  "Neuron",                        # Cluster 40
  "Astrocyte",                     # Cluster 41
  "Oligodendrocyte",               # Cluster 42
  "Neuron",                        # Cluster 43
  "Astrocyte",                     # Cluster 44
  "Oligodendrocyte precursor cell",# Cluster 45
  "Neuron",                        # Cluster 46
  "Oligodendrocyte",               # Cluster 47
  "Neuron",                        # Cluster 48
  "Astrocyte",                     # Cluster 49
  "Neuron",                        # Cluster 50
  "Neuron",                        # Cluster 51
  "Oligodendrocyte precursor cell",# Cluster 52
  "Oligodendrocyte",               # Cluster 53
  "Oligodendrocyte",               # Cluster 54
  "Neuron",                        # Cluster 55
  "Oligodendrocyte",               # Cluster 56
  "Oligodendrocyte",               # Cluster 57
  "Astrocyte",                     # Cluster 58
  "Oligodendrocyte precursor cell",# Cluster 59
  "Oligodendrocyte precursor cell",# Cluster 60
  "Oligodendrocyte"                # Cluster 61
)

# Verify label count matches cluster count before renaming (prevents silent misalignment)
stopifnot(length(new.cluster.ids) == length(levels(harmony_integrated)))

names(new.cluster.ids) <- levels(harmony_integrated)
harmony_integrated     <- RenameIdents(harmony_integrated, new.cluster.ids)

# Store cell type annotation in metadata for downstream use (e.g., subsetting, plotting)
harmony_integrated@meta.data$cells <- harmony_integrated@active.ident

# Visualize clusters after annotation
p_annotated <- DimPlot(
  harmony_integrated,
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
# STEP 11 — Assign Disease Group Labels to Cells
# ============================================================

# Sample ID vectors for each disease/phenotype group
# PD  = Parkinson's Disease samples
# CON = Healthy controls
# LBD = Lewy Body Dementia samples
# OLD = Aged controls (from GSE193688)

pd_samples <- c(
  "s.0096","s.0097","s.0098","s.0099","s.0100","s.0102","s.0103",
  "s.0104","s.0105","s.0107","s.0109","s.0110","s.0111","s.0116","s.0118","s.0119","s.0120",
  "pPDsHSrSNxi3873d200429PosA","pPDsHSrSNxi3873d200429PosE","pPDsHSrSNxi3873d200429PosF",
  "pPDsHSrSNxi3887d200429PosB","pPDsHSrSNxi3873d200429PosD","pPDsHSrSNxi4560d200429PosD",
  "pPDsHSrSNxi3873d200429DAPIA","pPDsHSrSNxi3887d200429DAPIB","pPDsHSrSNxi3873d200429PosB",
  "pPDsHSrSNxi3873d200429DAPIB","pPDsHSrSNxi4560d200429PosE","pPDsHSrSNxi3887d200429PosA",
  "pPDsHSrSNxi4560d200429DAPIB","pPDsHSrSNxi3873d200429PosC","pPDsHSrSNxi4560d200429DAPIA",
  "pPDsHSrSNxi4560d200429PosB","pPDsHSrSNxi4568d200429DAPIB","pPDsHSrSNxi2142d200429PosA",
  "pPDsHSrSNxi4560d200429PosC","pPDsHSrSNxi4568d200429PosB","pPDsHSrSNxi4560d200429PosA",
  "pPDsHSrSNxi1963d200429DAPIB","pPDsHSrSNxi3887d200429DAPIA","pPDsHSrSNxi4568d200429PosA",
  "pPDsHSrSNxi4568d200429DAPIA","pPDsHSrSNxi1963d200429DAPIC","pPDsHSrSNxi1963d200429DAPIA",
  "pPDsHSrSNxi2142d200429DAPIA","pPDsHSrSNxi2142d200429DAPIB","pPDsHSrSNxi1963d200429PosB",
  "pPDsHSrSNxi1963d200429PosA","pPDsHSrSNxi1963d200429PosC","pPDsHSrSNxi1963d200429PosD",
  "pPDsHSrSNxi1963d200429PosE",
  "GSM5818658_P02284_filtered_feature_bc_matrix.h5","GSM5818659_P4205_filtered_feature_bc_matrix.h5",
  "GSM5818660_P4560_filtered_feature_bc_matrix.h5","GSM5818661_P4653_filtered_feature_bc_matrix.h5",
  "GSM5818662_P4772_filtered_feature_bc_matrix.h5","GSM5818663_P4884_filtered_feature_bc_matrix.h5",
  "GSM5818664_P4919_filtered_feature_bc_matrix.h5","GSM5818665_P5318_filtered_feature_bc_matrix.h5",
  "GSM5818666_P5331_filtered_feature_bc_matrix.h5","GSM5818667_P5626_filtered_feature_bc_matrix.h5",
  "GSM5818668_P5662_filtered_feature_bc_matrix.h5","GSM5818669_P6259_filtered_feature_bc_matrix.h5",
  "GSM5818670_P6320_filtered_feature_bc_matrix.h5","GSM5818671_P6323_filtered_feature_bc_matrix.h5",
  "GSM5818672_P6326_filtered_feature_bc_matrix.h5"
)

con_samples <- c(
  "s.0127","s.0128","s.0129","s.0130","s.0131","s.0142","s.0147",
  "s.0151","s.0152","s.0153","s.0154","s.0158","s.0159","s.0165",
  "pPDCN4340DAPIA030419","pPDCN4340DAPIB030419",
  "pPDsHSrSNxi3482d200429PosB","pPDsHSrSNxi3346d200429PosB","pPDsHSrSNxi3482d200429PosA",
  "pPDsHSrSNxi3346d200429PosC","pPDsHSrSNxi5610d200429Pos","pPDsHSrSNxi3322d200429PosG",
  "pPDsHSrSNxi3322d200429PosE","pPDsHSrSNxi3346d200429PosA","pPDCN3839DAPIB030419",
  "pPDsHSrSNxi3345d200429DAPIB","pPDsHSrSNxi3482d200429DAPIB","pPDsHSrSNxi4956d200429PosA",
  "pPDsHSrSNxi3482d200429DAPIA","pPDsHSrSNxi3322d200429PosF","pPDsHSrSNxi6173d200429DAPIA",
  "pPDsHSrSNxi3345d200429DAPIA","pPDsHSrSNxi3346d200429DAPIB","pPDsHSrSNxi4956d200429DAPIB",
  "pPDsHSrSNxi3322d200429PosC","pPDsHSrSNxi3322d200429DAPIA","pPDsHSrSNxi4956d200429DAPIA",
  "pPDsHSrSNxi3345d200429DAPIC","pPDsHSrSNxi3322d200429PosD","pPDCN3839DAPIA030419",
  "pPDsHSrSNxi3345d200429PosA","pPDsHSrSNxi5610d200429Neg","pPDsHSrSNxi3322d200429DAPIB",
  "pPDCN5730NeuN22119","pPDsHSrSNxi4956d200429PosB","pPDsHSrSNxi6173d200429PosB",
  "pPDsHSrSNxi3322d200429PosA","pPDsHSrSNxi3298d200429DAPIC","pPDsHSrSNxi3322d200429PosB",
  "pPDsHSrSNxi3322d200429PosH","pPDsHSrSNxi3322d200429DAPIC","pPDsHSrSNxi6173d200429PosA",
  "pPDsHSrSNxi3298d200429DAPIA","pPDCN5730DAPI22119","pPDsHSrSNxi3298d200429DAPIB",
  "pPDCN3898DAPIB030419","pPDCN3898DAPIA030419","pPDsHSrSNxi3298d200429PosA",
  "pPDsHSrSNxi3298d200429PosB",
  "GSM5818673_YC4345_filtered_feature_bc_matrix.h5","GSM5818674_YC5276_filtered_feature_bc_matrix.h5",
  "GSM5818675_YC5346_filtered_feature_bc_matrix.h5","GSM5818676_YC5442_filtered_feature_bc_matrix.h5",
  "GSM5818677_YC5626_filtered_feature_bc_matrix.h5","GSM5818678_YC5751_filtered_feature_bc_matrix.h5",
  "GSM5818679_YC5813_filtered_feature_bc_matrix.h5","GSM5818680_YC5928_filtered_feature_bc_matrix.h5",
  "GSM5818681_YC6235_filtered_feature_bc_matrix.h5","GSM5818682_YC914_filtered_feature_bc_matrix.h5"
)

lbd_samples <- c(
  "pPDsHSrSNxi4775d200429PosA","pPDsHSrSNxi2544d200429DAPIB",
  "pPDsHSrSNxi4775d200429DAPIB","pPDsHSrSNxi4775d200429PosB","pPDsHSrSNxi4775d200429DAPIC",
  "pPDsHSrSNxi4775d200429DAPIA","pPDsHSrSNxi2561d200429PosC","pPDsHSrSNxi2544d200429PosB",
  "pPDsHSrSNxi2569d200429PosA","pPDsHSrSNxi2544d200429DAPIA","pPDsHSrSNxi2569d200429DAPIA",
  "pPDsHSrSNxi2561d200429PosA","pPDsHSrSNxi2561d200429DAPIB","pPDsHSrSNxi2544d200429PosA",
  "pPDsHSrSNxi2544d200429PosC","pPDsHSrSNxi2561d200429PosB","pPDsHSrSNxi2561d200429DAPIA",
  "pPDsHSrSNxi2569d200429DAPIB"
)

old_samples <- c(
  "GSM5818649_C16007_filtered_feature_bc_matrix.h5",
  "GSM5818650_C4176_filtered_feature_bc_matrix.h5",
  "GSM5818651_C4322_filtered_feature_bc_matrix.h5",
  "GSM5818652_C4361_filtered_feature_bc_matrix.h5",
  "GSM5818653_C4634_filtered_feature_bc_matrix.h5",
  "GSM5818654_C4660_filtered_feature_bc_matrix.h5",
  "GSM5818655_C4691_filtered_feature_bc_matrix.h5",
  "GSM5818656_C4814_filtered_feature_bc_matrix.h5",
  "GSM5818657_C5410_filtered_feature_bc_matrix.h5"
)

# Initialise group column; cells not matching any list will remain NA
harmony_integrated@meta.data$group <- NA

harmony_integrated@meta.data$group[
  harmony_integrated@meta.data$orig.ident %in% pd_samples]  <- "pd"
harmony_integrated@meta.data$group[
  harmony_integrated@meta.data$orig.ident %in% con_samples] <- "con"
harmony_integrated@meta.data$group[
  harmony_integrated@meta.data$orig.ident %in% lbd_samples] <- "lbd"
harmony_integrated@meta.data$group[
  harmony_integrated@meta.data$orig.ident %in% old_samples] <- "old"

# Validate that all cells received a group label
n_unassigned <- sum(is.na(harmony_integrated@meta.data$group))
if (n_unassigned > 0) {
  warning(n_unassigned, " cell(s) have no group label.",
          " Check that all orig.ident values are covered by the sample lists.")
} else {
  message("All cells successfully assigned to a group.")
}