# ============================================================
# Bulk Microarray Differential Expression Analysis
# Dataset: GSE7621 (Parkinson's Disease vs. Healthy Control)
# Platform: GPL570 (Affymetrix Human Genome U133 Plus 2.0)
# ============================================================

# ── Core Bioconductor / Analysis Packages ────────────────────
library(DESeq2)           
library(limma)            
library(ggplot2)         
library(ggrepel)          
library(pheatmap)         
library(RColorBrewer)     
library(dplyr)            
library(stringr)          
library(tidyverse)        
library(enrichplot)       
library(clusterProfiler)  
library(GEOquery)         

# ── Reserved for downstream analyses (not used in this script) ─
# library(qvalue)    
# library(cluster)   
# library(ggupset)   
# library(DOSE)     
# library(msigdbr)  
# library(fgsea)     
# library(Seurat)    
# library(tinyarray) 

# ============================================================
# STEP 1 — Download and Load GEO Dataset
# ============================================================

geo_number <- "GSE7621"
PATH       <- "D:/R/GEO/20250917-new/GSE7621"

# Download GSE7621 including platform annotation (GPL570)
geo_data  <- getGEO(geo_number, destdir = PATH, getGPL = TRUE)
geo_data0 <- geo_data[[1]]  # Extract the first (and only) ExpressionSet

# Extract raw expression matrix (probes × samples)
exp <- exprs(geo_data0)

# ============================================================
# STEP 2 — Load Platform Annotation and Map Probes to Genes
# ============================================================

# Read GPL570 annotation table; contains probe ID-to-gene symbol mappings
ids <- data.table::fread(
  "D:/R/GEO/20250917-new/GSE7621/GPL570-55999.txt",
  sep    = "\t",
  header = TRUE
)

# ============================================================
# STEP 3 — Build Annotated Expression Matrix
# ============================================================

# Convert expression matrix to data.frame and expose probe IDs as a join key
exp        <- as.data.frame(exp)
exp$ID     <- rownames(exp)

# Join probe-level expression with gene symbol annotation on probe ID.
# inner_join retains only probes present in both the annotation and expression data.
# Convert ids to data.frame first to avoid data.table single-bracket edge cases.
merged <- inner_join(
  as.data.frame(ids)[, c("ID", "Gene Symbol")],
  exp,
  by = "ID"
)

# Remove probes with missing or empty gene symbol annotations
merged <- na.omit(merged) %>%
  dplyr::filter(`Gene Symbol` != "")

# Deduplicate: when multiple probes map to the same gene symbol,
# retain only the probe with the highest average expression across all samples.
# NOTE: cur_data() is deprecated since dplyr 1.1.0; use pick() instead.
merged <- merged %>%
  group_by(`Gene Symbol`) %>%
  slice_max(
    order_by = rowMeans(pick(where(is.numeric))),
    n        = 1,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  na.omit()

# Set gene symbols as row names; drop the probe ID column
exp2 <- merged %>%
  tibble::column_to_rownames("Gene Symbol") %>%
  dplyr::select(-ID) %>%
  as.matrix()

# ============================================================
# STEP 4 — Log2 Transformation
# ============================================================

# Apply log2(x + 1) transformation to compress dynamic range and
# approximate normality. The pseudocount (+1) prevents log(0) errors.
# NOTE: This is log-scaling, NOT standardization (z-scoring).
# Verify that the raw GEO data has not already been log-transformed
# before applying this step to avoid double-transformation.
n_before <- nrow(exp2)
exp2     <- log2(exp2 + 1)
exp2     <- na.omit(exp2)
n_after  <- nrow(exp2)
if (n_before > n_after) {
  message(n_before - n_after, " gene(s) removed by na.omit() after log2 transformation.")
}

# ============================================================
# STEP 5 — Define Sample Groups and Verify Column Order
# ============================================================

# Hard-code expected sample count as a safety check
stopifnot(ncol(exp2) == 25)

# IMPORTANT: group labels are assigned by column position.
# Verify the actual sample names match the expected CON/PD order
# before constructing the design matrix to avoid silent mislabeling.
expected_con <- 9   # First 9 columns = healthy controls
expected_pd  <- 16  # Next 16 columns = Parkinson's Disease samples
message("First ", expected_con, " samples (CON): ",
        paste(colnames(exp2)[1:expected_con], collapse = ", "))
message("Next ", expected_pd, " samples (PD): ",
        paste(colnames(exp2)[(expected_con + 1):ncol(exp2)], collapse = ", "))

metadata <- data.frame(
  Sample = colnames(exp2),
  Group  = factor(
    c(rep("CON", expected_con), rep("PD", expected_pd)),
    levels = c("CON", "PD")   # CON is the reference level
  )
)

# ============================================================
# STEP 6 — Differential Expression Analysis (limma)
# ============================================================

# Build design matrix: intercept (CON baseline) + GroupPD coefficient
design <- model.matrix(~ Group, data = metadata)

# Fit linear model to each gene
fit <- lmFit(exp2, design)

# Apply empirical Bayes moderation with trend correction.
# trend = TRUE accounts for mean-variance relationship in log2 microarray data.
fit2 <- eBayes(fit, trend = TRUE)

# Extract full results table for the PD vs CON contrast
# adjust.method = "fdr": Benjamini-Hochberg FDR correction for multiple testing
DEG <- topTable(fit2, coef = "GroupPD", n = Inf, adjust.method = "fdr") %>%
  rownames_to_column("symbol")

# ============================================================
# STEP 7 — Classify Genes as Up / Down / Not Significant
# ============================================================

# NOTE: Two commonly used thresholds for significance:
#   - P.Value    < 0.05: raw p-value (more permissive, no multiple testing correction)
#   - adj.P.Val  < 0.05: FDR-adjusted p-value (recommended for multi-gene studies)
#
# This analysis uses raw P.Value for a broader exploratory gene set.
# Switch to adj.P.Val for publication-quality stringency.
DEG <- DEG %>%
  mutate(
    change = ifelse(
      P.Value < 0.05 & abs(logFC) >= 0.5,
      ifelse(logFC > 0, "Up", "Down"),
      "Not"
    )
  )

# Count and report DEG composition for quick diagnostic review
message("DEG summary — Up: ",   sum(DEG$change == "Up"),
        " | Down: ", sum(DEG$change == "Down"),
        " | Not significant: ", sum(DEG$change == "Not"))

# Retain only statistically significant differentially expressed genes
diff_gene <- DEG %>%
  dplyr::filter(P.Value < 0.05 & abs(logFC) >= 0.5)