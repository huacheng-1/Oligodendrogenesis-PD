# ============================================================
# Bulk Microarray Differential Expression Analysis
# Dataset: GSE49036 (Parkinson's Disease vs. Healthy Control)
# Platform: GPL570 (Affymetrix Human Genome U133 Plus 2.0)
# Note: GSE49036 contains multiple disease subgroups.
#       This script analyzes CON (controls) vs. PD Stage 3 only.
# ============================================================

# ── Core Analysis Packages ───────────────────────────────────
library(limma)            
library(ggplot2)          
library(pheatmap)        
library(RColorBrewer)     
library(dplyr)            
library(ggrepel)          
library(stringr)          
library(enrichplot)       
library(clusterProfiler)  
library(tidyverse)        
library(GEOquery)         

# ── Reserved for downstream analyses (not used in this script) ─
# library(DESeq2)     
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

geo_number <- "GSE49036"
PATH       <- "D:/R/GEO/20250917-NEW/NEW/GSE49036/"

# Download GSE49036 including GPL570 platform annotation
geo_data  <- getGEO(geo_number, destdir = PATH, getGPL = TRUE)
geo_data0 <- geo_data[[1]]  # Extract the first ExpressionSet object

# Extract the full probe-level expression matrix (probes × all samples)
exp <- exprs(geo_data0)

# ============================================================
# STEP 2 — Select Samples of Interest
# ============================================================

# GSE49036 contains multiple groups (controls, PD stages, LBD, etc.).
# This analysis compares:
#   CON     — healthy controls:        GSM1192691–GSM1192698 (n = 8)
#   PD Gr3  — Parkinson's Disease group 3: GSM1192711–GSM1192718 (n = 8)
# Only these 16 samples are retained; all other groups are excluded.
rowcon  <- paste0("GSM", 1192691:1192698)  # 8 control samples
rowpd3  <- paste0("GSM", 1192711:1192718)  # 8 PD Stage 3 samples

# Confirm that all expected sample IDs exist in the downloaded matrix
missing_ids <- setdiff(c(rowcon, rowpd3), colnames(exp))
if (length(missing_ids) > 0) {
  stop("The following sample IDs are missing from the expression matrix: ",
       paste(missing_ids, collapse = ", "))
}

# Subset to selected samples (CON first, then PD3 — order defines group labels)
exp <- exp %>%
  as.data.frame() %>%
  dplyr::select(all_of(c(rowcon, rowpd3))) %>%
  as.matrix()

# Verify that the expression matrix has the expected value range.
# Microarray data that has been log2-transformed typically has values in ~2–15.
# Raw (non-log) data typically has values in the thousands.
exp_range <- range(exp, na.rm = TRUE)
message("Expression value range: [", round(exp_range[1], 2),
        ", ", round(exp_range[2], 2), "]")
if (exp_range[2] > 100) {
  warning("Max expression value > 100. Data may NOT be log-transformed. ",
          "Verify before proceeding — applying log2(x+1) may be necessary.")
} else {
  message("Value range consistent with log2-transformed microarray data. ",
          "No further log transformation applied.")
}

# ============================================================
# STEP 3 — Load Platform Annotation
# ============================================================

# Read GPL570 annotation file; maps Affymetrix probe IDs to gene symbols
ids <- data.table::fread(
  "D:/R/GEO/20250917-new/GSE49036/GPL570-55999.txt",
  sep    = "\t",
  header = TRUE
)

# ============================================================
# STEP 4 — Build Annotated Expression Matrix
# ============================================================

# Convert expression matrix to data.frame and expose probe IDs as a join key
exp     <- as.data.frame(exp)
exp$ID  <- rownames(exp)

# Join probe-level expression with gene symbol annotation.
# inner_join retains only probes present in both tables.
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
# with_ties = FALSE ensures exactly one row per gene even when means are tied.
merged <- merged %>%
  group_by(`Gene Symbol`) %>%
  slice_max(
    order_by  = rowMeans(pick(where(is.numeric))),
    n         = 1,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  na.omit()

# Set gene symbols as row names and drop the probe ID column
exp2 <- merged %>%
  tibble::column_to_rownames("Gene Symbol") %>%
  dplyr::select(-ID) %>%
  as.matrix()

# Report gene loss from the second na.omit() for transparency
message("Genes retained after final na.omit(): ", nrow(exp2))

# ============================================================
# STEP 5 — Define Sample Groups and Verify Column Order
# ============================================================

# Hard-code expected total sample count as a safety check
stopifnot(ncol(exp2) == 16)

# Print the actual sample names assigned to each group so the
# CON/PD3 order can be visually confirmed before model fitting.
message("CON samples (columns 1–8):  ",
        paste(colnames(exp2)[1:8],  collapse = ", "))
message("PD3 samples (columns 9–16): ",
        paste(colnames(exp2)[9:16], collapse = ", "))

metadata <- data.frame(
  Sample = colnames(exp2),
  Group  = factor(
    c(rep("CON", 8), rep("PD", 8)),
    levels = c("CON", "PD")  # CON is the reference level for the contrast
  )
)

# ============================================================
# STEP 6 — Differential Expression Analysis (limma)
# ============================================================

# Build design matrix: intercept (CON baseline) + GroupPD coefficient
design <- model.matrix(~ Group, data = metadata)

# Fit a linear model to each gene across all samples
fit  <- lmFit(exp2, design)

# Apply empirical Bayes moderation.
# trend = TRUE: accounts for the mean-variance relationship in log2 microarray data.
fit2 <- eBayes(fit, trend = TRUE)

# Extract the full results table for the PD vs CON contrast.
# adjust.method = "none": no multiple testing correction is applied here.
# RATIONALE: With only n = 8 per group, FDR correction can be overly conservative
# and suppress biologically meaningful signals in exploratory analysis.
# For publication, consider switching to adjust.method = "BH" (Benjamini-Hochberg)
# and filtering on adj.P.Val < 0.05 instead of raw P.Value.
DEG <- topTable(fit2, coef = "GroupPD", n = Inf, adjust.method = "none") %>%
  rownames_to_column("symbol")

# ============================================================
# STEP 7 — Classify Genes as Up / Down / Not Significant
# ============================================================

# Significance criteria: raw P.Value < 0.05 AND |log2FC| >= 0.5
# (equivalent to ~1.41-fold change threshold)
DEG <- DEG %>%
  mutate(
    change = ifelse(
      P.Value < 0.05 & abs(logFC) >= 0.5,
      ifelse(logFC > 0, "Up", "Down"),
      "Not"
    )
  )

# Report DEG composition as a quick diagnostic check
message("DEG summary — Up: ",   sum(DEG$change == "Up"),
        " | Down: ", sum(DEG$change == "Down"),
        " | Not significant: ", sum(DEG$change == "Not"))

# Retain only the significant DEGs for downstream analysis
diff_gene <- DEG %>%
  dplyr::filter(P.Value < 0.05 & abs(logFC) >= 0.5)