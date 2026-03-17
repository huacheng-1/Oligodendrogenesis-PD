# ============================================================
# Target Gene Expression Heatmap
# Input: exp2 — normalized expression matrix from preprocessing
#        (GSE7621: 9 CON + 16 PD = 25 samples)
#        (GSE49036: 8 CON + 8 PD Braak stage 5-6 = 16 samples)
# ============================================================

# ── Dependencies ─────────────────────────────────────────────
library(pheatmap)  

# ============================================================
# STEP 1 — Select Active Dataset
# ============================================================

# Set this variable to switch all downstream steps between datasets.
# Valid values: "GSE7621" or "GSE49036"
# IMPORTANT: exp2 and diff_gene must already be loaded for the chosen dataset
#            before running this script.
active_dataset <- "GSE7621"  # <-- SET THIS: "GSE7621" or "GSE49036"
message("Active dataset: ", active_dataset)

# ============================================================
# STEP 2 — Define Target Gene Lists Per Dataset
# ============================================================

# Target genes for GSE7621:
# Stress response, apoptosis, transcription factors, and dopaminergic markers
target_genes_GSE7621 <- c(
  "ATF4",  "CAPN2", "BANF1", "PNPT1", "PTPRK",
  "FOXO3", "SIRPA", "STAU1", "ETV5",  "RELA",
  "KEAP1", "ATM",   "DHFR",  "ABL1"
)

# Target genes for GSE49036:
# Oxidative stress regulators and immune-related genes.
# NOTE: GSE49036 analysis is restricted to PD patients at Braak
# α-synuclein stages 5–6 and their matched healthy controls.
target_genes_GSE49036 <- c(
  "KEAP1", "ATM", "FOXO3", "RELA",
  "SETD2", "SIRPA", "BANF1"
)

# Select the gene list corresponding to the active dataset
target_genes <- switch(
  active_dataset,
  "GSE7621"  = target_genes_GSE7621,
  "GSE49036" = target_genes_GSE49036,
  stop("Unknown active_dataset value: '", active_dataset,
       "'. Must be 'GSE7621' or 'GSE49036'.")
)
message("Target genes selected: ", length(target_genes))

# ============================================================
# STEP 3 — Extract Target Gene Expression from Matrix
# ============================================================

# Subset exp2 to rows matching target gene symbols
heatmap_data <- exp2[rownames(exp2) %in% target_genes, , drop = FALSE]

# Report which target genes were found and which are absent in exp2
found_genes   <- intersect(target_genes, rownames(exp2))
missing_genes <- setdiff(target_genes, rownames(exp2))

message("Genes found in exp2:   ", length(found_genes),
        " — ", paste(found_genes, collapse = ", "))
if (length(missing_genes) > 0) {
  warning("Genes NOT found in exp2: ", length(missing_genes),
          " — ", paste(missing_genes, collapse = ", "),
          "\nCheck gene symbol case and spelling.")
}

# Stop early if no target genes were found — prevents cryptic pheatmap errors
if (nrow(heatmap_data) == 0) {
  stop("heatmap_data is empty: none of the target genes were found in exp2. ",
       "Check that exp2 row names use the same gene symbol format as target_genes.")
}

# ============================================================
# STEP 4 — Build Sample Metadata for Column Annotation
# ============================================================

# IMPORTANT: Only the active dataset's metadata frame is constructed here.
# Building both frames unconditionally would cause a length-mismatch error
# because each dataset has a different number of samples (25 vs 16),
# but both would try to use the same colnames(heatmap_data).

metadata <- switch(
  active_dataset,
  
  "GSE7621" = data.frame(
    # GSE7621: 9 healthy controls + 16 PD patients = 25 samples
    Sample = colnames(heatmap_data),
    Group  = factor(
      c(rep("CON", 9), rep("PD", 16)),
      levels = c("CON", "PD")
    ),
    stringsAsFactors = FALSE
  ),
  
  "GSE49036" = data.frame(
    # GSE49036 (Braak 5–6 subset): 8 controls + 8 PD patients = 16 samples
    Sample = colnames(heatmap_data),
    Group  = factor(
      c(rep("CON", 8), rep("PD", 8)),
      levels = c("CON", "PD")
    ),
    stringsAsFactors = FALSE
  )
)

# Verify metadata rows match heatmap columns before plotting
stopifnot(
  "Metadata row count does not match heatmap column count" =
    nrow(metadata) == ncol(heatmap_data)
)

# pheatmap requires annotation_col row names to match matrix column names
rownames(metadata) <- metadata$Sample

# ============================================================
# STEP 5 — Plot Heatmap with Group Annotation
# ============================================================

p <- pheatmap(
  heatmap_data,
  
  # Row-wise Z-score normalization: centers and scales each gene across samples
  # so that expression differences are comparable across genes
  scale        = "row",
  
  # Column clustering disabled to preserve CON / PD group order
  cluster_cols = FALSE,
  # Row clustering disabled to keep the predefined gene list order
  cluster_rows = FALSE,
  
  # 5-color diverging palette: black–navy–white–yellow–firebrick
  # 100 breaks ensures a smooth gradient (20 produces visible banding)
  color = colorRampPalette(
    c("black", "navy", "white", "yellow", "firebrick3")
  )(100),
  
  show_rownames = TRUE,
  show_colnames = TRUE,
  
  # Column annotation bar showing CON (blue) and PD (red) group membership
  annotation_col    = metadata["Group"],
  annotation_colors = list(Group = c(CON = "blue", PD = "red")),
  annotation_legend = TRUE,
  
  main = paste0("Target Gene Expression Heatmap — ", active_dataset),
  
  # Cell dimensions and label formatting
  cellwidth   = 8,
  cellheight  = 8,
  angle_col   = 45,   # Numeric angle for column labels (pheatmap docs recommend numeric)
  fontsize_row = 8,   # Row label font size (reduce if gene names overlap)
  fontsize_col = 6    # Column label font size
)