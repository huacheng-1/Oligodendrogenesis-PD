# ============================================================
# Volcano Plot Visualization for Differential Expression Results
# Input: DEG — output from limma topTable() (GSE7621 or GSE49036)
# Highlights biologically relevant target genes with gold points
# ============================================================

# ── Required Packages ────────────────────────────────────────
library(dplyr)      
library(ggplot2)    
library(ggrepel)    
library(tidyverse)  

# ── Reserved for other scripts (not used here) ───────────────
# library(DESeq2)        
# library(qvalue)       
# library(cluster)      
# library(pheatmap)      
# library(RColorBrewer)  
# library(stringr)       
# library(ggupset)       
# library(enrichplot)   
# library(DOSE)          
# library(clusterProfiler) 
# library(msigdbr)       
# library(fgsea)         
# library(Seurat)        
# library(GEOquery)      
# library(tinyarray)     
# library(limma)         

# ============================================================
# STEP 1 — Define Significance Thresholds as Named Constants
# ============================================================

# Using named constants avoids magic numbers scattered across the code.
# P_THRESH = -log10(0.05) ≈ 1.301: significance cutoff on the -log10 scale
# FC_THRESH = 0.5: minimum absolute log2 fold-change
P_THRESH <- -log10(0.05)  # ≈ 1.301
FC_THRESH <- 0.5

# ============================================================
# STEP 2 — Validate DEG Input and Build Volcano Data
# ============================================================

# Confirm that required columns are present in the DEG object
stopifnot(
  all(c("symbol", "P.Value", "logFC") %in% colnames(DEG)),
  !any(is.na(DEG$P.Value)),
  !any(is.na(DEG$logFC))
)

# Build volcano data: transform p-values to -log10 scale
vol_data <- data.frame(
  gene = DEG$symbol,
  pval = -log10(DEG$P.Value),   # -log10 scale: higher = more significant
  lfc  = DEG$logFC              # log2 fold-change: positive = up-regulated in PD
)

# ============================================================
# STEP 3 — Classify Each Gene by Direction and Significance
# ============================================================

# Genes are classified into three groups:
#   "Increased"      : logFC >  FC_THRESH AND -log10(p) >= P_THRESH
#   "Decreased"      : logFC < -FC_THRESH AND -log10(p) >= P_THRESH
#   "nonsignificant" : all other genes
# NOTE: Use bare column names inside mutate() — avoid vol_data$col syntax,
# which is redundant and fragile within dplyr pipelines.
vol_data <- vol_data %>%
  mutate(
    color = ifelse(
      lfc >  FC_THRESH & pval >= P_THRESH, "Increased",
      ifelse(
        lfc < -FC_THRESH & pval >= P_THRESH, "Decreased",
        "nonsignificant"
      )
    )
  )

# Report classification counts as a quick diagnostic
message("Volcano classification — ",
        "Increased: ",      sum(vol_data$color == "Increased"),
        " | Decreased: ",   sum(vol_data$color == "Decreased"),
        " | Not significant: ", sum(vol_data$color == "nonsignificant"))

# ============================================================
# STEP 4 — Define Target Gene Lists
# ============================================================

# Target genes represent candidates of biological interest for each dataset.
# These are highlighted with gold points and labeled on the volcano plot.

# Target genes for GSE7621 (stress response, dopaminergic, transcription factors)
target_genes_7621 <- c(
  "ATF4",  "CAPN2",  "BANF1",  "PNPT1",  "PTPRK",
  "FOXO3", "STAU1",  "ETV5",   "RELA",
  "KEAP1", "ATM",    "DHFR",   "ABL1",
  "TH",    "AGTR1",  "ALDH1A1","DLK1",   "SLC18A2","SLC6A3","PCSK1"
)

# Target genes for GSE49036 (oxidative stress, dopaminergic markers)
target_genes_49036 <- c(
  "NFE2L2", "HMOX1", "TXNIP", "CYBA", "S100A12",
  "TH", "AGTR1", "ALDH1A1", "DLK1", "SLC18A2", "SLC6A3", "PCSK1"
)

# ── IMPORTANT: Select the correct gene list for the current dataset ──────────
# Change this assignment when switching between GSE7621 and GSE49036.
# Mismatching the gene list with the dataset will produce incorrect highlights.
target_genes <- target_genes_7621  # <-- SET THIS: target_genes_7621 OR target_genes_49036
message("Active target gene list: ",
        ifelse(identical(target_genes, target_genes_7621), "GSE7621", "GSE49036"),
        " (", length(target_genes), " genes)")

# ============================================================
# STEP 5 — Annotate Target Genes for Highlighting and Labeling
# ============================================================

vol_data_labeled <- vol_data %>%
  mutate(
    # Label only target genes that also pass the significance threshold
    label = ifelse(gene %in% target_genes & pval >= P_THRESH, gene, NA),
    
    # Highlight column marks target genes above the threshold for a separate layer
    highlight = ifelse(gene %in% target_genes & pval >= P_THRESH, "Target", "Other")
  )

# Report how many target genes were not labeled (below significance threshold)
n_target_total    <- sum(vol_data_labeled$gene %in% target_genes)
n_target_labeled  <- sum(!is.na(vol_data_labeled$label))
n_target_unlabeled <- n_target_total - n_target_labeled
message("Target gene labeling — Labeled (significant): ", n_target_labeled,
        " | Unlabeled (below threshold): ", n_target_unlabeled,
        " | Not found in DEG: ",
        length(target_genes) - n_target_total)

# ============================================================
# STEP 6 — Draw Volcano Plot
# ============================================================

vol <- ggplot(vol_data_labeled, aes(x = lfc, y = pval, color = color)) +
  
  ggtitle(
    label    = "Volcano Plot",
    subtitle = "Colored by fold-change direction | Gold = target genes"
  ) +
  
  # ── Layer 1: All genes as background points ───────────────
  geom_point(size = 1, alpha = 0.6) +
  
  # ── Layer 2: Highlighted target genes (gold, larger, with border) ──
  # x and y are inherited from the global aes() — no need to re-specify.
  # color and fill are set outside aes() as fixed aesthetics (not mapped to a scale).
  geom_point(
    data   = vol_data_labeled %>% dplyr::filter(highlight == "Target"),
    color  = "black",   # Border color
    fill   = "gold",    # Fill color (shape 21 allows separate fill and border)
    shape  = 21,        # Fillable circle shape
    size   = 3,
    stroke = 0.5,       # Border thickness in mm
    alpha  = 0.8
  ) +
  
  # ── Color scale for background points ────────────────────
  scale_color_manual(
    name   = "Directionality",
    values = c(
      Increased      = "red",
      Decreased      = "blue",
      nonsignificant = "darkgray"
    )
  ) +
  
  # ── Layer 3: Labels for significant target genes ──────────
  # Labels are only shown for target genes passing the P and FC thresholds.
  # segment.linewidth requires ggrepel >= 0.9.2 (use segment.size on older versions)
  geom_text_repel(
    data             = vol_data_labeled %>% dplyr::filter(!is.na(label)),
    aes(label        = label),
    size             = 3.5,
    fontface         = "bold",
    color            = "black",
    box.padding      = 0.5,
    point.padding    = 0.2,
    arrow            = arrow(length = unit(0.01, "npc"), type = "closed"),
    segment.color    = "gold",
    segment.linewidth = 0.5,
    max.overlaps     = 30,
    direction        = "both"
  ) +
  
  # ── Significance threshold lines ─────────────────────────
  # Horizontal dashed line at -log10(0.05): marks p = 0.05 significance boundary
  geom_hline(
    yintercept = P_THRESH,
    colour     = "black",
    linetype   = "dashed",
    linewidth  = 0.3
  ) +
  # Vertical dashed lines at ±FC_THRESH: mark minimum fold-change boundaries
  geom_vline(
    xintercept = c(-FC_THRESH, FC_THRESH),
    colour     = "black",
    linetype   = "dashed",
    linewidth  = 0.3
  ) +
  
  # ── Axis labels and theme ─────────────────────────────────
  xlab(expression(log[2]("PD" / "CON"))) +
  ylab(expression(-log[10]("p-value"))) +
  theme_bw(base_size = 14) +
  theme(legend.position = "right")

print(vol)