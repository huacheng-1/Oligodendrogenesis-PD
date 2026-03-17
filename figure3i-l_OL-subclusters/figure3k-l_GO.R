# ============================================================
# GO Enrichment Analysis for Oligodendrocyte Subclusters
# Input: markers.RDS — FindAllMarkers() output from subclustered OL object
# Focus: Protein folding / heat stress response terms per cluster
# ============================================================

# ── Dependencies ─────────────────────────────────────────────
library(clusterProfiler) 
library(org.Hs.eg.db)    
# Replace with org.Mm.eg.db for mouse data
library(ggplot2)         
library(dplyr)           

# ============================================================
# STEP 1 — Load Subcluster Marker Data and Select Target Cluster
# ============================================================

markers <- readRDS("/data/h01005/hc/RE/markers.RDS")

# Set the cluster number to analyze.
# Valid values in this dataset: "3" or "6" (OL subclusters)
cluster <- "6"  # <-- SET THIS: "6" or "3"
message("Analyzing cluster: ", cluster)

# Filter to significant up-regulated markers for the selected cluster.
# NOTE: avg_log2FC > 0 is already required by the filter, so all retained
# genes will have change = "Up". The mutate() is kept for downstream
# compatibility but the "Down" branch is unreachable here.
cluster_markers_up <- markers %>%
  filter(cluster == !!cluster & p_val < 0.05 & avg_log2FC > 0)

message("Up-regulated markers in cluster ", cluster, ": ", nrow(cluster_markers_up))

if (nrow(cluster_markers_up) == 0) {
  stop("No significant up-regulated markers found for cluster ", cluster,
       ". Check the cluster ID and filtering thresholds.")
}

# ============================================================
# STEP 2 — Convert Gene Symbols to Entrez IDs
# ============================================================

# bitr() maps gene symbols to Entrez IDs required by clusterProfiler.
# Genes that cannot be mapped are silently dropped — report the loss.
gene_list <- bitr(
  cluster_markers_up$gene,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

n_input   <- nrow(cluster_markers_up)
n_mapped  <- nrow(gene_list)
message("Gene ID mapping — Input: ", n_input,
        " | Mapped: ",   n_mapped,
        " | Dropped: ",  n_input - n_mapped,
        " (", round(100 * (n_input - n_mapped) / n_input, 1), "%)")

gene_entrez <- gene_list$ENTREZID

# ============================================================
# STEP 3 — Run GO Enrichment Analysis
# ============================================================

# ont = "ALL": restricts to Biological Process terms only, consistent
# with the "GO Biological Process" axis label used in the plot.
# If MF or CC terms are also needed, change to ont = "ALL" and
# update the axis label and title accordingly.
#
# Two-step filtering: enrichGO uses pvalueCutoff = 0.1 as a wide net,
# then we manually tighten to p.adjust < 0.05 in STEP 4.
go_enrich <- enrichGO(
  gene          = gene_entrez,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "ALL",   # Biological Process only
  pAdjustMethod = "BH",   # Benjamini-Hochberg FDR correction
  pvalueCutoff  = 0.1,    # Wide initial filter (tightened below)
  qvalueCutoff  = 0.1
)

go_df <- as.data.frame(go_enrich)
message("GO-ALL terms returned (p.adjust < 0.1): ", nrow(go_df))

# ============================================================
# STEP 4 — Filter to Significant Terms and Compute Enrichment Score
# ============================================================

significant_go <- go_df[go_df$p.adjust < 0.05, ]
message("GO-ALL terms after p.adjust < 0.05 filter: ", nrow(significant_go))

if (nrow(significant_go) == 0) {
  stop("No GO terms passed p.adjust < 0.05. ",
       "Consider relaxing the cutoff or checking the input gene list.")
}

# Enrichment score: -log10(adjusted p-value); higher = more significant
significant_go$EnrichmentScore <- -log10(significant_go$p.adjust)

# Print available descriptions for keyword selection reference
message("Available GO descriptions:\n",
        paste(significant_go$Description, collapse = "\n"))

# ============================================================
# STEP 5 — Define Per-Cluster Keywords for Term Filtering
# ============================================================

# Keywords are matched exactly against GO term Description fields.
# Each cluster has a distinct set of relevant biological processes.

# Cluster 6: dominated by chaperone-mediated folding and heat stress
oxidative_stress_keywords_6 <- c(
  "chaperone-mediated protein folding",
  "protein folding",
  "response to heat",
  "cellular response to heat",
  "response to unfolded protein",
  "response to temperature stimulus",
  "response to topologically incorrect protein"
)

# Cluster 3: overlapping heat stress terms with de novo folding regulation
oxidative_stress_keywords_3 <- c(
  "chaperone-mediated protein folding",
  "'de novo' protein folding",
  "response to topologically incorrect protein",
  "response to unfolded protein",
  "response to heat",
  "cellular response to heat",
  "regulation of cellular response to heat"
)

# Select the keyword list matching the active cluster
# IMPORTANT: Change this when switching cluster values above
oxidative_stress_keywords <- switch(
  cluster,
  "6" = oxidative_stress_keywords_6,
  "3" = oxidative_stress_keywords_3,
  stop("No keyword list defined for cluster '", cluster, "'.")
)
message("Active keyword list: cluster ", cluster,
        " (", length(oxidative_stress_keywords), " keywords)")

# ============================================================
# STEP 6 — Filter GO Terms by Exact Keyword Matching
# ============================================================

# Exact matching: retain only terms whose Description is in the keyword list
significant_go_oxidative <- subset(
  significant_go,
  Description %in% oxidative_stress_keywords
)

message("GO terms matching keywords: ", nrow(significant_go_oxidative))

# Stop early if no terms matched — prevents cryptic ggplot errors
if (nrow(significant_go_oxidative) == 0) {
  stop("No GO terms matched the keyword list for cluster ", cluster, ". ",
       "Check that the keyword strings match the Description column exactly.")
}

# ============================================================
# STEP 7 — Visualize GO Enrichment Bar Chart
# ============================================================

# Upper x-axis limit with 15% padding to accommodate score labels
x_upper <- max(significant_go_oxidative$EnrichmentScore) * 1.15

ggplot(
  significant_go_oxidative,
  aes(
    x = EnrichmentScore,
    y = reorder(Description, EnrichmentScore)  # Sort bars by score
  )
) +
  
  # ── Bar layer ────────────────────────────────────────────
  geom_col(
    width = 0.6,
    fill  = "#2C7FB8",   # Blue — consistent with stress response themes
    color = "#1F5C87",   # Darker border for definition
    alpha = 0.8
  ) +
  
  # ── Score labels (right of each bar) ─────────────────────
  # Shows the numeric -log10(p.adjust) value to 2 decimal places.
  # NOTE: geom_label(Description) is intentionally omitted —
  # term names are already displayed on the y-axis via reorder().
  geom_text(
    aes(label = sprintf("%.2f", EnrichmentScore)),
    hjust = -0.2,
    size  = 3.5,
    color = "#333333"
  ) +
  
  # ── Row separator lines (between bars for readability) ───
  geom_hline(
    yintercept = seq(1.5, nrow(significant_go_oxidative), by = 1),
    color      = "#E0E0E0",
    linewidth  = 0.3
  ) +
  
  # ── Axis labels and titles ────────────────────────────────
  labs(
    x        = "Enrichment Score (-log10 adjusted p-value)",
    y        = "GO Biological Process",
    title    = paste0("Protein Folding / Heat Stress GO Enrichment — Cluster ", cluster),
    subtitle = "Significant terms (p.adjust < 0.05) | Exact keyword match"
  ) +
  
  # ── X-axis scale with right-side padding for labels ──────
  scale_x_continuous(
    limits = c(0, x_upper)
  ) +
  
  # ── Theme ────────────────────────────────────────────────
  # NOTE: axis.text.y face is set to a single fixed value ("plain").
  # Applying ifelse() to produce a per-row vector inside element_text()
  # causes an error because element_text(face =) expects a scalar, not
  # a vector. Use geom_text with fontface = "bold" for selective bolding.
  theme_minimal(base_size = 12) +
  theme(
    axis.title         = element_text(face = "bold", color = "#333333"),
    axis.text.x        = element_text(color = "#333333", size = 10),
    axis.text.y        = element_text(color = "#333333", size = 10,
                                      face = "plain"),
    plot.title         = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle      = element_text(size = 10, hjust = 0.5,
                                      margin = margin(b = 10)),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color = "#F5F5F5", linewidth = 0.2)
  )