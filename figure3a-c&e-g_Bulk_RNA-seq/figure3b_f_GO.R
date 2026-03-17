# ============================================================
# GO Enrichment Analysis and Visualization
# Input: diff_gene — DEG results from GSE7621 or GSE49036
# Focuses on up-regulated genes and oxidative-stress-related GO terms
# ============================================================

# ── Dependencies ─────────────────────────────────────────────
library(clusterProfiler) 
library(org.Hs.eg.db)    
# Replace with org.Mm.eg.db for mouse data
library(ggplot2)         
library(dplyr)          

# ============================================================
# STEP 1 — Extract Up-regulated Genes from DEG Results
# ============================================================

# Retain only genes classified as significantly up-regulated in PD vs CON
gene_up <- subset(diff_gene, change == "Up")
message("Up-regulated genes selected: ", nrow(gene_up))

# ============================================================
# STEP 2 — Convert Gene Symbols to Entrez IDs
# ============================================================

# bitr() maps gene symbols to Entrez IDs required by clusterProfiler.
# Genes that cannot be mapped are silently dropped — report the loss.
gene_list <- bitr(
  gene_up$symbol,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

n_input   <- nrow(gene_up)
n_mapped  <- nrow(gene_list)
n_dropped <- n_input - n_mapped
message("Gene ID mapping — Input: ", n_input,
        " | Mapped: ", n_mapped,
        " | Dropped (unmapped): ", n_dropped,
        " (", round(100 * n_dropped / n_input, 1), "%)")

# Extract the successfully mapped Entrez IDs
gene_entrez <- gene_list$ENTREZID

# ============================================================
# STEP 3 — Run GO Enrichment Analysis
# ============================================================

# ont = "BP": restricts analysis to Biological Process terms only.
# Using "ALL" would include MF and CC, which conflicts with the
# "GO Biological Process" axis label used in the plot below.
# If MF/CC are needed, change ont = "ALL" and update the plot label.
#
# Two-step filtering strategy:
#   enrichGO uses pvalueCutoff = 0.1 to cast a wide net,
#   then we manually filter to p.adjust < 0.05 in STEP 4.
#   This allows inspection of the 0.05–0.10 borderline terms if needed.
go_enrich <- enrichGO(
  gene          = gene_entrez,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",             # Biological Process only
  pAdjustMethod = "BH",             # Benjamini-Hochberg FDR correction
  pvalueCutoff  = 0.1,              # Wide initial net (filtered further below)
  qvalueCutoff  = 0.1
)

go_df <- as.data.frame(go_enrich)
message("Total GO-BP terms returned (p.adjust < 0.1): ", nrow(go_df))

# ============================================================
# STEP 4 — Filter to Significant Terms and Compute Enrichment Score
# ============================================================

# Tighten the significance filter from 0.1 to 0.05
significant_go <- go_df[go_df$p.adjust < 0.05, ]
message("GO-BP terms after p.adjust < 0.05 filter: ", nrow(significant_go))

if (nrow(significant_go) == 0) {
  stop("No GO terms passed the p.adjust < 0.05 threshold. ",
       "Consider relaxing the cutoff or checking input gene list.")
}

# Compute enrichment score as -log10(adjusted p-value):
# higher score = more significant enrichment
significant_go$EnrichmentScore <- -log10(significant_go$p.adjust)

# ============================================================
# STEP 5 — Define Oxidative Stress Keywords Per Dataset
# ============================================================

# Keywords are used for fuzzy matching against GO term descriptions.
# Each list is tailored to the biological signals observed in that dataset.

# Keywords for GSE7621 (hypoxia, apoptosis, NF-κB, HIF-1 signaling)
oxidative_stress_keywords_GSE7621 <- c(
  "response to hypoxia",
  "response to decreased oxygen levels",
  "response to oxygen levels",
  "canonical NF-kappaB signal transduction",
  "regulation of apoptotic signaling pathway",
  "HIF-1 signaling pathway",
  "intrinsic apoptotic signaling pathway"
)

# Keywords for GSE49036 (inflammation, phagocytosis, PI3K/Akt, TNF, xenobiotic)
oxidative_stress_keywords_GSE49036 <- c(
  "positive regulation of inflammatory response",
  "phagocytosis",
  "positive regulation of phosphatidylinositol 3-kinase/protein kinase B signal transduction",
  "tumor necrosis factor production",
  "canonical NF-kappaB signal transduction",
  "response to xenobiotic stimulus",
  "cellular response to abiotic stimulus",
  "cellular response to environmental stimulus"
)

# ── IMPORTANT: Select the keyword list matching the active dataset ────────────
# Change this assignment when switching between GSE7621 and GSE49036.
# Mismatching the keyword list with the dataset will silently produce
# irrelevant or empty results.
oxidative_stress_keywords <- oxidative_stress_keywords_GSE7621
# oxidative_stress_keywords <- oxidative_stress_keywords_GSE49036

# Confirm which keyword list is active at runtime
message("Active keyword list: ",
        ifelse(
          identical(oxidative_stress_keywords, oxidative_stress_keywords_GSE7621),
          "GSE7621", "GSE49036"
        ),
        " (", length(oxidative_stress_keywords), " keywords)")

# ============================================================
# STEP 6 — Filter GO Terms by Keyword Matching
# ============================================================

# Fuzzy matching: retain GO terms whose Description contains
# at least one of the oxidative stress keywords (case-insensitive).
# sapply() returns a logical matrix (rows = GO terms, cols = keywords).
# rowSums() > 0 flags terms matching at least one keyword.
keyword_match_matrix <- sapply(
  oxidative_stress_keywords,
  function(kw) grepl(kw, significant_go$Description, ignore.case = TRUE)
)
matched_rows <- rowSums(keyword_match_matrix) > 0

significant_go_oxidative <- significant_go[matched_rows, ]
message("GO terms matching oxidative stress keywords: ",
        nrow(significant_go_oxidative))

# Stop early with a clear error if no terms matched — prevents cryptic ggplot errors
if (nrow(significant_go_oxidative) == 0) {
  stop("No GO terms matched the oxidative stress keywords. ",
       "Check that the correct keyword list is selected and that ",
       "significant_go is non-empty.")
}

# Sort by enrichment score (descending) for consistent plot ordering
significant_go_oxidative <- significant_go_oxidative[
  order(-significant_go_oxidative$EnrichmentScore), 
]

# ============================================================
# STEP 7 — Plot GO Enrichment Bar Chart
# ============================================================

# The +1 padding on the x-axis upper limit prevents the rightmost bar
# from being clipped at the panel edge.
x_upper <- max(significant_go_oxidative$EnrichmentScore) + 1

ggplot(
  significant_go_oxidative,
  aes(
    x = EnrichmentScore,
    y = reorder(Description, EnrichmentScore)  # Sort bars by score
  )
) +
  # Bar layer: length encodes enrichment magnitude
  geom_bar(stat = "identity", fill = "steelblue", width = 0.7) +
  
  # NOTE: geom_text(label = Description) is intentionally omitted —
  # term descriptions are already displayed on the y-axis via reorder(),
  # so adding them as bar labels would create redundant, overlapping text.
  
  labs(
    x     = "Enrichment Score (-log10 adjusted p-value)",
    y     = "GO Biological Process",
    title = "GO Biological Process Enrichment",
    subtitle = paste0(
      "Up-regulated genes | Keywords: ",
      ifelse(
        identical(oxidative_stress_keywords, oxidative_stress_keywords_GSE7621),
        "GSE7621", "GSE49036"
      )
    )
  ) +
  # Prevent bar clipping with a small right-side margin
  coord_cartesian(xlim = c(0, x_upper)) +
  theme_minimal() +
  theme(
    axis.text.y  = element_text(size = 10, color = "black"),
    axis.text.x  = element_text(size = 10, color = "black"),
    axis.title   = element_text(size = 12, face = "bold"),
    plot.title   = element_text(size = 13, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray40")
  )