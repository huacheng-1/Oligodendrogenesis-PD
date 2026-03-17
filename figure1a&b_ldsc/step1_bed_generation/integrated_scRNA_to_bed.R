#===========================================================
# STEP 1: Generate cell-type-specific .bed files from single-cell RNA-seq data
#                  for use as LDSC (LD Score Regression) annotation input
# Input  : Integrated Seurat object (scRNA_combined) + NCBI hg19 gene coordinates
# Output : One .bed file per cell type (top 10% most specific genes, ±100 kb window)
#===========================================================

# --- Load required libraries ---
library(tidyverse)                 
library(snow)                        
library(Seurat)                      
library(org.Hs.eg.db)           
library(dplyr)                       
library(readr)                        

#===========================================================
# 1. Load and format gene coordinates from NCBI hg19 (GRCh37) gene location file
#       File format (no header) : ENTREZ | CHR | START | END | STRAND | SYMBOL
#             Window: ±100 kb flanking region added around each gene body
#===========================================================
gene_coordinates <- read_tsv(
	"/data/h01005/hc/magama/NCBI37.3.gene.loc",
	 col_names = FALSE, 
	 col_types = 'cciicc'	# col types: ENTREZ(c), CHR(c), START (i), END(i), ?(c), SYMBOL(c)
) %>%
  mutate(
# Extend gene boundaries by ±100 kb; clamp start at 0 to avoid negative coordinates
	start = ifelse(X3-100000 < 0, 0, X3-100000),
        end = X4+100000, 
	chr = X2, 
	ENTREZ = X1
) %>%
  dplyr::select(chr, start, end, ENTREZ)%>%  
# Add "chr" prefix to match standard BED chromosome naming (e.g., "1"-> "chr1")
  mutate(chr = paste0("chr", chr)) %>%
# Set Entrez ID as row names for downstream ID mapping
  column_to_rownames("ENTREZ") %>%
  dplyr::select(chr, start, end)

#===========================================================
# 2. Map Entrez IDs to HGNC gene symbols
#    Needed to join with expression data (which uses gene symbols, not Entrez IDs)
#    multiVals = "first": if one Entrez ID maps to multiple symbols, keep the first
#===========================================================

entrez_to_symbol <- mapIds(
	org.Hs.eg.db,
	keys = rownames(gene_coordinates),
	column = "SYMBOL",
	keytype = "ENTREZID",
	multiVals = "first"
)

# Attach gene symbols and restore Entrez IDs as a regular column
gene_coordinates$symbol <- entrez_to_symbol
gene_coordinates$ENTREZ <- rownames(gene_coordinates)


#===========================================================
# 3. Compute average gene expression per cell type from integrated scRNA-seq data
#
#	Input: scRNA_combined — integrated Seurat object (human single-cell dataset)
#	Group: metadata column "cells" (stores cell-type labels)
#===========================================================
avg_exp <- AverageExpression(
	scRNA_combined, 
	assays = "RNA",
	group.by = "cells", #  Metadata column containing cell-type annotations
	return.seurat = TRUE
)

#  Extract the normalized expression matrix (genes x cell types)
#  NOTE: Use layer = "data" for Seurat v5; replace with slot = "data" for Seurat v4
exp_matrix <- Seurat::GetAssayData(
	avg_exp,
	assay = "RNA",
	layer = "data" 
)

#===========================================================
# 4. Reshape expression matrix from wide to long format
#    Result columns:
#      Gene           — gene symbol (row identifier)
#      cell_type      — cell-type label (former column name)
#      Expr_sum_mean  — average normalized expression value
#===========================================================

exp_lvl5 <- exp_matrix %>%
  as.data.frame() %>%
  rownames_to_column("Gene") %>%
  pivot_longer(
    cols = -Gene,
    names_to = "cell_type",	# Cell-type label column
    values_to = "Expr_sum_mean"	 # Average normalized expression value
  )

#===========================================================
# 5. 	Remove genes with zero total expression across ALL cell types
#	These genes carry no specificity signal and would inflate noise downstream
#===========================================================

not_expressed <- exp_lvl5 %>%
  group_by(Gene) %>%
  summarise(total_sum = sum(Expr_sum_mean)) %>%
  filter(total_sum == 0) %>%
  pull(Gene)

exp_lvl5 <- exp_lvl5 %>%
  filter(!Gene %in% not_expressed)

# ===========================================================
# 6. Normalize expression within each cell type to counts per million (CPM)
#    Purpose: remove the effect of total expression level differences across
#             cell types, making values comparable for specificity calculation
#===========================================================


exp_lvl5 <- exp_lvl5 %>%
  group_by(cell_type) %>% 
  mutate(Expr_sum_mean = Expr_sum_mean * 1e6 / sum(Expr_sum_mean)) %>%
  ungroup()  # Always ungroup after grouped mutate to prevent unintended carry-over

#===========================================================
# 7. Calculate gene specificity score across cell types
#
#    Specificity = CPM in this cell type / sum of CPM across ALL cell types
#    Range: [0, 1]; value closer to 1 = expression more restricted to this cell type
#
#===========================================================
exp_lvl5 <- exp_lvl5 %>%
  group_by(Gene) %>%
  mutate(specificity = Expr_sum_mean / sum(Expr_sum_mean)) %>%
  ungroup()

#=================================================================
# 8. Select the top 10% most cell-type-specific genes per cell type for LDSC
#
#    Two-stage filter:
#      (a) CPM > 1  — gene must be detectably expressed in this cell type
#                     (removes lowly expressed genes that may have spuriously
#                      high specificity scores due to near-zero denominators)
#      (b) Top 10% by specificity score within each cell type
#
#    with_ties = FALSE: if tied values fall on the n boundary, drop them to
#                       return exactly n_genes_to_keep genes (ensures reproducibility)
#=================================================================

n_genes <- length(unique(exp_lvl5$Gene))
n_genes_to_keep <- round(n_genes * 0.1)

ldsc_gene_input <- exp_lvl5 %>%
  filter(Expr_sum_mean > 1) %>%		   # (a) expression filter
  group_by(cell_type) %>%
  slice_max(
      order_by = specificity, 
      n = n_genes_to_keep,
      with_ties = FALSE             # Return exactly n genes; drop boundary ties
  ) %>%
  ungroup() %>%
  dplyr::select(cell_type, Gene)  # dplyr:: prefix avoids masking by other packages

#=============================================================================
# 9. Join selected genes with genomic coordinates (matched by gene symbol)
#    inner_join: only genes present in BOTH the expression data AND the
#    coordinate file are retained; unmatched genes are silently dropped
#=============================================================================

ldsc_input <- ldsc_gene_input %>%
  inner_join(
    gene_coordinates %>% 
      dplyr::select(chr, start, end, symbol),
    by = c("Gene" = "symbol")  # Join on: expression "Gene" == coordinates "symbol"
  )
#=============================================================================
# 10. Write one BED file per cell type
#
#     BED format (tab-delimited, no header): chr | start | end | gene_name
#     make.names() converts cell-type labels to valid filenames:
#       spaces and special characters are replaced with "."
#       e.g., "Ex Neuron" -> "Ex.Neuron"
#=============================================================================

dir.create("/data/h01005/hc/ldsc/Bed", recursive = TRUE, showWarnings = FALSE)

ldsc_input %>%
  split(.$cell_type) %>%
  walk(~ write_tsv(
    .x %>% dplyr::select(chr, start, end, Gene),
    file = paste0(                          # file= explicit for readr >= 2.0
      "/data/h01005/hc/ldsc/Bed/",
      make.names(unique(.x$cell_type)),          # Cell type label -> valid filename
      ".bed"
    ),
    col_names = FALSE                            # BED format requires no header line
  ))
#####################################################################

