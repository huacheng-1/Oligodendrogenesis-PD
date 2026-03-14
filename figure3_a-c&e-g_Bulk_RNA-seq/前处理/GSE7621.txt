library(DESeq2) 
library(ggplot2)
library(qvalue)
library(cluster)
library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(ggrepel)
library(stringr)
library(ggupset)
library(enrichplot)
library(DOSE)
library(clusterProfiler)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(msigdbr)
library(fgsea)
library(Seurat)
library(tidyverse)
library(GEOquery)
library(tinyarray)
library(limma)
geo_number <- "GSE7621"
PATH <- 'D:/R/GEO/20250917-新文章/GSE7621'
geo_data <- getGEO(geo_number, destdir = PATH, getGPL = TRUE)
geo_data0 <- geo_data[[1]]
exp <- exprs(geo_data0)

# 2. 注释处理
ids <- data.table::fread("D:/R/GEO/20250917-新文章/GSE7621/GPL570-55999.txt", sep = "\t", header = T)


# 3. 表达矩阵整理
exp <- as.data.frame(exp)
exp$ID <- rownames(exp)
merged <- inner_join(ids[, c("ID", "Gene Symbol")], exp, by = "ID")
merged <- na.omit(merged) %>%
  dplyr::filter(`Gene Symbol` != "")

# 去重：保留高表达基因
merged <- merged %>% 
  group_by(`Gene Symbol`) %>% 
  slice_max(order_by = rowMeans(across(where(is.numeric)), na.rm = TRUE), n = 1) %>%
  ungroup() %>%
  na.omit()

exp2 <- merged %>% 
  tibble::column_to_rownames("Gene Symbol") %>% 
  dplyr::select(-ID) %>% 
  as.matrix()
exp2 <- log2(exp2 + 1) %>% na.omit()
#saveRDS(exp2,"D:/R/GEO/20250917-新文章/GSE7621/exp2.RDS")
# 4. 差异分析
# 统一使用 model.matrix 和 adj.P.Val
metadata <- data.frame(Sample = colnames(exp2), 
                       Group = factor(c(rep("CON", 9), rep("PD", 16)), levels = c("CON", "PD")))
design <- model.matrix(~ Group, data = metadata)  # 自动处理参考组

contrast.matrix <- makeContrasts(GroupPD = GroupPD, levels = design)
fit <- lmFit(exp2, design)
fit2 <- contrasts.fit(fit, contrast.matrix) %>% eBayes(trend = TRUE)

DEG <- topTable(fit2, coef = "GroupPD", n = Inf, adjust.method = "fdr") %>%
  rownames_to_column("symbol") %>%
  mutate(change = ifelse(P.Value < 0.05 & abs(logFC) >= 0.5, 
                         ifelse(logFC > 0.5, 'Up', 'Down'), 'Not'))

diff_gene <- DEG %>% filter(P.Value < 0.05 & abs(logFC) >= 0.5)#KEGG分析不使用abs(logFC) >= 0.5
# 5. 结果输出
write.csv(DEG, paste0(PATH, "/result/LIMMA_full_results.csv"), row.names = FALSE)
write.table(diff_gene, paste0(PATH, "/result/significant_genes.txt"), 
            sep = "\t", quote = FALSE, row.names = FALSE)