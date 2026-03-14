library(clusterProfiler)
library(org.Hs.eg.db)  # 人类基因数据库（其他物种需更换，如 org.Mm.eg.db 用于小鼠）
#读取差异基因数据命名为diff_gene
gene_up <- subset(diff_gene, subset = change %in% c("Up"))
gene_down <- subset(diff_gene, subset = change %in% c("Down"))

# 将基因符号转换为 Entrez ID（如果输入是 Gene Symbol）
gene_list <- bitr(gene_up$symbol, fromType = "SYMBOL", 
                  toType = "ENTREZID", 
                  OrgDb = org.Hs.eg.db)

# 提取 Entrez ID
gene_entrez <- gene_list$ENTREZID

# 运行 GO 富集分析（BP, MF, CC 分别分析）
go_enrich <- enrichGO(gene = gene_entrez,
                      OrgDb = org.Hs.eg.db,
                      keyType = "ENTREZID",
                      ont = "ALL",  # 可选 "BP", "MF", "CC" 或 "ALL"
                      pAdjustMethod = "none",
                      pvalueCutoff = 0.1,
                      qvalueCutoff = 0.1)  # 可选背景基因集 

go_df <- as.data.frame(go_enrich)
# 过滤显著结果（例如 p.adjust < 0.05）
significant_go <- go_df[go_df$p.adjust < 0.05, ]
significant_go$EnrichmentScore <- -log10(significant_go$p.adjust)
#significant_go$Description

# 筛选与氧化应激相关的通路
##这是两个数据集的目标基因，用于不同数据集的分析。因为处理步骤相同故放在一起
#GSE7621
oxidative_stress_keywords <- c("response to hypoxia", "response to decreased oxygen levels", "response to oxygen levels",
                               "canonical NF-kappaB signal transduction", "regulation of apoptotic signaling pathway",
                               "HIF-1 signaling pathway", "intrinsic apoptotic signaling pathway")
#GSE49036
oxidative_stress_keywords <- c("positive regulation of inflammatory response", "phagocytosis", 
                               "positive regulation of phosphatidylinositol 3-kinase/protein kinase B signal transduction",
                               "tumor necrosis factor production", 
                               "canonical NF-kappaB signal transduction",
                               "response to xenobiotic stimulus", 
                               "cellular response to abiotic stimulus","cellular response to environmental stimulus")

#模糊匹配
significant_go_oxidative <- significant_go[
  sapply(oxidative_stress_keywords, function(kw) 
    grepl(kw, significant_go$Description, ignore.case = TRUE)) %>% 
    rowSums() > 0,  # 至少匹配一个关键词
]

#精确匹配
significant_go_oxidative <- subset(significant_go, subset = Description %in% oxidative_stress_keywords)

# 按富集分数排序
significant_go_oxidative <- significant_go_oxidative[
  order(-significant_go_oxidative$EnrichmentScore), 
]


# 绘制条形图
ggplot(significant_go_oxidative, aes(x = EnrichmentScore, y = reorder(ID, EnrichmentScore))) +
  geom_bar(stat = "identity", fill = "steelblue", width = 0.7) +
  geom_text(aes(label = sprintf(Description)),  # 添加基因数量注释
            hjust = 1,vjust = 0, size = 5, color = "black") +
  labs(x = "Enrichment Score (-log10(p.adjust))", 
       y = "GO Biological Process",
       title = "Top GO Enrichment Terms") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 12, face = "bold")) +
  coord_cartesian(xlim = c(0, max(significant_go_oxidative$EnrichmentScore) + 1))  # 调整 x 轴范围
