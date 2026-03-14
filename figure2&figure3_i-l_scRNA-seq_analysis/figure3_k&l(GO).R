library(clusterProfiler)
library(org.Hs.eg.db)  # 人类基因数据库（其他物种需更换，如 org.Mm.eg.db 用于小鼠）
cluster <- "6" # 选择聚类编号
cluster_markers_up <- markers %>% 
  filter(cluster == !!cluster & p_val < 0.05 & avg_log2FC > 0) %>% 
  mutate(rank_pos = ifelse(avg_log2FC > 0, rank(-avg_log2FC), NA),
         rank_neg = ifelse(avg_log2FC < 0, rank(avg_log2FC), NA)) %>% 
  mutate(change = ifelse(abs(avg_log2FC) >= 0, 
                         ifelse(avg_log2FC > 0, 'Up', 'Down'), 'Not')) %>%
  dplyr::select(-rank_pos, -rank_neg)

cluster_markers_down <- markers %>% 
  filter(cluster == !!cluster & p_val < 0.05 & avg_log2FC < 0) %>% 
  mutate(rank_pos = ifelse(avg_log2FC > 0, rank(-avg_log2FC), NA),
         rank_neg = ifelse(avg_log2FC < 0, rank(avg_log2FC), NA)) %>% 
  mutate(change = ifelse(abs(avg_log2FC) >= 0, 
                         ifelse(avg_log2FC > 0, 'Up', 'Down'), 'Not')) %>%
  dplyr::select(-rank_pos, -rank_neg)

# 将基因符号转换为 Entrez ID（如果输入是 Gene Symbol）
gene_list <- bitr(cluster_markers_up$gene, fromType = "SYMBOL", 
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
significant_go$Description

#聚类6
oxidative_stress_keywords <- c("chaperone-mediated protein folding","protein folding","response to heat",
                               "cellular response to heat","response to unfolded protein","response to temperature stimulus",
                               "response to topologically incorrect protein")
#聚类3
oxidative_stress_keywords <- c("chaperone-mediated protein folding","'de novo' protein folding",
                               "response to topologically incorrect protein",
                               "response to unfolded protein", "response to heat",
                               "cellular response to heat","regulation of cellular response to heat")

#精确匹配
significant_go_oxidative <- subset(significant_go, subset = Description %in% oxidative_stress_keywords)

# 按富集分数排序
significant_go_oxidative <- significant_go_oxidative[
  order(-significant_go_oxidative$EnrichmentScore), 
]


# 绘制条形图
# 优化后的绘图代码
ggplot(significant_go_oxidative, 
       aes(x = EnrichmentScore, 
           y = reorder(ID, EnrichmentScore))) +
  # 使用更专业的柱状图样式
  geom_col(width = 0.6, 
           fill = "#2C7FB8",  # 蓝色系，适合氧化磷酸化主题
           color = "#1F5C87",  # 边框加深
           alpha = 0.8) +
  
  # 添加数值标签（优化显示）
  geom_text(aes(label = sprintf("%.2f", EnrichmentScore)),  # 显示两位小数
            hjust = -0.2,  # 标签在柱子右侧
            size = 3.5, 
            color = "#333333") +
  
  # 添加参考线（突出显著性阈值）
  geom_hline(yintercept = seq(1.5, nrow(significant_go_oxidative), by = 1), 
             color = "#E0E0E0", 
             size = 0.3) +
  
  # 坐标轴和标题优化
  labs(x = "Enrichment Score (-log10(p.adjust))", 
       y = "GO Biological Process",
       title = paste0("Oxidative Phosphorylation Pathway Enrichment in Cluster",cluster),
       subtitle = "Top significant GO terms (p.adjust < 0.05)") +
  
  # 主题优化（专业学术风格）
  theme_minimal(base_size = 12) +
  theme(
    axis.title = element_text(face = "bold", color = "#333333"),
    axis.text.x = element_text(color = "#333333", size = 10),
    axis.text.y = element_text(color = "#333333", 
                               face = ifelse(
                                 grepl("oxidative|mitochondrial|ATP", 
                                       significant_go_oxidative$ID),
                                 "bold", "plain")),  # 关键术语加粗
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, margin = margin(b = 10)),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color = "#F5F5F5", size = 0.2)
  ) +
  
  # 坐标轴范围优化
  scale_x_continuous(expand = expansion(mult = c(0, 0.1)),  # 右侧留白
                     limits = c(0, max(significant_go_oxidative$EnrichmentScore)*1.15)) +
  
  # 添加数据标签背景（可选）
  geom_label(aes(label = Description), 
             fill = "white", 
             alpha = 0.7,
             hjust = 1,
             vjust = 0.5,
             hjust = "inward")
