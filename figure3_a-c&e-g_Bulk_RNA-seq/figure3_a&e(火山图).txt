library(dplyr)
library(ggrepel)
library(DESeq2)
library(ggplot2)
library(qvalue)
library(cluster)
library(pheatmap)
library(RColorBrewer)
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

#DEG来源于前处理后的差异基因结果
#火山图
vol_data <- data.frame(gene=DEG$symbol,
                       pval=-log10(DEG$P.Value),
                       lfc=DEG$logFC)
# 设定上调与下调颜色
vol_data <- mutate(vol_data,
                   color=ifelse(vol_data$lfc > 0.5 & vol_data$pval > 1.3 ,"Increased",
                                ifelse(vol_data$lfc < -0.5 & vol_data$pval > 1.3 ,"Decreased",'nonsignificant'))
)

##这是两个数据集的目标基因，用于不同数据集的分析。因为处理步骤相同故放在一起
#GSE7621
target_genes <- c("ATF4", "CAPN2", "BANF1", "PNPT1", "PTPRK", 
                  "FOXO3", "STAU1", "ETV5", "RELA", 
                  "KEAP1", "ATM", "DHFR", "ABL1",
                  "TH", "AGTR1", "ALDH1A1", "DLK1", "SLC18A2", "SLCA63", "PCSK1")

#GSE49036
target_genes <- c("NFE2L2","HMOX1","TXNIP","CYBA","S100A12",
                  "TH", "AGTR1", "ALDH1A1", "DLK1", "SLC18A2", "SLCA63", "PCSK1")

# 标记目标基因，并新增一列用于高亮
vol_data_labeled <- vol_data %>%
  mutate(
    label = ifelse(gene %in% target_genes & pval >= 1.3, gene, NA),
    highlight = ifelse(gene %in% target_genes & pval >= 1.3, "Target", "Other")  # 新增高亮分组
  )

# 绘制火山图
vol <- ggplot(vol_data, aes(x = lfc, y = pval, color = color)) + 
  ggtitle(label = "Volcano Plot", subtitle = "Colored by fold-change direction") +
  
  # 普通点（背景）
  geom_point(size = 1, alpha = 0.6, na.rm = TRUE) +
  
  # 高亮目标基因点（增大尺寸 + 不同颜色 + 边框）
  geom_point(
    data = vol_data_labeled %>% filter(highlight == "Target"),
    aes(x = lfc, y = pval),
    color = "black",       # 黑色边框
    fill = "gold",         # 金色填充（可改其他颜色）
    shape = 21,            # 21表示可填充的点
    size = 3,              # 增大尺寸
    stroke = 0.5,          # 边框粗细
    alpha = 0.8            # 稍透明
  ) +
  
  # 自定义颜色映射
  scale_color_manual(
    name = "Directionality",
    values = c(Increased = "red", Decreased = "blue", nonsignificant = "darkgray")
  ) +
  
  # 添加标签（仅目标基因）
  geom_text_repel(
    data = vol_data_labeled %>% filter(!is.na(label)),
    aes(label = label),
    size = 3.5,                     # 增大字体
    fontface = "bold",              # 加粗
    color = "black",                # 黑色文字
    box.padding = 0.5,              # 标签与点的间距
    point.padding = 0.2,            # 避免标签覆盖其他点
    arrow = arrow(length = unit(0.01, "npc"), type = "closed"),
    segment.color = "gold",         # 箭头线颜色（与点一致）
    segment.size = 0.5,             # 箭头线粗细
    max.overlaps = 30,              # 允许更多重叠
    direction = "both"
  ) +
  
  # 添加阈值线
  geom_hline(yintercept = 1.3, colour = "black", linetype = "dashed", linewidth = 0.3) +
  geom_vline(xintercept = c(-0.5, 0.5), colour = "black", linetype = "dashed", linewidth = 0.3) +
  
  # 主题调整
  theme_bw(base_size = 14) +
  theme(legend.position = "right") +
  xlab(expression(log[2]("pd" / "con"))) +
  ylab(expression(-log[10]("adjusted p-value")))

# 显示图形
print(vol)