library(pheatmap)
#exp2为前处理中的表达矩阵数据
# 1. 提取目标基因的表达数据
##这是两个数据集的目标基因，用于不同数据集的分析。因为处理步骤相同故放在一起
#GSE7621
target_genes <- c("ATF4", "CAPN2", "BANF1", "PNPT1", "PTPRK", 
                  "FOXO3", "SIRPA", "STAU1", "ETV5", "RELA", 
                  "KEAP1", "ATM", "DHFR", "ABL1")

#GSE49036
target_genes <- c("KEAP1","ATM","FOXO3","RELA","SETD2","SIRPA","BANF1")

heatmap_data <- exp2[rownames(exp2) %in% target_genes, ]

# 2. 确保 metadata 的样本顺序与热图数据列顺序一致
##这是两个数据集的样本组成，用于不同数据集的分析。因为处理步骤相同故放在一起
#GSE7621
metadata <- data.frame(
  Sample = colnames(heatmap_data),  # 确保样本名匹配
  Group = factor(c(rep("CON", 9), rep("PD", 16)), levels = c("CON", "PD"))
)

#GSE49036
metadata <- data.frame(
  Sample = colnames(heatmap_data),  # 确保样本名匹配
  Group = factor(c(rep("CON",8),rep("PD",8)), levels = c("CON", "PD"))
)

rownames(metadata) <- metadata$Sample  # 设置行名为样本名（以便匹配）

# 3. 绘制热图，并添加列注释（分组信息）
p <- pheatmap(
  heatmap_data,
  scale = "row",  # 按行标准化（Z-score）
  annotation_legend = TRUE,
  cluster_cols = FALSE,
  cluster_rows = F,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  color = colorRampPalette(c("black","navy", "white","yellow", "firebrick3"))(20),
  show_rownames = TRUE,
  show_colnames = TRUE,
  annotation_col = metadata["Group"],
  annotation_colors = list(Group = c(CON = "blue", PD = "red")),
  main = "Expression Heatmap with Sample Groups in GSE49036", #"Expression Heatmap with Sample Groups in GSE7621"
  # 新增参数：调整行间距和列名倾斜
  cellwidth = 8,       # 调整单元格宽度（可选）
  cellheight = 8,       # 调整单元格高度（缩小行间距）
  angle_col = "45",     # 列名倾斜 45 度
  fontsize_row = 8,     # 调整行名字体大小（避免重叠）
  fontsize_col = 6      # 调整列名字体大小（可选）
)

z_score_data <- t(scale(t(heatmap_data)))  # 按行标准化
