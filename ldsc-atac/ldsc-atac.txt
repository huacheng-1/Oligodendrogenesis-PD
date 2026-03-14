library(data.table)

# 创建聚类与细胞类型映射关系
cell_type_map <- list(
  "Neurons" = c(1,2,3,4,5,6,7,11,12,18),
  "OPCs" = 8:10,
  "Astrocytes" = 13:17,
  "Oligodendrocytes" = 19:23,
  "Microglia" = 24
)

# 生成反向映射表（聚类编号 -> 细胞类型）
cluster_mapping <- data.table(
  Cluster = unlist(cell_type_map),
  Cell_Type = rep(names(cell_type_map), times = sapply(cell_type_map, length))
)

# 定义文件路径
data_dir <- "/data/h01005/hc/GSE147672/peak/"  # 请替换为实际路径
output_file <- "/data/h01005/hc/GSE147672/combined_peaks_with_celltypes.tsv.gz"

# 处理所有聚类文件
combined_data <- rbindlist(lapply(1:24, function(cluster_num) {
  # 读取文件
  file_path <- file.path(data_dir, paste0("Cluster", cluster_num, ".overlap.optimal.narrowPeak.gz"))
  dt <- fread(file_path, header = FALSE, sep = "\t", col.names = paste0("V", 1:10))
  
  # 添加细胞类型信息
  dt[, Cell_Type := cluster_mapping[Cluster == cluster_num, Cell_Type]]
  
  # 添加聚类编号列
  dt[, Cluster := cluster_num]
  
  return(dt)
}), fill = TRUE)

colname <- c("chrom","start","end","name","score","strand","signalValue","pValue","qValue","summit","cell_type","cluster")
colnames(combined_data) <- colname
# 保存结果
#saveRDS(combined_data, "/data/h01005/hc/GSE147672/原文聚类整合注释.rds")
#fwrite(combined_data, output_file, sep = "\t", compress = "gzip")

# 显示结果摘要
cat("整合完成！输出文件包含以下细胞类型分布：\n")
print(combined_data[, .N, by = .(Cluster, Cell_Type)][order(Cluster)])

significant_peaks <- combined_data %>%
  group_by(cell_type) %>%
  slice_min(order_by = qValue, prop = 0.2) %>%  # 取qValue最大的20%
  ungroup() %>%
  as.data.frame()

# 生成LDSC输入
ldsc_input <- significant_peaks %>%
    dplyr::select(chrom, start, end, name, cell_type)

dir.create("/data/h01005/hc/GSE147672/ldsc/bed/0.2", recursive = TRUE, showWarnings = FALSE)
ldsc_input %>%
  split(.$cell_type) %>%
  walk(~ write_tsv(.x %>% dplyr::select(chrom, start, end, name),
                   paste0("/data/h01005/hc/GSE147672/ldsc/bed/0.2/", make.names(unique(.x$cell_type)), ".bed"),
                   col_names = FALSE))
