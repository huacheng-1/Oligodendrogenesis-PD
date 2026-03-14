library(Seurat)
library(ggplot2)
library(tidyverse)
library(cowplot)
library(ggthemes)
# 用于数据的操作和处理,管道,过滤,排序等
library(dplyr)
# Load data
scRNA_combined <- readRDS("/data/h01005/hc/再整合/scRNA_combined.RDS")

# Check if SCTransform was used
DefaultAssay(scRNA_combined)  # Should be "SCT" if SCTransform was run

# Extract Oligodendrocyte subset
oligo_obj <- subset(oligo_obj, subset = cell == c("Oligodendrocyte"))
oligo_conpd <- subset(oligo_obj, subset = group == c("con","pd"))

# Re-run SCTransform (if needed) or ensure 'SCT' assay is present
oligo_obj <- SCTransform(oligo_conpd, verbose = TRUE)

# Set default assay to SCT
DefaultAssay(oligo_obj) <- "SCT"

# Dimensionality reduction & clustering
oligo_obj <- RunPCA(oligo_obj, verbose = TRUE)
ElbowPlot(oligo_obj) 
oligo_obj <- FindNeighbors(oligo_obj, dims = 1:10, verbose = T)
oligo_obj <- FindClusters(oligo_obj, resolution = 0.3, verbose = T)#1.0
oligo_obj <- RunUMAP(oligo_obj, dims = 1:10, verbose = T)
#########################3
oligo_obj <- oligo_conpd
# Visualize clusters
DimPlot(
  oligo_obj, 
  reduction = "umap", 
  group.by = "seurat_clusters", 
  label = FALSE, 
  pt.size = 0.5,
  cols = c("0" = "#f57c6e", "1" = "#f2b56f", "2" = "#fae69e", 
           "3" = "#84c3b7", "7" = "#88d8db", "4" = "#71b7ed", 
           "5" = "#b8aeeb", "6" = "#f2a7da")  # 如果存在聚类7
) + 
  ggtitle("Oligodendrocyte Subclusters")
oligo_obj <- oligo_conpd
meta_data <- oligo_obj@meta.data  
cell_types_to_compare <- c("Oligodendrocyte") 

filtered_meta <- meta_data[meta_data$cells %in% cell_types_to_compare, ]  
filtered_meta$group <- as.factor(filtered_meta$group)  
filtered_meta$cells <- as.factor(filtered_meta$cells)  

cell_counts <- filtered_meta %>%  
  group_by(group, seurat_clusters) %>%  
  summarise(count = n())
bar_per <- cell_counts %>% 
  group_by(group) %>%
  mutate(sum(count)) %>%
  mutate(percent = count / `sum(count)`)

##只留pd和con
# 在绘图代码前添加过滤步骤
bar_per <- bar_per %>% 
  filter(group %in% c("con", "pd"))  # 过滤保留 con 和 pd
#write_csv(bar_per,"/data/h01005/hc/再整合/20250902更新/堆积柱状图.csv")
###堆积柱状图
p1 <- ggplot(bar_per, aes(x = group, y = percent, fill = cells)) +  
  geom_bar(stat = "identity", position = "stack", width = 0.7) +  
  #geom_text(aes(label = count), position = position_stack(vjust = 0.5)) + # 在堆积条的中间显示数字
  theme(axis.ticks = element_line(linetype = "blank"),
        legend.position = "top",
        panel.grid.minor = element_line(colour = NA,linetype = "blank"), 
        panel.background = element_rect(fill = NA),
        plot.background = element_rect(colour = NA)) +
  labs(title = "Bar Chart",  
       x = "group",  
       y = "radio",  
       fill = "cell name") +  
  theme_minimal() + # 使用简洁的主题  
  scale_fill_manual(values = c("0" = "#f57c6e", "1" = "#f2b56f", "2" = "#fae69e", 
                               "3" = "#84c3b7", "7" = "#88d8db", "4" = "#71b7ed", 
                               "5" = "#b8aeeb", "6" = "#f2a7da"
  )) +
  theme_few()+
  theme(plot.title = element_text(size=12,hjust=0.5))
p1


# Prepare for marker detection (SCT-specific)
oligo_obj <- PrepSCTFindMarkers(oligo_obj, assay = "SCT")

# Find markers
markers <- FindAllMarkers(
  oligo_conpd,
  assay = "SCT",
  slot = "data",  # or "scale.data"
  group.by = "group",
  only.pos = F,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  verbose = TRUE
)
saveRDS(markers, "/data/h01005/hc/亚聚类/markersconpd_01_all.RDS")
# Save results
saveRDS(oligo_obj, "/data/h01005/hc/亚聚类/oligo_conpd.RDS")
#saveRDS(markers, "/data/h01005/hc/亚聚类/markers.RDS")
#write.csv(markers, "/data/h01005/hc/亚聚类/oligo_markers.csv")




