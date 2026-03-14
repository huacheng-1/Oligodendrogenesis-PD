# 主要用于做单细胞分析的各种函数
library(Seurat)
# 主要与批次效应有关
library(harmony)
# 广泛使用的数据分析和可视化的工具集
library(tidyverse)
# 是一个功能强大的R绘图包
library(ggplot2)
# 增强绘图和图形组合的包
library(cowplot)
# 用于数据的操作和处理,管道,过滤,排序等
library(dplyr)
# 用于组合多个图形
library(patchwork)
# 用于文件的处理
library(R.utils)
# 使用list.files函数列出文件夹中的所有文件
library(tidydr)
library(glmGamPoi)
library(dplyr)  
library(Seurat)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(clustree)
library(tidydr)
library(celldex)
library(SingleR)
library(tidyr)
library(reshape2)
library(gplots)
library(ggthemes)
library(tibble)
library(pheatmap)
library(networkD3)
library(clusterProfiler)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(data.table)
options(future.globals.maxSize = 5000 * 1024^2)  # 设置为5 GB
data <- fread("D://R//GEO//GSE243639//Filtered_count_table.csv")
data2 <- as.data.frame(data)  
rownames(data2) <- data$V1[1:nrow(data2)]
data2 <- data2[, -1]   
data2 <- data[rowSums(data2) > 0,]
data2 <- data[!duplicated(rownames(data)),]
se <- CreateSeuratObject(counts = data2 )#创建seurat[33537 x 83485]

dir_name <- list.files('D://R//GEO//GSE243639//GSE243639/')
level <- levels(se@meta.data$orig.ident)
filtered_dfs <- list()  
for (i in 1:length(level)){
  filtered_df <- data2 %>%  
    select(starts_with(level[i]))  
  filtered_dfs[[dir_name[i]]] <- filtered_df
}
se_list <- list()
for (i in 1:length(dir_name)){
  se_list[[i]] <- CreateSeuratObject(filtered_dfs[[i]],
                                     project = dir_name[i])  
}

for (i in 1:length(dir_name)){
  sc <- se_list[[i]] 
  sc <- SCTransform(sc)#标准化&高变基因&归一化&PCA
  sc <- FindVariableFeatures(sc, selection.method = "vst", nfeatures = 6000)
  sc <- RunPCA(sc, dims = 1:50)
  se_list[[i]] <- sc
  # 删除sc
  rm(sc)
}

reference_index <- which(sapply(se_list, function(x) "GSM7792123" %in% unique(x$orig.ident)))
combined_anchors <- FindIntegrationAnchors(object.list = se_list, reference = reference_index)
combin_integrated <- IntegrateData(anchorset =combined_anchors, dims = 1:20)
combin_integrated <- ScaleData(combin_integrated)
combin_integrated <- RunPCA(combin_integrated, dims = 1:50)
combin_integrated <- FindNeighbors(combin_integrated, dims = 1:10)#降维聚类umap
combin_integrated <- FindClusters(combin_integrated, resolution = 0.6) 
combin_integrated <- RunUMAP(combin_integrated, dims = 1:10) 
DimPlot(combin_integrated, reduction = "umap")
combin_integrated <- FindClusters(combin_integrated, resolution = 0.6) 
markers_all <- FindAllMarkers(combin_integrated, only.pos = TRUE, min.pct = 0.25,  logfc.threshold = 0.25, verbose = T)
top60 = markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)
g = unique(top60$gene)
tmp = split(top60$gene,top60$cluster)
cat(  unlist(lapply(names(tmp), function(x){
  paste0('cluster',x,':',paste(tmp[[x]],collapse=','))
})),sep = '\n')
new.cluster.ids <- c("Oligodendrocyte","Oligodendrocyte","Oligodendrocyte","Microglial cell",
                     "Astrocyte","Oligodendrocyte","Oligodendrocyte precursor cell","Astrocyte",
                     "Neurons","Oligodendrocyte",
                     "Microglial cell","Astrocyte","Microglial cell",
                     "Astrocyte","Oligodendrocyte","Neurons","Astrocyte",
                     "Vascular cell","Neurons","Oligodendrocyte precursor cell","Astrocyte",
                     "Microglial cell","Microglial cell","Oligodendrocyte precursor cell")
names(new.cluster.ids) <- levels(combin_integrated)
combin_integrated <- RenameIdents(combin_integrated, new.cluster.ids)

active_ident <- combin_integrated@active.ident  
combin_integrated@meta.data$cells <- active_ident 
