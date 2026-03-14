# 主要用于做单细胞分析的各种函数
library(Seurat)  
# 主要与批次效应有关
library(harmony)  
# 广泛使用的数据分析和可视化的工具集
library(tidyverse)  # 包含ggplot2/dplyr/tidyr等核心包
# 增强绘图和图形组合的包
library(cowplot)  
# 用于组合多个图形
library(patchwork)  
# 用于文件的处理
library(R.utils)  
# 用于细胞类型注释
library(celldex)  
library(SingleR)  
# 用于热图绘制
library(pheatmap)  
# 用于富集分析
library(clusterProfiler)  
library(org.Mm.eg.db)  # 小鼠基因注释
library(org.Hs.eg.db)  # 人类基因注释
organism = "org.Hs.eg.db"
options(future.globals.maxSize = 5000 * 1024^2)  # 设置为5 GB
dir_name <- list.files('D:/R/GEO/CNP0000892/CNS/mid')
dir_name

# 创建一个空的list对象来存储转化之后的seurat对象
scRNAlist <- list()
for (i in 1:length(dir_name)){
  counts <- Read10X(data.dir = paste("D:/R/GEO/CNP0000892/CNS/mid/",dir_name[i],"/filtered_feature_bc_matrix",sep = ''))
  scRNAlist[[i]] <- CreateSeuratObject(counts = counts,
                                       project = dir_name[i],
                                       min.cells = 3,#单个细胞至少表达三个基因
                                       min.features = 200)  #单个基因至少在200个细胞中表达
  
  
}
scRNAlist # 查看样本信息

for (i in 1:length(scRNAlist)){
  sc <- scRNAlist[[i]] # 获取scRNAlist中的第i个seurat对象
  # 计算线粒体的比例
  sc[['mt_percent']] <- PercentageFeatureSet(sc,pattern = '^mt-')
  # 假设sc是你的Seurat对象  
  # 计算每个细胞的转录本数量  
  transcript_counts <- Matrix::colSums(sc@assays$RNA$counts)  
  # 设置阈值  
  low_threshold <- 500  
  high_threshold <- 6000  
  # 找出要保留的细胞的索引  
  keep_cells <- (transcript_counts >= low_threshold) & (transcript_counts <= high_threshold)  
  # 过滤掉不符合条件的细胞  
  sc <- subset(sc, cells = names(sc@assays$RNA$counts)[keep_cells])
  scRNAlist[[i]] <- sc
  # 删除sc
  rm(sc)
}

scRNAlist <- lapply(scRNAlist, FUN = function(x){
  x <- subset(x,  mt_percent < 2)})


for (i in 1:length(dir_name)){
  sc <- scRNAlist[[i]] 
  sc <- SCTransform(sc)#标准化&高变基因&归一化&PCA
  sc <- FindVariableFeatures(sc, selection.method = "vst", nfeatures = 6000)
  sc <- RunPCA(sc, dims = 1:50)
  scRNAlist[[i]] <- sc
  # 删除sc
  rm(sc)
}
scrnalist2 <- scRNAlist

combined_anchors <- FindIntegrationAnchors(object.list = scRNAlist)
combin_integrated <- IntegrateData(anchorset =combined_anchors, dims = 1:20)
combin_integrated <- ScaleData(combin_integrated)
combin_integrated <- RunPCA(combin_integrated, dims = 1:20)
combin_integrated <- FindNeighbors(combin_integrated, dims = 1:20)#降维聚类umap
combin_integrated <- FindClusters(combin_integrated, resolution = 1) 
combin_integrated <- RunUMAP(combin_integrated, dims = 1:10, n.neighbors = 600, reduction = "pca") 
DimPlot(combin_integrated, reduction = "umap")
p1 <- DimPlot(combin_integrated, reduction = "umap", label = TRUE,  
              repel = T,pt.size = 0.8) + 
  theme_dr(xlength = 0.22, ylength = 0.22, arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed"))+
  theme(panel.grid = element_blank());p1


combin_integrated@meta.data$group <- NA 
CON_samples <- paste0("CNS01957", 44)  
PD_samples <- paste0("CNS01957", 47)  
combin_integrated@meta.data$group[combin_integrated@meta.data$orig.ident %in% PD_samples] <- "pd"  
combin_integrated@meta.data$group[combin_integrated@meta.data$orig.ident %in% CON_samples] <- "con"  

markers <- FindAllMarkers(combin_integrated, only.pos = TRUE, min.pct = 0.25,  logfc.threshold = 0.25, verbose = T)
top60 = markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)
g = unique(top60$gene)
tmp = split(top60$gene,top60$cluster)
cat(  unlist(lapply(names(tmp), function(x){
  paste0('cluster',x,':',paste(tmp[[x]],collapse=','))
})),sep = '\n')


new.cluster.ids <- c("Glutamatergic Neuron","Oligodendrocyte","GABAergic Neuron",
                     "GABAergic Neuron","Astrocyte","Glutamatergic Neuron","Oligodendrocyte",
                     "Microglial cell","Oligodendrocyte","GABAergic Neuron","Glutamatergic Neuron",
                     "Oligodendrocyte","Glutamatergic Neuron","GABAergic Neuron",
                     "Glutamatergic Neuron","Oligodendrocyte precursor cell","Dopaminergic Neuron","Oligodendrocyte precursor cell",
                     "Glutamatergic Neuron","Glutamatergic Neuron","Vascular cell","GABAergic Neuron",
                     "AP","GABAergic Neuron","Oligodendrocyte","Pericyte","Glutamatergic Neuron","Glutamatergic Neuron",
                     "Astrocyte","Serotonergic Neuron")
names(new.cluster.ids) <- levels(combin_integrated)
combin_integrated <- RenameIdents(combin_integrated, new.cluster.ids)
p1 <- DimPlot(combin_integrated, reduction = "umap", label = TRUE,  
              repel = T,pt.size = 0.8) + 
  theme_dr(xlength = 0.22, ylength = 0.22, arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed"))+
  theme(panel.grid = element_blank());p1


#将细胞注释与meta.data匹配
active_ident <- combin_integrated@active.ident  
combin_integrated@meta.data$cells <- active_ident 
meta_data <- combin_integrated@meta.data  
cell_types_to_compare <- c("Oligodendrocyte", "Oligodendrocyte precursor cell",
                           "Glutamatergic Neuron", "Astrocyte","Pericyte",
                           "GABAergic Neuron","Microglial cell","Serotonergic Neuron",
                           "Vascular cell","AP","Dopaminergic Neuron")  
filtered_meta <- meta_data[meta_data$cells %in% cell_types_to_compare, ]  
filtered_meta$group <- as.factor(filtered_meta$group)  
filtered_meta$cells <- as.factor(filtered_meta$cells)  
markers_all <- markers
markers_all$new_cluster <- new.cluster.ids[as.numeric(markers_all$cluster)]
markers <- markers_all
markers_re <- FindAllMarkers(combin_integrated, only.pos = TRUE, min.pct = 0.25,  logfc.threshold = 0.25, verbose = T)
cell_counts <- filtered_meta %>%  
  group_by(group, cells) %>%  
  summarise(count = n())
bar_per <- cell_counts %>% 
  group_by(group) %>%
  mutate(sum(count)) %>%
  mutate(percent = count / `sum(count)`)

radio <- bar_per$`sum(count)`[12]/(bar_per$`sum(count)`[1]+bar_per$`sum(count)`[12])

###堆积柱状图
p2 <- ggplot(bar_per, aes(x = group, y = percent, fill = cells)) +  
  geom_bar(stat = "identity", position = "stack", width = 0.7) +  
  #geom_text(aes(label = count), position = position_stack(vjust = 0.5)) + # 在堆积条的中间显示数字
  theme(axis.ticks = element_line(linetype = "blank"),
        legend.position = "top",
        panel.grid.minor = element_line(colour = NA,linetype = "blank"), 
        panel.background = element_rect(fill = NA),
        plot.background = element_rect(colour = NA)) +
  labs(title = "Stacked Bar Chart",  
       x = "group",  
       y = "radio",  
       fill = "cell name") +  
  theme_minimal() + # 使用简洁的主题  
  scale_fill_manual(values = c("Oligodendrocyte" = "red", "Astrocyte" = "#F78000", "Oligodendrocyte precursor cell" = "yellow", "Glutamatergic Neuron" ="green" ,
                               "Pericyte" = "#41BED1","GABAergic Neuron" = "pink", "Microglial cell" = "#3176B7", "Serotonergic Neuron" ="#e85b3f",
                               "Vascular cell" = "#d3eafa", "AP" = "brown", "Dopaminergic Neuron" ="purple")) +
  theme_few()+
  theme(plot.title = element_text(size=12,hjust=0.5))
p2
ggsave(p2,filename = "D:/R/GEO/CNP_RE/大聚类柱状堆积图2.pdf",width =7,height = 14)


#气泡图
top1 = markers_all %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.table(top1[,6:7],"D:/R/GEO/CNP_RE/new/top20gene.txt", row.names = F)
top1_dotplot <- DotPlot(combin_integrated, features = unique(top1$gene))+
  theme(axis.text.x = element_text(size = 6 ,angle = 45, vjust = 1, hjust = 1))
top1_dotplot


##气泡图联动热图
dotplot_data <- top1_dotplot$data
class(dotplot_data)
heatmap_data <- dotplot_data %>%
  dplyr::select(features.plot, id, avg.exp.scaled) %>%
  pivot_wider(names_from = features.plot, values_from = avg.exp.scaled)

#将细胞分群ID列设置为行名
heatmap_data = column_to_rownames(heatmap_data,var = "id")

#绘制热图
pheatmap(heatmap_data,
         cellwidth = 10,
         cellheight = 15,
         cluster_rows = F, 
         cluster_cols = F,
         legend_breaks = c( -2, -1, 0, 1, 2),
         color = colorRampPalette(c("#d3eafa", "white", "#e85b3f"))(100))