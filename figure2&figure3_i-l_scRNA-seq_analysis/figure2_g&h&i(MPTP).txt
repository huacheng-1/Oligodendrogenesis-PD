# 单细胞分析核心包
library(Seurat)      
library(harmony)        
# 数据科学工具链
library(tidyverse)      
# 增强可视化
library(cowplot)        
library(patchwork)      
library(pheatmap)     
library(ggthemes)      
library(SingleR)        
library(celldex)        
library(clusterProfiler) 
library(org.Mm.eg.db)  
library(org.Hs.eg.db)   
library(R.utils)        
library(data.table)     
library(readr)          
library(Matrix)        

options(future.globals.maxSize = 15 * 1024^2)  # 设置为15 GB
dir_name <- list.files('/data/h01005/hc/MPTP/rowdata/')
dir_name
scRNAlist <- list()
for (i in 1:length(dir_name)){
  counts <- Read10X(data.dir = paste('/data/h01005/hc/MPTP/rowdata/',dir_name[i],sep = ''))
  UMI_sums <- rowSums(counts)
  filtered_cells_matrix <- counts[UMI_sums <= 10000, ]
  # 使用createseuratobject来创建seurat对象,使用counts矩阵,设置样本名为目录名
  scRNAlist[[i]] <- CreateSeuratObject(counts = filtered_cells_matrix,
                                       project = dir_name[i],
                                       min.cells = 3,#单个细胞至少表达三个基因
                                       min.features = 200)  #单个基因至少在200个细胞中表达
  
  
}
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
  x <- subset(x,  mt_percent < 10 & nFeature_RNA> 200 & nFeature_RNA< 4000)})

for (i in 1:length(scRNAlist)){
  scRNA <- scRNAlist[[i]]
  scRNA <- SCTransform(scRNA)
  scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 3000)
  scRNA <- RunPCA(scRNA, dims = 1:30)
  scRNAlist[[i]] <- scRNA
}

integration_features <- SelectIntegrationFeatures(object.list = scRNAlist[1], nfeatures = 3000)
combined_anchors <- FindIntegrationAnchors(object.list = scRNAlist[1:6])
combin_integrated <- IntegrateData(anchorset =combined_anchors, dims = 1:20)
combin_integrated <- ScaleData(combin_integrated)
combin_integrated <- RunPCA(combin_integrated, dims = 1:50)
combin_integrated <- JackStraw(combin_integrated)
combin_integrated <- ScoreJackStraw(combin_integrated)
significant_pcs <- 1:30
combin_integrated <- combin_integrated[, , significant_pcs]
combin_integrated <- FindNeighbors(combin_integrated, dims = significant_pcs)#降维聚类umap
combin_integrated <- FindClusters(combin_integrated, resolution = 0.4) 
combin_integrated <- RunUMAP(combin_integrated, dims = 1:10) 
DimPlot(combin_integrated, reduction = "umap")
p1 <- DimPlot(combin_integrated, reduction = "umap", label = TRUE,  
              repel = T,pt.size = 0.8) + 
  theme_dr(xlength = 0.22, ylength = 0.22, arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed"))+
  theme(panel.grid = element_blank());p1

markers_tmp <- FindAllMarkers(combin_integrated, only.pos = TRUE, min.pct = 0.25,  logfc.threshold = 0.25, verbose = T)
top60 = markers_tmp %>% group_by(cluster) %>% top_n(n = 60, wt = avg_log2FC)
g = unique(top60$gene)
tmp = split(top60$gene,top60$cluster)
cat(  unlist(lapply(names(tmp), function(x){
  paste0('cluster',x,':',paste(tmp[[x]],collapse=','))
})),sep = '\n')
#saveRDS(markers,"/data/h01005/hc/MPTP/markers.rds")

#MSCs:间充质干细胞
new.cluster.ids <- c("Oligodendrocyte", "Neuron", "Astrocyte",
                     "Oligodendrocyte","Neuron", "Neuron",
                     "Microglial cell","Oligodendrocyte precursor cell","Neuron","Oligodendrocyte",
                     "Neuron","Neuron","Neuron",
                     "Neuron","Neuron","Astrocyte",
                     "Neuron","MSCs","Neuron","Neuron","Neuron","Astrocyte",
                     "MSCs","Microglial cell","Oligodendrocyte","Neuron","Neuron")
names(new.cluster.ids) <- levels(combin_integrated)
combin_integrated <- RenameIdents(combin_integrated, new.cluster.ids)

p1 <- DimPlot(combin_integrated, reduction = "umap", label = TRUE,  
              repel = T,pt.size = 0.8) + 
  theme_dr(xlength = 0.22, ylength = 0.22, arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed"))+
  theme(panel.grid = element_blank());p1
#ggsave(p1,filename = "/data/h01005/hc/MPTP/20250902更新/umap.pdf",width =7,height = 7)
active_ident <- combin_integrated@active.ident  
combin_integrated@meta.data$cells <- active_ident 


markers <- FindAllMarkers(combin_integrated, only.pos = TRUE, min.pct = 0.25,  logfc.threshold = 0.25, verbose = T)
markers <- na.omit(markers)

#在meta.data中创建一个名为group的空集
combin_integrated@meta.data$group <- NA 
con_samples <- c("blank_1","blank_2","blank_3")  
pd_samples <- c("model_1","model_2","model_3")  
combin_integrated@meta.data$group[combin_integrated@meta.data$orig.ident %in% pd_samples] <- "pd"  
combin_integrated@meta.data$group[combin_integrated@meta.data$orig.ident %in% con_samples] <- "con"  


meta_data <- combin_integrated@meta.data  
cell_types_to_compare <- c("Oligodendrocyte precursor cell",
                           "Astrocyte", "Microglial cell","Neuron",
                           "Oligodendrocyte","MSCs"
                           )  
filtered_meta <- meta_data[meta_data$cells %in% cell_types_to_compare, ]  
filtered_meta$group <- as.factor(filtered_meta$group)  
filtered_meta$cells <- as.factor(filtered_meta$cells) 
cell_counts <- filtered_meta %>%  
  group_by(group, cells) %>%  
  summarise(count = n())
bar_per <- cell_counts %>% 
  group_by(group) %>%
  mutate(sum(count)) %>%
  mutate(percent = count / `sum(count)`)

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
  scale_fill_manual(values = c("Oligodendrocyte" = "red", "Astrocyte" = "#F78000", "Oligodendrocyte precursor cell" = "yellow",
                               "Neuron" ="pink" ,"Microglial cell" = "#3176B7","MSCs" = "brown"
  )) +
  theme_few()+
  theme(plot.title = element_text(size=12,hjust=0.5))
p1
ggsave(p1,filename = "/data/h01005/hc/MPTP/20250902更新/柱状堆积图(竖).pdf",width =7,height = 10)


top10 = markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
genes <- top10$gene
top10_filt <- top10[top10$gene %in% genes, ]
top10_filt <- c("Opalin","Mal","Gm10754","Gm30382","Aqp4","Slc39a12","Dock8","P2ry12","Vcan","Vxn","Lama1","Cped1")
#气泡图
top1_dotplot <- DotPlot(combin_integrated, features = top10_filt)+
  theme(axis.text.x = element_text(size = 6 ,angle = 45, vjust = 1, hjust = 1))
top1_dotplot
#ODC GENE: Mog,Plekhh1,Mobp


###热图
top1 = markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top1 <- c("Opalin","Mal","Gm10754","Gm30382","Aqp4","Slc39a12","Dock8","P2ry12","Vcan","Vxn","Lama1","Cped1")
#write.table(top1[,6:7],"D:/R/GEO/CNP_RE/new/top20gene.txt", row.names = F)
top1_dotplot <- DotPlot(combin_integrated, features = top1)+
  theme(axis.text.x = element_text(size = 6 ,angle = 45, vjust = 1, hjust = 1))
top1_dotplot


dotplot_data <- top1_dotplot$data
class(dotplot_data)
heatmap_data <- dotplot_data %>%
  dplyr::select(features.plot, id, avg.exp.scaled) %>%
  pivot_wider(names_from = features.plot, values_from = avg.exp.scaled)


heatmap_data = column_to_rownames(heatmap_data,var = "id")


heatmap <- pheatmap(heatmap_data,
                    cellwidth = 10,
                    cellheight = 15,
                    cluster_rows = F, 
                    cluster_cols = F,
                    
                    legend_breaks = c( -2, -1, 0, 1, 2),
                    color = colorRampPalette(c("#d3eafa", "white", "#e85b3f"))(100))

# 转置矩阵并设置行列名
heatmap_data_transposed <- t(heatmap_data[, ]) # 去除id列后转置
colnames(heatmap_data_transposed) <- row.names(heatmap_data) # 原id列作为列名

# 绘制旋转后的热图
heatmap <- pheatmap(heatmap_data_transposed,
                    cellwidth = 15,  # 调整单元格宽度
                    cellheight = 10, # 调整单元格高度
                    cluster_rows = F,
                    cluster_cols = F,
                    angle_col = 45,  # 调整列标签角度
                    legend_breaks = c(-2, -1, 0, 1, 2),
                    color = colorRampPalette(c("#d3eafa", "white", "#e85b3f"))(100))

# 调整保存尺寸（宽高互换）
ggsave(heatmap,
       filename = "/data/h01005/hc/MPTP/20250902更新/heatmap_rotated.pdf",
       width = 9,      # 原高度
       height = 22)    # 原宽度

write.csv(as.data.frame(heatmap_data_transposed),"/data/h01005/hc/MPTP/20250902更新/heatmap_data_transposed.csv", row.names = T)



