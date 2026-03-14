第一步：环境准备与数据加载
options(future.globals.maxSize = 300000 * 1024^2)  # 设置系统最大调用内存为300 GB
#加载r包
# 主要用于做单细胞分析的各种函数
library(Seurat)
# 主要与批次效应有关
library(harmony)
# 广泛使用的数据分析和可视化的工具集
library(tidyverse)
# 用于数据的操作和处理,管道,过滤,排序等
library(dplyr)
# 用于文件的处理
library(R.utils)
# 使用list.files函数列出文件夹中的所有文件
library(tidydr)
library(glmGamPoi)
library(patchwork)
library(clustree)
library(SingleR)
library(tidyr)
library(networkD3)
library(clusterProfiler)
library(future)


第二步：读取需要整合的geo数据的原始count信息
#从头开始，com3965是一组已经把gse243639与gse178265整合完的seurat数据集
se_list <- list()
filtered_dfs <- list()  #预先创建之后会使用的列

#提取com3965的原始count信息
#将count信息根据orig.ident分类并输入filtered_dfs中
for (y in 1:length(com3965)){
  level <- levels(com3965[[y]]@meta.data$orig.ident)
  dir_name <- level
  data2 <- as.data.frame(com3965[[y]]@assays$RNA$counts)
  for (i in 1:length(level)){
    filtered_df <- data2 %>%  
      dplyr::select(starts_with(level[i]))  
    filtered_dfs[[dir_name[i]]] <- filtered_df  
    gc()
  }
  gc()
}
#根据orig.ident的count信息重新生成seurat信息，并输入se_list列中
#生成se_list列
se_list <- list()
for (i in 1:length(dir_name)){
  se_list[[i]] <- CreateSeuratObject(filtered_dfs[[i]],
                                     project = dir_name[i])  
}

#根据se_list中每一个样本进行循环
for (i in 1:length(se_list)){
  sc <- se_list[[i]] # 获取scRNAlist中的第i个seurat对象
  # 计算线粒体的比例
  sc[['mt_percent']] <- PercentageFeatureSet(sc,pattern = '^MT-')
  transcript_counts <- Matrix::colSums(sc@assays$RNA$counts)  
  # 设置阈值  
  low_threshold <- 200  
  high_threshold <- 6000  
  # 找出要保留的细胞的索引  
  keep_cells <- (transcript_counts >= low_threshold) & (transcript_counts <= high_threshold)  
  # 过滤掉不符合条件的细胞  
  sc <- subset(sc, cells = names(sc@assays$RNA$counts)[keep_cells])
  sc <- subset(sc, mt_percent < 10)
  sc <- NormalizeData(sc)
  gc()
  sc <- FindVariableFeatures(sc,nFeature_RNA = 4000) # 这里指定最后想要得到的基因数量
  gc()
  sc <- ScaleData(sc,vars.to.regress = c('mt_percent')) # 这里去除线粒体的影响
  sc <- RunPCA(sc) #运行降维
  gc()
  se_list[[i]] <- sc
  # 删除sc
  rm(sc)
  gc()
}
#将se_list列中的所有seurat数据集整合到一起
scRNA_combined <- merge(x = se_list[[1]],y = se_list[-1],add.cell.ids = names(se_list))

使用SCTransform函数进行降维，聚类和寻找高变基因分析
scRNA_combined2 <- SCTransform(
  scRNA_combined,
  vars.to.regress = c("mt_percent"), # 回归线粒体比例
  variable.features.n = 4000,   # 选择4000个高变基因
  verbose = FALSE
)#标准化&高变基因&归一化&PCA
#
gc()
scRNA_combined <- RunPCA(scRNA_combined2, dims = 1:30) # 计算前30个主成分
gc()
scRNA_combined <- FindNeighbors(scRNA_combined, dims = 1:10)#降维聚类umap
gc()
scRNA_combined <- FindClusters(scRNA_combined, resolution = 0.8) 
gc()
scRNA_combined <- RunUMAP(scRNA_combined, dims = 1:10, min.dist = 0.5, n.neighbors = 30) 
gc()

#加载GSE193688的seurat数据
load("~/hc/再整合/GSE193688-HOMO/注释前.RData")

#获取gse193688与scRNA_combined的共有基因
common_genes <- intersect(rownames(gse193688), rownames(scRNA_combined))

#筛选有共有基因的gse193688和scRNA_combined
gse193688_filtered <- gse193688[common_genes, ]
scRNA_combined_filtered <- scRNA_combined[common_genes, ]
gc()
rm(list = ls()[!ls() %in% c("gse193688_filtered","scRNA_combined_filtered")])

#重新对scRNA_combined和gse193688进行SCT，为后面整合做准备
scRNA__filtered_repro <- SCTransform(
  scRNA_combined_filtered,
  vars.to.regress = "mt_percent",
  variable.features.n = 3000,
  verbose = FALSE
)
gc()
gse193688_filtered_repro <- SCTransform(
  gse193688_filtered, 
  vars.to.regress = "mt_percent",  # 统一回归变量
  variable.features.n = 3000,      # 统一高变基因数量
  verbose = FALSE                  # 统一输出设置
)
gc()

#整合gse193688与scRNA_combined
combined <- merge(
  x = gse193688_filtered_repro,
  y = scRNA_filtered_repro,
  add.cell.ids = c("GSE193688", "scRNA_Combined"),  # 细胞名前缀
  project = "Integrated_Project"                     # 必须指定项目名
)
saveRDS(combined,"~/hc/combined.RDS")

# 运行Harmony（自动识别批次变量orig.ident）
# 使用更新后的对象进行合并
# 设置默认assay为SCT（因使用SCTransform）
DefaultAssay(combined) <- "SCT"

# 1. 访问模型特征属性 --------------------------------------------------------
# 查看可用的特征属性
names(combined@assays$SCT@SCTModel.list$model1@feature.attributes)

# 2. 提取高变基因索引 --------------------------------------------------------
# 根据方差选择高变基因（模拟FindVariableFeatures的vst方法）
var_genes_idx <- order(
  combined@assays$SCT@SCTModel.list$model1@feature.attributes$variance,
  decreasing = TRUE
)[1:3000]  # 提取方差最大的前3000个基因

# 3. 获取基因名称 -----------------------------------------------------------
var_genes <- rownames(combined)[var_genes_idx]
VariableFeatures(combined) <- var_genes
length(VariableFeatures(combined))  # 应返回3000
head(VariableFeatures(combined))    # 应显示基因名称


combined <- RunPCA(
  combined,
  assay = "SCT",
  npcs = 50,
  verbose = FALSE
)

# 运行Harmony（自动识别批次变量orig.ident）
combined <- combined %>%
  RunHarmony(
    group.by.vars = "orig.ident",  # 批次变量（默认使用合并时的orig.ident）
    dims.use = 1:30,               # 使用前30个PC
    theta = 2,                     # 聚类多样性参数（默认）
    lambda = 0.5,                  # 校正强度（可调）
    plot_convergence = TRUE        # 显示收敛曲线
  )



# 运行UMAP
scRNA_combined <- RunUMAP(
  combined,
  reduction = "harmony",
  dims = 1:30,
  n.neighbors = 30,
  min.dist = 0.3
)

#生成未注释的聚类umap图
p1 <- DimPlot(scRNA_combined, reduction = "umap", label = TRUE,  
              repel = T,pt.size = 0.8) + 
  theme_dr(xlength = 0.22, ylength = 0.22, arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed"))+
  theme(panel.grid = element_blank());p1

#寻找每个聚类具有显著差异的基因
scRNA_combined <- PrepSCTFindMarkers(scRNA_combined)
markers <- FindAllMarkers(scRNA_combined, only.pos = TRUE, min.pct = 0.25,  logfc.threshold = 0.25, verbose = FALSE)

# 根据avg_log2FC由大到小排序生成每个聚类前300的基因
top300 = markers %>% group_by(cluster) %>% top_n(n = 300, wt = avg_log2FC)
g = unique(top300$gene)

tmp = split(top300$gene,top300$cluster)
cat(  unlist(lapply(names(tmp), function(x){
  paste0('cluster',x,':',paste(tmp[[x]],collapse=','))
})),sep = '\n')
#打印每个聚类的信息，在ACT基因注释网站手动注释细胞类型。

#聚类信息，储存在new.cluster.ids中
new.cluster.ids <-
 c("Oligodendrocyte","Oligodendrocyte","Microglial cell","Oligodendrocyte","Neuron","Neuron","Astrocyte","Oligodendrocyte","Oligodendrocyte","Astrocyte","Oligodendrocyte precursor cell",
"Oligodendrocyte","Fibroblast","Oligodendrocyte","Microglial cell",
"Neuron","Neuron","Microglial cell", "Oligodendrocyte","Oligodendrocyte",
"Oligodendrocyte","Neuron","Oligodendrocyte","Astrocyte",         "Neuron","Neuron","Neuron","Astrocyte","Oligodendrocyte","Neuron",         "Oligodendrocyte","Neuron","Neuron","Neuron","Neuron","Neuron","Oligodendrocyte","Neuron","Fibroblast","Neuron","Neuron","Astrocyte","Oligodendrocyte","Neuron","Astrocyte","Oligodendrocyte precursor cell",
"Neuron","Oligodendrocyte", "Neuron","Astrocyte","Neuron","Neuron",
"Oligodendrocyte precursor cell","Oligodendrocyte",
"Oligodendrocyte","Neuron","Oligodendrocyte","Oligodendrocyte","Astrocyte","Oligodendrocyte precursor cell","Oligodendrocyte precursor cell","Oligodendrocyte")

#将聚类信息返回到scRNA_combined中
names(new.cluster.ids) <- levels(scRNA_combined)
scRNA_combined <- RenameIdents(scRNA_combined, new.cluster.ids)

#生成聚类注释后的umap图
p1 <- DimPlot(scRNA_combined, reduction = "umap", label = TRUE,  
              repel = T,pt.size = 0.8) + 
  theme_dr(xlength = 0.22, ylength = 0.22, arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed"))+
  theme(panel.grid = element_blank());p1
#将ident信息返回到scRNA_combined中
active_ident <- scRNA_combined@active.ident  
#将细胞注释与meta.data匹配
scRNA_combined@meta.data$cells <- active_ident 

#添加分组信息
scRNA_combined@meta.data$group <- NA 
pd_samples <- c("s.0096","s.0098","s.0100","s.0098","s.0100","s.0102","s.0103",
"s.0104","s.0105","s.0107","s.0109","s.0110","s.0111","s.0116","s.0118","s.0119","s.0120","pPDsHSrSNxi3873d200429PosA","pPDsHSrSNxi3873d200429PosE","pPDsHSrSNxi3873d200429PosF","pPDsHSrSNxi3887d200429PosB","pPDsHSrSNxi3873d200429PosD","pPDsHSrSNxi4560d200429PosD","pPDsHSrSNxi3873d200429DAPIA","pPDsHSrSNxi3887d200429DAPIB","pPDsHSrSNxi3873d200429PosB","pPDsHSrSNxi3873d200429DAPIB","pPDsHSrSNxi4560d200429PosE","pPDsHSrSNxi3887d200429PosA","pPDsHSrSNxi4560d200429DAPIB","pPDsHSrSNxi3873d200429PosC","pPDsHSrSNxi4560d200429DAPIA","pPDsHSrSNxi4560d200429PosB","pPDsHSrSNxi4568d200429DAPIB","pPDsHSrSNxi2142d200429PosA","pPDsHSrSNxi4560d200429PosC","pPDsHSrSNxi4568d200429PosB","pPDsHSrSNxi4560d200429PosA","pPDsHSrSNxi1963d200429DAPIB","pPDsHSrSNxi3887d200429DAPIA","pPDsHSrSNxi4568d200429PosA","pPDsHSrSNxi4568d200429DAPIA","pPDsHSrSNxi1963d200429DAPIC","pPDsHSrSNxi1963d200429DAPIA","pPDsHSrSNxi2142d200429DAPIA","pPDsHSrSNxi2142d200429DAPIB","pPDsHSrSNxi1963d200429PosB","pPDsHSrSNxi1963d200429PosA","pPDsHSrSNxi1963d200429PosC","pPDsHSrSNxi1963d200429PosD","pPDsHSrSNxi1963d200429PosE" ,"GSM5818658_P02284_filtered_feature_bc_matrix.h5","GSM5818659_P4205_filtered_feature_bc_matrix.h5","GSM5818660_P4560_filtered_feature_bc_matrix.h5" ,"GSM5818661_P4653_filtered_feature_bc_matrix.h5","GSM5818662_P4772_filtered_feature_bc_matrix.h5" ,"GSM5818663_P4884_filtered_feature_bc_matrix.h5" ,"GSM5818664_P4919_filtered_feature_bc_matrix.h5" ,"GSM5818665_P5318_filtered_feature_bc_matrix.h5" ,"GSM5818666_P5331_filtered_feature_bc_matrix.h5" ,"GSM5818667_P5626_filtered_feature_bc_matrix.h5"  ,"GSM5818668_P5662_filtered_feature_bc_matrix.h5" ,"GSM5818669_P6259_filtered_feature_bc_matrix.h5" ,"GSM5818670_P6320_filtered_feature_bc_matrix.h5"  ,"GSM5818671_P6323_filtered_feature_bc_matrix.h5" ,"GSM5818672_P6326_filtered_feature_bc_matrix.h5" )  

con_samples <- c("s.0127","s.0128","s.0129","s.0130","s.0131","s.0142","s.0147",
"s.0151","s.0152","s.0153","s.0154","s.0158","s.0159","s.0165","pPDCN4340DAPIA030419","pPDCN4340DAPIB030419","pPDsHSrSNxi3482d200429PosB","pPDsHSrSNxi3346d200429PosB","pPDsHSrSNxi3482d200429PosA","pPDsHSrSNxi3346d200429PosC",  "pPDsHSrSNxi5610d200429Pos","pPDsHSrSNxi3322d200429PosG","pPDsHSrSNxi3322d200429PosE","pPDsHSrSNxi3346d200429PosA","pPDCN3839DAPIB030419", "pPDsHSrSNxi3345d200429DAPIB","pPDsHSrSNxi3482d200429DAPIB", "pPDsHSrSNxi4956d200429PosA","pPDsHSrSNxi3482d200429DAPIA","pPDsHSrSNxi3322d200429PosF","pPDsHSrSNxi6173d200429DAPIA","pPDsHSrSNxi3345d200429DAPIA","pPDsHSrSNxi3346d200429DAPIB","pPDsHSrSNxi4956d200429DAPIB", "pPDsHSrSNxi3322d200429PosC","pPDsHSrSNxi3322d200429DAPIA","pPDsHSrSNxi4956d200429DAPIA","pPDsHSrSNxi3345d200429DAPIC","pPDsHSrSNxi3322d200429PosD","pPDCN3839DAPIA030419","pPDsHSrSNxi3345d200429PosA","pPDsHSrSNxi5610d200429Neg","pPDsHSrSNxi3322d200429DAPIB","pPDCN5730NeuN22119","pPDsHSrSNxi4956d200429PosB","pPDsHSrSNxi6173d200429PosB",pPDsHSrSNxi3322d200429PosA","pPDsHSrSNxi3298d200429DAPIC","pPDsHSrSNxi3322d200429PosB","pPDsHSrSNxi3322d200429PosH","pPDsHSrSNxi3322d200429DAPIC","pPDsHSrSNxi6173d200429PosA","pPDsHSrSNxi3298d200429DAPIA","pPDCN5730DAPI22119","pPDsHSrSNxi3298d200429DAPIB","pPDCN3898DAPIB030419","pPDCN3898DAPIA030419","pPDsHSrSNxi3298d200429PosA","pPDsHSrSNxi3298d200429PosB","GSM5818673_YC4345_filtered_feature_bc_matrix.h5","GSM5818674_YC5276_filtered_feature_bc_matrix.h5","GSM5818675_YC5346_filtered_feature_bc_matrix.h5","GSM5818676_YC5442_filtered_feature_bc_matrix.h5","GSM5818677_YC5626_filtered_feature_bc_matrix.h5","GSM5818678_YC5751_filtered_feature_bc_matrix.h5" ,"GSM5818679_YC5813_filtered_feature_bc_matrix.h5","GSM5818680_YC5928_filtered_feature_bc_matrix.h5","GSM5818681_YC6235_filtered_feature_bc_matrix.h5" ,"GSM5818682_YC914_filtered_feature_bc_matrix.h5")  

lbd_samples <- c("pPDsHSrSNxi4775d200429PosA","pPDsHSrSNxi2544d200429DAPIB",
"pPDsHSrSNxi4775d200429DAPIB","pPDsHSrSNxi4775d200429PosB","pPDsHSrSNxi4775d200429DAPIC","pPDsHSrSNxi4775d200429DAPIA","pPDsHSrSNxi2561d200429PosC","pPDsHSrSNxi2544d200429PosB","pPDsHSrSNxi2569d200429PosA","pPDsHSrSNxi2544d200429DAPIA","pPDsHSrSNxi2569d200429DAPIA","pPDsHSrSNxi2561d200429PosA","pPDsHSrSNxi2561d200429DAPIB","pPDsHSrSNxi2544d200429PosA","pPDsHSrSNxi2544d200429PosC","pPDsHSrSNxi2561d200429PosB","pPDsHSrSNxi2561d200429DAPIA","pPDsHSrSNxi2569d200429DAPIB")

old_samples <- c("GSM5818649_C16007_filtered_feature_bc_matrix.h5"
                 ,"GSM5818650_C4176_filtered_feature_bc_matrix.h5" 
                 ,"GSM5818651_C4322_filtered_feature_bc_matrix.h5" 
                 ,"GSM5818652_C4361_filtered_feature_bc_matrix.h5" 
                 ,"GSM5818653_C4634_filtered_feature_bc_matrix.h5" 
                 ,"GSM5818654_C4660_filtered_feature_bc_matrix.h5" 
                 ,"GSM5818655_C4691_filtered_feature_bc_matrix.h5" 
                 ,"GSM5818656_C4814_filtered_feature_bc_matrix.h5" 
                 ,"GSM5818657_C5410_filtered_feature_bc_matrix.h5" )

scRNA_combined@meta.data$group[scRNA_combined@meta.data$orig.ident %in% pd_samples] <- "pd"  
scRNA_combined@meta.data$group[scRNA_combined@meta.data$orig.ident %in% con_samples] <- "con"  
scRNA_combined@meta.data$group[scRNA_combined@meta.data$orig.ident %in% lbd_samples] <- "lbd"
scRNA_combined@meta.data$group[scRNA_combined@meta.data$orig.ident %in% old_samples] <- "old"
meta_data <- scRNA_combined@meta.data  
cell_types_to_compare <- c("Oligodendrocyte", "Oligodendrocyte precursor cell",
                           "Astrocyte", "Microglial cell","Neuron", "Fibroblast") 

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

##只留pd和con
bar_per <- bar_per %>% 
  filter(group %in% c("con", "pd"))  # 过滤保留 con 和 pd

##堆积柱状图
p2 <- ggplot(bar_per, aes(x = group, y = percent, fill = cells)) +  
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
  scale_fill_manual(values = c("Oligodendrocyte" = "red", "Astrocyte" = "#F78000", "Oligodendrocyte precursor cell" = "yellow",  "Neuron" ="pink" ,"Microglial cell" = "#3176B7","Fibroblast" = "brown"
  )) +
  theme_few()+
  theme(plot.title = element_text(size=12,hjust=0.5))
p2


##热图
top1 = markers_after %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top1 <- c("OPALIN","MOG","DOCK8","P2RY12","MAP2","DCX","AQP4","SLC39A12","VCAN","OLIG2","PDGFRB","FOXC1")

top1_dotplot <- DotPlot(scRNA_combined, features = top1)+
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
ggsave(heatmap,filename = "D:/R/GEO/CNP_RE/20250902更新/heatmap.pdf",width =22,height = 9)

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
write.csv(as.data.frame(heatmap_data_transposed),"/data/h01005/hc/再整合/20250911/heatmap_data_transposed.csv", row.names = T)
# 调整保存尺寸（宽高互换）
ggsave(heatmap,
       filename = "/data/h01005/hc/再整合/20250911/heatmap_rotated.pdf",
       width = 9,      # 原高度
       height = 22)    # 原宽度
