# install.packages('Seurat')      https://satijalab.org/seurat/articles/install_v5



# 提高Seurat速度和性能
# setRepositories(ind = 1:3, addURLs = c('https://satijalab.r-universe.dev', 
#                                        'https://bnprks.r-universe.dev/'))
# options()$repos
# install.packages(c("BPCells", "presto"))
# BiocManager::install("glmGamPoi")

# 支持使用其他整合和差异表达方法
# remotes::install_github("satijalab/seurat-wrappers", quiet = TRUE)
# Myeloid_data_ob_v3 <- readRDS("D:/.GitHub/Myeloid_data_ob_v3.rds")
# Myeloid_data_ob_v3
# UpdateSeuratObject(Myeloid_data_ob_v3) ->mye

# 
# mye@assays$RNA@var.features
# mye@assays$RNA@meta.features
# mye@assays$RNA@counts
# mye@meta.data


# mye   
# # An object of class Seurat 
# # 32285 features across 24503 samples within 1 assay 
# # Active assay: RNA (32285 features, 2000 variable features)
# # 3 layers present: counts, data, scale.data
# # 3 dimensional reductions calculated: mnn, pca, umap


# myeloid_expr_counts <- mye@assays$RNA@counts

# saveRDS(myeloid_expr_counts,file = "myeloid_expr_counts.rds")

mye_counts <- readRDS("E:/_GitHub/_data/myeloid_expr_counts.rds")


# 自动注释 BiocManager::install("celldex") ####

library(celldex)
# options(timeout = 500)
# surveyReferences()
ref1 <- fetchReference("mouse_rnaseq", "2024-02-26")
ref2 <- celldex::ImmGenData()


# library(SummarizedExperiment)
# assay(x = ref,"logcounts")

table(ref1$label.fine)
table(ref2$label.fine)

library(SingleCellExperiment)

sce <- SingleCellExperiment(assays = list(counts = mye_counts))
sce






library(SingleR)

# 单个参考集
pred <- SingleR(test = sce, ref = ref, 
                labels = ref$label.fine, assay.type.test=1,
                de.method="wilcox")
colnames(pred)

pred %>% as.data.frame() %>% group_by(pruned.labels) %>% 
    summarise(n=n())

pred$scores
plotScoreHeatmap(pred)

colData(sce) <- pred
sce


to.remove <- is.na(pred$pruned.labels)
table(Label=pred$labels, Removed=to.remove)

to.remove <- pruneScores(pred, min.diff.med=0.2)
table(Label=pred$labels, Removed=to.remove)
plotDeltaDistribution(pred)

# 多个参考集

com.res2 <- SingleR(test = sce, assay.type.test=1,
                    ref = list(mouse_rnaseq=ref1, Immune=ref2), 
                    labels = list(ref1$label.main, ref2$label.main),
                    #BPPARAM=SnowParam(5)
                    )

# Check the final label from the combined assignment.
table(com.res2$labels) 




# Seurat   ####
library(Seurat)
mye <- CreateSeuratObject(counts = mye_counts, 
                   project = "myeloid", 
                   # min.cells = 3, 
                   # min.features = 200
                   )
mye@meta.data %>% head()
Layers(mye)
Seurat::Assays(mye)
mye[["RNA"]]
LayerData(mye, assay = "RNA", layer = "counts")

colnames(mye[[]])
VariableFeatures(mye)

levels(mye)
table(Idents(mye))


# 细胞
Cells(mye)
ncol(mye)
str_match(Cells(mye),pattern = "...\\_\\d\\_") %>% table()

# 基因
Features(mye)
nrow(mye)


mye_counts[c("Cd3d", "Cd14", "Ms4a1"),]

#           预处理  小鼠线粒体、血红蛋白基因####
# 质量控制 细胞选择和过滤
str_view(Features(mye),pattern = "^mt-")
mye[["percent_mito"]] <- PercentageFeatureSet(mye, pattern = "^mt-", )
summary(mye@meta.data$percent_mito)

str_view(Features(mye),pattern = "^Hb[ab]")
mye[["percent_Hb"]] <- PercentageFeatureSet(mye, 
                                            pattern = "^Hb[ab]")
summary(mye@meta.data$percent_Hb)

VlnPlot(mye, 
        features = c("nFeature_RNA", "nCount_RNA", 
                     "percent_mito", "percent_Hb"),
        ncol = 2,layer = "counts",
        group.by = "orig.ident",pt.size = 0.1)



FeatureScatter(mye, feature1 = "nCount_RNA", feature2 = "percent_mito") +
    FeatureScatter(mye, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

mye_QC <- subset(mye, subset = nFeature_RNA < 6500 & nCount_RNA <37500
               & percent_mito < 5 & percent_Hb<1)

# 标准化、高度可变基因、缩放
mye_normalized <- NormalizeData(mye_QC, normalization.method = "LogNormalize",
                     scale.factor = 10000)
mye_HVG <- FindVariableFeatures(mye_normalized, selection.method = "vst",
                            nfeatures = 2000)
top10 <- VariableFeatures(mye_HVG) %>% head(.,n=10)

p1 <- VariableFeaturePlot(mye_HVG)
LabelPoints(plot = p1,points = top10,repel = T)

# 默认情况下，仅缩放可变特征。
mye_scaled <- ScaleData(mye_HVG,features = Features(mye_HVG))

mye <- mye_scaled 
# PCA降维
mye <- RunPCA(object = mye,npcs = 50, verbose = FALSE)

ElbowPlot(mye, ndims = 50) # 找到“肘部点”，即变化开始趋于平缓的位置
# VizDimLoadings(mye, dims = 1:3, reduction = "pca")
# DimHeatmap(mye, dims = 1:3, cells = 500, balanced = TRUE)

# 相互最近邻 MNN  ?IntegrateLayers()    BiocManager::install("batchelor")
mye <- IntegrateLayers(object = mye, method = FastMNNIntegration,
                       batch = mye$orig.ident,k = 30,
                       assay = "RNA",layers = "data",scale.layer="scale.data",
                       orig.reduction = "pca",new.reduction = 'integrated.mnn', 
                       verbose = FALSE)

DimPlot(mye, reduction = "integrated.mnn",group.by = "orig.ident")
# HarmonyIntegration

# 聚类
dims_to_use <- 1:30 
mye <- FindNeighbors(mye, dims = dims_to_use,reduction = "integrated.mnn")
mye <- FindClusters(mye, resolution = seq(0.1,2,by = 0.1), verbose = FALSE)
head(Idents(mye),5)
DimPlot(mye, reduction = "integrated.mnn")


# tsne

mye <- RunTSNE(mye, reduction = "pca")
DimPlot(mye,reduction = "tsne", 
        label = TRUE, repel = TRUE,
        label.size = 4) 
# umap
mye <- FindClusters(mye, resolution = 1, verbose = FALSE)

mye <- RunUMAP(mye, dims = dims_to_use,umap.method = "uwot")

mye@meta.data %>% colnames()

mye <- JoinLayers(mye)
DimPlot(mye,reduction = "umap", 
        label = TRUE, repel = TRUE,
        label.size = 4) 

# 细胞注释
mye$celltype <- com.res2$pruned.labels
DimPlot(mye, group.by = "celltype", 
        label = TRUE, repel = TRUE, seed = 123, 
        label.size = 4) 



# 单核细胞marker
str_view(Features(mye),pattern = "(Cd14)|(Fcgr3)")
FeaturePlot(mye, features = c("Cd14", "Fcgr3"), 
            reduction = "umap", cols = c("lightgrey", "red"))

# 巨噬细胞marker
str_view(Features(mye),pattern = "(Cd68)|(Csf1r)|(Itgax)")

FeaturePlot(mye, features = c("Cd68", "Csf1r"),  
            reduction = "umap", cols = c("lightgrey", "blue"))

FeaturePlot(mye, features = c("Itgax"),  
            reduction = "umap", cols = c("lightgrey", "blue"))


# T细胞marker
str_view(Features(mye),pattern = "(Cd3[a-z]$)|(Cd4$)|(Cd8a)")
FeaturePlot(
    mye,
    features = c("Cd3g", "Cd3d", "Cd3e", "Cd4", "Cd8a"),
    reduction = "umap",
    cols = c("lightgrey", "blue"),
) 

# NK细胞marker
str_view(Features(mye),pattern = "(Ncam1)|(Klrd1)|(Nkg7)")
FeaturePlot(mye, features = c("Klrd1", "Nkg7", "Ncam1"), 
            reduction = "umap", cols = c("lightgrey", "red"))

# B细胞marker
str_view(Features(mye),pattern = "(Cd19)|(Ms4a1$)")
FeaturePlot(mye, features = c("Cd19", "Ms4a1"),  
            reduction = "umap", cols = c("lightgrey", "blue"))





# 手动注释 #### 

mks <- FindAllMarkers(mye, 
                      min.pct = 0.25, 
                      logfc.threshold = 0.25, 
                      assay = "RNA")

mks <- mks %>% group_by(cluster) %>% slice_max(avg_log2FC,n = 10) %>% 
    select(cluster,gene) %>% 
    summarise(gene = str_flatten_comma(gene))

mks




# 初始化细胞类型# 初始化细胞gene类型
mye$celltype <- "Unknown"

# 注释T细胞
mye$celltype[WhichCells(mye, expression = Cd3e > 1)] <- "T cell"

# 注释NK细胞
mye$celltype[WhichCells(mye, expression = Nkg7 > 1 | Klrd1 > 1)] <- "NK cell"

# 注释B细胞
mye$celltype[WhichCells(mye, expression = Ms4a1 > 1)] <- "B cell"







