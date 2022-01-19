
# MESSAGE -----------------------------------------------------------------
#
# author: Yulin Lyu, Cheng Li Lab, Peking University
# email: lvyulin@pku.edu.cn
#
# ---

# load package ------------------------------------------------------------

library(tidyverse)
library(magrittr)
library(data.table)

library(Seurat)

`%ni%` <- function(a, b) return(! a %in% b)

# WOT ---------------------------------------------------------------------

allSo <- readRDS("middata/hADSC_0618_scRNA/allSo_qc.rds")

# prepare input

usedSample <- c(
  "ADSC", "S1D0", "S1D0.5", "S1D1", "S1D2", "S1D4", "S1D8",
  "S2D4", "S2D8", "S2D12", "S2D20",
  "S3D4", "S3D8", "S3D12",
  "S4D1", "S4D2", "S4D4", "S4D10")

sampleTime <- c(
  -1, c(0, 0.5, 1, 2, 4, 8),
  -1 + 8 + c(4, 8, 12, 20),
  -1 + 8 + 20 + c(4, 8, 12),
  -1 + 8 + 20 + 12 + c(1, 2, 4, 10)) + 1

names(sampleTime) <- usedSample

allSo <- allSo[, allSo$sample %in% usedSample]

allSo$time <- sampleTime[allSo$sample]
table(allSo@meta.data[, c("time", "sample")])

wotExp <- allSo[["RNA"]]@data[VariableFeatures(allSo), ] %>% as.matrix %>% t
wotExp <- cbind(rownames(wotExp), wotExp)
colnames(wotExp)[1] <- "id"

write.table(wotExp, "middata/hADSC_0618_scRNA/wot/wotExp.txt", sep = "\t", row.names = F, quote = F)

wotDay <- as.matrix(allSo$time)
rownames(wotDay) <- colnames(allSo)
wotDay <- cbind(rownames(wotDay), wotDay)
colnames(wotDay) <- c("id", "day")
write.table(wotDay, "middata/hADSC_0618_scRNA/wot/wotDay.txt", sep = "\t", quote = F, row.names = F)

# perform WOT analysis according to developer instructions (https://broadinstitute.github.io/wot/)

# identify S4D10 celltype

end <- allSo[, allSo$sample == "S4D10"]
end %<>% FindVariableFeatures
end %<>% ScaleData %>% RunPCA
end %>% ElbowPlot(ndims = 50)
end %<>% FindNeighbors(dims = 1:20)
end %<>% FindClusters(resolution = 1)
end %<>% RunUMAP(dims = 1:20)

DimPlot(end, group.by = "seurat_clusters", label = T)

FeaturePlot(end, c("POU5F1", "SOX2", "NANOG", "KLF4"))
FeaturePlot(end, c("COL1A2", "DCN"))

end$celltype <- "other"
end$celltype[end$seurat_clusters %in% c()] <- "hCiPS"
end$celltype[end$seurat_clusters %in% c()] <- "Stromal"

# the celltypes were used to calculate right and wrong probability

allSo$CiPS <- readRDS("middata/hADSC_0618_scRNA/CiPS_score.rds")[colnames(allSo)]
allSo$Stromal <- readRDS("middata/hADSC_0618_scRNA/Stromal_score.rds")[colnames(allSo)]

saveRDS(allSo, "middata/hADSC_0618_scRNA/allSo_traj.rds")

# identify right trajectory -----------------------------------------------

allSo <- readRDS("middata/hADSC_0618_scRNA/allSo_traj.rds")

allSo$traj <- "F"
allSo$traj[allSo$sample == "prime"] <- "prime"
allSo$traj[allSo$sample == "H1"] <- "ES"
allSo$traj[allSo$sample %in% c("ADSC", "S1D0", "S1D0.5", "S1D1", "S1D2")] <- structure(
  names = c("ADSC", "S1D0", "S1D0.5", "S1D1", "S1D2"),
  str_c("R_", 0:4)
)[allSo$sample[allSo$sample %in% c("ADSC", "S1D0", "S1D0.5", "S1D1", "S1D2")]]

s1 <- allSo[, allSo$sample %in% c("S1D4", "S1D8")]
s1 %<>% FindVariableFeatures
s1 %<>% ScaleData %>% RunPCA
s1 %>% ElbowPlot(ndims = 50)
s1 %<>% FindNeighbors(dims = 1:20)
s1 %<>% FindClusters(resolution = 1)
s1 %<>% RunUMAP(dims = 1:20)

DimPlot(s1, group.by = "sample", label = T)
DimPlot(s1, group.by = "seurat_clusters", label = T)

FeaturePlot(s1, c("CiPS", "Stromal"))
VlnPlot(s1, c("CiPS", "Stromal"), group.by = "seurat_clusters", pt.size = 0)

s1$traj[s1$seurat_clusters %in% c()] <- "R5"
s1$traj[s1$seurat_clusters %in% c()] <- "R6"

allSo@meta.data[colnames(s1), "traj"] <- s1$traj

s2 <- allSo[, allSo$sample %in% c("S2D4", "S2D8", "S2D12", "S2D20")]
s2 %<>% FindVariableFeatures
s2 %<>% ScaleData %>% RunPCA
s2 %>% ElbowPlot(ndims = 50)
s2 %<>% FindNeighbors(dims = 1:20)
s2 %<>% FindClusters(resolution = 1)
s2 %<>% RunUMAP(dims = 1:20)

DimPlot(s2, group.by = "sample", label = T)
DimPlot(s2, group.by = "seurat_clusters", label = T)

FeaturePlot(s2, c("CiPS", "Stromal"))
VlnPlot(s2, c("CiPS", "Stromal"), group.by = "seurat_clusters", pt.size = 0)

s2$traj[s2$seurat_clusters %in% c()] <- "R_7"
s2$traj[s2$seurat_clusters %in% c()] <- "R_8"
s2$traj[s2$seurat_clusters %in% c()] <- "R_9"

allSo@meta.data[colnames(s2), "traj"] <- s2$traj

s34 <- allSo[, allSo$sample %in% c("S3D4", "S3D8", "S3D12", "S4D1", "S4D2", "S4D4", "S4D10")]
s34 %<>% FindVariableFeatures
s34 %<>% ScaleData %>% RunPCA
s34 %>% ElbowPlot(ndims = 50)
s34 %<>% FindNeighbors(dims = 1:20)
s34 %<>% FindClusters(resolution = 1)
s34 %<>% RunUMAP(dims = 1:20)

DimPlot(s34, group.by = "sample", label = T)
DimPlot(s34, group.by = "seurat_clusters", label = T)

FeaturePlot(s34, c("CiPS", "Stromal"))
VlnPlot(s34, c("CiPS", "Stromal"), group.by = "seurat_clusters", pt.size = 0)

s34$branch <- "other"
s34$branch[s34$seurat_clusters %in% c()] <- "Stromal-branch"

s34r <- s34[, s34$branch != "Stromal-branch"]
s34r %<>% FindVariableFeatures
s34r %<>% ScaleData %>% RunPCA
s34r %>% ElbowPlot(ndims = 50)
s34r %<>% FindNeighbors(dims = 1:20)
s34r %<>% FindClusters(resolution = 1)
s34r %<>% RunUMAP(dims = 1:20)

DimPlot(s34r, group.by = "sample", label = T)
DimPlot(s34r, group.by = "seurat_clusters", label = T)

FeaturePlot(s34r, c("CiPS", "Stromal"))
VlnPlot(s34r, c("CiPS", "Stromal"), group.by = "seurat_clusters", pt.size = 0)

FeaturePlot(s34r, c("POU5F1", "SOX2", "NANOG", "KLF4"))
FeaturePlot(s34r, c("APOA4", "SERPINA1", "RBP2", "CPA2"))

s34r$branch[s34r$seurat_clusters %in% c()] <- "Stromal-bias"
s34r$branch[s34r$seurat_clusters %in% c()] <- "CiPS-tip"
s34r$branch[s34r$seurat_clusters %in% c()] <- "BranchX-tip"
s34r$branch[s34r$seurat_clusters %in% c()] <- "Start-tip"

s34rr <- s34r[, s34r$branch != "Stromal-bias"]
s34rr %<>% FindVariableFeatures
s34rr %<>% ScaleData %>% RunPCA
s34rr %>% ElbowPlot(ndims = 50)
s34rr %<>% FindNeighbors(dims = 1:20)
s34rr %<>% FindClusters(resolution = 1)
s34rr %<>% RunUMAP(dims = 1:20)

# perform FateID analysis (https://github.com/dgrun/FateID) to get trajectory bias probability for Start, CiPS, and BranchX

FeaturePlot(s34rr, c("p_Start", "p_CiPS", "p_BranchX"))

s34rr$traj[s34rr$p_BranchX > 0.5] <- "BranchX"

s34rrr <- s34rr[, s34rr$traj != "BranchX"]
s34rrr %<>% FindVariableFeatures
s34rrr %<>% ScaleData %>% RunPCA
s34rrr %>% ElbowPlot(ndims = 50)
s34rrr %<>% FindNeighbors(dims = 1:20)
s34rrr %<>% FindClusters(resolution = 1)
s34rrr %<>% RunUMAP(dims = 1:20)

DimPlot(s34rrr, group.by = "sample", label = T)
DimPlot(s34rrr, group.by = "seurat_clusters", label = T)

FeaturePlot(s34rrr, c("CiPS", "p_CiPS"))
VlnPlot(s34rrr, c("CiPS", "p_CiPS"), group.by = "seurat_clusters", pt.size = 0)

FeaturePlot(s34rrr, c("POU5F1", "SOX2", "NANOG", "KLF4"))

s34rrr$traj[s34rrr$seurat_clusters %in% c()] <- "R_10"
s34rrr$traj[s34rrr$seurat_clusters %in% c()] <- "R_11"
s34rrr$traj[s34rrr$seurat_clusters %in% c()] <- "R_12"
s34rrr$traj[s34rrr$seurat_clusters %in% c()] <- "R_13"

allSo@meta.data[colnames(s34rrr), "traj"] <- s34rrr$traj

saveRDS(allSo, "middata/hADSC_0618_scRNA/allSo_traj.rds")

# integrate H1 ------------------------------------------------------------

# perform LSI projection (https://github.com/sajuukLyu/projectLSI) to get LSI coordiates lsiData

allSo[["lsi"]] <- CreateDimReducObject(
  lsiData[colnames(allSo), ], assay = "RNA", key = "LSI_"
)

umapRes <- uwot::umap(
  lsiData[, 1:20],
  n_neighbors = 40,
  min_dist = 0.3,
  metric = "euclidean",
  ret_model = T,
  verbose = T)

umapData <- umapRes$embedding
rownames(umapData) <- colnames(allSo)
colnames(umapData) <- str_c("UMAP_", 1:2)

allSo[["umap"]] <- CreateDimReducObject(
  umapData[colnames(allSo), ], assay = "RNA", key = "UMAP_"
)

saveRDS(allSo, "middata/hADSC_0618_scRNA/allSo_traj.rds")

# PAGA --------------------------------------------------------------------

fso <- allSo[, allSo$traj == "F"]
fso %<>% ScaleData %>% RunPCA
fso %<>% FindNeighbors(dims = 1:20)
fso %<>% FindClusters(resolution = .4)

allSo$part <- allSo$traj %>% as.character
allSo@meta.data[colnames(fso), "part"] <- fso$seurat_clusters %>% as.vector %>% str_c("F_", .)

cellMeta <- allSo@meta.data[, c("part", "sample", "traj")]
geneMeta <- GetAssay(allSo)[[]]
usedGene <- VariableFeatures(allSo)

library(reticulate)
sc <- import("scanpy")
adata <- sc$AnnData(
  X = allSo[["RNA"]]@data[usedGene, ] %>% Matrix::t(),
  obs = cellMeta,
  var = geneMeta[usedGene, ]
)
adata$obsm$update(X_pca = Embeddings(allSo, "lsi"))
adata$obsm$update(X_umap = Embeddings(allSo, "umap"))

adata$write_h5ad("middata/hADSC_0618_scRNA/allAnn.h5ad")

# perform PAGA analysis (https://github.com/theislab/paga) to get FA coordinates of trajectory

adata <- sc$read_h5ad("middata/hADSC_0618_scRNA/allAnn_paga.h5ad")

fa_emb <- adata$obsm['X_draw_graph_fa']
colnames(fa_emb) <- c("FA_1", "FA_2")
row.names(fa_emb) <- adata$obs_names$to_list()

allSo[["fa"]] <- CreateDimReducObject(embeddings = fa_emb, key = "FA_", assay = "RNA")

DimPlot(allSo, group.by = "sample", reduction = "fa", label = T)

saveRDS(allSo, "middata/hADSC_0618_scRNA/allSo_traj.rds")

# trajectory for S3 and S4 ------------------------------------------------

s34 <- allSo[, allSo$sample %in% c("S3D4", "S3D8", "S3D12", "S4D1", "S4D2", "S4D4", "S4D10")]
s34 %<>% FindVariableFeatures
s34 %<>% ScaleData %>% RunPCA
s34 %>% ElbowPlot(ndims = 50)
s34 %<>% FindNeighbors(dims = 1:20)
s34 %<>% FindClusters(resolution = 1)
s34 %<>% RunUMAP(dims = 1:20)

DimPlot(s34, group.by = "sample", label = T)

usedData <- allSo[["RNA"]]@data[VariableFeatures(allSo), ] %>% as.matrix

library(phateR)

dataPhate <- phate(usedData, ndim = 2, knn = 10, npca = 100)
s34[["phate"]] <- CreateDimReducObject(
  embeddings = dataPhate$embedding * 100,
  assay = "RNA",
  key = "PHATE_"
)

DimPlot(s34, group.by = "sample", reduction = "phate", dims = c(1, 3), label = T)

s34rt <- s34[, s34$traj %in% c("R_10", "R_11", "R_12", "R_13")]

# perform PAGA analysis to get pseudotime of s34rt

s34$ptime <- NA
s34@meta.data[colnames(s34rt), "ptime"] <- s34rt$ptime

saveRDS(s34, "middata/hADSC_0618_scRNA/s34_traj.rds")

