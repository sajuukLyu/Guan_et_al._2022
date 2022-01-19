
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
library(SingleCellExperiment)
library(scran)

`%ni%` <- function(a, b) return(! a %in% b)

# load data ---------------------------------------------------------------

usedSample <- c(
  "ASF", "S1D16", "S2D24", "prime_ASF")

dataList <- map(str_c("data/hASF_scRNA/", usedSample), Read10X)

names(dataList) <- usedSample
dataList <- imap(dataList, ~ {set_colnames(.x, str_c(.y, "_", colnames(.x)))})
map_int(dataList, ncol)

countData <- do.call(cbind, dataList)
allSo <- CreateSeuratObject(countData)

saveRDS(allSo, "middata/hASF_scRNA/allSo_raw.rds")

# quality control ---------------------------------------------------------

allSo <- readRDS("middata/hASF_scRNA/allSo_raw.rds")
allSo$sample <- str_extract(colnames(allSo), ".*_") %>% str_remove("_$")
table(allSo$sample)

mtGene <- rownames(allSo) %>% str_subset("^MT-")
allSo$mt <- PercentageFeatureSet(allSo, features = mtGene)

allSo <- allSo[, allSo$nFeature_RNA > 500 & allSo$nCount_RNA > 1000 & allSo$mt < 20]
allSo %<>% NormalizeData

# check doublet
soList <- SplitObject(allSo, "sample")
soList <- map(soList, FindVariableFeatures)
SCElist <- map(soList, ~ SingleCellExperiment(assays = list(counts = .x[["RNA"]]@counts, logcounts = .x[["RNA"]]@data)))
varGeneList <- map(soList, VariableFeatures)
dblDen <- map2(SCElist, varGeneList, ~ doubletCells(.x, subset.row = intersect(rownames(.x), .y), d = 20))
dblDen <- map(dblDen, ~ log10(.x + 1))
allSo$dbl <- unlist(dblDen) %>% unname

allSo %<>% FindVariableFeatures
allSo %<>% ScaleData %>% RunPCA
allSo %>% ElbowPlot(ndims = 50)
allSo %<>% FindNeighbors(dims = 1:20)
allSo %<>% FindClusters(resolution = 1)
allSo %<>% RunUMAP(dims = 1:20)

DimPlot(allSo, group.by = "sample", label = T)

VlnPlot(allSo, "nCount_RNA", group.by = "seurat_clusters", pt.size = 0)
VlnPlot(allSo, "nFeature_RNA", group.by = "seurat_clusters", pt.size = 0)
VlnPlot(allSo, "dbl", group.by = "seurat_clusters", pt.size = 0)

badClu <- c() # low quality clusters
allSo <- allSo[, allSo$seurat_clusters %ni% badClu]

saveRDS(allSo, "middata/hASF_scRNA/allSo_qc.rds")
