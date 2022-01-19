
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

# public human limb bud ---------------------------------------------------

library(loomR)

rawData <- connect(filename = "data/pub_human/GSM4274188_CS13_limbbud_rawdata.loom", mode = "r+")
rawData <- connect(filename = "data/pub_human/GSM4274189_CS15_limbbud_rawdata.loom", mode = "r+")
rawData <- connect(filename = "data/pub_human/GSM4274190_CS15_2__limbbud_rawdata.loom", mode = "r+")

countData <- rawData[["matrix"]][, ] %>% t %>% as("dgCMatrix")
rownames(countData) <- rawData[["row_attrs"]][["Gene"]][]
colnames(countData) <- rawData[["col_attrs"]][["CellID"]][]

dataList <- list()
dataList$emb1 <- countData
dataList$emb2 <- countData
dataList$emb3 <- countData

countData <- do.call(cbind, dataList)
hg <- CreateSeuratObject(countData)
hg$sample <- rep(names(dataList), map_int(dataList, ncol))

hg %<>% NormalizeData()
hg %<>% FindVariableFeatures %>% ScaleData %>% RunPCA
hg %>% ElbowPlot(ndims = 50)
hg %<>% RunUMAP(dims = 1:20)
hg %<>% FindNeighbors(dims = 1:20) %>% FindClusters(resolution = 0.6)

DimPlot(hg, group.by = "sample", label = T)
DimPlot(hg, group.by = "seurat_clusters", label = T)

mtGene <- str_subset(rownames(hg), "^MT-")
hg$mt <- PercentageFeatureSet(hg, features = mtGene)

hg$nCount_RNA
FeaturePlot(hg, "nCount_RNA")
FeaturePlot(hg, "mt")

VlnPlot(hg, "PDGFRA", group.by = "seurat_clusters", pt.size = 0)

usedClu <- c()
hg <- hg[, hg$seurat_clusters %in% usedClu]
hg$celltype <- hg$sample

saveRDS(hg, "middata/pub/hg.rds")

# pub axolotl data --------------------------------------------------------

dataRaw <- fread("data/pub_axolotl/GSE106269_Table_S7.csv")

dataRe <- dataRaw[, `0610011F06RIK`:ZZZ3] %>% as.matrix %>% t
colnames(dataRe) <- dataRaw[, cell_id]

axo <- CreateSeuratObject(dataRe)

axo$celltype <- dataRaw[, celltype]
axo$time <- dataRaw[, timepoint]

table(axo@meta.data[, c("time", "celltype")])

axo %<>% FindVariableFeatures %>% NormalizeData %>% ScaleData %>% RunPCA
axo %>% ElbowPlot(ndims = 50)
axo %<>% RunUMAP(dims = 1:15)

DimPlot(axo, group.by = "celltype", label = T)
DimPlot(axo, group.by = "time", label = T)

saveRDS(axo, "middata/pub/axo.rds")

# pub frog data -----------------------------------------------------------

dataPath <- "data/pub_frog/GSE165901_RAW"

dataFile <- list.files(dataPath, ".gz")
dataSample <- str_remove_all(dataFile, "_feature.*|_barcode.*|_matrix.*")
dataBase <- str_remove_all(dataFile, ".*_")

str_c(dataPath, "/", unique(dataSample)) %>% map(dir.create)
file.copy(list.files(dataPath, ".gz", full.names = T), str_c(dataPath, "/", dataSample, "/", dataBase))

samplePath <- list.dirs(dataPath)[-1]
sampleData <- map(samplePath, Read10X)

names(sampleData) <- str_remove_all(samplePath, ".*/")
sampleData <- imap(sampleData, ~ {colnames(.x) <- str_c(.y, "_", colnames(.x)); .x})

ctGene <- head(rownames(sampleData[[1]]), 37)

xenlaList <- list()
map_int(sampleData, nrow)
map_int(sampleData, ncol)

# preprocess according to the original article

ppSeurat <- function(
  mtx, sampleName, ctGene = ctGene,
  UMIhigh = 1e5, UMIlow = 0, MThigh = 100) {
  require(Seurat)
  
  so <- CreateSeuratObject(mtx)
  so$sample <- sampleName
  so$MTpct <- PercentageFeatureSet(so, features = ctGene)
  
  so <- subset(so, nCount_RNA >= UMIlow & nCount_RNA <= UMIhigh & MTpct <= MThigh)
  
  so <- NormalizeData(so) %>%
    FindVariableFeatures() %>%
    ScaleData() %>%
    RunPCA()
  so
}

xenlaList$`0dpa` <- sampleData$GSM5057655_107606_Xen_FACS_BL0dpa %>%
  ppSeurat("0dpa", ctGene = ctGene, 30000, 2000, 10)
xenlaList$`3dpa` <- sampleData$GSM5057656_107606_Xen_FACS_BL3dpa %>%
  ppSeurat("3dpa", ctGene = ctGene, 60000, 2000, 10)
xenlaList$EarlyBL <- sampleData$GSM5057661_107606_Xen_Pool_BL7_10_14dpa %>%
  ppSeurat("EarlyBL", ctGene = ctGene, 40000, 3000, 20)
xenlaList$LateBL <- sampleData$GSM5057660_107606_Xen_Pool_BL14_20_52dpa %>%
  ppSeurat("LateBL", ctGene = ctGene, 40000, 2000, 10)
xenlaList$`14dpa` <- sampleData$GSM5057665_88527_Xen_BL14dpa %>%
  ppSeurat("14dpa", ctGene = ctGene, 40000, 3000, 10)
xenlaList$s50 <- sampleData$GSM5057657_107606_Xen_LBst50 %>%
  ppSeurat("s50", ctGene = ctGene, 40000, 2000, 20)
xenlaList$s51 <- sampleData$GSM5057658_107606_Xen_LBst51 %>%
  ppSeurat("s51", ctGene = ctGene, 40000, 2000, 20)
xenlaList$s52 <- sampleData$GSM5057659_107606_Xen_LBst52 %>%
  ppSeurat("s52", ctGene = ctGene, 40000, 2000, 10)
xenlaList$s54_0dpa <- sampleData$GSM5057666_88528_Xen_Pool_LBst54_BL0dpa %>%
  ppSeurat("s54_0dpa", ctGene = ctGene, 40000, 3000, 10)

map_int(xenlaList, ncol)

# split mixed samples

usedSo <- xenlaList$EarlyBL
usedSo %>% ElbowPlot(ndims = 50)
usedSo <- RunUMAP(usedSo, dims = 1:20)

FeaturePlot(usedSo, "Venus", label = F)
FeaturePlot(usedSo, "mCherryB51", label = F)
FeaturePlot(usedSo, "mTFP1", label = F)

usedSo$seem_7dpa <- usedSo[["RNA"]]@data["Venus", ] > 0
usedSo$seem_10dpa <- usedSo[["RNA"]]@data["mCherryB51", ] > 0
usedSo$seem_14dpa <- usedSo[["RNA"]]@data["mTFP1", ] > 0

table(usedSo@meta.data[, c("seem_7dpa", "seem_14dpa")])

usedSo$sample <- "10dpa"
usedSo$sample[usedSo$seem_7dpa & !usedSo$seem_14dpa] <- "7dpa"
usedSo$sample[!usedSo$seem_7dpa & usedSo$seem_14dpa] <- "14dpa"
xenlaList$EarlyBL <- usedSo


usedSo <- xenlaList$LateBL
usedSo %>% ElbowPlot(ndims = 50)
usedSo <- RunUMAP(usedSo, dims = 1:20)

FeaturePlot(usedSo, "mTFP1", label = F)
FeaturePlot(usedSo, "mCherryB51", label = F)
FeaturePlot(usedSo, "Venus", label = F)

usedSo$seem_14dpa <- usedSo[["RNA"]]@data["mTFP1", ] > 0
usedSo$seem_20dpa <- usedSo[["RNA"]]@data["mCherryB51", ] > 0
usedSo$seem_52dpa <- usedSo[["RNA"]]@data["Venus", ] > 0

table(usedSo@meta.data[, c("seem_14dpa", "seem_52dpa")])

usedSo$sample <- "20dpa"
usedSo$sample[usedSo$seem_14dpa & !usedSo$seem_52dpa] <- "14dpa"
usedSo$sample[!usedSo$seem_14dpa & usedSo$seem_52dpa] <- "52dpa"
xenlaList$LateBL <- usedSo


usedSo <- xenlaList$s54_0dpa

FeaturePlot(usedSo, "mTFP1", label = F)
FeaturePlot(usedSo, "tdTomato", label = F)

usedSo$seem_s54 <- usedSo[["RNA"]]@data["mTFP1", ] > 0
usedSo$seem_0dpa <- usedSo[["RNA"]]@data["tdTomato", ] > 0

table(usedSo@meta.data[, c("seem_s54", "seem_0dpa")])

usedSo$sample <- "unknown"
usedSo$sample[usedSo$seem_s54 & !usedSo$seem_0dpa] <- "s54"
usedSo$sample[!usedSo$seem_s54 & usedSo$seem_0dpa] <- "0dpa"
xenlaList$s54_0dpa <- usedSo

# combine data

names(xenlaList)
map_int(xenlaList, nrow)

countData <- xenlaList %>% map(~ {.x[["RNA"]]@counts})
countData$`0dpa` %<>% {.[-nrow(.), ]}
countData$`3dpa` %<>% {.[-nrow(.), ]}
map_int(countData, nrow)

countData %<>% purrr::reduce(cbind)
xenlaAll <- CreateSeuratObject(countData)
xenlaAll$sample <- map(xenlaList, ~ {.x$sample}) %>% do.call(c, .)

table(xenlaAll$sample)

xenla <- subset(xenlaAll, sample %in% c("0dpa", "3dpa", "7dpa", "10dpa", "14dpa", "20dpa", "s50", "s51", "s52", "s54"))

xenla %<>% NormalizeData
xenla %<>% FindVariableFeatures %>% ScaleData %>% RunPCA
xenla %>% ElbowPlot(ndims = 50)
xenla %<>% RunUMAP(dims = 1:20)
xenla %<>% FindNeighbors(dims = 1:20) %>% FindClusters(resolution = .6)

DimPlot(xenla, group.by = "sample", label = T)
DimPlot(xenla, group.by = "seurat_clusters", label = T)

FeaturePlot(xenla, "prrx1.L")
FeaturePlot(xenla, "dcn.L")
FeaturePlot(xenla, "col1a2.L")

VlnPlot(xenla, "prrx1.L", group.by = "seurat_clusters", pt.size = 0)

usedClu <- c()
xenla <- xenla[, xenla$seurat_clusters %in% usedClu]

saveRDS(xenla, "middata/pub/xenlaRaw.rds")

keepGene <- rownames(xenla[["RNA"]]@data) %>% str_subset("\\.L$")

xenlaHomoData <- xenla[["RNA"]]@data[keepGene, ]
rownames(xenlaHomoData) %<>% str_remove_all("\\.L$") %>% toupper()
rownames(xenlaHomoData) %<>% str_remove_all("-LIKE$")

xenlaHomoData %<>% apply(2, function(x) tapply(x, rownames(xenlaHomoData), mean))

xenlaHomo <- CreateSeuratObject(xenlaHomoData)
xenlaHomo$sample <- xenla$sample

saveRDS(xenlaHomo, "middata/pub/xenlaHomo.rds")

# pub mouse data ----------------------------------------------------------

# mouse blastema

mmRaw <- readRDS("data/pub_mouse/Mouse_Limb_BGI_Celltype.RDS")

countData <- mmRaw[["RNA"]]@counts
metaData <- mmRaw@meta.data

mm <- CreateSeuratObject(countData)

mm$sample <- metaData[colnames(mm), "stim"]
mm$celltype <- metaData[colnames(mm), "Celltype"]

mm %<>% NormalizeData
mm %<>% FindVariableFeatures %>% ScaleData %>% RunPCA
mm %>% ElbowPlot(ndims = 50)
mm %<>% RunUMAP(dims = 1:20)
mm %<>% FindNeighbors(dims = 1:20) %>% FindClusters(resolution = 0.6)

DimPlot(mm, group.by = "sample")
DimPlot(mm, group.by = "celltype")

FeaturePlot(mm, "ENSMUSG00000045394-Epcam")
FeaturePlot(mm, "ENSMUSG00000029661-Col1a2")
FeaturePlot(mm, "ENSMUSG00000026586-Prrx1")

table(mm$celltype)
mm <- mm[, mm$celltype %in% c("CT_cell", "Epithelial_cell")]

usedGene <- str_remove(rownames(countData), "ENSMUSG[0-9]+-")
homoRes <- homologene::mouse2human(usedGene)

dupMouseGene <- homoRes$mouseGene %>% {.[duplicated(.)]} %>% unique()

for(i in dupMouseGene) {
  message(i)
  id <- homoRes[homoRes$mouseGene == i, "humanID"]
  idRemove <- id %>% setdiff(min(.))
  
  homoRes <- homoRes[!(homoRes$mouseGene == i & homoRes$humanID %in% idRemove), ]
}

mmHomoData <- mm[["RNA"]]@counts[match(homoRes$mouseGene, usedGene), ]
rownames(mmHomoData) <- homoRes$humanGene

mmHomoData <- apply(as.matrix(mmHomoData), 2, function(x) tapply(x, homoRes$humanGene, sum))
mmHomoData <- as(mmHomoData, "dgCMatrix")

mm <- CreateSeuratObject(mmHomoData)

# mouse limb bud

dataPath <- "data/pub_mouse/GSE135985_RAW"
samplePath <- list.dirs(dataPath)[-1]
sampleData <- map(samplePath, Read10X)
names(sampleData) <- str_remove_all(samplePath, ".*/")
sampleData <- imap(sampleData, ~ {colnames(.x) <- str_c(.y, "_", colnames(.x)); .x})

ctGene <- rownames(sampleData[[1]]) %>% str_subset("mt-")
mmList <- imap(sampleData, ~ ppSeurat(.x, .y, ctGene, 200, 20))
map_int(mmList, ncol)

countData <- mmList %>% map(~ {.x[["RNA"]]@counts})
usedGene <- map(countData, rownames) %>% reduce(intersect)
countData <- map(countData, ~ .x[usedGene, ]) %>% reduce(cbind)

homoRes <- homologene::mouse2human(usedGene)

dupMouseGene <- homoRes$mouseGene %>% {.[duplicated(.)]} %>% unique()

for(i in dupMouseGene) {
  message(i)
  id <- homoRes[homoRes$mouseGene == i, "humanID"]
  idRemove <- id %>% setdiff(min(.))
  
  homoRes <- homoRes[!(homoRes$mouseGene == i & homoRes$humanID %in% idRemove), ]
}

mmHomoData <- countData[homoRes$mouseGene, ]
rownames(mmHomoData) <- homoRes$humanGene

mmHomoData <- apply(as.matrix(mmHomoData), 2, function(x) tapply(x, homoRes$humanGene, sum))
mmHomoData <- as(mmHomoData, "dgCMatrix")

mmLimb <- CreateSeuratObject(mmHomoData)
mmLimb$sample <- map(mmList, ~ {.x$sample}) %>% reduce(c)

mmLimb <- mmLimb[, mmLimb$sample %in% c(
  "12_Processed_E11", "13_Processed_E14")]
mmLimb$state <- structure(
  c("E11", "E14"),
  names = c("12_Processed_E11", "13_Processed_E14")
)[mmLimb$sample]

mmLimb %<>% NormalizeData
mmLimb %<>% FindVariableFeatures %>% ScaleData %>% RunPCA
mmLimb %>% ElbowPlot(ndims = 50)
mmLimb %<>% RunUMAP(dims = 1:20)
mmLimb %<>% FindNeighbors(dims = 1:20) %>% FindClusters(resolution = 0.6)

DimPlot(mmLimb, group.by = "sample", label = T)
DimPlot(mmLimb, group.by = "seurat_clusters", label = T)

FeaturePlot(mmLimb, "nCount_RNA")
FeaturePlot(mmLimb, "PDGFRA")

VlnPlot(mmLimb, "PDGFRA", group.by = "seurat_clusters", pt.size = 0)

usedClu <- c()
mmLimb <- mmLimb[, mmLimb$seurat_clusters %in% usedClu]

# combine mouse data

usedGene <- intersect(rownames(mm), rownames(mmLimb))

normData <- cbind(
  mmLimb[["RNA"]]@data[usedGene, ],
  mm[["RNA"]]@data[usedGene, ])
countData <- cbind(
  mmLimb[["RNA"]]@counts[usedGene, ],
  mm[["RNA"]]@counts[usedGene, ])

mmAll <- CreateSeuratObject(countData)
mmAll[["RNA"]]@data <- normData

mmAll$time <- c(mmLimb$state, str_c(mm$sample, "pa"))

mmAll %<>% FindVariableFeatures %>% ScaleData %>% RunPCA
mmAll %>% ElbowPlot(ndims = 50)
mmAll %<>% RunUMAP(dims = 1:20)
mmAll %<>% FindNeighbors(dims = 1:20) %>% FindClusters(resolution = 0.6)

DimPlot(mmAll, "time")

saveRDS(mmAll, "middata/pub/mmHomo.rds")

# human integrate ---------------------------------------------------------

adsc <- readRDS("middata/hADSC_0618_scRNA/allSo_traj.rds")
asf <- readRDS("middata/hASF_scRNA/allSo_qc.rds")
hef <- readRDS("middata/HEF_scRNA/allSo_qc.rds")
limb <- readRDS("middata/pub/hg.rds")

table(adsc$sample)
adsc <- adsc[, adsc$sample %in% c(
  "ADSC", "S1D0", "S1D0.5", "S1D1", "S1D2", "S1D4", "S1D8",
  "S2D4", "S2D8", "S2D12", "S2D20",
  "prime", "H1"
)]
adsc$batch <- "adsc"
adsc$batch[adsc$sample == "H1"] <- "H1"

table(asf$sample)
asf$batch <- "asf"

table(hef$sample)
hef <- hef[, hef$sample %in% c("HEF", "S2D20", "prime_HEF")]
hef$sample <- str_replace(hef$sample, "S2D20", "S2D20_HEF")
hef$batch <- "hef"

table(limb$sample)
limb$batch <- limb$sample

usedGene <- intersect(rownames(adsc), rownames(limb))

countData <- list(
  adsc = adsc[["RNA"]]@counts[usedGene, ],
  hef = hef[["RNA"]]@counts[usedGene, ],
  asf = asf[["RNA"]]@counts[usedGene, ],
  limb = limb[["RNA"]]@counts[usedGene, ]
)
colnames(countData$hef) %<>% str_replace("^S2D20", "S2D20_HEF")
countData %<>% do.call(cbind, .)

hgAll <- CreateSeuratObject(countData)
hgAll$sample <- c(adsc$sample %>% as.vector, hef$sample %>% as.vector, asf$sample %>% as.vector, limb$sample %>% as.vector
hgAll$batch <- c(adsc$batch, hef$batch, asf$batch, limb$batch)

table(hgAll$sample)
table(hgAll$batch)

hgAll %<>% NormalizeData()
hgAll %<>% FindVariableFeatures()
hgAll %<>% ScaleData() %>% RunPCA()

library(harmony)

hgAll %<>% RunHarmony("batch", sigma = 0.075)
hgAll %>% ElbowPlot(ndims = 50, reduction = "harmony")
hgAll %<>% RunUMAP(dims = 1:30, reduction = "harmony", n.neighbors = 50, min.dist = .3)
hgAll %<>% FindNeighbors(reduction = "harmony", dims = 1:30) %>% FindClusters(resolution = .6)

DimPlot(hgAll, group.by = "sample")
DimPlot(hgAll, group.by = "seurat_clusters")

hgAll$celltype <- "other"
hgAll$celltype[hgAll$sample %in% c("emb1", "emb2", "emb3")] <- "limb"
hgAll$celltype[hgAll$sample %in% c("ADSC")] <- "ADSC"
hgAll$celltype[hgAll$sample %in% c("HEF")] <- "HEF"
hgAll$celltype[hgAll$sample %in% c("ASF")] <- "ASF"
hgAll@meta.data[colnames(adsc)[adsc$traj == "R_6"], "celltype"] <- "ADSC_S1"
hgAll@meta.data[colnames(adsc)[adsc$traj == "R_9"], "celltype"] <- "ADSC_S2"

# celltypes of hef and asf are identified separately according to similarity with adsc

table(hgAll$celltype)

saveRDS(hgAll, "middata/early/hgAll.rds")

# dediff identity score ---------------------------------------------------

# human

hgAll <- readRDS("middata/early/hgAll.rds")
hgAll <- hgAll[, hgAll$celltype != "other"]

table(hgAll$celltype)
Idents(hgAll) <- hgAll$celltype
hgDiff <- FindMarkers(
  hgAll, ident.1 = c("limb"), ident.2 = c("ADSC", "HEF", "ASF"),
  min.pct = .25, only.pos = F)

hgUpGene <- hgDiff %>% {rownames(.)[.$p_val_adj < 0.05 & .$avg_logFC > 0.25 & .$pct.1 > .$pct.2]}
hgDownGene <- hgDiff %>% {rownames(.)[.$p_val_adj < 0.05 & .$avg_logFC < -0.25 & .$pct.1 < .$pct.2]}
hgDiffGene <- c(hgUpGene, hgDownGene)

an <- hgAll[["RNA"]]@data[hgDiffGene, hgAll$celltype == "limb"] %>% as.matrix %>% rowMeans
bn <- hgAll[["RNA"]]@data[hgDiffGene, hgAll$celltype %in% c("ADSC", "HEF", "ASF")] %>% as.matrix %>% rowMeans

Dmat <- 2 * matrix(c(sum(an^2), sum(an*bn), sum(an*bn), sum(bn^2)), nrow = 2)

Amat <- matrix(c(1, 0, 0, 1, -1, -1), nrow = 2)
bvec <- c(0, 0, -1)

dData <- apply(as.matrix(hgAll[["RNA"]]@data[hgDiffGene, ]), 2, function(x) 2 * c(sum(an*x), sum(bn*x)))
dim(dData)

library(quadprog)

sovData <- apply(dData, 2, function(x) (solve.QP(Dmat, x, Amat, bvec))$solution)
dim(sovData)

hgAll$score <- sovData[1, ] / colSums(sovData)
VlnPlot(hgAll, "score", group.by = "celltype", pt.size = 0)

saveRDS(hgAll, "middata/early/hgSub.rds")

# axolotl

axo <- readRDS("middata/pub/axo.rds")

table(axo$time)
Idents(axo) <- axo$time
axoDiff <- FindMarkers(axo, ident.1 = c("Stage28", "Stage40", "Stage44"), ident.2 = "0dpa", min.pct = 0.1, only.pos = F)

axoUpGene <- axoDiff %>% {rownames(.)[.$p_val_adj < 0.05 & .$avg_logFC > 0.25 & .$pct.1 > .$pct.2]}
axoDownGene <- axoDiff %>% {rownames(.)[.$p_val_adj < 0.05 & .$avg_logFC < -0.25 & .$pct.1 < .$pct.2]}
axoDiffGene <- c(axoUpGene, axoDownGene)

an <- axo[["RNA"]]@data[axoDiffGene, axo$time %in% c("Stage28", "Stage40", "Stage44")] %>% as.matrix %>% rowMeans
bn <- axo[["RNA"]]@data[axoDiffGene, axo$time %in% c("0dpa")] %>% as.matrix %>% rowMeans

Dmat <- 2 * matrix(c(sum(an^2), sum(an*bn), sum(an*bn), sum(bn^2)), nrow = 2)

Amat <- matrix(c(1, 0, 0, 1, -1, -1), nrow = 2)
bvec <- c(0, 0, -1)

dData <- apply(as.matrix(axo[["RNA"]]@data[axoDiffGene, ]), 2, function(x) 2 * c(sum(an*x), sum(bn*x)))
sovData <- apply(dData, 2, function(x) (solve.QP(Dmat, x, Amat, bvec))$solution)

axo$score <- sovData[1, ] / colSums(sovData)
VlnPlot(axo, "score", group.by = "time", pt.size = 0)

saveRDS(axo, "middata/pub/axo.rds")

# frog

xenla <- readRDS("middata/pub/xenlaHomo.rds")

Idents(xenla) <- xenla$sample
xenlaDiff <- FindMarkers(xenla, ident.1 = c("s50", "s51", "s52", "s54"), ident.2 = "0dpa", min.pct = 0.1, only.pos = F)

xenlaUpGene <- xenlaDiff %>% {rownames(.)[.$p_val_adj < 0.05 & .$avg_logFC > 0.25 & .$pct.1 > .$pct.2]}
xenlaDownGene <- xenlaDiff %>% {rownames(.)[.$p_val_adj < 0.05 & .$avg_logFC < -0.25 & .$pct.1 < .$pct.2]}
xenlaDiffGene <- c(xenlaUpGene, xenlaDownGene)

an <- xenla[["RNA"]]@data[xenlaDiffGene, xenla$sample %in% c("s50", "s51", "s52", "s54")] %>% as.matrix %>% rowMeans
bn <- xenla[["RNA"]]@data[xenlaDiffGene, xenla$sample %in% c("0dpa")] %>% as.matrix %>% rowMeans

Dmat <- 2 * matrix(c(sum(an^2), sum(an*bn), sum(an*bn), sum(bn^2)), nrow = 2)

Amat <- matrix(c(1, 0, 0, 1, -1, -1), nrow = 2)
bvec <- c(0, 0, -1)

dData <- apply(as.matrix(xenla[["RNA"]]@data[xenlaDiffGene, ]), 2, function(x) 2 * c(sum(an*x), sum(bn*x)))
sovData <- apply(dData, 2, function(x) (solve.QP(Dmat, x, Amat, bvec))$solution)

xenla$score <- sovData[1, ] / colSums(sovData)
VlnPlot(xenla, "score", group.by = "sample", pt.size = 0)

saveRDS(xenla, "middata/pub/xenlaHomo.rds")

# mouse

mm <- readRDS("middata/pub/mmHomo.rds")

Idents(mm) <- mm$time
mmDiff <- FindMarkers(mm, ident.1 = c("E11", "E14"), ident.2 = "0dpa", min.pct = 0.1, only.pos = F)

mmUpGene <- mmDiff %>% {rownames(.)[.$p_val_adj < 0.05 & .$avg_logFC > 0.25 & .$pct.1 > .$pct.2]}
mmDownGene <- mmDiff %>% {rownames(.)[.$p_val_adj < 0.05 & .$avg_logFC < -0.25 & .$pct.1 < .$pct.2]}
mmDiffGene <- c(mmUpGene, mmDownGene)

an <- mm[["RNA"]]@data[mmDiffGene, mm$time %in% c("E11", "E14")] %>% as.matrix %>% rowMeans
bn <- mm[["RNA"]]@data[mmDiffGene, mm$time %in% c("0dpa")] %>% as.matrix %>% rowMeans

Dmat <- 2 * matrix(c(sum(an^2), sum(an*bn), sum(an*bn), sum(bn^2)), nrow = 2)

Amat <- matrix(c(1, 0, 0, 1, -1, -1), nrow = 2)
bvec <- c(0, 0, -1)

dData <- apply(as.matrix(mm23[["RNA"]]@data[mmDiffGene, ]), 2, function(x) 2 * c(sum(an*x), sum(bn*x)))
sovData <- apply(dData, 2, function(x) (solve.QP(Dmat, x, Amat, bvec))$solution)

mm$score <- sovData[1, ] / colSums(sovData)
VlnPlot(mm, "score", group.by = "time", pt.size = 0)

saveRDS(mm, "middata/pub/mmHomo.rds")
