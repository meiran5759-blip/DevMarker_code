library(Seurat)
library(data.table)
library(dplyr)
library(tidyr)
library(ggpubr)
library(stringr)
library(SeuratDisk)  
library(Matrix)
library(DoubletFinder)
library(harmony)
library(clustree)
library(tibble)
library(RColorBrewer)
library(writexl)
library(parallel)
library(pbapply)
library(future)
library(future.apply)

Convert("/gpfs/hpc/home/chenchao/ranm/SCvsSN/SN_Datasets/merge_sn/merge_sn_count.h5ad",dest = "h5seurat", overwrite = TRUE, assay = "RNA")
merge_sn_count <- LoadH5Seurat("/gpfs/hpc/home/chenchao/ranm/SCvsSN/SN_Datasets/merge_sn/merge_sn_count.h5seurat", assays ="RNA",slots='counts',meta.data = T, misc = T)
metadata <- read.csv("/gpfs/hpc/home/chenchao/ranm/SCvsSN/SN_Datasets/merge_sn/merge_sn_count_obs.csv", row.names = 1)
merge_sn_count@meta.data <- metadata 
merge_sn_count <- NormalizeData(merge_sn_count)
merge_sn_count <- FindVariableFeatures(merge_sn_count, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(merge_sn_count)
merge_sn_count <- ScaleData(merge_sn_count, features = all.genes)
merge_sn_count <- RunPCA(merge_sn_count, features = VariableFeatures(object = merge_sn_count))
####DoubletFinder
dim.usage <- 50
merge_sn_count <- RunUMAP(merge_sn_count, dims = 1:50)
merge_sn_count <- FindNeighbors(merge_sn_count, dims = 1:50)
merge_sn_count <- FindClusters(merge_sn_count, resolution = 0.8)
sweep.res.list <- paramSweep(merge_sn_count, PCs = 1:dim.usage, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pK_bcmvn <- as.numeric(bcmvn$pK[which.max(bcmvn$BCmetric)]) 
homotypic.prop <- modelHomotypic(merge_sn_count$seurat_clusters) 
nExp_poi <- round(0.08*ncol(merge_sn_count))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
p <-as.numeric(as.vector(bcmvn[bcmvn$MeanBC==max(bcmvn$MeanBC),]$pK))
merge_sn_count <- doubletFinder(merge_sn_count, PCs = 1:dim.usage, pN = 0.25, pK = p, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
colnames(merge_sn_count@meta.data)[ncol(merge_sn_count@meta.data)] = "doublet_info"
seu_list <- SplitObject(merge_sn_count, split.by = "Dataset")  
for (dataset_name in names(seu_list)) {
    assign(paste0("seu_", dataset_name), seu_list[[dataset_name]])
}

plan(multicore, workers = 18)
task <- function(seu) {
  library(Seurat)
  seu <- RunUMAP(seu, dims = 1:50)
  seu <- FindNeighbors(seu, dims = 1:50)
  seu <- FindClusters(seu, resolution = 0.8)
  sweep.res.list <- paramSweep(seu, PCs = 1:dim.usage, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pK_bcmvn <- as.numeric(bcmvn$pK[which.max(bcmvn$BCmetric)]) 
  homotypic.prop <- modelHomotypic(seu$seurat_clusters) 
  nExp_poi <- round(0.08*ncol(seu))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop)) 
  p <-as.numeric(as.vector(bcmvn[bcmvn$MeanBC==max(bcmvn$MeanBC),]$pK))
  seu <- doubletFinder(seu, PCs = 1:dim.usage, pN = 0.25, pK = p, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  colnames(seu@meta.data)[ncol(seu@meta.data)] = "doublet_info"
  cat("Task 1 (seu) has completed.\n")
  return(seu)
}
seurat_objects <- list(
  seu_cameron = seu_cameron,
  seu_Velmeshev = seu_Velmeshev,
  seu_Trevino = seu_Trevino,
  seu_Herring = seu_Herring,
  seu_Ramos = seu_Ramos,
  seu_Zhu = seu_Zhu
)
results <- future_lapply(seurat_objects, task, future.seed = TRUE)
seu_cameron <- results[[1]]
seu_Velmeshev <- results[[2]]
seu_Trevino <- results[[3]]
seu_Herring <- results[[4]]
seu_Ramos <- results[[5]]
seu_Zhu <- results[[6]]
merge_sn_count <- Reduce(function(x, y) merge(x, y),results)
plan(sequential)
doublet_cells <- rownames(merge_sn_count@meta.data)[which(merge_sn_count@meta.data$doublet_info == 'Doublet')]
write.table(doublet_cells, file = "/gpfs/hpc/home/chenchao/ranm/SCvsSN/SN_Datasets/merge_sn/Split_sn_Doublet_cells.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

merge_sn_count <-subset(merge_sn_count, subset=doublet_info=="Singlet")
del_columns <- grep("pANN_",colnames(merge_sn_count@meta.data))
merge_sn_count@meta.data <- merge_sn_count@meta.data[,-del_columns]
merge_sn_count <- RunHarmony(merge_sn_count, [group.by](http://group.by/) = "Dataset", plot_convergence = TRUE)
merge_sn_count <- RunUMAP(merge_sn_count, reduction = "harmony", dims = 1:50)
merge_sn_count <- FindNeighbors(merge_sn_count, reduction = "harmony", dims = 1:50)
merge_sn_count <- FindClusters(merge_sn_count, resolution = 0.8)



                         
