library(Seurat)
library(SeuratDisk)
library(dplyr)
library(ggplot2)
library(tidyr)
library(scRNAseq)
source('/gpfs/hpc/home/chenchao/ranm/SCvsSN/SN_Datasets/cv_sc.r')


Convert("/gpfs/hpc/home/chenchao/ranm/SCvsSN/SN_Datasets/merge_sn_count.h5ad",dest = "h5seurat", overwrite = TRUE, assay = "RNA")
merge_sn_count <- LoadH5Seurat("/gpfs/hpc/home/chenchao/ranm/SCvsSN/SN_Datasets/merge_sn_count.h5seurat", assays ="RNA",slots='counts',meta.data = T, misc = T)
meta <- read.csv("/gpfs/hpc/home/chenchao/ranm/SCvsSN/SN_Datasets/merge_sn/merge_sn_count_obs.csv")
rownames(meta) <- meta$X
merge_sn_count@meta.data <- meta
Ramos_meta <- meta[which(meta$Dataset== 'Ramos'),]
counts_sn <- GetAssayData(object = merge_sn, assay = "RNA", slot = "counts")
Ramos_counts <- counts_sn[, colnames(counts_sn) %in% rownames(Ramos_meta)]
Ramos_meta_class <- Ramos_meta[,c('individual','Celltype')] # extract sample id and cell identity from metadata
colnames(Ramos_meta_class) <- c('sample','celltype')
Ramos_meta_class$sample <- "Group"
Ramos_cvlist <- get_cv_for_replicates(Ramos_meta_class,Ramos_counts)#calculate CV value for each gene, each cell type, in each sample
Ramos_stable <- get_sample(Ramos_meta_class)#get the sample id which has the largest cell numbers
Ramos_perm <- get_perm(Ramos_meta_class,Ramos_counts,Ramos_cvlist,100)
Ramos_celltype_cv <- list()
for(i in 1:length(Ramos_perm)){
    gene_cv <- na.omit(rowMeans(Ramos_perm[[i]]))
    Ramos_celltype_cv[[i]] <- gene_cv
}
names(Ramos_celltype_cv) <- names(Ramos_perm)
Ramos_Astcv <- names(Ramos_celltype_cv[[1]])[which(Ramos_celltype_cv[[1]] < 0.1)]
Ramos_ENcv <- names(Ramos_celltype_cv[[2]])[which(Ramos_celltype_cv[[2]] < 0.1)]
Ramos_Glioblastcv <- names(Ramos_celltype_cv[[3]])[which(Ramos_celltype_cv[[3]] < 0.1)]
Ramos_INcv <- names(Ramos_celltype_cv[[4]])[which(Ramos_celltype_cv[[4]] < 0.1)]
Ramos_IPCcv <- names(Ramos_celltype_cv[[5]])[which(Ramos_celltype_cv[[5]] < 0.1)]
Ramos_Miccv <- names(Ramos_celltype_cv[[6]])[which(Ramos_celltype_cv[[6]] < 0.1)]
Ramos_OPCcv <- names(Ramos_celltype_cv[[7]])[which(Ramos_celltype_cv[[7]] < 0.1)]
Ramos_RGcv <- names(Ramos_celltype_cv[[8]])[which(Ramos_celltype_cv[[8]] < 0.1)]

