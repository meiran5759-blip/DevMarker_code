library(Seurat)
library(SeuratDisk)
library(MetaMarkers)
library(dplyr)
library(ggplot2)
library(tidyr)
library(singleCellNet)


calculate_celltype_expression <- function(seurat_obj, gene_list, target_celltype) {
  if (!inherits(seurat_obj, "Seurat")) stop("Input must be a Seurat object")
  if (!"RNA" %in% names(seurat_obj@assays)) stop("RNA assay not found")
  if (!target_celltype %in% seurat_obj$Celltype) {
    stop("Target celltype '", target_celltype, "' not found in metadata")
  }
  
  counts <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")
  meta <- seurat_obj@meta.data

  if (!identical(colnames(counts), rownames(meta))) {
    warning("Aligning counts and metadata - some cells may be dropped")
    common_cells <- intersect(colnames(counts), rownames(meta))
    counts <- counts[, common_cells]
    meta <- meta[common_cells, ]
  }
  
  gene_list <- intersect(gene_list, rownames(counts))
  if (length(gene_list) == 0) stop("No valid genes found")
  
  valid_cells <- !is.na(meta$Celltype)
  counts <- counts[, valid_cells]
  meta <- meta[valid_cells, ]
  cell_indices <- split(seq_len(ncol(counts)), meta$Celltype)
  
  result_df <- do.call(cbind, lapply(cell_indices, function(idx) {
    rowMeans(counts[gene_list, idx, drop = FALSE] > 0)
  }))
  
  final_result <- lapply(gene_list, function(gene) {
    target_prop <- result_df[gene, target_celltype]
    
    other_celltypes <- setdiff(colnames(result_df), target_celltype)
    second_celltype <- other_celltypes[which.max(result_df[gene, other_celltypes])]
    second_prop <- result_df[gene, second_celltype]
    
    target_mean <- mean(counts[gene, cell_indices[[target_celltype]]])
    second_mean <- mean(counts[gene, cell_indices[[second_celltype]]])
    fc_second <- target_mean / max(second_mean, 1e-10)  # 避免除以0
    
    data.frame(
      Gene = gene,
      Target_Celltype = target_celltype,
      Second_Celltype = second_celltype,
      Target_prop = target_prop,
      Second_prop = second_prop,
      FC_second = fc_second,
      stringsAsFactors = FALSE
    )
  }) %>% bind_rows() 
  return(final_result)
}
Ast_FCsecond <- calculate_celltype_expression(seurat_obj = merge_sn,
                                             gene_list = Ast_COSG_CV,
                                             target_celltype = "Ast")
RG_FCsecond <- calculate_celltype_expression(seurat_obj = merge_sn,
                                             gene_list = RG_COSG_CV,
                                             target_celltype = "RG")
Glioblast_FCsecond <- calculate_celltype_expression(seurat_obj = merge_sn,
                                             gene_list = Glioblast_COSG_CV,
                                             target_celltype = "Glioblast")
OPC_FCsecond <- calculate_celltype_expression(seurat_obj = merge_sn,
                                             gene_list = OPC_COSG_CV,
                                             target_celltype = "OPC")
IPC_FCsecond <- calculate_celltype_expression(seurat_obj = merge_sn,
                                             gene_list = IPC_COSG_CV,
                                             target_celltype = "IPC")                                             
Mic_FCsecond <- calculate_celltype_expression(seurat_obj = merge_sn,
                                             gene_list = Mic_COSG_CV,
                                             target_celltype = "Mic") 
EN_FCsecond <- calculate_celltype_expression(seurat_obj = merge_sn,
                                             gene_list = EN_COSG_CV,
                                             target_celltype = "EN")                                                                                             
IN_FCsecond <- calculate_celltype_expression(seurat_obj = merge_sn,
                                             gene_list = IN_COSG_CV,
                                             target_celltype = "IN")                                                                                         
merge_FC_second <- rbind(RG_FCsecond,Glioblast_FCsecond,Ast_FCsecond,OPC_FCsecond,
                         IPC_FCsecond,Mic_FCsecond,EN_FCsecond,IN_FCsecond)
colnames(merge_FC_second)[4] <- "Target_prop" 
merge_FC_second_filter <- merge_FC_second %>%
                   filter(Target_prop > 0.4, Second_prop < 0.25, FC_second > 2)
