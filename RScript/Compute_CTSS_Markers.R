library(Seurat)
library(SeuratDisk)
library(MetaMarkers)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(scRNAseq)
library(broom)

counts <- Read10X(data.dir = '/gpfs/hpc/home/chenchao/ranm/SCvsSN/SN_Datasets/2019_ASD_Velmeshev/raw_matrix')
ASD_2019 <- CreateSeuratObject(counts, 
                               project = 'ASD_2019',
                               min.cells = 3,
                               min.features =500)
ASD_2019@meta.data <- cbind(ASD_2019@meta.data,meta_data)							   
ASD_2019 <- subset(ASD_2019, subset = region == "PFC" & diagnosis == "Control")
cluster_mapping <- c(
  "AST-FB" = "Ast",
  "AST-PP" = "Ast",
  "Endothelial" = "Endo",
  "IN-PV" = "IN",
  "IN-SST" = "IN",
  "IN-SV2C" = "IN",
  "IN-VIP" = "IN",
  "L2/3" = "EN",
  "L4" = "EN",
  "L5/6" = "EN",
  "L5/6-CC" = "EN",
  "Microglia" = "Mic",
  "Neu-mat" = "EN",
  "Neu-NRGN-I" = "EN",
  "Neu-NRGN-II" = "EN",
  "Oligodendrocytes" = "Olig",
  "OPC" = "OPC"
)
ASD_2019@meta.data <- ASD_2019@meta.data %>%
                      mutate(Celltype = recode(cluster, !!!cluster_mapping))
ASD_2019$Stage <- ifelse(ASD_2019$age < 20, "Adolescence", "Adult")
ASD_2019@meta.data <- ASD_2019@meta.data %>%
  select(-orig.ident, -nCount_RNA, -nFeature_RNA, -cell, -cluster, -individual, -region, -diagnosis)
colnames(ASD_2019@meta.data)[colnames(ASD_2019@meta.data) == "sample"] <- "individual"
SaveH5Seurat(ASD_2019, filename = "/gpfs/hpc/home/chenchao/ranm/SCvsSN/SN_Datasets/2019_ASD_Velmeshev/ASD_2019.h5Seurat")
ASD_2019 <- LoadH5Seurat("/gpfs/hpc/home/chenchao/ranm/SCvsSN/SN_Datasets/2019_ASD_Velmeshev/ASD_2019.h5Seurat", assays ="RNA",slots='counts')
velmeshev_count <- LoadH5Seurat("/gpfs/hpc/home/chenchao/ranm/2nd_3rd_exp/Velmeshev_count.h5seurat", assays ="RNA",slots='counts',meta.data = F, misc = F)
velmeshev_cells_obs <- read.csv("/gpfs/hpc/home/chenchao/ranm/2nd_3rd_exp/velmeshev_cells_obs.csv")
velmeshev_count@meta.data <- velmeshev_cells_obs
rownames(velmeshev_count@meta.data) <- velmeshev_count$X
velmeshev_count <- subset(velmeshev_count, subset = Stage != "Fetus")
GSE157827 <- readRDS("/gpfs/hpc/home/chenchao/pengss/scAD/GSE157827/results/gse157_scCDC_fil_db_anno.rds")
GSE157827 <- subset(GSE157827, subset = diagnosis == 'NC')
GSE157827@meta.data <- GSE157827@meta.data[, c("celltype", "sampleid")]
colnames(GSE157827@meta.data) <- c('Celltype','individual')
GSE157827$Stage <- 'Adult'				
celltype_mapping <- c(
  "Exc" = "EN",
  "Oli" = "Olig",
  "Ast" = "Ast",
  "Inh" = "IN",
  "OPC" = "OPC",
  "Mic" = "Mic",
  "End" = "Endo")
GSE157827@meta.data <- GSE157827@meta.data %>%
                      mutate(Celltype = recode(Celltype, !!!celltype_mapping))
GSE157827_v4 <- CreateSeuratObject(LayerData(GSE157827, layer = "counts"), meta.data = GSE157827@meta.data)
del_cell <- filter(velmeshev_cells_obs, Stage != "Fetus")$X
Stage_count <- Stage_count[,!colnames(Stage_count) %in% del_cell]

common_genes <- reduce(list(velmeshev_count, Wang_2024_count, ASD_2019, GSE157827_v4), 
                      ~ intersect(.x, rownames(.y)), .init = rownames(velmeshev_count))
common_meta <- reduce(list(velmeshev_count, Wang_2024_count, ASD_2019, GSE157827_v4), 
                     ~ intersect(.x, colnames(.y@meta.data)), .init = colnames(velmeshev_count@meta.data))
velmeshev_count <- velmeshev_count[common_genes,]
Wang_2024_count <- Wang_2024_count[common_genes,]
ASD_2019 <- ASD_2019[common_genes,]
GSE157827_v4 <- GSE157827_v4[common_genes,]
Stage_count <- merge(velmeshev_count, y = list(Wang_2024_count, ASD_2019, GSE157827_v4))

Stage_count <- subset(Stage_count, 
  subset = !((Stage == "Adolescent" & Celltype %in% c("IPC-EN", "IPC-Glia", "RG")) |  
             (Stage == "Fetus" & Celltype == "Olig") | 
             (Stage == "Infancy" & Celltype %in% c("IPC-EN", "IPC-Glia"))))
Stage_count <- SCTransform(Stage_count, verbose = FALSE)
Stage_count@meta.data <- Stage_count@meta.data %>%
  select(where(~ !any(is.na(.))))
Stage_count$Stage <- case_when(
  Stage_count$Stage %in% c("adole", "Adolescence", "Adolescent") ~ "Adolescent",
  Stage_count$Stage %in% c("Adult") ~ "Adult",
  Stage_count$Stage %in% c("fetus", "Fetus") ~ "Prenatal",
  Stage_count$Stage %in% c("infancy", "Infancy") ~ "Infancy")  

Adult_marker_gene <- read.csv("/gpfs/hpc/home/chenchao/ranm/SCvsSN/SN_Datasets/2023_Setti/Adult_Stable_Marker_FCsecond.csv", 
                 header = TRUE,
                 row.names = NULL)
Fetus_marker_gene <- read.csv("/gpfs/hpc/home/chenchao/ranm/SCvsSN/SN_Datasets/marker_identify/merge_FC_second_filter.csv", 
                 header = TRUE,
                 row.names = NULL)	

Calculate_Gene_Exp_Prop <- function(seurat_obj, gene, target_celltype) {
  if (!inherits(seurat_obj, "Seurat")) stop("input must be Seurat object")
  gene = intersect(gene, rownames(seurat_obj))
  if (!gene %in% rownames(seurat_obj)) stop("gene is not in Seurat object")
  if (!target_celltype %in% seurat_obj$Celltype) {
    stop("target_celltype is not in Seurat object")
  }
  SLA_Celltype_prop <- seurat_obj@meta.data %>%
    mutate(
      SLA_expr = GetAssayData(seurat_obj, assay = "SCT", slot = "data")[gene, ],
      SLA_count = GetAssayData(seurat_obj, assay = "SCT", slot = "counts")[gene, ],
      Cell_count = 1
    ) %>%
    group_by(individual, Stage, Celltype) %>%
    summarise(
      Gene = gene,
      Cell_count = sum(Cell_count),
      Mean_expression = mean(SLA_expr),
      Expression_proportion = mean(SLA_count > 0),
      .groups = "drop"
    )
  valid_individuals <- SLA_Celltype_prop %>%
    filter(Celltype == target_celltype) %>%
    pull(individual) %>%
    unique()
  result <- SLA_Celltype_prop %>%
    filter(individual %in% valid_individuals) %>%
    group_by(Stage, Celltype) %>%
    mutate(group_mean = mean(Mean_expression)) %>%
    ungroup() %>%
    arrange(individual, ifelse(Celltype == target_celltype, 0, 1), -group_mean) %>%
    group_by(individual) %>%
    slice_head(n = 2) %>%
    mutate(
      Rankcelltype = case_when(
        row_number() == 1 ~ as.character(Celltype),
        row_number() == 2 ~ "Rank-2",
        TRUE ~ NA_character_
      )
    ) %>% as.data.frame() %>%
    ungroup()
  
  attr(result, "gene") <- gene
  attr(result, "target_celltype") <- target_celltype
  
  return(result)
}

analyze_expression_trend <- function(data) {
  data <- data %>%
    mutate(Stage = factor(Stage, 
                         levels = c("Prenatal", "Infancy", "Adolescent", "Adult"),
                         ordered = TRUE))  
  results <- data %>%
    group_by(Gene) %>%
    do({
      model <- lm(Mean_expression ~ as.numeric(Stage), data = .)
      tidy(model) %>%
        filter(term == "as.numeric(Stage)")
    }) %>%
    ungroup() 
  results <- results %>%
    mutate(p.adj = p.adjust(p.value, method = "fdr")) %>% 
    dplyr::select(Gene, estimate, std.error, statistic, p.value, p.adj) %>% # p.adj
    mutate(
      direction = ifelse(estimate > 0, "Increasing", "Decreasing"),
      significance = ifelse(p.adj < 0.05, "Significant", "Non-significant"),
      trend = paste(direction, significance, sep = ", ")
    )
  return(results)
}

analyze_gene_trends <- function(seurat_obj, gene_list, target_celltype) {
  exp_prop <- lapply(gene_list, function(g) {
    Calculate_Gene_Exp_Prop(seurat_obj, g, target_celltype)
  }) %>% bind_rows()
    exp_prop <- exp_prop %>%
    mutate(Stage = factor(Stage, 
                         levels = c("Prenatal", "Infancy", "Adolescent", "Adult"), 
                         ordered = TRUE))
    celltype_data <- exp_prop %>% 
    filter(Rankcelltype == target_celltype) %>%
    dplyr::select(Gene, Stage, individual, Mean_expression)
    trend_results <- analyze_expression_trend(celltype_data)
  return(trend_results)
}
EN_specific <- setdiff(union(Fetus_marker_gene$Gene[Fetus_marker_gene$Celltype=='EN'],
                             Adult_marker_gene$Gene[Adult_marker_gene$Celltype=='EN']),
                       intersect(Fetus_marker_gene$Gene[Fetus_marker_gene$Celltype=='EN'],
                                 Adult_marker_gene$Gene[Adult_marker_gene$Celltype=='EN'])) %>%
                       intersect(rownames(Stage_count))
Ast_specific <- setdiff(union(Fetus_marker_gene$Gene[Fetus_marker_gene$Celltype=='Ast'],
                              Adult_marker_gene$Gene[Adult_marker_gene$Celltype=='Ast']),
                        intersect(Fetus_marker_gene$Gene[Fetus_marker_gene$Celltype=='Ast'],
                                  Adult_marker_gene$Gene[Adult_marker_gene$Celltype=='Ast'])) %>%
                intersect(rownames(Stage_count))
Mic_specific <- setdiff(union(Fetus_marker_gene$Gene[Fetus_marker_gene$Celltype=='Mic'],
                              Adult_marker_gene$Gene[Adult_marker_gene$Celltype=='Mic']),
                        intersect(Fetus_marker_gene$Gene[Fetus_marker_gene$Celltype=='Mic'],
                                  Adult_marker_gene$Gene[Adult_marker_gene$Celltype=='Mic'])) %>%
                intersect(rownames(Stage_count))

OPC_specific <- setdiff(union(Fetus_marker_gene$Gene[Fetus_marker_gene$Celltype=='OPC'],
                              Adult_marker_gene$Gene[Adult_marker_gene$Celltype=='OPC']),
                        intersect(Fetus_marker_gene$Gene[Fetus_marker_gene$Celltype=='OPC'],
                                  Adult_marker_gene$Gene[Adult_marker_gene$Celltype=='OPC'])) %>%
                intersect(rownames(Stage_count))
IN_specific <- setdiff(union(Fetus_marker_gene$Gene[Fetus_marker_gene$Celltype=='IN'],
                             Adult_marker_gene$Gene[Adult_marker_gene$Celltype=='IN']),
                       intersect(Fetus_marker_gene$Gene[Fetus_marker_gene$Celltype=='IN'],
                                 Adult_marker_gene$Gene[Adult_marker_gene$Celltype=='IN'])) %>%
               intersect(rownames(Stage_count))

Ast_trend_45 <- analyze_gene_trends(
  seurat_obj = Stage_count,
  gene_list = Ast_specific,
  target_celltype = "Ast"
)
EN_trend_32 <- analyze_gene_trends(
  seurat_obj = Stage_count,
  gene_list = EN_specific,
  target_celltype = "EN"
)
OPC_trend_22 <- analyze_gene_trends(
  seurat_obj = Stage_count,
  gene_list = OPC_specific,
  target_celltype = "OPC"
)
IN_trend_12 <- analyze_gene_trends(
  seurat_obj = Stage_count,
  gene_list = IN_specific,
  target_celltype = "IN"
)
Mic_trend_53 <- analyze_gene_trends(
  seurat_obj = Stage_count,
  gene_list = Mic_specific,
  target_celltype = "Mic"
)
Trend_164 <- rbind(Ast_trend_45, EN_trend_32, OPC_trend_22, IN_trend_12, Mic_trend_53)



