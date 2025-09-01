#Functions script to source from
#packages
library(tibble)
library(dplyr)
library(DESeq2)
#Pseudobulk info found here
#https://hbctraining.github.io/Pseudobulk-for-scRNAseq/schedule/self-learning.html 

createPseudoBulk <- function(data, variables){
  
  meta_columns <- variables
  meta <- data[[]] %>%
    select(all_of(meta_columns)) %>%
    unique() %>%
    remove_rownames()
  
  bulk <- AggregateExpression(
    data,
    return.seurat = TRUE,
    assays = "RNA",
    group.by = meta_columns
  )
  
  n_cells <- data[[]] %>% 
    dplyr::group_by(across(all_of(meta_columns))) %>% 
    dplyr::summarise('n_cells' = n())
  
  meta_bulk <- left_join(bulk[[]], n_cells)
  rownames(meta_bulk) <- meta_bulk$orig.ident
  bulk[[]] <- meta_bulk
  
  # Get count matrix
  cluster_counts <- FetchData(bulk, layer="counts", vars=rownames(bulk))
  
  designForm <- reformulate(variables)
  dds <- DESeqDataSetFromMatrix(t(cluster_counts),
                                colData = bulk@meta.data,
                                design = designForm)
}

#Reprocess a subset seurat object (rerun normalizing, scaling etc...)
prepSeuratObj <- function(obj){
  #Rerun through data processing and visualization
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj)
  obj <- ScaleData(obj)
  obj <- RunPCA(obj)
  obj
}

prepUmapSeuratObj <- function(obj, nDims, reductionName){
  obj <- FindNeighbors(obj, dims = 1:nDims, reduction = "pca")
  obj <- FindClusters(obj, resolution = 2)  
  obj <- RunUMAP(obj, dims = 1:nDims, reduction = "pca", reduction.name = reductionName)
  obj
}



