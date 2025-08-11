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
  
  dds <- DESeqDataSetFromMatrix(t(cluster_counts),
                                colData = bulk@meta.data,
                                design = ~ Genotype + Treatment + Timepoint + isInfected)
}


