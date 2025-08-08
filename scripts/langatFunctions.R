#Functions script to source from
#packages
library(tibble)

#Pseudobulk info found here
#https://hbctraining.github.io/Pseudobulk-for-scRNAseq/schedule/self-learning.html 

createPseudoBulk <- function(data, variables){
  
}


meta_columns <- c("Genotype", "Treatment", "Timepoint", "Organ", "isInfected")
meta <- cerebellumObj[[]] %>%
  select(all_of(meta_columns)) %>%
  unique() %>%
  remove_rownames()

bulk <- AggregateExpression(
  cerebellumObj,
  return.seurat = TRUE,
  assays = "RNA",
  group.by = c("Genotype", "Treatment", "Timepoint", "isInfected")
)

n_cells <- cerebellumObj[[]] %>% 
  dplyr::count(Genotype, Treatment, Timepoint, isInfected) %>% 
  rename("n_cells"="n")

bulk$n_cells <- n_cells$n_cells
meta_bulk <- left_join(bulk[[]], n_cells)
rownames(meta_bulk) <- meta_bulk$orig.ident
bulk@meta.data <- meta_bulk
