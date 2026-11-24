library(Seurat)
library(dplyr)
library(ggplot2)

#Load in data from 01_data_preprocessing.R script
ParseSeuratObj_int <- LoadSeuratRds("~/Documents/ÖverbyLab//data/FilteredRpcaIntegratedDat.rds")

#Look at groups of doublets
DimPlot(ParseSeuratObj_int, group.by = 'scDblFinderLabel', reduction = "umap.integrated")

#Look at clusters with highest scDblFinder percentage
ParseSeuratObj_int[[]] %>% group_by(seurat_clusters, scDblFinderLabel) %>% dplyr::summarise(n = n()) %>%
  dplyr::mutate(freq = n / sum(n)) %>% 
  dplyr::filter(scDblFinderLabel == 'doublet') %>% 
  dplyr::arrange(desc(freq)) %>% 
  ggplot(aes(x = seurat_clusters, y = freq))+
  geom_bar(stat = 'identity')

#Manual removal of clusters with a) high proprtion of doublets from scDblFinder or
# b) expression of marker genes from distinct cell types

#16 have high proportion of doublets and shows both astrocyte + microglia markers
ParseSeuratObj_int <- subset(ParseSeuratObj_int, seurat_clusters != 16)

#32 high proportion of doublets and shows astrocyte + endothelial markers
ParseSeuratObj_int <- subset(ParseSeuratObj_int, seurat_clusters != 32)

#19 high proportion of doublets and shows astrocyte + oligo + microglia markers
ParseSeuratObj_int <- subset(ParseSeuratObj_int, seurat_clusters != 19)

#Getting rid of 27 for having both endothelial and pericyte markers
ParseSeuratObj_int <- subset(ParseSeuratObj_int, seurat_clusters != 27)

#Remove cluster 38 - includes markers that should be distinct to macrophages or microglia
ParseSeuratObj_int <- subset(ParseSeuratObj_int, seurat_clusters != 38)

#SaveSeuratRds(ParseSeuratObj_int, "./data/FilteredRpcaIntegratedDatNoDoublets.rds")