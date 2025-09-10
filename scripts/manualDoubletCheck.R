#Manual check of labelled doublets
library(Seurat)

#Source useful functions
source('./scripts/langatFunctions.R')

#Load in data. Have now removed genes w both Xist and Eif2s3y, so must reload
#data from beginning to see which were identified
ParseSeuratObj_int <- LoadSeuratRds("./data/FilteredRpcaIntegratedDat.rds")
DimPlot(ParseSeuratObj_int, group.by = 'scDblFinderLabel', reduction = "umap.integrated")

DimPlot(ParseSeuratObj_int, group.by = 'singleR_labels', reduction = "umap.integrated",
        label = TRUE)

#Add manual annotations from end of LangatCellAnnotations script

#Think cluster 14 is doublets, identified by scdblfindr and many cells with both sex chromosomes

#Look at clusters with highest scDblFinder percentage
ParseSeuratObj_int[[]] %>% group_by(seurat_clusters, scDblFinderLabel) %>% dplyr::summarise(n = n()) %>%
  dplyr::mutate(freq = n / sum(n)) %>% 
  dplyr::filter(scDblFinderLabel == 'doublet') %>% 
  dplyr::arrange(desc(freq))

DimPlot(ParseSeuratObj_int, reduction = "umap.integrated",
        label = TRUE)

#16 have high proportion of doublets and shows both astrocyte + microglia markers. Removing
ParseSeuratObj_int <- subset(ParseSeuratObj_int, seurat_clusters != 16)

#32 high proportion of doublets and shows astrocyte + endothelial markers
ParseSeuratObj_int <- subset(ParseSeuratObj_int, seurat_clusters != 32)

#19 high proportion of doublets and shows astrocyte + oligo + microglia markers
ParseSeuratObj_int <- subset(ParseSeuratObj_int, seurat_clusters != 19)

#Getting rid of 27 for having both endothelial and pericyte markers
ParseSeuratObj_int <- subset(ParseSeuratObj_int, seurat_clusters != 27)

#Remove cluster 38 - includes genes that Allen cell atlas says should only be expressed
#in either microglia or macrophages, not both
ParseSeuratObj_int <- subset(ParseSeuratObj_int, seurat_clusters != 38)

#Check feature counts across data
FeaturePlot(ParseSeuratObj_int, 'nFeature_RNA', reduction = 'umap.integrated')

#SaveSeuratRds(ParseSeuratObj_int, "./data/FilteredRpcaIntegratedDatNoDoublets.rds")
