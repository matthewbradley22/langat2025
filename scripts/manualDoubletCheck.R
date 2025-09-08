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
ParseSeuratObj_int[[]] %>% group_by(seurat_clusters, scDblFinderLabel) %>% summarise(n = n()) %>%
  mutate(freq = n / sum(n)) %>% 
  filter(scDblFinderLabel == 'doublet') %>% 
  arrange(desc(freq))

DimPlot(ParseSeuratObj_int, reduction = "umap.integrated",
        label = TRUE)

#16 have high proportion of doublets and shows both astrocyte + microglia markers. Removing
ParseSeuratObj_int <- subset(ParseSeuratObj_int, seurat_clusters != 16)

#32 high proportion of doublets and shows astrocyte + endothelial markers
ParseSeuratObj_int <- subset(ParseSeuratObj_int, seurat_clusters != 32)


#Evaluate microglia, macrophages
micro_mac <- subset(ParseSeuratObj_int, manualAnnotation == 'Macrophage/Monocytes' |
                      manualAnnotation == 'Microglia' | seurat_clusters == 34)
table(micro_mac$Organ, micro_mac$scDblFinderLabel)
micro_mac <- prepSeuratObj(micro_mac)
ElbowPlot(micro_mac, ndims = 30)
micro_mac <- prepUmapSeuratObj(micro_mac, 25, reductionName = 'subsetUMAP_25')

DimPlot(micro_mac, reduction = 'subsetUMAP_25', group.by = 'scDblFinderLabel')
DimPlot(micro_mac, reduction = 'subsetUMAP_25')
DimPlot(micro_mac, reduction = 'subsetUMAP_25', group.by = 'manualAnnotation')

FeaturePlot(micro_mac, 'Hexb', reduction = 'subsetUMAP_25')
FeaturePlot(micro_mac, 'Inpp5d', reduction = 'subsetUMAP_25')
FeaturePlot(micro_mac, 'Ms4a4a', reduction = 'subsetUMAP_25')
FeaturePlot(micro_mac, 'Cd74', reduction = 'subsetUMAP_25')
FeaturePlot(micro_mac, 'Cd209a', reduction = 'subsetUMAP_25')
FeaturePlot(micro_mac, 'Tmem119', reduction = 'subsetUMAP_25')
