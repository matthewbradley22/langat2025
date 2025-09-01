#Manual check of labelled doublets
library(Seurat)

#Source useful functions
source('./scripts/langatFunctions.R')

#Load in data
ParseSeuratObj_int <- LoadSeuratRds("./data/FilteredRpcaIntegratedDat.rds")
DimPlot(ParseSeuratObj_int, group.by = 'scDblFinderLabel', reduction = "umap.integrated")

DimPlot(ParseSeuratObj_int, group.by = 'singleR_labels', reduction = "umap.integrated",
        label = TRUE)

#Add manual annotations from end of LangatCellAnnotations script

#Start with microglia, macrophages
micro_mac <- subset(ParseSeuratObj_int, manualAnnotation == 'Macrophage/Monocytes' |
                      manualAnnotation == 'Microglia')
table(micro_mac$Organ, micro_mac$scDblFinderLabel)
micro_mac <- prepSeuratObj(micro_mac)
ElbowPlot(micro_mac, ndims = 30)
micro_mac <- prepUmapSeuratObj(micro_mac, 20, reductionName = 'subsetUMAP_20')

DimPlot(micro_mac, reduction = 'subsetUMAP_20', group.by = 'scDblFinderLabel')
DimPlot(micro_mac, reduction = 'subsetUMAP_20')
DimPlot(micro_mac, reduction = 'subsetUMAP_20', group.by = 'manualAnnotation')


