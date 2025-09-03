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
#Check for any cells expressing male and female genes
sexGenes <- ParseSeuratObj_int[['RNA']]$data[c('Xist', 'Eif2s3y'),]
sexGenePresence <- colSums(sexGenes > 0)
ParseSeuratObj_int$sexGenePresence <- case_when(sexGenePresence == 0 ~ 'None',
                                                sexGenePresence == 1 ~ 'One',
                                                sexGenePresence == 2 ~ 'Two')

table(subset(ParseSeuratObj_int, sexGenePresence == 'Two')$seurat_clusters) %>% sort()

#Think cluster 14 is doublets, identified by scdblfindr and many cells with both sex chromosomes


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


