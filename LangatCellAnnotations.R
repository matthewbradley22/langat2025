#CELL ANNOTATION 
#Load in libraries
library(Seurat)
library(SingleR)
library(celldex)
library(pheatmap)

#Load in data
ParseSeuratObj_int <- LoadSeuratRds('./ccaIntegratedDat.rds')

####Annotate cell types ####
#Will try automatic annotation
#Mouse rna seq reference
ref <- fetchReference("mouse_rnaseq", "2024-02-26")
pred <- SingleR(test=ParseSeuratObj_int[['RNA']]$data, ref=ref, labels=ref$label.main)
ParseSeuratObj_int$singleR_labels <- pred$labels
plotScoreHeatmap(pred)

#UMAP with cell assignments
DimPlot(ParseSeuratObj_int, group.by = 'singleR_labels', label = TRUE)
clusterAssignments <- table(Assigned=pred$pruned.labels, Cluster= ParseSeuratObj_int$seurat_clusters)
pheatmap(log2(clusterAssignments+10), color=colorRampPalette(c("white", "blue"))(101))

#### Look at canonical gene markers for annotation ####
#Paper has some canonical markers for cell types
#https://umu.diva-portal.org/smash/get/diva2:1897514/FULLTEXT01.pdf

#This paper https://link.springer.com/article/10.1007/s10571-021-01159-3
#Heterogeneity and Molecular Markers for CNS Glial Cells Revealed by Single-Cell Transcriptomics also 
#has lists of markers

#Some neuron markers
FeaturePlot(ParseSeuratObj_int, 'Snap25')
FeaturePlot(ParseSeuratObj_int, 'Pcp2')
FeaturePlot(ParseSeuratObj_int, 'Rbfox3')
FeaturePlot(ParseSeuratObj_int, 'Dpp10')
FeaturePlot(ParseSeuratObj_int, 'Syt1')


#Endothelial cells
FeaturePlot(ParseSeuratObj_int, 'Flt1')

#Astrocytes
FeaturePlot(ParseSeuratObj_int, 'Gfap')
FeaturePlot(ParseSeuratObj_int, 'Aqp4')
FeaturePlot(ParseSeuratObj_int, 'Fgfr3')
FeaturePlot(ParseSeuratObj_int, 'Aldh1l1')
FeaturePlot(ParseSeuratObj_int, 'Slc1a3')

#Microglia
FeaturePlot(ParseSeuratObj_int, 'Ctss')
FeaturePlot(ParseSeuratObj_int, 'Csf1r')
FeaturePlot(ParseSeuratObj_int, 'Cx3cr1')
FeaturePlot(ParseSeuratObj_int, 'C1qa')
FeaturePlot(ParseSeuratObj_int, 'Tmem119')
FeaturePlot(ParseSeuratObj_int, 'P2ry12')

#Oligo
FeaturePlot(ParseSeuratObj_int, 'Mag')
FeaturePlot(ParseSeuratObj_int, 'Mog')

#OPC
FeaturePlot(ParseSeuratObj_int, 'Pdgfra')
FeaturePlot(ParseSeuratObj_int, 'Cspg4')

#Macrophage markers
FeaturePlot(ParseSeuratObj_int, 'Ptprc')
FeaturePlot(ParseSeuratObj_int, 'Ccr2')
FeaturePlot(ParseSeuratObj_int, 'Fn1')

#T cells
FeaturePlot(ParseSeuratObj_int, 'Cd3g')

#### Manual annotation ####
#Look at specific clusters top markers to confirm cell types
#Cluster 0 upregulated w Chil3, Ms4a8a, Saa3, GM15056. Looks like macrophages
possibleMacrophages <- FindMarkers(ParseSeuratObj_int, ident.1 = 0, group.by = 'seurat_clusters', only.pos = TRUE)

#Cluster 7: Adgre4, MS4a genes (monocyte?)... This one more confusing
possibleMacrophages2 <- FindMarkers(ParseSeuratObj_int, ident.1 = 7, group.by = 'seurat_clusters', only.pos = TRUE)

#Cluster 26 (granulocytes?): Retnlg, s100a9, s100a8, Mrgpra2a
possibleGranulocytes <- FindMarkers(ParseSeuratObj_int, ident.1 = 26, 
                                    group.by = 'seurat_clusters', only.pos = TRUE)

#Interesting cluster w Foxb1 near the top, but not a ton of astrocyte markers afaik
possibleAstrocytes <- FindMarkers(ParseSeuratObj_int, ident.1 = 8, 
                                    group.by = 'seurat_clusters', only.pos = TRUE)

possibleAstrocytes2 <- FindMarkers(ParseSeuratObj_int, ident.1 = 19, 
                                  group.by = 'seurat_clusters', only.pos = TRUE)
#Manual annotation
clusters <- ParseSeuratObj_int$seurat_clusters
ParseSeuratObj_int$markerBasedAnnotation <- case_when(clusters == 0 ~ 'macrophages',
                                                      cluster == 26 ~ 'granulocytes')




