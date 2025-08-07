#CELL ANNOTATION 
#Load in libraries
library(Seurat)
library(SingleR)
library(celldex)
library(pheatmap)

#Load in data
#ParseSeuratObj_int <- LoadSeuratRds('./ccaIntegratedDat.rds')
ParseSeuratObj_int <- LoadSeuratRds("./FilteredRpcaIntegratedDat.rds")
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

####Data is loaded in with doublets, so can start here ####
#Look at doublets in integrated data
DimPlot(ParseSeuratObj_int, group.by = 'scDblFinderLabel')
ParseSeuratObj_int <- subset(ParseSeuratObj_int, scDblFinderLabel == 'singlet')

#Plot viral counts vs well/treatment. 
#removing doublets doesn't get rid of viral contamination in PBS samples
ggplot(ParseSeuratObj_int[[]], aes(x = orig.ident, y = virusCount, col = Treatment))+
  geom_point() +
  scale_x_discrete(labels= wellMap$well)+
  theme(axis.text.x = element_text(angle = 90))

ParseSeuratObj_int[[]] %>% mutate(virusPresent = ifelse(virusCount>0, 'yes', 'no')) %>% 
  group_by(Treatment, virusPresent) %>% dplyr::summarise(count = n()) %>% 
  ggplot(aes(x = Treatment, y = count, fill = virusPresent))+
  geom_bar(stat = 'identity', position = 'dodge')

#### Manual annotation ####

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


#T cells
FeaturePlot(ParseSeuratObj_int, 'Cd3g')

DimPlot(ParseSeuratObj_int, label = TRUE)

#Look at specific clusters top markers to confirm cell types

#Manual annotation
clusters <- ParseSeuratObj_int$seurat_clusters
ParseSeuratObj_int$markerBasedAnnotation <- case_when(clusters == 27 ~ 'neurons',
                                                      )




