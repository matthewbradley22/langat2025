#CELL ANNOTATION 
#Load in libraries
library(Seurat)
library(SingleR)
library(celldex)
library(pheatmap)
library(dplyr)

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
DimPlot(ParseSeuratObj_int, group.by = 'singleR_labels', reduction = "umap.integrated",
        label = TRUE)
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

#Look at neo1 across genotypes
ParseSeuratObj_int$neo_exp <- ParseSeuratObj_int[['RNA']]$data['Neo1',]
ggplot(ParseSeuratObj_int[[]], aes(x = Well, y = neo_exp))+
  geom_point() +
  theme(axis.text.x = element_text(angle = 90))


#Count of virus presence
ParseSeuratObj_int[[]] %>% mutate(virusPresent = ifelse(virusCount>0, 'yes', 'no')) %>% 
  group_by(Treatment, virusPresent) %>% dplyr::summarise(count = n()) %>% 
  ggplot(aes(x = Treatment, y = count, fill = virusPresent))+
  geom_bar(stat = 'identity', position = 'dodge')

#Plot of virus presence over small amount
ParseSeuratObj_int[[]] %>% mutate(virusHight = ifelse(virusCount>3, 'yes', 'no')) %>% 
  group_by(Treatment, virusHight) %>% dplyr::summarise(count = n()) %>% 
  ggplot(aes(x = Treatment, y = count, fill = virusHight))+
  geom_bar(stat = 'identity', position = 'dodge')

#Look at virus presence across celltypes
ParseSeuratObj_int[[]] <- ParseSeuratObj_int[[]] %>% mutate(virusPresence = ifelse(virusCount > 0, 'yes', 'no')) 
  DimPlot(ParseSeuratObj_int, group.by = 'virusPresence')

ParseSeuratObj_int[[]] %>% group_by(singleR_labels, virusPresence) %>% dplyr::summarise(count = n()) %>% 
  ggplot(aes(x = singleR_labels, y = count, fill = virusPresence))+
  geom_bar(stat = 'identity', position = 'dodge')+
  theme(axis.text.x = element_text(angle = 90))

#Look at cell counts per well
table(ParseSeuratObj_int$Well) %>% barplot()

#Most pbs cells with virus have very low amount, could just filter by amount
pbsWithVirus <- subset(ParseSeuratObj_int, Treatment == 'PBS' & virusCountPAdj > 0)
summary(pbsWithVirus$virusCountPAdj)
hist(pbsWithVirus$virusCountPAdj)

#### Manual annotatTimepoint#### Manual annotation ####

#Paper has some canonical markers for cell types
#https://umu.diva-portal.org/smash/get/diva2:1897514/FULLTEXT01.pdf

#This paper https://link.springer.com/article/10.1007/s10571-021-01159-3
#Heterogeneity and Molecular Markers for CNS Glial Cells Revealed by Single-Cell Transcriptomics also 
#has lists of markers

#Some neuron markers
FeaturePlot(ParseSeuratObj_int, 'Snap25', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Pcp2', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Rbfox3', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Dpp10', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Syt1', reduction = 'umap.integrated')


#Endothelial cells
FeaturePlot(ParseSeuratObj_int, 'Flt1', reduction = 'umap.integrated')

#Astrocytes
FeaturePlot(ParseSeuratObj_int, 'Gfap', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Aqp4', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Fgfr3', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Aldh1l1', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Slc1a3', reduction = 'umap.integrated')

#Microglia
FeaturePlot(ParseSeuratObj_int, 'Ctss', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Csf1r', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Cx3cr1', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'C1qa', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Tmem119', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'P2ry12', reduction = 'umap.integrated')

#Oligo
FeaturePlot(ParseSeuratObj_int, 'Mag',reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Mog', reduction = 'umap.integrated')

#OPC
FeaturePlot(ParseSeuratObj_int, 'Pdgfra', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Cspg4', reduction = 'umap.integrated')

#Macrophage markers
FeaturePlot(ParseSeuratObj_int, 'Ptprc', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Ccr2', reduction = 'umap.integrated')


#T cells
FeaturePlot(ParseSeuratObj_int, 'Cd3g', reduction = 'umap.integrated')

#Pericytes
#Along with PangloaDB, this paper has pericyte markers: 
#https://www.sciencedirect.com/science/article/pii/S1537189124001605
FeaturePlot(ParseSeuratObj_int, 'Pdgfrb', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Acta2', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Atp13a5', reduction = 'umap.integrated')


#Dimplots for convenience
DimPlot(ParseSeuratObj_int, label = TRUE, reduction = 'umap.integrated')
DimPlot(ParseSeuratObj_int, label = TRUE, group.by = 'singleR_labels', reduction = 'umap.integrated')

