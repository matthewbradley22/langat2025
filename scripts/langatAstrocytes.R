#Explore astrocytes in data

#Load packages
library(Seurat)
library(ggplot2)
library(dplyr)

#Source functions
source('./scripts/langatFunctions.R')

#Load in data
ParseSeuratObj_int <- LoadSeuratRds("./data/FilteredRpcaIntegratedDat.rds")

#Not filtering to singlets yet, first look at data and make sure it makes sense (ie doublets cluster together)

#Plot astrocytes
DimPlot(ParseSeuratObj_int, reduction = 'umap.integrated', label = TRUE) +
  NoLegend()

DimPlot(ParseSeuratObj_int, reduction = 'umap.integrated', group.by = 'singleR_labels',
        cols = c('gray', 'red', rep('gray', 16))) 

#Astrocytes appear to make up clusters 4, 5, 12, 15, 18, 21, 35, 40
#Look more closely at 31, 32, 34, 36 - 39
#33 and 14 are labelled epithelial cells by singleR, but don't express the canonical markers, need to check
#if they are astrocytes
table(subset(ParseSeuratObj_int, singleR_labels == 'Astrocytes')$seurat_clusters)

#Conservative approach to labelling astrocytes
astrocytes <- subset(ParseSeuratObj_int, singleR_labels == 'Astrocytes' & 
                       seurat_clusters %in% c('4', '5', '12', '15', '18', '21', '35', '40'))

#Initial look at variable distributions
table(astrocytes$Genotype)
table(astrocytes$Treatment) %>% barplot(main = 'Astrocytes across treatments') 
table(astrocytes$Timepoint)%>% barplot(main = 'Astrocytes across times') 

#Many more astrocytes in cerebrum
table(astrocytes$Organ)%>% barplot(main = 'Astrocytes across organs') 

#Check virus presence. Make this slightly stricter later on (need > 1 virus) for closer analysis
astrocytes$virusPresent <- ifelse(astrocytes$virusCountPAdj > 0, 1, 0)
table(astrocytes$virusPresent) %>% barplot(main = 'Astrocytes with viral reads') 
table(astrocytes$virusPresent, astrocytes$Organ)

#Rerun through data processing and visualization
astrocytes <- NormalizeData(astrocytes)
astrocytes <- FindVariableFeatures(astrocytes)
astrocytes <- ScaleData(astrocytes)
astrocytes <- RunPCA(astrocytes)
ElbowPlot(astrocytes, ndims = 40)
astrocytes <- FindNeighbors(astrocytes, dims = 1:30, reduction = "pca")
astrocytes <- FindClusters(astrocytes, resolution = 2, cluster.name = "neuron_clusters")  
astrocytes <- RunUMAP(astrocytes, dims = 1:30, reduction = "pca", reduction.name = "umap")

DimPlot(astrocytes, reduction = 'umap', label = TRUE)

#fetal/ proliferative markers from https://pmc.ncbi.nlm.nih.gov/articles/PMC5890820/
FeaturePlot(astrocytes, features = 'Top2a', reduction = 'umap')
FeaturePlot(astrocytes, features = 'Mki67', reduction = 'umap')




