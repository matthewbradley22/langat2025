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

#Astrocytes appear to make up clusters 4, 5, 12, 15, 21, 35
#Look more closely at 31, 32, 34, 36 - 40
#33 and 14 are labelled epithelial cells by singleR, but don't express the canonical markers, need to check
#if they are astrocytes
table(subset(ParseSeuratObj_int, singleR_labels == 'Astrocytes')$seurat_clusters)





