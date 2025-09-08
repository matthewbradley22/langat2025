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

# are labelled epithelial cells by singleR, but don't express the canonical markers, need to check
#if they are astrocytes
table(subset(ParseSeuratObj_int, singleR_labels == 'Astrocytes')$seurat_clusters)


#Rerun through data processing and visualization
astrocytes <- prepSeuratObj(astrocytes)
ElbowPlot(astrocytes, ndims = 40)
astrocytes <- FindNeighbors(astrocytes, dims = 1:30, reduction = "pca")
astrocytes <- FindClusters(astrocytes, resolution = 2, cluster.name = "astrocyte_clusters")  
astrocytes <- RunUMAP(astrocytes, dims = 1:30, reduction = "pca", reduction.name = "umap")

DimPlot(astrocytes, reduction = 'umap', label = TRUE, group.by = 'astrocyte_clusters')

#Look at doublets
DimPlot(astrocytes, reduction = 'umap', label = TRUE, group.by = 'scDblFinderLabel')

#Remove clusters comprised of mostly doublets
astrocytes_singlets <- subset(astrocytes, !(astrocyte_clusters %in% c('17', '20', '22',
                                                                      '24', '25', '8', '34',
                                                                      '18', '30', '31', '32',
                                                                      '21') ) &
                                scDblFinderLabel == 'singlet')

DimPlot(astrocytes_singlets, reduction = 'umap', label = TRUE, group.by = 'scDblFinderLabel')

astrocytes_singlets <- NormalizeData(astrocytes_singlets)
astrocytes_singlets <- FindVariableFeatures(astrocytes_singlets)
astrocytes_singlets <- ScaleData(astrocytes_singlets)
astrocytes_singlets <- RunPCA(astrocytes_singlets)
ElbowPlot(astrocytes_singlets, ndims = 40)
astrocytes_singlets <- FindNeighbors(astrocytes_singlets, dims = 1:30, reduction = "pca")
astrocytes_singlets <- FindClusters(astrocytes_singlets, resolution = 2, cluster.name = "astrocyte_clusters")  
astrocytes_singlets <- RunUMAP(astrocytes_singlets, dims = 1:30, reduction = "pca", reduction.name = "umap")


DimPlot(astrocytes_singlets, reduction = 'umap', label = TRUE, group.by = 'astrocyte_clusters')
DimPlot(astrocytes_singlets, reduction = 'umap', label = TRUE, group.by = 'astrocyte_clusters')

#### Descriptive bar plots ####
#Check virus presence. Make this slightly stricter later on (need > 1 virus) for closer analysis
astrocytes_singlets$virusPresent <- ifelse(astrocytes_singlets$virusCountPAdj > 0, 1, 0)
table(astrocytes_singlets$virusPresent) %>% barplot(main = 'Astrocytes with viral reads') 
table(astrocytes_singlets$virusPresent, astrocytes_singlets$Organ)

astrocytes_singlets[[]] %>% group_by(Timepoint) %>% summarise(virusProportion = mean(virusPresent)) %>% 
ggplot(aes(x = Timepoint, y = virusProportion, fill = Timepoint))+
  geom_bar(stat = 'identity')

sum(astrocytes_singlets$virusCountPAdj > 15)

DimPlot(astrocytes_singlets, reduction = 'umap', group.by = 'Organ')

#Look at variable distributions
table(astrocytes_singlets$Genotype)
table(astrocytes_singlets$Treatment) %>% barplot(main = 'Astrocytes across treatments') 
table(astrocytes_singlets$Timepoint)%>% barplot(main = 'Astrocytes across times') 

#Many more astrocytes in cerebrum
table(astrocytes_singlets$Organ)%>% barplot(main = 'Astrocytes across organs') 

#### Begin analysis of astrocytes ####



