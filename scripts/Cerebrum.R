#Script analyzing cerebellum

#Load packages
library(Seurat)

#Load in data
ParseSeuratObj_int <- LoadSeuratRds("./data/FilteredRpcaIntegratedDat.rds")
ParseSeuratObj_int <- subset(ParseSeuratObj_int, scDblFinderLabel == 'singlet')
#Good to load in manual annotations from LangatCellAnnotations.R as well

#Subset to cerebellum data
cerebrumObj <- subset(ParseSeuratObj_int, Organ == 'Cerebrum')
cerebrumObj$virusPresence = ifelse(cerebrumObj$virusCountPAdj>0, 'yes', 'no')

#Rerun through data processing and visualization
cerebrumObj <- NormalizeData(cerebrumObj)
cerebrumObj <- FindVariableFeatures(cerebrumObj)
cerebrumObj <- ScaleData(cerebrumObj)
cerebrumObj <- RunPCA(cerebrumObj)
ElbowPlot(cerebrumObj, ndims = 40)
cerebrumObj <- FindNeighbors(cerebrumObj, dims = 1:30, reduction = "pca")
cerebrumObj <- FindClusters(cerebrumObj, resolution = 2, cluster.name = "cerebellum_clusters")  
cerebrumObj <- RunUMAP(cerebrumObj, dims = 1:30, reduction = "pca", reduction.name = "umap")

#Look at data
DimPlot(cerebrumObj, reduction = 'umap', group.by = 'manualAnnotation', label = TRUE) + NoLegend()
DimPlot(cerebrumObj, reduction = 'umap', group.by = 'virusPresence', label = TRUE)
DimPlot(cerebrumObj, reduction = 'umap', group.by = 'Timepoint')
DimPlot(cerebrumObj, reduction = 'umap', group.by = 'Genotype')
DimPlot(cerebrumObj, reduction = 'umap', group.by = 'cerebellum_clusters', label = TRUE)

#Timepoint differences
cerebrumObj[[]] <- cerebrumObj[[]] %>%  mutate(virusPresence = ifelse(virusCountPAdj > 2, 'yes', 'no'))

cerebrumObj[[]]  %>% filter(Treatment == 'rChLGTV') %>% group_by(Timepoint, Genotype) %>% 
  dplyr::summarise(virusPresenceProp = mean(virusPresence == 'yes')) %>% 
  ggplot(aes(x = Genotype, y = virusPresenceProp, fill = Timepoint))+
  geom_bar(stat = 'identity', position = 'dodge')+
  ggtitle('rChLGTV')


#Plot proportion infected per cell type
cerebrumObj[[]] %>% group_by(manualAnnotation) %>% dplyr::count(virusPresence) %>% 
  ggplot(aes(x = manualAnnotation, y = n, fill = virusPresence))+
  geom_bar(position = 'fill', stat = 'identity')+
  theme(axis.text.x = element_text(angle = 90))



