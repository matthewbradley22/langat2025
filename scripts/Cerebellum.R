#Script analyzing cerebellum

#Load packages
library(Seurat)

#Load in data
ParseSeuratObj_int <- LoadSeuratRds("./data/FilteredRpcaIntegratedDat.rds")

#Subset to cerebellum data
cerebellumObj <- subset(ParseSeuratObj_int, Organ == 'Cerebellum' & scDblFinderLabel == 'singlet')
cerebellumObj$isInfected = ifelse(cerebellumObj$virusCountPAdj>0, 'yes', 'no')

#Rerun through data processing and visualization
cerebellumObj <- NormalizeData(cerebellumObj)
cerebellumObj <- FindVariableFeatures(cerebellumObj)
cerebellumObj <- ScaleData(cerebellumObj)
cerebellumObj <- RunPCA(cerebellumObj)
ElbowPlot(cerebellumObj, ndims = 40)
cerebellumObj <- FindNeighbors(cerebellumObj, dims = 1:30, reduction = "pca")
cerebellumObj <- FindClusters(cerebellumObj, resolution = 2, cluster.name = "cerebellum_clusters")  
cerebellumObj <- RunUMAP(cerebellumObj, dims = 1:30, reduction = "pca", reduction.name = "umap")

#Look at doublets
DimPlot(cerebellumObj, reduction = 'umap', group.by = 'singleR_labels', label = TRUE)
DimPlot(cerebellumObj, reduction = 'umap', group.by = 'Timepoint')
DimPlot(cerebellumObj, reduction = 'umap', group.by = 'Genotype')
DimPlot(cerebellumObj, reduction = 'umap', group.by = 'cerebellum_clusters', label = TRUE)


#Timepoint differences
cerebellumObj[[]] %>%  mutate(virusPresence = ifelse(virusCountPAdj > 0, 'yes', 'no')) %>% 
  group_by(Timepoint, Genotype) %>% 
  dplyr::summarise(virusPresenceProp = mean(virusPresence == 'yes')) %>% 
  ggplot(aes(x = Genotype, y = virusPresenceProp, fill = Timepoint))+
  geom_bar(stat = 'identity', position = 'dodge')


table(cerebellumObj$Treatment)


