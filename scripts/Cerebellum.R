#Script analyzing cerebellum

#Load packages
library(Seurat)
library(RColorBrewer)
library(tidyr)

#Load in data
ParseSeuratObj_int <- LoadSeuratRds("./data/seuratSingletsAnnotated.rds")

#Good to load in manual annotations from LangatCellAnnotations.R as well

#Subset to cerebellum data
cerebellumObj <- subset(ParseSeuratObj_int, Organ == 'Cerebellum' & scDblFinderLabel == 'singlet')
cerebellumObj$virusPresence = ifelse(cerebellumObj$virusCountPAdj>0, 'yes', 'no')

#Rerun through data processing and visualization
cerebellumObj <- NormalizeData(cerebellumObj)
cerebellumObj <- FindVariableFeatures(cerebellumObj)
cerebellumObj <- ScaleData(cerebellumObj)
cerebellumObj <- RunPCA(cerebellumObj)
ElbowPlot(cerebellumObj, ndims = 40)
cerebellumObj <- FindNeighbors(cerebellumObj, dims = 1:30, reduction = "pca")
cerebellumObj <- FindClusters(cerebellumObj, resolution = 2, cluster.name = "cerebellum_clusters")  
cerebellumObj <- RunUMAP(cerebellumObj, dims = 1:30, reduction = "pca", reduction.name = "umap")

#Look at data
DimPlot(cerebellumObj, reduction = 'umap', group.by = 'manualAnnotation', label = TRUE) + NoLegend()
DimPlot(cerebellumObj, reduction = 'umap', group.by = 'Timepoint')
DimPlot(cerebellumObj, reduction = 'umap', group.by = 'Genotype')
DimPlot(cerebellumObj, reduction = 'umap', group.by = 'cerebellum_clusters', label = TRUE)


#Timepoint differences
cerebellumObj[[]] <- cerebellumObj[[]] %>%  mutate(virusPresence = ifelse(virusCountPAdj > 2, 'yes', 'no'))
cerebellumObj[[]]  %>% group_by(Timepoint, Genotype) %>% 
  dplyr::summarise(virusPresenceProp = mean(virusPresence == 'yes')) %>% 
  ggplot(aes(x = Genotype, y = virusPresenceProp, fill = Timepoint))+
  geom_bar(stat = 'identity', position = 'dodge')

#Plot proportion infected per cell type
cerebellumObj[[]] %>% group_by(manualAnnotation) %>% dplyr::count(virusPresence) %>% 
 ggplot(aes(x = manualAnnotation, y = n, fill = virusPresence))+
  geom_bar(position = 'fill', stat = 'identity')+
  theme(axis.text.x = element_text(angle = 90))


#Cell type proportions by timepoint
newCols <-  c(brewer.pal(12, 'Paired'), '#99FFE6', '#CE99FF', '#18662E', '#737272')
newCols[2] = "#CAB2D6"
newCols[9] = "#1F78B4"

cerebellumObj[[]]  %>% filter(Treatment == 'rChLGTV') %>% group_by(Timepoint, Genotype) %>% 
  dplyr::count(manualAnnotation) %>%  ggplot(aes(x = Timepoint, y = n, fill = manualAnnotation))+
  geom_bar(stat = 'identity', position = 'fill')+
  scale_fill_manual(values=newCols)+
  ggtitle('Cell type proportions in cerebellum')
