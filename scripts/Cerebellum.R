#Script analyzing cerebellum

#Load packages
library(Seurat)
library(RColorBrewer)
library(tidyr)

source('./scripts/langatFunctions.R')

#Load in data
ParseSeuratObj_int <- LoadSeuratRds("./data/FilteredRpcaIntegratedDat.rds")

#Load in manual annotations from the end of LangatCellAnnotations.R (will save manual labels with data once finalized)

#Subset to cerebellum data
cerebellumObj <- subset(ParseSeuratObj_int, Organ == 'Cerebellum')
cerebellumObj$virusPresence = ifelse(cerebellumObj$virusCountPAdj>4, 'yes', 'no')

#Rerun through data processing and visualization
cerebellumObj <- prepSeuratObj(cerebellumObj)
ElbowPlot(cerebellumObj, ndims = 40)
cerebellumObj <- prepUmapSeuratObj(cerebellumObj, nDims = 30, reductionName = 'umap')

#Look at data
DimPlot(cerebellumObj, reduction = 'umap', group.by = 'manualAnnotation', label = TRUE) + NoLegend()
DimPlot(cerebellumObj, reduction = 'umap', group.by = 'Timepoint')
DimPlot(cerebellumObj, reduction = 'umap', group.by = 'Genotype')
DimPlot(cerebellumObj, reduction = 'umap', label = TRUE)

#Look at opc and astrocyte markers
#Astrocytes
FeaturePlot(cerebellumObj, 'Gfap', reduction = 'umap')
FeaturePlot(cerebellumObj, 'Aqp4', reduction = 'umap')
FeaturePlot(cerebellumObj, 'Fgfr3', reduction = 'umap')
FeaturePlot(cerebellumObj, 'Aldh1l1', reduction = 'umap')
FeaturePlot(cerebellumObj, 'Slc1a3', reduction = 'umap')
FeaturePlot(cerebellumObj, 'Gli3', reduction = 'umap')
FeaturePlot(cerebellumObj, 'Slc39a12', reduction = 'umap')

#Oligo
FeaturePlot(cerebellumObj, 'Mag',reduction = 'umap')
FeaturePlot(cerebellumObj, 'Mog', reduction = 'umap')
FeaturePlot(cerebellumObj, 'Gjc3',reduction = 'umap')
FeaturePlot(cerebellumObj, 'Mbp', reduction = 'umap')

#OPC
FeaturePlot(cerebellumObj, 'Stk32a', reduction = 'umap')

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
