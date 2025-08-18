#Script analyzing cerebellum

#Load packages
library(Seurat)
library(RColorBrewer)

#Load in data
ParseSeuratObj_int <- LoadSeuratRds("./data/seuratSingletsAnnotated.rds")
DimPlot(ParseSeuratObj_int, group.by = 'scDblFinderLabel', reduction = "umap.integrated")
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

#Viral infection over time
cerebrumObj[[]]  %>% filter(Treatment == 'rChLGTV') %>% group_by(Timepoint) %>% 
  dplyr::summarise(virusPresenceProp = mean(virusPresence == 'yes')) %>% 
  ggplot(aes(x = Timepoint, y = virusPresenceProp, fill = Timepoint))+
  geom_bar(stat = 'identity', position = 'dodge')+
  ggtitle('rChLGTV over time - cerebrum')+
  ylab('Proportion infected cells')

cerebrumObj[[]]  %>% filter(Treatment == 'rLGTV') %>% group_by(Timepoint) %>% 
  dplyr::summarise(virusPresenceProp = mean(virusPresence == 'yes')) %>% 
  ggplot(aes(x = Timepoint, y = virusPresenceProp, fill = Timepoint))+
  geom_bar(stat = 'identity', position = 'dodge')+
  ggtitle('rLGTV over time - cerebrum')+
  ylab('Proportion infected cells')

#Cell type proportions by timepoint
newCols <-  c(brewer.pal(12, 'Paired'), '#99FFE6', '#CE99FF')
cerebrumObj[[]]  %>% filter(Treatment == 'rChLGTV') %>% group_by(Timepoint, Genotype) %>% 
  dplyr::count(manualAnnotation) %>%  ggplot(aes(x = Timepoint, y = n, fill = manualAnnotation))+
  geom_bar(stat = 'identity', position = 'fill')+
  scale_fill_manual(values=newCols)+
  ggtitle('Cell type proportions in cerebrum')

#Line chart of celltype change over time
cerebrumObj[[]]  %>% filter(Treatment == 'rChLGTV') %>% group_by(Timepoint, Genotype) %>% 
  dplyr::count(manualAnnotation) %>% mutate(day = as.numeric(substr(Timepoint, 5, 5))) %>% 
  ggplot(aes(x = day, y = n, color = manualAnnotation))+
  geom_smooth(se = FALSE)+
  scale_color_manual(values=newCols)+
  scale_x_continuous(breaks=c(3, 4, 5))

#Plot proportion infected per cell type

cerebrumObj[[]] %>% filter(Treatment == 'rChLGTV') %>%  group_by(manualAnnotation) %>% dplyr::count(virusPresence) %>% 
  ggplot(aes(x = manualAnnotation, y = n, fill = virusPresence))+
  geom_bar(position = 'fill', stat = 'identity')+
  theme(axis.text.x = element_text(angle = 90))

cerebrumObj[[]] %>% filter(Treatment == 'rChLGTV') %>%  group_by(manualAnnotation, Timepoint) %>% dplyr::count(virusPresence) %>% 
  mutate(ratio = n / sum(n)) %>% filter(virusPresence == 'yes') %>% 
  mutate(day = as.numeric(substr(Timepoint, 5, 5))) %>% 
  ggplot(aes(x = day, y = ratio, color = manualAnnotation))+
  geom_line(size = 1)+
  scale_color_manual(values=newCols)

#NA group looks high, but this could make sense if doublets (higher chance of having virus present
#in one of the two cells?)
cerebrumObj[[]] %>% filter(Treatment == 'rLGTV') %>%  group_by(manualAnnotation) %>% 
  dplyr::count(virusPresence) %>% 
  ggplot(aes(x = manualAnnotation, y = n, fill = virusPresence))+
  geom_bar(position = 'fill', stat = 'identity')+
  theme(axis.text.x = element_text(angle = 90))

cerebrumObj[[]] %>% filter(Treatment == 'rLGTV') %>%  group_by(manualAnnotation, Timepoint) %>% 
  dplyr::count(virusPresence) %>% 
  mutate(ratio = n / sum(n)) %>% filter(virusPresence == 'yes') %>% 
  mutate(day = as.numeric(substr(Timepoint, 5, 5))) %>% 
  ggplot(aes(x = day, y = ratio, color = manualAnnotation))+
  geom_line(se = FALSE)+
  scale_color_manual(values=newCols)+
  ylab('Proportion of infected cells')


