library(Seurat)

source('./scripts/langatFunctions.R')
#Load in singlet data
ParseSeuratObj_int <- LoadSeuratRds("./data/FilteredRpcaIntegratedDatNoDoublets.rds") 

#Confirm main ependymal markers expressed compared to rest of data
markersEpendymal <- FindMarkers(ParseSeuratObj_int, group.by = 'manualAnnotation', ident.1 = 'Ependymal',
                         only.pos = TRUE)

#looks good
markersEpendymal['Foxj1',]
markersEpendymal['Pifo',]
markersEpendymal['Dynlrb2',]

#Top 6 genes all involved in flagella/cilia formation. Weird but 
#looks like ependymal cells cilia are important https://pmc.ncbi.nlm.nih.gov/articles/PMC9315228/
#Rgs22 also ependymal cilia related https://link.springer.com/article/10.1007/s11427-024-2720-8
head(markersEpendymal, n = 40)
write.csv(rownames(markersEpendymal)[1:50], 
          '~/Documents/ÖverbyLab/data/geneOntologyData/ependymalMarkers.csv',
          quote = FALSE,
          row.names = FALSE)

#Subset to ependymal cells
ependymal <- subset(ParseSeuratObj_int, manualAnnotation == 'Ependymal')
hist(log(ependymal$virusCountPAdj))

#How many infected cells?
sum(ependymal$virusCountPAdj >= 10)
sum(ependymal$virusCountPAdj < 10)

ependymal$virusPresence <- ifelse(ependymal$virusCountPAdj >= 10, 1, 0)

ependymal <- prepSeuratObj(ependymal)
ElbowPlot(ependymal, ndims = 40)

ependymal <- prepUmapSeuratObj(ependymal, nDims = 25, reductionName = 'ependymalUMAP')
DimPlot(ependymal, reduction = 'ependymalUMAP', group.by = 'Treatment')
DimPlot(ependymal, reduction = 'ependymalUMAP', group.by = 'virusPresence')
DimPlot(ependymal, reduction = 'ependymalUMAP', label = TRUE) + theme(legend.position = 'None')

#What is cluster 22?
#Top marker is a neuropeptide receptor Npsr1? Neurons? Cdh8 also seems to be present in neurons mostly
#MYH2 mostly in muscle fibers... Doublets?
markers22 <- FindMarkers(ependymal, group.by = 'seurat_clusters', ident.1 = 22,
                                only.pos = TRUE)

#Doesn't seem to highly express canonical neuron markers
markers22['Snap25',]
head(markers22, n = 20)

#virus across variables
ependymal[[]] %>% group_by(Treatment, Genotype) %>% 
  dplyr::summarise(virusProp = mean(virusPresence), cellCount = n()) %>% 
  ggplot(aes(x = Treatment, y = virusProp, fill = Genotype))+
  geom_bar(stat = 'identity', position = 'dodge')+
  ylab('proportion infected cells')

#Not a ton of Day 4 IPS cells, still enough that I think it's meaningful
ependymal[[]] %>% filter(Treatment == 'rChLGTV') %>% group_by(Timepoint, Genotype) %>% 
  dplyr::summarise(virusProp = mean(virusPresence), cellCount = n()) %>% 
  ggplot(aes(x = Timepoint, y = virusProp, fill = Genotype))+
  geom_bar(stat = 'identity', position = 'dodge')+
  ylab('proportion infected cells')

ependymal[[]] %>% group_by(Organ, Genotype) %>% 
  dplyr::summarise(virusProp = mean(virusPresence), cellCount = n()) %>% 
  ggplot(aes(x = Organ, y = virusProp, fill = Genotype))+
  geom_bar(stat = 'identity', position = 'dodge')

#Check if any timepoints have different distributions across organs, not concerningly so
table(ependymal$Organ, ependymal$Timepoint)
table(ependymal$Treatment, ependymal$Timepoint)

#DEGs in infected cells
#Hyou1 upregulated in many diseases https://www.tandfonline.com/doi/full/10.2147/OTT.S297332
#Several top genes upregulated in tumors, involved in general transcription etc...
markersInfected <- FindMarkers(ependymal, group.by = 'virusPresence', ident.1 = '1',
                                only.pos = TRUE)
head(markersInfected, n = 20)

write.csv(rownames(markersInfected)[1:50], 
          '~/Documents/ÖverbyLab/data/geneOntologyData/infectedEpendymalMarkers.csv',
          quote = FALSE,
          row.names = FALSE)
