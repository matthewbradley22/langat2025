library(Seurat)

source('./scripts/langatFunctions.R')
#Load in singlet data
ParseSeuratObj_int <- LoadSeuratRds("./data/FilteredRpcaIntegratedDatNoDoublets.rds") 
ParseSeuratObj_int$hasVirus = ifelse(ParseSeuratObj_int$virusCountPAdj >= 10, 1, 0)

#Check data
newCols <-  c(brewer.pal(12, 'Paired'), '#99FFE6', '#CE99FF', '#18662E',  '#FF8AEF','#737272')
DimPlot(ParseSeuratObj_int, label = FALSE, group.by = 'manualAnnotation', reduction = 'umap.integrated',
        cols = newCols)

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

#Plot infection over time
ependymal[[]] %>% dplyr::filter(Treatment == 'rLGTV' & Genotype == 'IPS1') %>% 
  dplyr::group_by(Timepoint, hasVirus) %>% 
  dplyr::summarise(count = n()) %>% 
  dplyr::mutate(freq = count/sum(count)) %>%  
  ggplot(aes(x = Timepoint, y = count, fill = factor(hasVirus)))+
  geom_bar(position="dodge", stat="identity")+
  ggtitle('IPS1 Ependymal infection LGTV') +
  guides(fill=guide_legend(title="Has virus")) +
  ylab('Cell count')+
  geom_text(aes(label=round(freq, digits = 2)), vjust=0,
            position = position_dodge(width = .9))

ependymal[[]] %>% dplyr::filter(Treatment == 'rChLGTV'  & Genotype == 'WT') %>% 
  dplyr::group_by(Timepoint, hasVirus) %>% 
  dplyr::summarise(count = n()) %>% 
  dplyr::mutate(freq = count/sum(count)) %>%  
  ggplot(aes(x = Timepoint, y = count, fill = factor(hasVirus)))+
  geom_bar(position="dodge", stat="identity")+
  ggtitle('WT Ependymal infection ChLGTV')+
  guides(fill=guide_legend(title="Has virus"))+
  geom_text(aes(label=round(freq, digits = 2)), vjust=0,
            position = position_dodge(width = .9))

#How many infected cells?
table(ependymal$hasVirus)

ependymal$virusPresence <- ifelse(ependymal$virusCountPAdj >= 10, 1, 0)

ependymal <- prepSeuratObj(ependymal)
ElbowPlot(ependymal, ndims = 40)

ependymal <- prepUmapSeuratObj(ependymal, nDims = 25, reductionName = 'ependymalUMAP')
DimPlot(ependymal, reduction = 'ependymalUMAP', group.by = 'Treatment')
DimPlot(ependymal, reduction = 'ependymalUMAP', group.by = 'hasVirus')
DimPlot(ependymal, reduction = 'ependymalUMAP', group.by = 'Genotype')
DimPlot(ependymal, reduction = 'ependymalUMAP', label = TRUE) + theme(legend.position = 'None')

#What is cluster 22?
#Top marker is a neuropeptide receptor Npsr1? Neurons? Cdh8 also seems to be present in neurons mostly
#MYH2 mostly in muscle fibers... Doublets?
markers22 <- FindMarkers(ependymal, group.by = 'seurat_clusters', ident.1 = 22,
                                only.pos = TRUE)

#Doesn't seem to highly express canonical neuron markers
markers22['Snap25',]
head(markers22, n = 20)

#Cerebrum and cerebellum glance
ependymal[[]] %>% dplyr::filter(Treatment == 'rLGTV'  & Organ == 'Cerebrum') %>% 
  dplyr::group_by(Organ, hasVirus) %>% 
  dplyr::summarise(count = n()) %>% 
  dplyr::mutate(freq = count/sum(count)) %>%  
  ggplot(aes(x = Organ, y = count, fill = factor(hasVirus)))+
  geom_bar(position="dodge", stat="identity")+
  ggtitle('Cerebellum Ependymal infection ChLGTV')+
  guides(fill=guide_legend(title="Has virus"))+
  geom_text(aes(label=round(freq, digits = 2)), vjust=0,
            position = position_dodge(width = .9))


#Check if any timepoints have different distributions across organs, not concerningly so
table(ependymal$Organ, ependymal$Timepoint)
table(ependymal$Treatment, ependymal$Timepoint)

#DEGs in infected cells
#Hyou1 upregulated in many diseases https://www.tandfonline.com/doi/full/10.2147/OTT.S297332
#Several top genes upregulated in tumors, involved in general transcription etc...
markersInfected <- FindMarkers(ependymal, group.by = 'hasVirus', ident.1 = '1',
                                only.pos = TRUE)
head(markersInfected, n = 20)

write.csv(rownames(markersInfected)[1:50], 
          '~/Documents/ÖverbyLab/data/geneOntologyData/infectedEpendymalMarkers.csv',
          quote = FALSE,
          row.names = FALSE)

#What is separating clusters
DimPlot(ependymal, reduction = 'ependymalUMAP', label = TRUE)

ependymal[[]] <-  ependymal[[]] %>% mutate(umapGroup = case_when(seurat_clusters %in% c(14, 8, 13, 27, 21, 23, 25, 26,
                                                                          5, 19) ~ 'top',
                                               seurat_clusters == 22 ~ 'outlier',
                                               .default = 'bottom'))
DimPlot(ependymal, reduction = 'ependymalUMAP', group.by = 'umapGroup')
bottomMarkers <- FindMarkers(ependymal, group.by = 'umapGroup', ident.1 = 'bottom',
                               only.pos = TRUE)
bottomMarkers$gene = rownames(bottomMarkers)
bottomMarkers$pct_dif <- bottomMarkers$pct.1 - bottomMarkers$pct.2
head(bottomMarkers, n = 20)
bottomMarkers %>% arrange(desc(pct_dif)) %>% head(n = 20)
FeaturePlot(ependymal, reduction = 'ependymalUMAP', features = 'Dnah3')

