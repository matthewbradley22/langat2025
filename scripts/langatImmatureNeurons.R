library(Seurat)

source('./scripts/langatFunctions.R')
#Load in singlet data
ParseSeuratObj_int <- LoadSeuratRds("./data/FilteredRpcaIntegratedDatNoDoublets.rds") 
ParseSeuratObj_int$hasVirus = ifelse(ParseSeuratObj_int$virusCountPAdj >= 10, 1, 0)

#Check data
newCols <-  c(brewer.pal(12, 'Paired'), '#99FFE6', '#CE99FF', '#18662E',  '#FF8AEF','#737272')
DimPlot(ParseSeuratObj_int, label = FALSE, group.by = 'manualAnnotation', reduction = 'umap.integrated',
        cols = newCols)

#Confirm main immatureNeu markers expressed compared to rest of data
markersImmatureNeurons <- FindMarkers(ParseSeuratObj_int, group.by = 'manualAnnotation', 
                                      ident.1 = 'Immature Neurons',
                                only.pos = TRUE)

#Check known markers
markersImmatureNeurons['Dcx',]
markersImmatureNeurons['Sbk1',]

#Subset to immature neurons 
immatureNeu <- subset(ParseSeuratObj_int, manualAnnotation == 'Immature Neurons')
hist(log(immatureNeu$virusCountPAdj))

#Plot infection over time
immatureNeu[[]] %>% dplyr::filter(Treatment == 'rLGTV' & Genotype == 'WT') %>% 
  dplyr::group_by(Timepoint, hasVirus) %>% 
  dplyr::summarise(count = n()) %>% 
  dplyr::mutate(freq = count/sum(count)) %>%  
  ggplot(aes(x = Timepoint, y = count, fill = factor(hasVirus)))+
  geom_bar(position="dodge", stat="identity")+
  ggtitle('WT Immature neu infection LGTV') +
  guides(fill=guide_legend(title="Has virus")) +
  ylab('Cell count')+
  geom_text(aes(label=round(freq, digits = 2)), vjust=0,
            position = position_dodge(width = .9))

immatureNeu[[]] %>% dplyr::filter(Treatment == 'rChLGTV'  & Genotype == 'IPS1') %>% 
  dplyr::group_by(Timepoint, hasVirus) %>% 
  dplyr::summarise(count = n()) %>% 
  dplyr::mutate(freq = count/sum(count)) %>%  
  ggplot(aes(x = Timepoint, y = count, fill = factor(hasVirus)))+
  geom_bar(position="dodge", stat="identity")+
  ggtitle('IPS1 Immature neu infection ChLGTV')+
  guides(fill=guide_legend(title="Has virus"))+
  geom_text(aes(label=round(freq, digits = 2)), vjust=0,
            position = position_dodge(width = .9))

#Very few day 4 captured
table(immatureNeu$Timepoint)

immatureNeu <- prepSeuratObj(immatureNeu)
ElbowPlot(immatureNeu, ndims = 40)

immatureNeu <- prepUmapSeuratObj(immatureNeu, nDims = 20, reductionName = 'immatureNeuUMAP')
DimPlot(immatureNeu, reduction = 'immatureNeuUMAP', group.by = 'Treatment')
DimPlot(immatureNeu, reduction = 'immatureNeuUMAP', group.by = 'hasVirus')
DimPlot(immatureNeu, reduction = 'immatureNeuUMAP', group.by = 'Timepoint')
DimPlot(immatureNeu, reduction = 'immatureNeuUMAP', label = TRUE) + theme(legend.position = 'None')


#DEGs in infected cells
#Hyou1 upregulated in many diseases https://www.tandfonline.com/doi/full/10.2147/OTT.S297332
#Several top genes upregulated in tumors, involved in general transcription etc...
markersInfected <- FindMarkers(immatureNeu, group.by = 'hasVirus', ident.1 = '1',
                               only.pos = TRUE)
head(markersInfected, n = 20)

write.csv(rownames(markersInfected)[1:50], 
          '~/Documents/ÖverbyLab/data/geneOntologyData/infectedimmatureNeuMarkers.csv',
          quote = FALSE,
          row.names = FALSE)

#What is separating clusters
DimPlot(immatureNeu, reduction = 'immatureNeuUMAP', label = TRUE)

immatureNeu[[]] <-  immatureNeu[[]] %>% mutate(umapGroup = case_when(seurat_clusters %in% c(14, 8, 13, 27, 21, 23, 25, 26,
                                                                                        5, 19) ~ 'top',
                                                                 seurat_clusters == 22 ~ 'outlier',
                                                                 .default = 'bottom'))
DimPlot(immatureNeu, reduction = 'immatureNeuUMAP', group.by = 'umapGroup')
bottomMarkers <- FindMarkers(immatureNeu, group.by = 'umapGroup', ident.1 = 'bottom',
                             only.pos = TRUE)
bottomMarkers$gene = rownames(bottomMarkers)
bottomMarkers$pct_dif <- bottomMarkers$pct.1 - bottomMarkers$pct.2
head(bottomMarkers, n = 20)
bottomMarkers %>% arrange(desc(pct_dif)) %>% head(n = 20)
FeaturePlot(immatureNeu, reduction = 'immatureNeuUMAP', features = 'Dnah3')
