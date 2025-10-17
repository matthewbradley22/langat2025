#Packages and functions
library(gridExtra)
source('./scripts/langatFunctions.R')

#Load data
ParseSeuratObj_int <- LoadSeuratRds("./data/FilteredRpcaIntegratedDatNoDoublets.rds") 
ParseSeuratObj_int$hasVirus = ifelse(ParseSeuratObj_int$virusCountPAdj >= 10, 1, 0)

#Check data
newCols <-  c(brewer.pal(12, 'Paired'), '#99FFE6', '#8A5730', '#18662E',  '#FF8AEF','#737272')
DimPlot(ParseSeuratObj_int, label = FALSE, group.by = 'manualAnnotation', reduction = 'umap.integrated',
        cols = newCols)

#ReUMAP
wt_cerebrum <- subset(ParseSeuratObj_int, Treatment %in% c('PBS', 'rLGTV') & Organ == 'Cerebrum' & Genotype == 'WT')

wt_cerebrum <- prepSeuratObj(wt_cerebrum)
ElbowPlot(wt_cerebrum, ndims = 40)
wt_cerebrum <- prepUmapSeuratObj(wt_cerebrum, nDims = 20, reductionName = 'wt.cerebrum.umap')

DimPlot(wt_cerebrum, reduction = 'wt.cerebrum.umap', label = TRUE)
DimPlot(wt_cerebrum, reduction = 'wt.cerebrum.umap', group.by = 'manualAnnotation', cols = newCols)+
  ggtitle('WT Cerebrum UMAP')

#Barplot of proportions
barDat <- wt_cerebrum[[]] %>% dplyr::group_by(Treatment, manualAnnotation) %>% 
  dplyr::summarise(total = n()) %>% dplyr::mutate(prop = total/sum(total)) %>% 
  arrange(desc(total)) 

order <- subset(barDat, Treatment == 'PBS')
barDat$manualAnnotation = factor(barDat$manualAnnotation, levels = order$manualAnnotation)
barDat %>% ggplot(aes(x = Treatment, y = prop, fill = manualAnnotation))+
  geom_bar(stat = 'identity', position = 'stack')+
  scale_fill_manual(values = newCols)

#Immune cell subset
#Need to finish labelling last mo/microglia cluster
immune <- subset(wt_cerebrum, manualAnnotation %in% c('Microglia', 'Macrophage/Monocytes', 'B Cells',
                                                             'T cells', 'Granulocytes', 'Nk cells'))

#ReUMAP
immune <- prepSeuratObj(immune)
ElbowPlot(immune, ndims = 40)
immune <- prepUmapSeuratObj(immune, nDims = 25, reductionName = 'immune.umap')

#Plug granulocytes into allan brain atlas and check
DimPlot(immune, reduction = 'immune.umap', label = TRUE)
DimPlot(immune, reduction = 'immune.umap', group.by = 'manualAnnotation')

#Dotplot
DotPlot(immune, features = c('Lgals3', 'Adgre1', 'Ptprc', 'Cd68', 'Ccr1', 'Ccr2', 'Ccr3', 'Ccr5', 'Tmem119'),
        group.by = 'manualAnnotation')+
  scale_color_gradient2(low = 'blue', mid = 'white', high = 'red')+
  scale_size(range = c(2, 10))+
   theme(axis.text.x = element_text(angle = 45, vjust = 0.8))

#Featureplots
geneList = c('Lgals3', 'Adgre1', 'Ptprc', 'Ccr1',
             'Ccr2', 'Ccr3', 'Ccr5', 'Cd68',
             'Cd86', 'Tmem119', 'Tspo', 'Csf1r')
FeaturePlot(immune, geneList, reduction = 'immune.umap')

featurePlotLight <- function(gene){
  FeaturePlot(immune, gene, reduction = 'immune.umap') +  
    theme(line = element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_blank())
}

plotList <- list(featurePlotLight('Lgals3'),
                 featurePlotLight('Adgre1'),
                 featurePlotLight('Ptprc'),
                 featurePlotLight('Ccr1'),
                 featurePlotLight('Ccr2'),
                 featurePlotLight('Ccr3'),
                 featurePlotLight('Ccr5'),
                 featurePlotLight('Cd68'),
                 featurePlotLight('Cd86'),
                 featurePlotLight('Tmem119'),
                 featurePlotLight('Tspo'),
                 featurePlotLight('Csf1r'))
do.call(grid.arrange, plotList)


#Look at mock and infected separately
wt_cerebrum_mock <- subset(wt_cerebrum, Treatment == ('PBS'))
wt_cerebrum_infected <- subset(wt_cerebrum, Treatment == ('rLGTV'))

wt_cerebrum_mock <- prepSeuratObj(wt_cerebrum_mock)
ElbowPlot(wt_cerebrum_mock, ndims = 40)
wt_cerebrum_mock <- prepUmapSeuratObj(wt_cerebrum_mock, nDims = 20, reductionName = 'wt.cerebrum.mock.umap')

DimPlot(wt_cerebrum_mock, reduction = 'wt.cerebrum.mock.umap', label = TRUE)
DimPlot(wt_cerebrum_mock, reduction = 'wt.cerebrum.mock.umap', group.by = 'manualAnnotation',
        cols = newCols)

immune_wt_mock <- subset(wt_cerebrum_mock, manualAnnotation %in% c('Microglia', 'Macrophage/Monocytes', 'B Cells',
                                                      'T cells', 'Granulocytes', 'Nk cells'))

immune_wt_mock <- prepSeuratObj(immune_wt_mock)
ElbowPlot(immune_wt_mock, ndims = 40)
immune_wt_mock <- prepUmapSeuratObj(immune_wt_mock, nDims = 20, reductionName = 'wt.immune.mock.umap')

DimPlot(immune_wt_mock, reduction = 'wt.immune.mock.umap', label = TRUE)
DimPlot(immune_wt_mock, reduction = 'wt.immune.mock.umap', group.by = 'manualAnnotation',
        cols = newCols)

#Infected now
wt_cerebrum_infected <- prepSeuratObj(wt_cerebrum_infected)
ElbowPlot(wt_cerebrum_infected, ndims = 40)
wt_cerebrum_infected <- prepUmapSeuratObj(wt_cerebrum_infected, nDims = 20, reductionName = 'wt.cerebrum.infected.umap')

DimPlot(wt_cerebrum_infected, reduction = 'wt.cerebrum.infected.umap', label = TRUE)
DimPlot(wt_cerebrum_infected, reduction = 'wt.cerebrum.infected.umap', group.by = 'manualAnnotation',
        cols = newCols)

immune_wt_infected <- subset(wt_cerebrum_infected, manualAnnotation %in% c('Microglia', 'Macrophage/Monocytes', 'B Cells',
                                                                   'T cells', 'Granulocytes', 'Nk cells'))

immune_wt_infected <- prepSeuratObj(immune_wt_infected)
ElbowPlot(immune_wt_infected, ndims = 40)
immune_wt_infected <- prepUmapSeuratObj(immune_wt_infected, nDims = 20, reductionName = 'wt.immune.infected.umap')

DimPlot(immune_wt_infected, reduction = 'wt.immune.infected.umap', label = TRUE)
DimPlot(immune_wt_infected, reduction = 'wt.immune.infected.umap', group.by = 'manualAnnotation',
        cols = newCols)

#Now look at macrophage subsets
#plug into allen brain + look at markers + literature scan
macrophages_wt_mock <- subset(wt_cerebrum_mock, manualAnnotation %in% c('Macrophage/Monocytes'))

macrophages_wt_mock <- prepSeuratObj(macrophages_wt_mock)
ElbowPlot(macrophages_wt_mock, ndims = 40)
macrophages_wt_mock <- prepUmapSeuratObj(macrophages_wt_mock, nDims = 15, reductionName = 'wt.immune.mac.umap')

DimPlot(macrophages_wt_mock, reduction = 'wt.immune.mac.umap', label = TRUE, group.by = 'seurat_clusters')
macMarkers <- FindAllMarkers(macrophages_wt_mock, only.pos = TRUE, assay = 'RNA',
                             test.use = 'MAST')
View(macMarkers %>% dplyr::filter(p_val_adj < 0.01))

#Type 2 macrophage group on left?
#Dab2 https://pubmed.ncbi.nlm.nih.gov/26927671/. Mrc1 known marker
FeaturePlot(macrophages_wt_mock, 'Mrc1', reduction = 'wt.immune.mac.umap')
FeaturePlot(macrophages_wt_mock, 'Csf1r', reduction = 'wt.immune.mac.umap')
FeaturePlot(macrophages_wt_mock, 'Dab2', reduction = 'wt.immune.mac.umap')

#APOE inflammation control marker https://pmc.ncbi.nlm.nih.gov/articles/PMC10436255/
FeaturePlot(macrophages_wt_mock, 'Apoe', reduction = 'wt.immune.mac.umap')

#Type 1 markers?
#Cd74 https://www.nature.com/articles/s41598-024-58899-7
FeaturePlot(macrophages_wt_mock, 'Cd74', reduction = 'wt.immune.mac.umap')
FeaturePlot(macrophages_wt_mock, 'Ccr2', reduction = 'wt.immune.mac.umap')

#infected macs
macrophages_wt_infected <- subset(wt_cerebrum_infected, manualAnnotation %in% c('Macrophage/Monocytes'))

macrophages_wt_infected <- prepSeuratObj(macrophages_wt_infected)
ElbowPlot(macrophages_wt_infected, ndims = 40)
macrophages_wt_infected <- prepUmapSeuratObj(macrophages_wt_infected, nDims = 20, reductionName = 'wt.infected.mac.umap')

DimPlot(macrophages_wt_infected, reduction = 'wt.infected.mac.umap', label = TRUE, group.by = 'seurat_clusters')
macMarkers <- FindAllMarkers(macrophages_wt_infected, only.pos = TRUE, assay = 'RNA',
                             test.use = 'MAST')
View(macMarkers %>% dplyr::filter(p_val_adj < 0.01))

#interesting genes
FeaturePlot(macrophages_wt_infected, reduction = 'wt.infected.mac.umap', features = 'Ccr7')
FeaturePlot(macrophages_wt_infected, reduction = 'wt.infected.mac.umap', features = 'Cxcl9') #m1 marker
FeaturePlot(macrophages_wt_infected, reduction = 'wt.infected.mac.umap', features = 'Adgre1')
FeaturePlot(macrophages_wt_infected, reduction = 'wt.infected.mac.umap', features = 'Lgals3')

#Some type 1 markers per https://www.nature.com/articles/s41598-020-73624-w
#https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0119751
#Poppovich paper seems like good source https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0145342
FeaturePlot(macrophages_wt_infected, reduction = 'wt.infected.mac.umap', features = 'Fcgr2b')
FeaturePlot(macrophages_wt_infected, reduction = 'wt.infected.mac.umap', features = 'Fcgr1')
FeaturePlot(macrophages_wt_infected, reduction = 'wt.infected.mac.umap', features = 'Cd80')
FeaturePlot(macrophages_wt_infected, reduction = 'wt.infected.mac.umap', features = 'Il6')
FeaturePlot(macrophages_wt_infected, reduction = 'wt.infected.mac.umap', features = 'Il12a')
FeaturePlot(macrophages_wt_infected, reduction = 'wt.infected.mac.umap', features = 'Il12b')
FeaturePlot(macrophages_wt_infected, reduction = 'wt.infected.mac.umap', features = 'Tnf')
FeaturePlot(macrophages_wt_infected, reduction = 'wt.infected.mac.umap', features = 'Slc7a2')
FeaturePlot(macrophages_wt_infected, reduction = 'wt.infected.mac.umap', features = 'Cxcl9')
FeaturePlot(macrophages_wt_infected, reduction = 'wt.infected.mac.umap', features = 'Il1b')
FeaturePlot(macrophages_wt_infected, reduction = 'wt.infected.mac.umap', features = 'Il1a')
FeaturePlot(macrophages_wt_infected, reduction = 'wt.infected.mac.umap', features = 'Ptgs2')

#Type 2
FeaturePlot(macrophages_wt_infected, reduction = 'wt.infected.mac.umap', features = 'Cd163')
FeaturePlot(macrophages_wt_infected, reduction = 'wt.infected.mac.umap', features = 'Msr1')
FeaturePlot(macrophages_wt_infected, reduction = 'wt.infected.mac.umap', features = 'Mrc1')

#M2c macs
FeaturePlot(macrophages_wt_infected, reduction = 'wt.infected.mac.umap', features = 'Il10')
FeaturePlot(macrophages_wt_infected, reduction = 'wt.infected.mac.umap', features = 'Tgfb1')

#Cluster 15 seems interesting, no idea which type
#Look at how many macros coexpress f480 and lgals3 

