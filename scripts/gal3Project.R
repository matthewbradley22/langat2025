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
wt_cerebrum <- subset(ParseSeuratObj_int, Treatment %in% c('PBS', 'rLGTV') & Organ == 'Cerebrum')

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





