#Packages and functions
library(gridExtra)
library(ggpubr)
library(Seurat)
source('./scripts/langatFunctions.R')

#Load data
ParseSeuratObj_int <- LoadSeuratRds("./data/FilteredRpcaIntegratedDatNoDoublets.rds") 
ParseSeuratObj_int$hasVirus = ifelse(ParseSeuratObj_int$virusCountPAdj >= 10, 1, 0)

#Check data
newCols <-  c(brewer.pal(12, 'Paired'), '#99FFE6', '#CE99FF', '#18662E','#737272',  '#FF8AEF')
newCols[11] =  '#FF8AEF'
DimPlot(ParseSeuratObj_int, label = FALSE, group.by = 'manualAnnotation', reduction = 'umap.integrated',
        cols = newCols)

false_macrophages_toremove <- colnames(subset(immune_wt_infected, seurat_clusters == 11))
ParseSeuratObj_int <- subset(ParseSeuratObj_int, cells = false_macrophages_toremove, invert = TRUE)
#ReUMAP
wt_cerebrum <- subset(ParseSeuratObj_int, Treatment %in% c('PBS', 'rLGTV') & Organ == 'Cerebrum' & Genotype == 'WT')

wt_cerebrum <- prepSeuratObj(wt_cerebrum)
ElbowPlot(wt_cerebrum, ndims = 40)
wt_cerebrum <- prepUmapSeuratObj(wt_cerebrum, nDims = 20, reductionName = 'wt.cerebrum.umap')

DimPlot(wt_cerebrum, reduction = 'wt.cerebrum.umap', label = TRUE)

pdf(file = '~/Documents/ÖverbyLab/scPlots/galectin3_proj/wt_cerebrum_heatmap.pdf',
    width = 7, height = 5)
DimPlot(wt_cerebrum, reduction = 'wt.cerebrum.umap', group.by = 'manualAnnotation', cols = newCols)+
  ggtitle('WT Cerebrum UMAP')+
  xlab('Umap 1')+
  ylab('Umap 2')+  
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank())
dev.off()

#Barplot of proportions
barDat <- wt_cerebrum[[]] %>% dplyr::group_by(Treatment, manualAnnotation) %>% 
  dplyr::summarise(total = n()) %>% dplyr::mutate(prop = total/sum(total)) %>% 
  arrange(desc(total)) 

barDat %>% dplyr::filter(manualAnnotation %in% c('Microglia', 'Macrophage/Monocytes'))

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

pdf(file = '~/Documents/ÖverbyLab/scPlots/galectin3_proj/immune_cell_umap.pdf',
    width = 7, height = 5)
DimPlot(immune, reduction = 'immune.umap', group.by = 'manualAnnotation')+
  ggtitle('Immune cells')+
  xlab('Umap 1')+
  ylab('Umap 2')+  
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank())
dev.off()

#Dotplot
DotPlot(immune, features = c('Lgals3', 'Adgre1', 'Ptprc', 'Cd68', 'Cd86', 'Ccr1', 'Ccr2', 
                             'Ccr3', 'Ccr5', 'Tmem119', 'Tspo', 'Csf1r'),
        group.by = 'manualAnnotation')+
  scale_color_gradient2(low = 'blue', mid = 'white', high = 'red')+
  scale_size(range = c(2, 10))+
   theme(axis.text.x = element_text(angle = 45, vjust = 0.7))

#Featureplots
geneList = c('Lgals3', 'Adgre1', 'Ptprc', 'Ccr1',
             'Ccr2', 'Ccr3', 'Ccr5', 'Cd68',
             'Cd86', 'Tmem119', 'Tspo', 'Csf1r')
FeaturePlot(immune, geneList, reduction = 'immune.umap')

featurePlotLight <- function(gene, data, reduction_choice, scale = FALSE){
  dat = FeaturePlot(data, gene, reduction = reduction_choice)$data
  colnames(dat) = c('umap1', 'umap2', 'ident', 'expression')
  ggplot(dat, aes(x = umap1, y = umap2, color = expression))+
    geom_point(size = 0.1)+  
    theme(line = element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_blank(),
          panel.background = element_rect(fill = '#F2F2F2', color = '#F2F2F2'))+
    scale_color_gradient(low = 'lightgrey', high = 'blue', limits = c(0,6))+
    ggtitle(gene)
    # annotate(x = min(dat$umap1) - 1, xend = min(dat$umap1) - 1, 
    #          y = min(dat$umap2) - 1, yend = min(dat$umap2) + 1, geom = 'segment')+
    # annotate(x = min(dat$umap1)-1, xend = min(dat$umap1) + 1, 
    #          y = min(dat$umap2) - 1, yend = min(dat$umap2) - 1, geom = 'segment')
}



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

pdf(file = '~/Documents/ÖverbyLab/scPlots/galectin3_proj/wt_immune_mock_cerebrum.pdf',
    width = 7, height = 5)
DimPlot(immune_wt_mock, reduction = 'wt.immune.mock.umap', group.by = 'manualAnnotation',
        cols = newCols)
dev.off()

plotList <- list(featurePlotLight('Lgals3', data = immune_wt_mock, reduction_choice = 'wt.immune.mock.umap'),
                 featurePlotLight('Adgre1', data = immune_wt_mock, reduction_choice = 'wt.immune.mock.umap'),
                 featurePlotLight('Ptprc', data = immune_wt_mock, reduction_choice = 'wt.immune.mock.umap'),
                 featurePlotLight('Ccr1', data = immune_wt_mock, reduction_choice = 'wt.immune.mock.umap'),
                 featurePlotLight('Ccr2', data = immune_wt_mock, reduction_choice = 'wt.immune.mock.umap'),
                 featurePlotLight('Ccr3', data = immune_wt_mock, reduction_choice = 'wt.immune.mock.umap'),
                 featurePlotLight('Ccr5', data = immune_wt_mock, reduction_choice = 'wt.immune.mock.umap'),
                 featurePlotLight('Cd68', data = immune_wt_mock, reduction_choice = 'wt.immune.mock.umap'),
                 featurePlotLight('Cd86', data = immune_wt_mock, reduction_choice = 'wt.immune.mock.umap'),
                 featurePlotLight('Tmem119', data = immune_wt_mock, reduction_choice = 'wt.immune.mock.umap'),
                 featurePlotLight('Tspo', data = immune_wt_mock, reduction_choice = 'wt.immune.mock.umap'),
                 featurePlotLight('Csf1r', data = immune_wt_mock, reduction_choice = 'wt.immune.mock.umap'))

pdf(file = '~/Documents/ÖverbyLab/scPlots/galectin3_proj/wt_immune_mock_features.pdf',
    width = 9, height = 6)
do.call(ggarrange, c(plotList, common.legend = TRUE, legend = 'right'))
dev.off()

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

#Microglia
FeaturePlot(immune_wt_infected, 'Ctss', reduction = 'wt.immune.infected.umap')
FeaturePlot(immune_wt_infected, 'Csf1r', reduction = 'wt.immune.infected.umap')
FeaturePlot(immune_wt_infected, 'Tmem119', reduction = 'wt.immune.infected.umap')
FeaturePlot(immune_wt_infected, 'P2ry12', reduction = 'wt.immune.infected.umap')
FeaturePlot(immune_wt_infected, 'Cx3cr1', reduction = 'wt.immune.infected.umap')
FeaturePlot(immune_wt_infected, 'Itgam', reduction = 'wt.immune.infected.umap')

#https://actaneurocomms.biomedcentral.com/articles/10.1186/s40478-019-0665-y has some markers
FeaturePlot(immune_wt_infected, 'Gda', reduction = 'wt.immune.infected.umap')
FeaturePlot(immune_wt_infected, 'Sell', reduction = 'wt.immune.infected.umap')

#Micro/macrophage markers from Allen atlas https://knowledge.brain-map.org/celltypes/CCN202002013
FeaturePlot(immune_wt_infected, 'Hexb', reduction = 'wt.immune.infected.umap')
FeaturePlot(immune_wt_infected, 'Inpp5d', reduction = 'wt.immune.infected.umap')
FeaturePlot(immune_wt_infected, 'Ms4a4a', reduction = 'wt.immune.infected.umap')
FeaturePlot(immune_wt_infected, 'Cd74', reduction = 'wt.immune.infected.umap')
FeaturePlot(immune_wt_infected, 'Klra2', reduction = 'wt.immune.infected.umap')
FeaturePlot(immune_wt_infected, 'Cd209a', reduction = 'wt.immune.infected.umap')

#Macrophage/monocyte markers
FeaturePlot(immune_wt_infected, 'Ptprc', reduction = 'wt.immune.infected.umap')
FeaturePlot(immune_wt_infected, 'Ccr2', reduction = 'wt.immune.infected.umap')
FeaturePlot(immune_wt_infected, 'Lyz2', reduction = 'wt.immune.infected.umap')

FeaturePlot(immune_wt_infected, 'Adgre1', reduction = 'wt.immune.infected.umap')

#Possible bam cluster? genes from https://journals.lww.com/nrronline/fulltext/2026/01000/changes_in_border_associated_macrophages_after.38.aspx
#and https://link.springer.com/article/10.1186/s12974-024-03059-x
FeaturePlot(immune_wt_infected, features =  'Mrc1', reduction = 'wt.immune.infected.umap')
FeaturePlot(immune_wt_infected, features =  'Pf4', reduction = 'wt.immune.infected.umap')
FeaturePlot(immune_wt_infected, features =  'Cd163', reduction = 'wt.immune.infected.umap')
FeaturePlot(immune_wt_infected, features =  'Apoe', reduction = 'wt.immune.infected.umap')
FeaturePlot(immune_wt_infected, features =  'Ms4a7', reduction = 'wt.immune.infected.umap')

#Look for clust 11 markers against microglia and macrophages
clust11_markers_vsMicro <- FindMarkers(immune_wt_infected, group.by = 'seurat_clusters', ident.1 = 11, 
                               ident.2 = 5,
            test.use = 'MAST') 
clust11_markers_vsMac<- FindMarkers(immune_wt_infected, group.by = 'seurat_clusters', ident.1 = 11, 
                                       ident.2 = 1,
                                       test.use = 'MAST') 
clust11_markers_vsMicro
topGenesVsMicro <- head(rownames(clust11_markers_vsMicro), n = 20)
FeaturePlot(immune_wt_infected, topGenesVsMicro, reduction = 'wt.immune.infected.umap')

clust11_markers_vsMac
topGenesVsMac <- head(rownames(clust11_markers_vsMac), n = 20)
FeaturePlot(immune_wt_infected, features =  'Nav3', reduction = 'wt.immune.infected.umap')
clust11_markers_vsMac['Tmem119',]
#Plot important markers
plotList_infected <- list(featurePlotLight('Lgals3', data = immune_wt_infected, reduction_choice = 'wt.immune.infected.umap'),
                 featurePlotLight('Adgre1', data = immune_wt_infected, reduction_choice = 'wt.immune.infected.umap'),
                 featurePlotLight('Ptprc', data = immune_wt_infected, reduction_choice = 'wt.immune.infected.umap'),
                 featurePlotLight('Ccr1', data = immune_wt_infected, reduction_choice = 'wt.immune.infected.umap'),
                 featurePlotLight('Ccr2', data = immune_wt_infected, reduction_choice = 'wt.immune.infected.umap'),
                 featurePlotLight('Ccr3', data = immune_wt_infected, reduction_choice = 'wt.immune.infected.umap'),
                 featurePlotLight('Ccr5', data = immune_wt_infected, reduction_choice = 'wt.immune.infected.umap'),
                 featurePlotLight('Cd68', data = immune_wt_infected, reduction_choice = 'wt.immune.infected.umap'),
                 featurePlotLight('Cd86', data = immune_wt_infected, reduction_choice = 'wt.immune.infected.umap'),
                 featurePlotLight('Tmem119', data = immune_wt_infected, reduction_choice = 'wt.immune.infected.umap'),
                 featurePlotLight('Tspo', data = immune_wt_infected, reduction_choice = 'wt.immune.infected.umap'),
                 featurePlotLight('Csf1r', data = immune_wt_infected, reduction_choice = 'wt.immune.infected.umap'))

lapply(plotList_infected, FUN = function(x){
  #Make sure 6 is high enough scale for all plots
  dat = x$data
  print(max(dat[4]))
})
pdf(file = '~/Documents/ÖverbyLab/scPlots/galectin3_proj/wt_immune_infected_features.pdf',
    width = 9, height = 6)
do.call(ggarrange, c(plotList_infected, common.legend = TRUE, legend = 'right'))
dev.off()

#Nonclustering mac group seems to fit in inflammatory monocyte signature from here https://www.nature.com/articles/s41467-021-21407-w
#Ly6c2high, Ccr2high, and Tgfbilow but also probably in some microglia groups they mention, check further

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

DimPlot(macrophages_wt_infected, reduction = 'wt.infected.mac.umap', label = TRUE, group.by = 'seurat_clusters')+
  ggtitle('WT Infected Macrophages')

macMarkers <- FindAllMarkers(macrophages_wt_infected, only.pos = TRUE, assay = 'RNA',
                             test.use = 'MAST')

macMarkers$pct_dif = macMarkers$pct.1 - macMarkers$pct.2
macMarkers %>% arrange(desc(pct_dif)) %>% 
  group_by(cluster) %>% dplyr::slice_head(n = 5) %>% 
  arrange(desc(pct_dif)) %>% View()

macMarkers %>% dplyr::filter(p_val_adj < 0.01 & cluster == 14) %>% head(n = 150) %>%  dplyr::select(gene) %>% 
  remove_rownames() %>% write.csv(file = "~/Documents/ÖverbyLab/macrophage_markers_clust9.csv", quote = FALSE, row.names = FALSE)

#Cluster 4 - lots about response to virus
FeaturePlot(macrophages_wt_infected, 'Nos2', reduction = 'wt.infected.mac.umap', slot = 'data') #proinflammatory
FeaturePlot(macrophages_wt_infected, 'Bnip3', reduction = 'wt.infected.mac.umap', slot = 'data') #Apoptosis
FeaturePlot(macrophages_wt_infected, 'Ccr1', reduction = 'wt.infected.mac.umap', slot = 'data') #recruitment
FeaturePlot(macrophages_wt_infected, 'Cxcl2', reduction = 'wt.infected.mac.umap', slot = 'data') #recruitment
FeaturePlot(macrophages_wt_infected, 'Ccl5', reduction = 'wt.infected.mac.umap', slot = 'data')  #recruitment

#Cluster 9 markers
FeaturePlot(macrophages_wt_infected, 'Gpnmb', reduction = 'wt.infected.mac.umap', slot = 'data') #Often anti inflammatory https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2021.674739/full
FeaturePlot(macrophages_wt_infected, 'Spp1', reduction = 'wt.infected.mac.umap', slot = 'data') #

#Clust 11
FeaturePlot(macrophages_wt_infected, 'Cd74', reduction = 'wt.infected.mac.umap', slot = 'data') #histocompatibility chaperone
FeaturePlot(macrophages_wt_infected, 'H2-Ab1', reduction = 'wt.infected.mac.umap', slot = 'data') 
FeaturePlot(macrophages_wt_infected, 'H2-Eb1', reduction = 'wt.infected.mac.umap', slot = 'data') 
FeaturePlot(macrophages_wt_infected, 'Ciita', reduction = 'wt.infected.mac.umap', slot = 'data') 

#Clust 14
FeaturePlot(macrophages_wt_infected, 'Vcan', reduction = 'wt.infected.mac.umap', slot = 'data') 

#Write out data for allen mapmycells
macrophages_wt_infected[['RNA']]$counts %>% t() %>% write.csv(file = './data/wt_infected_macrophages.csv', row.names = TRUE)

#I think this doesn't work well becaues mapMy has almost no monocytes to map to
wt_infected_macro_mapMY <- read_csv("/Users/matthewbradley/Documents/ÖverbyLab/data/wt_infected_macrophagescsv_mapMy/wt_infected_macrophagescsv_10xWholeMouseBrain(CCN20230722)_HierarchicalMapping_UTC_1761034702698.csv", 
                           skip = 4)

#interesting genes
FeaturePlot(macrophages_wt_infected, reduction = 'wt.infected.mac.umap', features = 'Cxcl9') #m1 marker

pdf(file = '~/Documents/ÖverbyLab/scPlots/galectin3_proj/infected_macs_Adgre1.pdf',
    width = 7, height = 5)
FeaturePlot(macrophages_wt_infected, reduction = 'wt.infected.mac.umap', features = 'rna_Adgre1')+
  theme(line = element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        panel.background = element_rect(fill = '#F2F2F2', color = '#F2F2F2'))+
  ggtitle('Adgre1')
dev.off()

pdf(file = '~/Documents/ÖverbyLab/scPlots/galectin3_proj/infected_macs_Lgals3.pdf',
    width = 7, height = 5)
FeaturePlot(macrophages_wt_infected, reduction = 'wt.infected.mac.umap', features = 'rna_Lgals3')+
  theme(line = element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        panel.background = element_rect(fill = '#F2F2F2', color = '#F2F2F2'))+
  ggtitle('Lgals3')
dev.off()

FeaturePlot(macrophages_wt_infected, reduction = 'wt.infected.mac.umap', features = 'Slc7a2')

#Some type 1 markers per https://www.nature.com/articles/s41598-020-73624-w
#https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0119751
#Poppovich paper seems like good source https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0145342
#Canonical m1 markers
m1Canonical <- list(featurePlotLight('Tnf', data = macrophages_wt_infected, reduction_choice = 'wt.infected.mac.umap'),
                 featurePlotLight('Nos2', data = macrophages_wt_infected, reduction_choice = 'wt.infected.mac.umap'),
                 featurePlotLight('Il1b', data = macrophages_wt_infected, reduction_choice = 'wt.infected.mac.umap'),
                 featurePlotLight('Il6', data = macrophages_wt_infected, reduction_choice = 'wt.infected.mac.umap'),
                 featurePlotLight('Il12b', data = macrophages_wt_infected, reduction_choice = 'wt.infected.mac.umap'),
                 featurePlotLight('Ccr7', data = macrophages_wt_infected, reduction_choice = 'wt.infected.mac.umap'),
                 featurePlotLight('Inhba', data = macrophages_wt_infected, reduction_choice = 'wt.infected.mac.umap'))
do.call(grid.arrange, m1Canonical)

FeaturePlot(macrophages_wt_infected, reduction = 'wt.infected.mac.umap', features = 'Tnf')
FeaturePlot(macrophages_wt_infected, reduction = 'wt.infected.mac.umap', features = 'Nos2')
FeaturePlot(macrophages_wt_infected, reduction = 'wt.infected.mac.umap', features = 'Il1b')
FeaturePlot(macrophages_wt_infected, reduction = 'wt.infected.mac.umap', features = 'Il6')
FeaturePlot(macrophages_wt_infected, reduction = 'wt.infected.mac.umap', features = 'Il12b')
FeaturePlot(macrophages_wt_infected, reduction = 'wt.infected.mac.umap', features = 'Ccr7')

#Other m1
FeaturePlot(macrophages_wt_infected, reduction = 'wt.infected.mac.umap', features = 'Fcgr2b')
FeaturePlot(macrophages_wt_infected, reduction = 'wt.infected.mac.umap', features = 'Fcgr1')
FeaturePlot(macrophages_wt_infected, reduction = 'wt.infected.mac.umap', features = 'Cd80')
FeaturePlot(macrophages_wt_infected, reduction = 'wt.infected.mac.umap', features = 'Il12a')
FeaturePlot(macrophages_wt_infected, reduction = 'wt.infected.mac.umap', features = 'Slc7a2')
FeaturePlot(macrophages_wt_infected, reduction = 'wt.infected.mac.umap', features = 'Cxcl9')
FeaturePlot(macrophages_wt_infected, reduction = 'wt.infected.mac.umap', features = 'Il1a')
FeaturePlot(macrophages_wt_infected, reduction = 'wt.infected.mac.umap', features = 'Ptgs2')

#Type 2
FeaturePlot(macrophages_wt_infected, reduction = 'wt.infected.mac.umap', features = 'Cd163')
FeaturePlot(macrophages_wt_infected, reduction = 'wt.infected.mac.umap', features = 'Msr1')
FeaturePlot(macrophages_wt_infected, reduction = 'wt.infected.mac.umap', features = 'Mrc1')

#M2 canonical
m2Canonical <- list(featurePlotLight('Arg1', data = macrophages_wt_infected, reduction_choice = 'wt.infected.mac.umap'),
                    featurePlotLight('Chil3', data = macrophages_wt_infected, reduction_choice = 'wt.infected.mac.umap'),
                    featurePlotLight('Egr2', data = macrophages_wt_infected, reduction_choice = 'wt.infected.mac.umap'),
                    featurePlotLight('Fn1', data = macrophages_wt_infected, reduction_choice = 'wt.infected.mac.umap'),
                    featurePlotLight('Mrc1', data = macrophages_wt_infected, reduction_choice = 'wt.infected.mac.umap'))
do.call(grid.arrange, m2Canonical)
#Other m2 markers
FeaturePlot(macrophages_wt_infected, reduction = 'wt.infected.mac.umap', features = 'Il10')
FeaturePlot(macrophages_wt_infected, reduction = 'wt.infected.mac.umap', features = 'Tgfb1')

#Plot cluster markers
clusterMarkerList = list()
for(i in 1:length(unique(macMarkers$cluster))){
  clusterMarkers <- head(macMarkers[macMarkers$cluster == i-1,], n = 10)
  clusterMarkers <- clusterMarkers$gene
  clusterPlots <- lapply(clusterMarkers, FUN = featurePlotLight, 
                          data = macrophages_wt_infected, reduction_choice = 'wt.infected.mac.umap')
  clusterMarkerList[[i]] = clusterPlots
}

#Number from list indicates cluster +1, so plotting [[4]] is cluster 3
do.call(grid.arrange, clusterMarkerList[[16]])

#Cluster 15 seems interesting, no idea which type
#Look at how many macros coexpress f480 and lgals3 
macrophages_wt_infected$Lgals3 = FetchData(macrophages_wt_infected, vars = 'Lgals3', layer = 'counts')
macrophages_wt_infected$Adgre1 = FetchData(macrophages_wt_infected, vars = 'Adgre1', layer = 'counts')
table(macrophages_wt_infected$Lgals3>0, macrophages_wt_infected$Adgre1>0 )
macrophages_wt_infected[[]] <- macrophages_wt_infected[[]] %>% mutate(Lgals_Adgre_both = 
                                         case_when(Lgals3 > 0 & Adgre1> 0 ~ 'Both',
                                                   Lgals3 > 0 & Adgre1 == 0 ~ 'Lgals3',
                                                   Lgals3 == 0 & Adgre1 > 0 ~ 'Adgre1',
                                                   Lgals3 == 0 & Adgre1 == 0 ~ 'Neither'))

pdf(file = '~/Documents/ÖverbyLab/scPlots/galectin3_proj/infected_macs_Lgals3_Adgre1_overlay.pdf',
    width = 7, height = 5)
DimPlot(macrophages_wt_infected, reduction = 'wt.infected.mac.umap',
        group.by = 'Lgals_Adgre_both')+
  ggtitle('Lgals3 Adgre1 Expression WT Infected Macrophages')+
  theme(line = element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        panel.background = element_rect(fill = '#F2F2F2', color = '#F2F2F2'))+
  scale_color_manual(values = c('#0FA7FF', '#FC5656', '#65BD40', 'gray'))
dev.off()

#Look at unknown cells in micro/macro clusters
unknown <- subset(wt_cerebrum, manualAnnotation == 'unknown' & seurat_clusters %in% c(0,1,2,15,31))


unknown <- prepSeuratObj(unknown)
ElbowPlot(unknown, ndims = 40)
unknown <- prepUmapSeuratObj(unknown, nDims = 20, reductionName = 'wt.cerebrum.unknown.umap')

DimPlot(unknown, reduction = 'wt.cerebrum.unknown.umap', label = TRUE)

plotList_unknown <- list(featurePlotLight('Lgals3', data = unknown, reduction_choice = 'wt.cerebrum.unknown.umap'),
                 featurePlotLight('Adgre1', data = unknown, reduction_choice = 'wt.cerebrum.unknown.umap'),
                 featurePlotLight('Ptprc', data = unknown, reduction_choice = 'wt.cerebrum.unknown.umap'),
                 featurePlotLight('Ccr1', data = unknown, reduction_choice = 'wt.cerebrum.unknown.umap'),
                 featurePlotLight('Ccr2', data = unknown, reduction_choice = 'wt.cerebrum.unknown.umap'),
                 featurePlotLight('Ccr3', data = unknown, reduction_choice = 'wt.cerebrum.unknown.umap'),
                 featurePlotLight('Ccr5', data = unknown, reduction_choice = 'wt.cerebrum.unknown.umap'),
                 featurePlotLight('Cd68', data = unknown, reduction_choice = 'wt.cerebrum.unknown.umap'),
                 featurePlotLight('Cd86', data = unknown, reduction_choice = 'wt.cerebrum.unknown.umap'),
                 featurePlotLight('Tmem119', data = unknown, reduction_choice = 'wt.cerebrum.unknown.umap'),
                 featurePlotLight('Tspo', data = unknown, reduction_choice = 'wt.cerebrum.unknown.umap'),
                 featurePlotLight('Csf1r', data = unknown, reduction_choice = 'wt.cerebrum.unknown.umap'),
                 featurePlotLight('Gda', data = unknown, reduction_choice = 'wt.cerebrum.unknown.umap'),
                 featurePlotLight('Sell', data = unknown, reduction_choice = 'wt.cerebrum.unknown.umap'))


do.call(ggarrange, c(plotList_unknown, common.legend = TRUE, legend = 'right'))




