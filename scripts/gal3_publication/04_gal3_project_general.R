library(gridExtra)
library(ggpubr)
library(Seurat)
library(gprofiler2)
source('./langatFunctions.R')

#Load data
ParseSeuratObj_int <- LoadSeuratRds("~/Documents/ÖverbyLab/data/FilteredRpcaIntegratedDatNoDoublets.rds") 

# - - - - - - - - - - - - - - - - 
#### Further doublet removal ####
# - - - - - - - - - - - - - - - - 

#Have spotted one more group of doublets to remove for future analyses
immune_wt_infected_all <- subset(ParseSeuratObj_int, Organ == 'Cerebrum' & 
                                   Genotype == 'WT' &
                                   manualAnnotation %in% c('Microglia', 'Macrophage/Monocytes', 
                                                                      'T cells', 'Granulocytes', 'Nk cells') &
                                   Treatment == ('rLGTV')) 

#Function from langatFunctions.R
immune_wt_infected_all <- prepSeuratObj(immune_wt_infected_all) 
ElbowPlot(immune_wt_infected_all, ndims = 40)

#Function from langatFunctions.R
immune_wt_infected_all <- prepUmapSeuratObj(immune_wt_infected_all, nDims = 20, 
                                            reductionName = 'wt.immune.infected.umap', resolution_value = 0.8) 

#Large cluster of labelled doublets that have markers of both macrophages and microglia
DimPlot(immune_wt_infected_all, reduction = 'wt.immune.infected.umap', label = TRUE)
DimPlot(immune_wt_infected_all, reduction = 'wt.immune.infected.umap', label = TRUE, group.by = 'scDblFinderLabel')

false_macs_to_remove <- colnames(subset(immune_wt_infected_all, seurat_clusters == 8))

ParseSeuratObj_int <- subset(ParseSeuratObj_int, cells = false_macs_to_remove, invert = TRUE)

#SaveSeuratRds(ParseSeuratObj_int, "./data/seuratSingletsAnnotated.rds")

# - - - - - - - - - - - - - - - - - - 
#### WT Cerebrum day 5 analysis ####
# - - - - - - - - - - - - - - - - - - 

#Keep day 5 infected cells, and all mock cells
wt_cerebrum_day5 <-  subset(ParseSeuratObj_int, Treatment %in% c('PBS', 'rLGTV') & Organ == 'Cerebrum' & 
                              Genotype == 'WT' & (Timepoint == 'Day 5' | Treatment == 'PBS'))

#Prepare for umap of subset cells (functions from langatFunctions.R)
wt_cerebrum_day5 <- prepSeuratObj(wt_cerebrum_day5)
ElbowPlot(wt_cerebrum_day5, ndims = 40)
wt_cerebrum_day5 <- prepUmapSeuratObj(wt_cerebrum_day5, nDims = 20, reductionName = 'wt.cerebrum.umap',  num_neighbors = 30L)

#Prepare colors and order for UMAP
umap_color_list <- c(   "#7047A1", "#6D92F8", "#FF96A2","#F76363", "#D6644B",
                        "#8a0000", "#074F00", "#208d1f","#7bcd79", "#fdc087","#F08C3A", "#B370AE","#6DC3F8", "#166DF0", "#292270",
                        "gray")

wt_cerebrum_day5$manualAnnotation <- factor(wt_cerebrum_day5$manualAnnotation, 
                                            levels = c('Astrocytes', 'Choroid Plexus', 'Endothelial',
                                                       'Ependymal', 'Immature Neurons', 'Microglia', 'Muscle cells',
                                                       'Neurons', 'Oligodendrocytes', 'Pericytes', 'B Cells',
                                                       'Granulocytes', 'Macrophage/Monocytes', 'Nk cells', 'T cells',
                                                       'unknown'))

#UMAP
DimPlot(wt_cerebrum_day5, reduction = 'wt.cerebrum.umap', group.by = 'manualAnnotation', cols = umap_color_list)+
  ggtitle('WT Cerebrum Day 5 + PBS')+
  xlab('Umap 1')+
  ylab('Umap 2')+  
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank())

#Proportion of cells barplot
barDat <- wt_cerebrum_day5[[]] %>% dplyr::group_by(Treatment, manualAnnotation) %>% 
  dplyr::summarise(total = n()) %>% dplyr::mutate(prop = total/sum(total)) %>% 
  arrange(desc(total)) 

barDat %>% ggplot(aes(x = Treatment, y = prop, fill = manualAnnotation))+
  geom_bar(stat = 'identity', position = 'stack')+
  theme(text = element_text(size = 23))+
  scale_fill_manual(values = umap_color_list)+
  guides(fill=guide_legend(title="Cell Type"))

#Porportions of WT cerebrum cells at each timepoint
wt_cerebrum <-  subset(ParseSeuratObj_int, Treatment %in% c('PBS', 'rLGTV') & Organ == 'Cerebrum' & 
                              Genotype == 'WT')

barDat_time <- wt_cerebrum[[]] %>% dplyr::group_by(Timepoint, Treatment, manualAnnotation) %>% 
  dplyr::summarise(total = n()) %>% dplyr::mutate(prop = total/sum(total)) %>% 
  arrange(desc(total)) 

barDat_time$manualAnnotation <- factor(barDat_time$manualAnnotation, 
                                            levels = c('Astrocytes', 'Choroid Plexus', 'Endothelial',
                                                       'Ependymal', 'Immature Neurons', 'Microglia', 'Muscle cells',
                                                       'Neurons', 'Oligodendrocytes', 'Pericytes', 'B Cells',
                                                       'Granulocytes', 'Macrophage/Monocytes', 'Nk cells', 'T cells',
                                                       'unknown'))

ggplot(barDat_time, aes(x = Timepoint, y = prop, fill = manualAnnotation))+
  geom_bar(stat = 'identity', position = 'stack')+
  facet_grid(~Treatment, scales = 'free_x') +
  scale_fill_manual(values = umap_color_list)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 16))+
  xlab('')+
  ylab('')

#Dotplot of cell marker gene expression for known cell types
wt_cerebrum_day5_lgtv <- subset(wt_cerebrum_day5, manualAnnotation != 'unknown' &  Treatment == 'rLGTV')
wt_cerebrum_day5_mock <- subset(wt_cerebrum_day5, manualAnnotation != 'unknown' &  Treatment == 'PBS')

#Change order of cell types for plot
wt_cerebrum_day5_lgtv$manualAnnotation <- factor(wt_cerebrum_day5_lgtv$manualAnnotation, 
                                                  levels = c( 'unknown',  'T cells',   'Nk cells', 'Macrophage/Monocytes', 
                                                              'Granulocytes', 'B Cells',  'Pericytes', 'Oligodendrocytes','Neurons',
                                                              'Muscle cells', 'Microglia', 'Immature Neurons', 'Ependymal','Endothelial', 
                                                              'Choroid Plexus', 'Astrocytes'))
wt_cerebrum_day5_mock$manualAnnotation <- factor(wt_cerebrum_day5_mock$manualAnnotation, 
                                                 levels = c( 'unknown',  'T cells',   'Nk cells', 'Macrophage/Monocytes', 
                                                             'Granulocytes', 'B Cells',  'Pericytes', 'Oligodendrocytes','Neurons',
                                                             'Muscle cells', 'Microglia', 'Immature Neurons', 'Ependymal','Endothelial', 
                                                             'Choroid Plexus', 'Astrocytes'))

DotPlot(wt_cerebrum_day5_lgtv, features = c('Tspo', 'Lgals3', 'Adgre1', 'Cd68',  'Cd86', 'Ptprc',  'Ccr1', 'Ccr2', 
                          'Ccr3', 'Ccr5', 'Tmem119',  'Csf1r', 'Aif1'),
        group.by = 'manualAnnotation', scale = TRUE)+
  scale_color_gradient2(low = 'blue', mid = 'white', high = 'red', labels = c(-1, 0, 1,2),
                        breaks = c(-1, 0, 1,2), , limits = c(-2,2.5))+
  scale_size(range = c(2, 10), limits = c(0,100))+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.7))+
  theme(legend.position = "bottom",
        legend.justification = "center",
        legend.direction = "horizontal",
        legend.title = element_text(hjust = 0.3),
        legend.spacing.x = unit(2, "cm"))+
  guides(size = guide_legend(title.position = "top", title = 'Percent Expressed', order = 2),
         color = guide_colorbar(title.position = "top", title = 'Average Scaled Expression'))

DotPlot(wt_cerebrum_day5_mock, features = c('Tspo', 'Lgals3', 'Adgre1', 'Cd68',  'Cd86', 'Ptprc',  'Ccr1', 'Ccr2', 
                                            'Ccr3', 'Ccr5', 'Tmem119',  'Csf1r', 'Aif1'),
        group.by = 'manualAnnotation', scale = TRUE)+
  scale_color_gradient2(low = 'blue', mid = 'white', high = 'red', labels = c(-1, 0, 1,2),
                        breaks = c(-1, 0, 1,2), , limits = c(-2,2.5))+
  scale_size(range = c(2, 10), limits = c(0,100))+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.7))+
  theme(legend.position = "bottom",
        legend.justification = "center",
        legend.direction = "horizontal",
        legend.title = element_text(hjust = 0.3),
        legend.spacing.x = unit(2, "cm"))+
  guides(size = guide_legend(title.position = "top", title = 'Percent Expressed', order = 2),
         color = guide_colorbar(title.position = "top", title = 'Average Scaled Expression'))

# - - - - - - - - - - - - - - - 
#### Only immune cell types ####
# - - - - - - - - - - - - - - - 

immune <- subset(wt_cerebrum_day5, manualAnnotation %in% c('Microglia', 'Macrophage/Monocytes', 'B Cells',
                                                           'T cells', 'Granulocytes', 'Nk cells'))

#Prepare for umap of subset cells (functions from langatFunctions.R)
immune <- prepSeuratObj(immune)
ElbowPlot(immune, ndims = 40)
immune <- prepUmapSeuratObj(immune, nDims = 25, reductionName = 'immune.umap', resolution_value = 0.8)

DimPlot(immune, reduction = 'immune.umap', group.by = 'manualAnnotation', 
        cols = c("#8a0000", "#F08C3A", "#B370AE", "#6DC3F8","#166DF0", "#292270"))+
  ggtitle('Immune cells')+
  xlab('Umap 1')+
  ylab('Umap 2')+  
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank())

# - - - - - - - - - - - - - - - - - - 
#### Only infected immune cells ####
# - - - - - - - - - - - - - - - - - - 

immune_wt_infected <- subset(wt_cerebrum_day5, manualAnnotation %in% c('Microglia', 'Macrophage/Monocytes', 
                                                                                'T cells', 'Granulocytes', 'Nk cells') &
                               Treatment == 'rLGTV') 

#Prepare for umap of subset cells (functions from langatFunctions.R)
immune_wt_infected <- prepSeuratObj(immune_wt_infected)
ElbowPlot(immune_wt_infected, ndims = 40)
immune_wt_infected <- prepUmapSeuratObj(immune_wt_infected, nDims = 20, reductionName = 'wt.immune.infected.umap', 
                                        num_neighbors = 20L)

DimPlot(immune_wt_infected, reduction = 'wt.immune.infected.umap', group.by = 'manualAnnotation',
        cols = c("#B370AE", "#166DF0", "#292270","#6DC3F8", "#8a0000"))+
  ggtitle('Immune WT Infected UMAP')+
  xlab('Umap 1')+
  ylab('Umap 2')+  
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank())

# - - - - - - - - - - - - - - - - - - 
#### Only PBS mock immune cells ####
# - - - - - - - - - - - - - - - - - - 

#Don't include extreme low low cell counts here
immune_wt_mock <- subset(wt_cerebrum_day5, manualAnnotation %in% c('Microglia', 'Macrophage/Monocytes', 
                                                                   'T cells', 'B Cells') &
                           Treatment == 'PBS')

#Prepare for umap of subset cells (functions from langatFunctions.R)
immune_wt_mock <- prepSeuratObj(immune_wt_mock)
ElbowPlot(immune_wt_mock, ndims = 40)
immune_wt_mock <- prepUmapSeuratObj(immune_wt_mock, nDims = 20, reductionName = 'wt.immune.mock.umap')

DimPlot(immune_wt_mock, reduction = 'wt.immune.mock.umap', group.by = 'manualAnnotation',
        cols = c("#F08C3A", "#292270", "#6DC3F8", "#8a0000" ))+
  ggtitle('Day 5 + pbs immune mock')

#TSPO dotplot to compare cell type expression levels
wt_cerebrum_day5 <- subset(wt_cerebrum_day5, manualAnnotation != 'unknown')
create_dot_plot(dat = wt_cerebrum_day5, gene = 'Tspo', main_title = 'Tspo', x_var = 'treatment',
                combine_pbs = TRUE, flip_coords = TRUE)

#Feature plot example
FeaturePlot(immune_wt_mock, features = 'Tspo', reduction = 'wt.immune.mock.umap')
FeaturePlot(immune_wt_infected, features = 'Tspo', reduction = 'wt.immune.infected.umap')



