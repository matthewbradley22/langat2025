#Need to organize this and the cell_death_pathways.R script
#there is not an isg regressed umap created in cell_death_pathways.R that should just be done here for organization
#Packages and functions
library(gridExtra)
library(ggpubr)
library(Seurat)
library(gprofiler2)
library(RColorBrewer)
source('~/Documents/ÖverbyLab//scripts/langatFunctions.R')

#Load data
ParseSeuratObj_int <- LoadSeuratRds("~/Documents/ÖverbyLab/data/FilteredRpcaIntegratedDatNoDoublets.rds") 
ParseSeuratObj_int$hasVirus = ifelse(ParseSeuratObj_int$virusCountPAdj >= 10, 1, 0)

#Check data
newCols <-  c(brewer.pal(12, 'Paired'), '#99FFE6', '#CE99FF', '#18662E','#737272',  '#FF8AEF')
newCols[11] =  '#FF8AEF'
DimPlot(ParseSeuratObj_int, label = FALSE, group.by = 'manualAnnotation', reduction = 'umap.integrated',
        cols = newCols)

#Can run this to remove odd macrophage group, but need to create immune_wt_infected below first
#write.csv(false_macrophages_toremove, "~/Documents/ÖverbyLab/false_macs_to_remove.csv")

false_macrophages_toremove <- colnames(subset(immune_wt_infected, seurat_clusters == 11))
#or load in from csv
false_macs_to_remove <- read.csv("~/Documents/ÖverbyLab/false_macs_to_remove.csv")[[2]]
#should be 409 cells
length(false_macs_to_remove)

ParseSeuratObj_int <- subset(ParseSeuratObj_int, cells = false_macs_to_remove, invert = TRUE)

#Check why t cells don't show up much for day 5 wt cerebrum
wt <- subset(ParseSeuratObj_int, Genotype == 'WT')
table(wt$Treatment, wt$manualAnnotation, wt$Timepoint)

#Create all subsets that will be used
#Subset to same cells as in gal3 project, but keep all timepoints for now, rather than just day 5
wt_cerebrum <-  subset(ParseSeuratObj_int, Treatment %in% c('PBS', 'rLGTV') & Organ == 'Cerebrum' & 
                              Genotype == 'WT')

wt_cerebrum <- prepSeuratObj(wt_cerebrum)
ElbowPlot(wt_cerebrum, ndims = 30)
wt_cerebrum <- prepUmapSeuratObj(wt_cerebrum, nDims = 20, reductionName = 'wt_cerebrum', resolution_value = 1)

DimPlot(wt_cerebrum, reduction = 'wt_cerebrum', group.by = 'manualAnnotation', cols = newCols)

#Polarization stimulus genes from https://www.mdpi.com/1422-0067/25/22/12078 / https://doi.org/10.3390/ijms252212078
polarization_stimulus_genes <- c('Ifng', 'Csf2', 'Il4', 'Il13', 'Il1b', 'Il10',
                                 'Tgfb1', 'Csf3', 'Pf4')

DotPlot(wt_cerebrum, features = polarization_stimulus_genes, scale = FALSE, group.by = 'manualAnnotation')+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))

#Split and look by timepoint
wt_cerebrum_day3 <- subset(wt_cerebrum, Treatment == 'rLGTV' & Timepoint == 'Day 3')
wt_cerebrum_day4 <- subset(wt_cerebrum, Treatment == 'rLGTV' & Timepoint == 'Day 4')
wt_cerebrum_day5 <- subset(wt_cerebrum, Treatment == 'rLGTV' & Timepoint == 'Day 5')

#Only keep cell types with at least 50 cells
day3_cells_with_enough <- c(table(wt_cerebrum_day3$manualAnnotation) %>% as.data.frame() %>% 
  dplyr::filter(Freq > 50) %>% dplyr::select(Var1))$Var1
wt_cerebrum_day3_dot_dat <- subset(wt_cerebrum_day3, manualAnnotation %in% day3_cells_with_enough)

DotPlot(wt_cerebrum_day3_dot_dat, features = polarization_stimulus_genes, scale = FALSE, group.by = 'manualAnnotation')+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))+
  scale_size_continuous(range = c(0.5,6), limits = c(0,100))+
  scale_color_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white"), 
                        values = c(1.0,0.7,0.4,0),
                        limits = c(0,5))+
  ggtitle('Day 3 LGTV')

day4_cells_with_enough <- c(table(wt_cerebrum_day4$manualAnnotation) %>% as.data.frame() %>% 
                              dplyr::filter(Freq > 50) %>% dplyr::select(Var1))$Var1
wt_cerebrum_day4_dot_dat <- subset(wt_cerebrum_day4, manualAnnotation %in% day4_cells_with_enough)

DotPlot(wt_cerebrum_day4_dot_dat, features = polarization_stimulus_genes, scale = FALSE, group.by = 'manualAnnotation')+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))+
  scale_size_continuous(range = c(0.5,6), limits = c(0,100))+
  scale_color_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white"), 
                        values = c(1.0,0.7,0.4,0),
                        limits = c(0,5))+
  ggtitle('Day 4')

day5_cells_with_enough <- c(table(wt_cerebrum_day5$manualAnnotation) %>% as.data.frame() %>% 
                              dplyr::filter(Freq > 50) %>% dplyr::select(Var1))$Var1
wt_cerebrum_day5_dot_dat <- subset(wt_cerebrum_day5, manualAnnotation %in% day5_cells_with_enough)

DotPlot(wt_cerebrum_day5_dot_dat, features = polarization_stimulus_genes, scale = FALSE, group.by = 'manualAnnotation')+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))+
  scale_size_continuous(range = c(0.5,6), limits = c(0,100))+
  scale_color_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white"), 
                        values = c(1.0,0.7,0.4,0),
                        limits = c(0,5))+
  ggtitle('Day 5 LGTV')

#Also look at pbs
pbs <- subset(wt_cerebrum, Treatment == 'PBS')

pbs_cells_with_enough <- c(table(pbs$manualAnnotation) %>% as.data.frame() %>% 
                              dplyr::filter(Freq > 50) %>% dplyr::select(Var1))$Var1
pbs_dot_dat <- subset(pbs, manualAnnotation %in% pbs_cells_with_enough)

DotPlot(pbs_dot_dat, features = polarization_stimulus_genes, scale = FALSE, group.by = 'manualAnnotation')+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))+
  scale_size_continuous(range = c(0.5,6), limits = c(0,100))+
  scale_color_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white"), 
                        values = c(1.0,0.7,0.4,0),
                        limits = c(0,5))+
  ggtitle('PBS')


#Which celltypes have enough cells for dotplots
table(wt_cerebrum_day3$manualAnnotation)
table(wt_cerebrum_day4$manualAnnotation)
table(wt_cerebrum_day5$manualAnnotation)

#Look at arrest/stop vs transmigration markers from https://academic.oup.com/view-large/3854487
arrest_only_markers <- c('Vcam1', 'Ccr1')
transmigration_only_markers <- c('Mcam', 'F11r', 'Jam2', 'Pecam1', 'Pvr', 'Cd99', 'Cd99l2', 'Cdh5', 'Epha1', 'Ephb1')
all_arrest_markers <- c('Itgal', 'Itgb2', 'Itgam', 'Vcam1', 'Amica1', 'Selp', 'Ccr1')

#Plot arrest and transmigration markers
DotPlot(wt_cerebrum_day3_dot_dat, features = arrest_only_markers, scale = FALSE, group.by = 'manualAnnotation')+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))+
  scale_size_continuous(range = c(0.5,6), limits = c(0,100))+
  scale_color_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white"), 
                        values = c(1.0,0.7,0.4,0),
                        limits = c(0,2.5))+
  ggtitle('Day 3 Arrest Markers')

DotPlot(wt_cerebrum_day3_dot_dat, features = transmigration_only_markers, scale = FALSE, group.by = 'manualAnnotation')+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))+
  scale_size_continuous(range = c(0.5,6), limits = c(0,100))+
  scale_color_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white"), 
                        values = c(1.0,0.7,0.4,0),
                        limits = c(0,2.5))+
  ggtitle('Day 3 Transmigration')

DotPlot(wt_cerebrum_day4_dot_dat, features = arrest_only_markers, scale = FALSE, group.by = 'manualAnnotation')+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))+
  scale_size_continuous(range = c(0.5,6), limits = c(0,100))+
  scale_color_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white"), 
                        values = c(1.0,0.7,0.4,0),
                        limits = c(0,2.5))+
  ggtitle('Day 4 Arrest Markers')

DotPlot(wt_cerebrum_day4_dot_dat, features = transmigration_only_markers, scale = FALSE, group.by = 'manualAnnotation')+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))+
  scale_size_continuous(range = c(0.5,6), limits = c(0,100))+
  scale_color_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white"), 
                        values = c(1.0,0.7,0.4,0),
                        limits = c(0,2.5))+
  ggtitle('Day 4 Transmigration')

DotPlot(wt_cerebrum_day5_dot_dat, features = arrest_only_markers, scale = FALSE, group.by = 'manualAnnotation')+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))+
  scale_size_continuous(range = c(0.5,6), limits = c(0,100))+
  scale_color_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white"), 
                        values = c(1.0,0.7,0.4,0),
                        limits = c(0,2.5))+
  ggtitle('Day 5 Arrest Markers')

DotPlot(wt_cerebrum_day5_dot_dat, features = transmigration_only_markers, scale = FALSE, group.by = 'manualAnnotation')+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))+
  scale_size_continuous(range = c(0.5,6), limits = c(0,100))+
  scale_color_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white"), 
                        values = c(1.0,0.7,0.4,0),
                        limits = c(0,2.5))+
  ggtitle('Day 5 Transmigration')



#Look at main adhesion/integrin pairs from https://link.springer.com/article/10.1186/1742-2094-9-270
#"During viral infection of the brain, we have found that recruitment of monocytes to the CNS is also VLA-4-dependent"
adhesion_pairs <- c('Itga4', 'Itgal', 'Vcam1', 'Icam1')

DotPlot(wt_cerebrum_day3_dot_dat, features = adhesion_pairs, scale = FALSE, group.by = 'manualAnnotation')+
  ggtitle('Day 3')

DotPlot(wt_cerebrum_day4_dot_dat, features = adhesion_pairs, scale = FALSE, group.by = 'manualAnnotation')+
  ggtitle('Day 4')

DotPlot(wt_cerebrum_day5_dot_dat, features = adhesion_pairs, scale = FALSE, group.by = 'manualAnnotation')+
  ggtitle('Day 5')

#Chemokine plots
ccl_chemokines <- rownames(wt_cerebrum@assays$RNA$data)[grep('Ccl', rownames(wt_cerebrum@assays$RNA$data))]
cxcl_chemokines <- rownames(wt_cerebrum@assays$RNA$data)[grep('Cxcl', rownames(wt_cerebrum@assays$RNA$data))]
wt_cerebrum$treatment_celltype <- paste(wt_cerebrum$Treatment, wt_cerebrum$manualAnnotation)

#Subset each timepoint of data to cell types with enough cells
wt_cerebrum_day3_all <- subset(wt_cerebrum, Timepoint == 'Day 3')
day3_all_cells_with_enough <- c(table(wt_cerebrum_day3_all$manualAnnotation) %>% as.data.frame() %>% 
                              dplyr::filter(Freq > 75) %>% dplyr::select(Var1))$Var1
wt_cerebrum_day3_all_dot_dat <- subset(wt_cerebrum_day3_all, manualAnnotation %in% day3_all_cells_with_enough)

wt_cerebrum_day4_all <- subset(wt_cerebrum, Timepoint == 'Day 4')
day4_all_cells_with_enough <- c(table(wt_cerebrum_day4_all$manualAnnotation) %>% as.data.frame() %>% 
                                  dplyr::filter(Freq > 75) %>% dplyr::select(Var1))$Var1
wt_cerebrum_day4_all_dot_dat <- subset(wt_cerebrum_day4_all, manualAnnotation %in% day4_all_cells_with_enough)

wt_cerebrum_day5_all <- subset(wt_cerebrum, Timepoint == 'Day 5')
day5_all_cells_with_enough <- c(table(wt_cerebrum_day5_all$manualAnnotation) %>% as.data.frame() %>% 
                                  dplyr::filter(Freq > 75) %>% dplyr::select(Var1))$Var1
wt_cerebrum_day5_all_dot_dat <- subset(wt_cerebrum_day5_all, manualAnnotation %in% day3_all_cells_with_enough)


DotPlot(wt_cerebrum_day3_all_dot_dat, features = ccl_chemokines, group.by = 'treatment_celltype')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
  ylab('')+
  xlab('')

DotPlot(wt_cerebrum_day3_all_dot_dat, features = cxcl_chemokines, group.by = 'treatment_celltype')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
  ylab('')+
  xlab('')

pdf('~/Documents/ÖverbyLab/scPlots/galectin3_proj/day3_allgal3_chemokine_dot.pdf', width = 8, height = 6)
DotPlot(wt_cerebrum_day3_all_dot_dat, features = c('Ccl1', 'Ccl2', 'Ccl3', 'Ccl4', 'Ccl5',
                                                   'Ccl6', 'Ccl7', 'Ccl8', 'Ccl9', 'Ccl10', 'Ccl11',
                                                   'Cxcl10'), group.by = 'treatment_celltype', scale = FALSE)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
  ylab('')+
  xlab('')+
  ggtitle('Day 3 Chemokines')+
  scale_size_continuous(limits = c(0,100), range = c(0.01, 6))+
  scale_color_gradient2(limits = c(0, 3.8), low = 'white', high = 'blue')
dev.off()


pdf('~/Documents/ÖverbyLab/scPlots/galectin3_proj/day5_allgal3_chemokine_dot.pdf', width = 8, height = 6)
DotPlot(wt_cerebrum_day5_all_dot_dat, features = c('Ccl1', 'Ccl2', 'Ccl3', 'Ccl4', 'Ccl5',
                                                   'Ccl6', 'Ccl7', 'Ccl8', 'Ccl9', 'Ccl10', 'Ccl11',
                                                   'Cxcl10'), group.by = 'treatment_celltype', scale = FALSE)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
  ylab('')+
  xlab('')+
  ggtitle('Day 5 Chemokines')+
  scale_size_continuous(limits = c(0,100), range = c(0.01, 5))+
  scale_color_gradient2(limits = c(0, 3.8), low = 'white', high = 'blue')
dev.off()

#Split by celltype instead
wt_endothelial <- subset(wt_cerebrum, Treatment == 'rLGTV' & manualAnnotation == 'Endothelial')
wt_mac_mono <- subset(wt_cerebrum, Treatment == 'rLGTV' & manualAnnotation == 'Macrophage/Monocytes')

DotPlot(wt_endothelial, features = c('Vcam1', 'Icam1', 'Icam2', 'Sele', 'Selp', 'Mcam', 'F11r', 'Jam2', 'Pecam1', 'Pvr', 'Cd99l2',
                                     'Cdh5'), group.by = 'Timepoint', scale = FALSE)+
  coord_flip()+
  scale_size_continuous(range = c(0.5,6), limits = c(0,100))+
  scale_color_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white"), 
                        values = c(1.0,0.7,0.4,0),
                        limits = c(0,2.5))+
  ggtitle('Endothelial trafficking genes')

DotPlot(wt_mac_mono, features = c('Pecam1', 'Pvr', 'Cd99l2', 'Epha1', 'Ephb1', 'Itga4', 'Itgb1', 'Ccr1', 'Ccr2', 'Ccr5'), group.by = 'Timepoint', scale = FALSE)+
  coord_flip()+
  scale_size_continuous(range = c(0.5,6), limits = c(0,100))+
  scale_color_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white"), 
                        values = c(1.0,0.7,0.4,0),
                        limits = c(0,2.5))+
  ggtitle('Monocyte trafficking')

#Same plots but pbs samples
wt_endothelial_pbs <- subset(wt_cerebrum, Treatment == 'PBS' & manualAnnotation == 'Endothelial')
wt_mac_mono_pbs <- subset(wt_cerebrum, Treatment == 'PBS' & manualAnnotation == 'Macrophage/Monocytes')

DotPlot(wt_endothelial_pbs, features = c('Vcam1', 'Icam1', 'Icam2', 'Sele', 'Selp', 'Mcam', 'F11r', 'Jam2', 'Pecam1', 'Pvr', 'Cd99l2',
                                     'Cdh5'), group.by = 'Timepoint', scale = FALSE)+
  coord_flip()+
  scale_size_continuous(range = c(0.5,6), limits = c(0,100))+
  scale_color_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white"), 
                        values = c(1.0,0.7,0.4,0),
                        limits = c(0,2.5))+
  ggtitle('PBS Endothelial trafficking')


DotPlot(wt_mac_mono_pbs, features = c('Pecam1', 'Pvr', 'Cd99l2', 'Epha1', 'Ephb1', 'Itga4', 'Itgb1', 'Ccr1', 'Ccr2', 'Ccr5'), group.by = 'Timepoint', scale = FALSE)+
  coord_flip()+
  scale_size_continuous(range = c(0.5,6), limits = c(0,100))+
  scale_color_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white"), 
                        values = c(1.0,0.7,0.4,0),
                        limits = c(0,2.5))+
  ggtitle('PBS Monocyte trafficking')
