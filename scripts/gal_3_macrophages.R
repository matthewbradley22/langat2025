#Need to organize this and the cell_death_pathways.R script
#there is not an isg regressed umap created in cell_death_pathways.R that should just be done here for organization
#Packages and functions
library(gridExtra)
library(ggpubr)
library(Seurat)
library(gprofiler2)
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

#Create all subsets that will be used
#Subset to same cells as in gal3 project
wt_cerebrum_day5 <-  subset(ParseSeuratObj_int, Treatment %in% c('PBS', 'rLGTV') & Organ == 'Cerebrum' & 
                              Genotype == 'WT' & (Timepoint == 'Day 5' | Treatment == 'PBS'))

#For figure after imaging figures - going to see what separates day 5 from day 3 and 4 macs 
#Earlier figures only focus on day 5 LGTV
wt_cerebrum_macrophages <-  subset(ParseSeuratObj_int, Treatment %in% c('PBS', 'rLGTV') & Organ == 'Cerebrum' & 
                                     Genotype == 'WT' & manualAnnotation == 'Macrophage/Monocytes')

#Prepare umap for infected samples
macrophages_wt_infected <- subset(wt_cerebrum_macrophages, Treatment == 'rLGTV')  
macrophages_wt_infected <- prepSeuratObj(macrophages_wt_infected)
ElbowPlot(macrophages_wt_infected, ndims = 40)
macrophages_wt_infected <- prepUmapSeuratObj(macrophages_wt_infected, nDims = 20, reductionName = 'wt.infected.mac.umap',
                                             resolution_value = 0.8)

#Prepare UMAP for macrophages
#Switch to macrophages_wt_infected to get plot in previous powerpoints/meetings
wt_cerebrum_macrophages <- prepSeuratObj(wt_cerebrum_macrophages)
ElbowPlot(wt_cerebrum_macrophages, ndims = 40)

#Higher num neighbors for fewer clusters
wt_cerebrum_macrophages <- prepUmapSeuratObj(wt_cerebrum_macrophages, nDims = 20, reductionName = 'wt.cerebrum.mac.umap',
                                             resolution_value = 0.8)

DimPlot(wt_cerebrum_macrophages, reduction = 'wt.cerebrum.mac.umap', label = TRUE, group.by = 'seurat_clusters',
        label.size = 6)+
  ggtitle('WT Infected Macrophages')

DimPlot(wt_cerebrum_macrophages, reduction = 'wt.cerebrum.mac.umap', label = FALSE, group.by = 'Timepoint',
        label.size = 6)+
  ggtitle('WT Infected Macrophages')

DimPlot(wt_cerebrum_macrophages, reduction = 'wt.cerebrum.mac.umap', label = FALSE, group.by = 'Treatment',
        label.size = 6)+
  ggtitle('WT Infected Macrophages')

#Examine what separates day 5 from other timepoints
day5_macro_markers <- FindMarkers(macrophages_wt_infected, group.by = 'Timepoint', ident.1 = 'Day 5',
                                  test.use = 'MAST')
day5_up_markers <- dplyr::filter(day5_macro_markers, (avg_log2FC) > 1 & p_val_adj < 0.01)

#Gene ontology 
day5_paths <- gprofiler2::gost(query = rownames(day5_up_markers), organism = 'mmusculus', evcodes = TRUE)
day5_paths$result[8,]

#Function for nicer featureplots
featurePlotLight <- function(gene, data, reduction_choice, scale = FALSE, minLim = 0, maxLim = 5){
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
    scale_color_gradient(low = 'lightgrey', high = 'blue', limits = c(minLim,maxLim))+
    ggtitle(gene)
}

#Monocyte/macrophage markers from nat com paper https://www.nature.com/articles/s41467-023-37698-0#MOESM5
FeaturePlot(wt_cerebrum_macrophages, features = c('Slfn4', 'Ms4a8a', 'Clec4e', 'Itga4'), reduction = 'wt.cerebrum.mac.umap')#Macrophage marker
FeaturePlot(wt_cerebrum_macrophages, features = c('Ccl3', 'Ccl4', 'Il12b', 'Tnf'), reduction = 'wt.cerebrum.mac.umap')

plotList <- lapply(c('Slfn4', 'Ms4a8a', 'Clec4e', 'Itga4'), featurePlotLight, data = wt_cerebrum_macrophages, 
                   reduction_choice = 'wt.cerebrum.mac.umap', maxLim = 6)
do.call(ggarrange, c(plotList, common.legend = TRUE, legend = 'right'))

plotList <- list(featurePlotLight('Ccl3', data = wt_cerebrum_macrophages, reduction_choice = 'wt.cerebrum.mac.umap', maxLim = 6),
                 featurePlotLight('Ccl4', data = wt_cerebrum_macrophages, reduction_choice = 'wt.cerebrum.mac.umap', maxLim = 6),
                 featurePlotLight('Il12b', data = wt_cerebrum_macrophages, reduction_choice = 'wt.cerebrum.mac.umap', maxLim = 6),
                 featurePlotLight('Tnf', data = wt_cerebrum_macrophages, reduction_choice = 'wt.cerebrum.mac.umap', maxLim = 6))

do.call(ggarrange, c(plotList, common.legend = TRUE, legend = 'right'))
#Marker paper: https://pmc.ncbi.nlm.nih.gov/articles/PMC12362712/
#Gene lists from supplementary table 1: Gene set scores
FeaturePlot(wt_cerebrum_macrophages, features = 'Nos2', reduction = 'wt.cerebrum.mac.umap')
FeaturePlot(wt_cerebrum_macrophages, features = 'Apoe', reduction = 'wt.cerebrum.mac.umap') #Phagocytosis
FeaturePlot(wt_cerebrum_macrophages, features = 'Mrc1', reduction = 'wt.cerebrum.mac.umap') #bam marker/alternate activation

FeaturePlot(wt_cerebrum_macrophages, features = 'Cxcl2', reduction = 'wt.cerebrum.mac.umap') #Oxidative stress
FeaturePlot(wt_cerebrum_macrophages, features = 'Prdx5', reduction = 'wt.cerebrum.mac.umap')
FeaturePlot(wt_cerebrum_macrophages, features = 'Txn1', reduction = 'wt.cerebrum.mac.umap')
FeaturePlot(wt_cerebrum_macrophages, features = 'Gsr', reduction = 'wt.cerebrum.mac.umap')

FeaturePlot(wt_cerebrum_macrophages, features = 'H2-Aa', reduction = 'wt.cerebrum.mac.umap') #Antigen presentation
FeaturePlot(wt_cerebrum_macrophages, features = 'H2-Eb1', reduction = 'wt.cerebrum.mac.umap')
FeaturePlot(wt_cerebrum_macrophages, features = 'H2-Ab1', reduction = 'wt.cerebrum.mac.umap')

#Complement cascade
FeaturePlot(wt_cerebrum_macrophages, features = 'C1qc', reduction = 'wt.cerebrum.mac.umap')
FeaturePlot(wt_cerebrum_macrophages, features = 'C1qb', reduction = 'wt.cerebrum.mac.umap')

#Other interesting genes
FeaturePlot(wt_cerebrum_macrophages, features = 'F10', reduction = 'wt.cerebrum.mac.umap')
FeaturePlot(wt_cerebrum_macrophages, features = 'Cd38', reduction = 'wt.cerebrum.mac.umap')
FeaturePlot(wt_cerebrum_macrophages, features = 'Ccr1', reduction = 'wt.cerebrum.mac.umap')
FeaturePlot(wt_cerebrum_macrophages, features = 'Socs3', reduction = 'wt.cerebrum.mac.umap')

#Monocyte score
monocyte_score_genes <- c('S100a4', 'Itgb7', 'Napsa', 'Cd300lg', 'Adora2b', 'Emb', 'Ly6c2', 'Ms4a4c', 'Fn1', 'Sell', 'Padi2', 'Lilra6', 
                          'Ccnb2', 'Galnt9', 'Upb1', 'Lmo1', 'F13a1', 'Ccr2', 'Gm15987', 'AI839979')
select_monocyte_genes <- c('S100a4', 'Itgb7', 'Ly6c2', 'Sell', 'Ccr2')
wt_cerebrum_macrophages <- AddModuleScore(wt_cerebrum_macrophages, features = list(monocyte_score_genes), name = 'monocyte_score')
FeaturePlot(wt_cerebrum_macrophages, features = 'monocyte_score1', reduction = 'wt.cerebrum.mac.umap')
FeaturePlot(wt_cerebrum_macrophages, features = select_monocyte_genes, reduction = 'wt.cerebrum.mac.umap')

plotList <- list(featurePlotLight('S100a4', data = wt_cerebrum_macrophages, reduction_choice = 'wt.cerebrum.mac.umap', maxLim = 5.5),
                 featurePlotLight('Itgb7', data = wt_cerebrum_macrophages, reduction_choice = 'wt.cerebrum.mac.umap', maxLim = 5.5),
                 featurePlotLight('Ly6c2', data = wt_cerebrum_macrophages, reduction_choice = 'wt.cerebrum.mac.umap', maxLim = 5.5),
                 featurePlotLight('Sell', data = wt_cerebrum_macrophages, reduction_choice = 'wt.cerebrum.mac.umap', maxLim = 5.5),
                 featurePlotLight('Ccr2', data = wt_cerebrum_macrophages, reduction_choice = 'wt.cerebrum.mac.umap', maxLim = 5.5))

do.call(ggarrange, c(plotList, common.legend = TRUE, legend = 'right'))


#Macrophage score
macrophage_score_genes <- c('Adgre1', 'Csf1r', 'H2-Ab1', 'Cd68', 'Lyz2', 'Itgam', 'Mertk')
wt_cerebrum_macrophages <- AddModuleScore(wt_cerebrum_macrophages, features = list(macrophage_score_genes), name = 'macrophage_score')
FeaturePlot(wt_cerebrum_macrophages, features = 'macrophage_score1', reduction = 'wt.cerebrum.mac.umap')
FeaturePlot(wt_cerebrum_macrophages, features = macrophage_score_genes, reduction = 'wt.cerebrum.mac.umap')

plotList <- list(featurePlotLight('Adgre1', data = wt_cerebrum_macrophages, reduction_choice = 'wt.cerebrum.mac.umap', maxLim = 5.5),
                 featurePlotLight('Csf1r', data = wt_cerebrum_macrophages, reduction_choice = 'wt.cerebrum.mac.umap', maxLim = 5.5),
                 featurePlotLight('H2-Ab1', data = wt_cerebrum_macrophages, reduction_choice = 'wt.cerebrum.mac.umap', maxLim = 5.5),
                 featurePlotLight('Cd68', data = wt_cerebrum_macrophages, reduction_choice = 'wt.cerebrum.mac.umap', maxLim = 5.5),
                 featurePlotLight('Lyz2', data = wt_cerebrum_macrophages, reduction_choice = 'wt.cerebrum.mac.umap', maxLim = 5.5),
                 featurePlotLight('Itgam', data = wt_cerebrum_macrophages, reduction_choice = 'wt.cerebrum.mac.umap', maxLim = 5.5),
                 featurePlotLight('Mertk', data = wt_cerebrum_macrophages, reduction_choice = 'wt.cerebrum.mac.umap', maxLim = 5.5))
do.call(ggarrange, c(plotList, common.legend = TRUE, legend = 'right'))

#Other macrophage scores
#Phagocytosis
phagocytosis_score <- c('Aif1', 'Cd36', 'Coro1a', 'Ccr2', 'Cnn2', 'Cyba', 'Fcer1g', 'Fcgr1', 'Fcgr3', 'Gas6', 'Irf8', 'Il2rg', 'Itgb1', 
                        'Itgb2', 'Lrp1', 'Ncf4', 'Slc11a1', 'Pros1', 'Msr1', 'Thbs1', 'Tyrobp', 'Pycard', 'Cd209b', 'Trem2', 'Gsn', 'P2ry6', 'Myo1g')
macrophages_wt_infected <- AddModuleScore(macrophages_wt_infected, features = list(phagocytosis_score), name = 'phagocytosis_score')
FeaturePlot(macrophages_wt_infected, features = 'phagocytosis_score1', reduction = 'wt.cerebrum.mac.umap')

#Inflammatory response
inflammatory_score <- c('Adam8', 'Aif1', 'Alox5', 'Alox5ap', 'Fabp4', 'Apoe', 'App', 'C1qa', 'Ciita', 'C3ar1', 'C5ar1', 'Cd14', 'Cd44', 'Cd68', 'Ccr2', 
'Csf1r', 'Ctsc', 'Cyba', 'Cybb', 'Ednrb', 'Fcer1g', 'Fcgr1', 'Fcgr3', 'Fn1', 'B4galt1', 'Grn', 'Cxcl1', 'Hp', 'Igf1', 'Il6', 'Itgam', 'Itgb2', 'Jun', 'Lpl',
'Lrp1', 'Mif', 'Nfkb1', 'Slc11a1', 'Ptger4', 'Rel', 'Ccl6', 'Cxcl2', 'Thbs1', 'Tnfrsf1b', 'Tyrobp', 
'Pik3cg', 'Cxcl13', 'Ccl24', 'Pf4', 'Pycard', 'Trem2', 'Cd163', 'Gpsm3', 'Tlr7', 'Stab1')
macrophages_wt_infected <- AddModuleScore(macrophages_wt_infected, features = list(inflammatory_score), name = 'inflammatory_score')
FeaturePlot(macrophages_wt_infected, features = 'inflammatory_score1', reduction = 'wt.cerebrum.mac.umap')

#Oxidative stress
oxidative_stress_score <- c('Prdx5', 'Txn1', 'Gsr', 'Ptgs2', 'Ccs', 'Prdx6', 'Gpx4', 'Sesn1', 'Sod3', 'Ltc4s')
macrophages_wt_infected <- AddModuleScore(macrophages_wt_infected, features = list(oxidative_stress_score), name = 'oxidative_stress_score')
FeaturePlot(macrophages_wt_infected, features = 'oxidative_stress_score1', reduction = 'wt.cerebrum.mac.umap')

#ECM Organization
ecm_score <- c('Col1a1', 'Nid1', 'Dpt', 'B4galt1', 'Lum', 'Col3a1', 'Ccdc80', 'Ramp2', 'Serpinh1', 'Ddr2')
macrophages_wt_infected <- AddModuleScore(macrophages_wt_infected, features = list(ecm_score), name = 'ecm_score')
FeaturePlot(macrophages_wt_infected, features = 'ecm_score1', reduction = 'wt.cerebrum.mac.umap') +
  ggtitle('ECM Organization score')

#Cytokine cytokine receptor interacitons
cytokine_score <- c('Ccl6', 'Ccl24', 'Cxcl1', 'Cxcl2', 'Pf4', 'Cxcl13', 'Cxcl12', 'Cxcl16', 'Il6', 'Ccr2', 'Il2rg', 'Csf1r', 'Tnfrsf1b')
macrophages_wt_infected <- AddModuleScore(macrophages_wt_infected, features = list(cytokine_score), name = 'cytokine_score')
FeaturePlot(macrophages_wt_infected, features = 'cytokine_score1', reduction = 'wt.cerebrum.mac.umap') 

#Migratory markers
ccl_chemokines <- rownames(wt_cerebrum_macrophages[['RNA']]$data)[grepl('Ccl', rownames(wt_cerebrum_macrophages[['RNA']]$data))]
DotPlot(wt_cerebrum_macrophages, features = ccl_chemokines, group.by = 'Timepoint', scale = FALSE)+
  theme(axis.text.x = element_text(angle = 90))

DotPlot(wt_cerebrum_macrophages, features = c('Sell', 'Itgam', 'Tnf', 'Ccr1', 'Ccr2', 'Ccr9', 'Cx3cr1', 'Il10', 'Il18'), 
        group.by = 'Timepoint', scale = FALSE)+
  theme(axis.text.x = element_text(angle = 90))

#Markers from https://www.mdpi.com/1422-0067/25/22/12078#The_Characteristics_of_M%CF%86s
#Only look at day 5
wt_cerebrum_macrophages_day5 <-  subset(wt_cerebrum_macrophages, Timepoint == 'Day 5')
wt_cerebrum_macrophages_day4 <-  subset(wt_cerebrum_macrophages, Timepoint == 'Day 4')
wt_cerebrum_macrophages_day3 <-  subset(wt_cerebrum_macrophages, Timepoint == 'Day 3')

#Prepare UMAP for macrophages at all time points
wt_cerebrum_macrophages_day5 <- prepSeuratObj(wt_cerebrum_macrophages_day5)
ElbowPlot(wt_cerebrum_macrophages_day5, ndims = 40)
wt_cerebrum_macrophages_day5 <- prepUmapSeuratObj(wt_cerebrum_macrophages_day5, nDims = 20, reductionName = 'day5_macs_wt',
                                             resolution_value = 0.8)

wt_cerebrum_macrophages_day4 <- prepSeuratObj(wt_cerebrum_macrophages_day4)
ElbowPlot(wt_cerebrum_macrophages_day4, ndims = 40)
wt_cerebrum_macrophages_day4 <- prepUmapSeuratObj(wt_cerebrum_macrophages_day4, nDims = 20, reductionName = 'day4_macs_wt',
                                                  resolution_value = 0.8)

wt_cerebrum_macrophages_day3 <- prepSeuratObj(wt_cerebrum_macrophages_day3)
ElbowPlot(wt_cerebrum_macrophages_day3, ndims = 40)
wt_cerebrum_macrophages_day3 <- prepUmapSeuratObj(wt_cerebrum_macrophages_day3, nDims = 15, reductionName = 'day3_macs_wt',
                                                  resolution_value = 0.8)

DimPlot(wt_cerebrum_macrophages_day5, reduction = 'day5_macs_wt')
DimPlot(wt_cerebrum_macrophages_day5, reduction = 'day5_macs_wt', group.by = 'Treatment')

DimPlot(wt_cerebrum_macrophages_day4, reduction = 'day4_macs_wt', group.by = 'Treatment' )
DimPlot(wt_cerebrum_macrophages_day3, reduction = 'day3_macs_wt', group.by = 'Treatment')

#M1
m1_genes = c('Tnf', 'Il1b', 'Il12a', 'Il23a', 'Nos2', 'Cd68', 'Cd80','Cd86', 'Fcgr1', 'Fcgr2b', 'Fcgr3')

#M1 day 5
plotList_m1_day5 <- lapply(m1_genes, featurePlotLight, data = wt_cerebrum_macrophages_day5, 
                   reduction_choice = 'day5_macs_wt', maxLim = 6)
do.call(ggarrange, c(plotList_m1_day5, common.legend = TRUE, legend = 'right'))

#M1 day 4
plotList_m1_day4 <- lapply(m1_genes, featurePlotLight, data = wt_cerebrum_macrophages_day4, 
                           reduction_choice = 'day4_macs_wt', maxLim = 6)
do.call(ggarrange, c(plotList_m1_day4, common.legend = TRUE, legend = 'right'))

#Could subset to infected cells too for these
DotPlot(wt_cerebrum_macrophages, features = c('Tnf', 'Il1b', 'Il12a', 'Il23a', 'Nos2', 'Cd68', 'Cd80',
                                              'Cd86', 'Fcgr1', 'Fcgr2b', 'Fcgr3'), group.by = 'Timepoint', scale = FALSE) +
  theme(axis.text.x = element_text(angle = 90))

#Mature macrophage signature per https://www.cell.com/article/S1074-7613(14)00235-0/fulltext . Maybe
#transcripts expressed day 3 indicate divergence into macrophages leading to expression of genes like Nos2 above
#Got to read this too https://www.nature.com/articles/ni.2419
DotPlot(wt_cerebrum_macrophages, features = c('Mertk'), group.by = 'Timepoint', scale = FALSE) +
       theme(axis.text.x = element_text(angle = 90))

#M2a
FeaturePlot(wt_cerebrum_macrophages_day5, features = c('Mrc1', 'Cd68', 'Cd163', 'Arg1', 'Il10', 'Tgfb1', 'Il1b') , reduction = 'day5_macs_wt')

DotPlot(wt_cerebrum_macrophages, features = c('Mrc1', 'Cd68', 'Cd163', 'Arg1', 'Il10', 'Tgfb1', 'Il1b'), group.by = 'Timepoint', scale = FALSE) +
  theme(axis.text.x = element_text(angle = 90))


#M2b
FeaturePlot(wt_cerebrum_macrophages_day5, features = c('Cd86', 'Cd68', 'Mrc1', 'Il10',
                                                       'Il1b', 'Il6', 'Tnf') , reduction = 'day5_macs_wt')

DotPlot(wt_cerebrum_macrophages, features = c('Cd86', 'Cd68', 'Mrc1', 'Il10',
                                              'Il1b', 'Il6', 'Tnf'), group.by = 'Timepoint', scale = FALSE) +
  theme(axis.text.x = element_text(angle = 90))

#M2c
FeaturePlot(wt_cerebrum_macrophages_day5, features = c('Cd163', 'Cd68', 'Mrc1', 'Arg1',
                                                       'Il10', 'Tgfb1') , reduction = 'day5_macs_wt')

DotPlot(wt_cerebrum_macrophages, features = c('Cd163', 'Cd68', 'Mrc1', 'Arg1',
                                              'Il10', 'Tgfb1'), group.by = 'Timepoint', scale = FALSE) +
  theme(axis.text.x = element_text(angle = 90))

#M2d
FeaturePlot(wt_cerebrum_macrophages_day5, features = c('Cd68', 'Mrc1', 'Il10') , reduction = 'day5_macs_wt')

DotPlot(wt_cerebrum_macrophages, features = c('Cd68', 'Mrc1', 'Il10'), group.by = 'Timepoint', scale = FALSE) +
  theme(axis.text.x = element_text(angle = 90))

#Try dotplot with all above markers

all_markers_sorted <- c('Il12a', 'Il23a', 'Nos2', 'Cd80', 'Fcgr1', 'Fcgr2b', 'Tnf', 'Il1b', 'Cd86',
                        'Cd68', 'Il6', 'Mrc1', 'Cd163', 'Arg1', 'Il10', 'Tgfb1')

#Match markers with which mac type they identify
mac_sorted_type <- c('M1', 'M1', 'M1', 'M1', 'M1', 'M1', 'M1/M2', 'M1/M2', 'M1/M2', 'M1/M2', 'M1/M2', 'M2 - multiple', 'M2 - multiple',
                     'M2 - multiple', 'M2 - multiple', 'M2 - multiple')

all_surface_markers_sorted <-  c('Cd80', 'Fcgr2b', 'Fcgr1', 'Nos2', 'Il1r1', 'Tlr2', 'Tlr4', 'Cd86',
                                 'Cd68', 'Mrc1', 'Cd163', 'Arg1', 'Cd81','Siglec1', 'Itgam', 'Vcam1',
                                 'H2-Ab1', 'H2-Aa', 'H2-Eb1')

all_cytokines_sorted <- c('Il12a', 'Il18', 'Il23a', 'Il1a', 'Nos2', 'Arg1', 'Tgfb1', 'Il1b', 'Tnf', 'Il6',
                          'Il10', 'Ccl2', 'Ccl5', 'Ccl22', 'Vegfa')


mac_identifier_type_colors <- case_when(mac_sorted_type == 'M1'~ 'orange',
                                        mac_sorted_type == 'M1/M2'~ 'green',
                                        mac_sorted_type == 'M2 - multiple'~ 'darkred')

DotPlot(wt_cerebrum_macrophages, features = all_markers_sorted, group.by = 'Timepoint', scale = FALSE) +
  theme(axis.text.x = element_text(angle = 90))+
  geom_point(aes(size = pct.exp), color = rep(mac_identifier_type_colors, 3),  pch = 21) #pch > 20 lets you outline points

DotPlot(wt_cerebrum_macrophages, features = all_surface_markers_sorted, group.by = 'Timepoint', scale = FALSE) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))+
  ggtitle('Surface markers')+
  coord_flip()

DotPlot(wt_cerebrum_macrophages, features = all_cytokines_sorted, group.by = 'Timepoint', scale = FALSE) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))+
  ggtitle('Cytokines')+
  coord_flip()

#Markers from feb 4 meeting w anna
monocyte_markers_feb4 <- c('Lyz2', 'Ctss', 'Fcgr1', 'Lgals3', 'Tyrobp', 'Aif1',
                           'Ly6c2', 'Ccr2', 'S100a8', 'S100a9', 'Plac8')
nonclassic_monocyte_markers_feb4 <- c('Cx3cr1', 'Nr4a1', 'Ifitm3', 'Ltb')

DotPlot(wt_cerebrum_macrophages, features = monocyte_markers_feb4, group.by = 'Timepoint', scale = FALSE) +
  theme(axis.text.x = element_text(angle = 90))+
  ggtitle('Monocyte markers')+
  coord_flip()

plotList_classic_mono <- lapply(monocyte_markers_feb4, featurePlotLight, data = wt_cerebrum_macrophages, 
                       reduction_choice = 'wt.infected.mac.umap', maxLim = 6)
do.call(ggarrange, c(plotList_classic_mono, common.legend = TRUE, legend = 'right'))


DotPlot(wt_cerebrum_macrophages, features = nonclassic_monocyte_markers_feb4, group.by = 'Timepoint', scale = FALSE) +
  theme(axis.text.x = element_text(angle = 90))+
  ggtitle('Nonclassical monocyte markers')+
  coord_flip()

plotList_nonclassic_mono <- lapply(nonclassic_monocyte_markers_feb4, featurePlotLight, data = wt_cerebrum_macrophages, 
                                reduction_choice = 'wt.cerebrum.mac.umap', maxLim = 6)
do.call(ggarrange, c(plotList_nonclassic_mono, common.legend = TRUE, legend = 'right'))


#From Li et al Fig 1 https://pmc.ncbi.nlm.nih.gov/articles/PMC6542613/
#M1 markers
li_type1_markers <- lapply(c('Nos2', 'Cd86', 'Tnf'), featurePlotLight, data = wt_cerebrum_macrophages, 
                       reduction_choice = 'wt.cerebrum.mac.umap', maxLim = 6)
do.call(ggarrange, c(li_type1_markers, common.legend = TRUE, legend = 'right'))
#M2 markers
li_type2_markers <- lapply(c('Chil3', 'Arg1', 'Mgl2', 'Retnla'), featurePlotLight, data = wt_cerebrum_macrophages, 
                           reduction_choice = 'wt.cerebrum.mac.umap', maxLim = 6)
do.call(ggarrange, c(li_type2_markers, common.legend = TRUE, legend = 'right'))


#macropahge linaege markers
FeaturePlot(wt_cerebrum_macrophages, features = c('Ptprc') , reduction = 'wt.cerebrum.mac.umap')
FeaturePlot(wt_cerebrum_macrophages, features = c('Lyz2') , reduction = 'wt.cerebrum.mac.umap')
FeaturePlot(wt_cerebrum_macrophages, features = c('Adgre1') , reduction = 'wt.cerebrum.mac.umap')

#Other markers from paper
FeaturePlot(wt_cerebrum_macrophages_day5, features = c('Siglec1') , reduction = 'day5_macs_wt')
FeaturePlot(wt_cerebrum_macrophages_day5, features = c('H2-Ab1', 'H2-Aa', 'H2-Eb1') , reduction = 'day5_macs_wt')
FeaturePlot(wt_cerebrum_macrophages_day5, features = c('Vegfa') , reduction = 'day5_macs_wt')
FeaturePlot(wt_cerebrum_macrophages_day5, features = c('Jak2') , reduction = 'day5_macs_wt')

#Markers from https://pmc.ncbi.nlm.nih.gov/articles/PMC3666862/ (which looked at monocytes vs macrophages in humans)
#Mac markers
FeaturePlot(wt_cerebrum_macrophages_day5, features = c('Selenop') , reduction = 'day5_macs_wt')
FeaturePlot(wt_cerebrum_macrophages_day5, features = c('C1qc') , reduction = 'day5_macs_wt')
FeaturePlot(wt_cerebrum_macrophages_day5, features = c('Apoe') , reduction = 'day5_macs_wt')

#Monocyte markers
FeaturePlot(wt_cerebrum_macrophages_day5, features = c('Dnm2') , reduction = 'day5_macs_wt')
FeaturePlot(wt_cerebrum_macrophages_day5, features = c('Ccr2') , reduction = 'day5_macs_wt')
FeaturePlot(wt_cerebrum_macrophages_day5, features = c('Itgb2') , reduction = 'day5_macs_wt')

#More markers... From https://link.springer.com/article/10.1186/s12974-016-0581-z
M_il4_macs <- c('Arg1', 'Mrc1', 'Chil3', 'Tgm2', 'Il1rn', 'Msr1', 'Pdcd1lg2', 'Cd209a', 'Clec10a', 'Retnla', 'Alox15',
                'Socs2', 'Irf4', 'Ppard', 'Pparg', 'Ccl17', 'Ccl24')
mhc_2_macs <- c('H2-Aa', 'H2-Ab1', 'H2-DMa', 'H2-DMb1', 'H2-DMb2', 'H2-Eb1')
M_lps_ifng_macs <- c('Tnf', 'Il1b', 'Il6', 'Il12a', 'Il23a', 'Il27', 'Nos2', 'Ido1', 'Irf5', 'Nfkblz',
                     'Socs1', 'Marco', 'Cd80', 'Cd86', 'Cd274', 'C1d', 'C1qa', 'C1qb', 'C1qbp', 'C1qc',
                     'C1qtnf1', 'C3', 'C3ar1', 'C5ar1', 'Itgb2', 'Ccl2', 'Ccl3', 'Ccl4', 'Ccl5', 'Cxcl9',
                     'Cxcl10', 'Cxcl11', 'Cxcl16', 'Ccr7')
M_il10_macs <- c('Il10', 'Il4ra', 'Tgfb1', 'Nfil3', 'Sbno2', 'Socs3', 'Fcgr1', 'Fcgr2b', 'Fcgr3', 'Marco',
                 'Cxcl13')
M_ic_macs <- c('Il10', 'Il6', 'Nos2', 'Ccl1', 'Ccl20', 'Cxcl3', 'Cxcl13')


DotPlot(wt_cerebrum_macrophages, features = M_il4_macs, group.by = 'Timepoint', scale = FALSE) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))+
  ggtitle('M IL4 (alternative activation) mac markers')+
  coord_flip()

DotPlot(wt_cerebrum_macrophages, features = mhc_2_macs, group.by = 'Timepoint', scale = FALSE) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))+
  ggtitle('MHC-2 mac markers')+
  coord_flip()

DotPlot(wt_cerebrum_macrophages, features = M_lps_ifng_macs, group.by = 'Timepoint', scale = FALSE) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))+
  ggtitle('LPS Ifng (classical activation) mac markers')+
  coord_flip()

DotPlot(wt_cerebrum_macrophages, features = M_il10_macs, group.by = 'Timepoint', scale = FALSE) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))+
  ggtitle('Il10 (M2c, tissue remodeling) mac markers')+
  coord_flip()

DotPlot(wt_cerebrum_macrophages, features = M_ic_macs, group.by = 'Timepoint', scale = FALSE) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))+
  ggtitle('Ic mac markers')+
  coord_flip()

#Feature plot genes
#Removing genes with no expression
plotList_il4 <- lapply(M_il4_macs[!M_il4_macs %in% c('Retnla', 'Alox15', 'Ccl17')], featurePlotLight, data = wt_cerebrum_macrophages, 
                   reduction_choice = 'wt.cerebrum.mac.umap', maxLim = 6)
do.call(ggarrange, c(plotList_il4, common.legend = TRUE, legend = 'right'))

plotList_lps <- lapply(M_lps_ifng_macs[!M_lps_ifng_macs %in% c('Nfkblz', 'Marco', 'Cxcl11')], featurePlotLight, data = wt_cerebrum_macrophages, 
                       reduction_choice = 'wt.cerebrum.mac.umap', maxLim = 6)
do.call(ggarrange, c(plotList_lps, common.legend = TRUE, legend = 'right'))

plotList_mhc <- lapply(mhc_2_macs, featurePlotLight, data = wt_cerebrum_macrophages, 
                       reduction_choice = 'wt.cerebrum.mac.umap', maxLim = 6)
do.call(ggarrange, c(plotList_mhc, common.legend = TRUE, legend = 'right'))

plotList_il10 <- lapply(M_il10_macs[!M_il10_macs %in% c('Marco')], featurePlotLight, data = wt_cerebrum_macrophages, 
                       reduction_choice = 'wt.cerebrum.mac.umap', maxLim = 6)
do.call(ggarrange, c(plotList_il10, common.legend = TRUE, legend = 'right'))


plotList_ic <- lapply(M_ic_macs, featurePlotLight, data = wt_cerebrum_macrophages, 
                        reduction_choice = 'wt.cerebrum.mac.umap', maxLim = 6)
do.call(ggarrange, c(plotList_ic, common.legend = TRUE, legend = 'right'))

#Look at functional differences between seurat clusters
#Examine what separates day 5 from other timepoints
macro_markers <- FindAllMarkers(macrophages_wt_infected, group.by = 'seurat_clusters',
                                  test.use = 'MAST')
macro_markers_up <- dplyr::filter(macro_markers, (avg_log2FC) > 1 & p_val_adj < 0.01)

cluster0_paths <- gprofiler2::gost(query = macro_markers_up[macro_markers_up$cluster == 3,]$gene, organism = 'mmusculus', evcodes = TRUE)
cluster0_paths$result[1,]

FeaturePlot(wt_cerebrum_macrophages, features = c('Lyz2') , reduction = 'wt.cerebrum.mac.umap')

