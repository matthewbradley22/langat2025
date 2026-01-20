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
macrophages_wt_infected <- subset(wt_cerebrum_macrophages, Treatment == 'rLGTV')  

#Prepare UMAP for macrophages
macrophages_wt_infected <- prepSeuratObj(macrophages_wt_infected)
ElbowPlot(macrophages_wt_infected, ndims = 40)

#Higher num neighbors for fewer clusters
macrophages_wt_infected <- prepUmapSeuratObj(macrophages_wt_infected, nDims = 20, reductionName = 'wt.infected.mac.umap',
                                             resolution_value = 0.8)

DimPlot(macrophages_wt_infected, reduction = 'wt.infected.mac.umap', label = TRUE, group.by = 'seurat_clusters',
        label.size = 6)+
  ggtitle('WT Infected Macrophages')

DimPlot(macrophages_wt_infected, reduction = 'wt.infected.mac.umap', label = FALSE, group.by = 'Timepoint',
        label.size = 6)+
  ggtitle('WT Infected Macrophages')

#Examine what separates day 5 from other timepoints
day5_macro_markers <- FindMarkers(macrophages_wt_infected, group.by = 'Timepoint', ident.1 = 'Day 5',
                                  test.use = 'MAST')
day5_up_markers <- dplyr::filter(day5_macro_markers, (avg_log2FC) > 1 & p_val_adj < 0.01)

#Gene ontology 
day5_paths <- gprofiler2::gost(query = rownames(day5_up_markers), organism = 'mmusculus', evcodes = TRUE)
day5_paths$result

#Marker paper: https://pmc.ncbi.nlm.nih.gov/articles/PMC12362712/
FeaturePlot(macrophages_wt_infected, features = 'Nos2', reduction = 'wt.infected.mac.umap')
FeaturePlot(macrophages_wt_infected, features = 'Apoe', reduction = 'wt.infected.mac.umap') #Phagocytosis
FeaturePlot(macrophages_wt_infected, features = 'Mrc1', reduction = 'wt.infected.mac.umap') #bam marker/alternate activation

FeaturePlot(macrophages_wt_infected, features = 'Cxcl2', reduction = 'wt.infected.mac.umap') #Oxidative stress
FeaturePlot(macrophages_wt_infected, features = 'Prdx5', reduction = 'wt.infected.mac.umap')
FeaturePlot(macrophages_wt_infected, features = 'Txn1', reduction = 'wt.infected.mac.umap')
FeaturePlot(macrophages_wt_infected, features = 'Gsr', reduction = 'wt.infected.mac.umap')

FeaturePlot(macrophages_wt_infected, features = 'H2-Aa', reduction = 'wt.infected.mac.umap') #Antigen presentation
FeaturePlot(macrophages_wt_infected, features = 'H2-Eb1', reduction = 'wt.infected.mac.umap')
FeaturePlot(macrophages_wt_infected, features = 'H2-Ab1', reduction = 'wt.infected.mac.umap')

#Complement cascade
FeaturePlot(macrophages_wt_infected, features = 'C1qc', reduction = 'wt.infected.mac.umap')
FeaturePlot(macrophages_wt_infected, features = 'C1qb', reduction = 'wt.infected.mac.umap')

#Other interesting genes
FeaturePlot(macrophages_wt_infected, features = 'F10', reduction = 'wt.infected.mac.umap')
FeaturePlot(macrophages_wt_infected, features = 'Cd38', reduction = 'wt.infected.mac.umap')
FeaturePlot(macrophages_wt_infected, features = 'Ccr1', reduction = 'wt.infected.mac.umap')
FeaturePlot(macrophages_wt_infected, features = 'Socs3', reduction = 'wt.infected.mac.umap')

#Monocyte score
monocyte_score_genes <- c('S100a4', 'Itgb7', 'Napsa', 'Cd300lg', 'Adora2b', 'Emb', 'Ly6c2', 'Ms4a4c', 'Fn1', 'Sell', 'Padi2', 'Lilra6', 
                          'Ccnb2', 'Galnt9', 'Upb1', 'Lmo1', 'F13a1', 'Ccr2', 'Gm15987', 'AI839979')
macrophages_wt_infected <- AddModuleScore(macrophages_wt_infected, features = list(monocyte_score_genes), name = 'monocyte_score')
FeaturePlot(macrophages_wt_infected, features = 'monocyte_score1', reduction = 'wt.infected.mac.umap')

#Macrophage score
macrophage_score_genes <- c('Adgre1', 'Csf1r', 'H2-Ab1', 'Cd68', 'Lyz2', 'Itgam', 'Mertk')
macrophages_wt_infected <- AddModuleScore(macrophages_wt_infected, features = list(macrophage_score_genes), name = 'macrophage_score')
FeaturePlot(macrophages_wt_infected, features = 'macrophage_score1', reduction = 'wt.infected.mac.umap')

#Other macrophage scores
#Phagocytosis
phagocytosis_score <- c('Aif1', 'Cd36', 'Coro1a', 'Ccr2', 'Cnn2', 'Cyba', 'Fcer1g', 'Fcgr1', 'Fcgr3', 'Gas6', 'Irf8', 'Il2rg', 'Itgb1', 
                        'Itgb2', 'Lrp1', 'Ncf4', 'Slc11a1', 'Pros1', 'Msr1', 'Thbs1', 'Tyrobp', 'Pycard', 'Cd209b', 'Trem2', 'Gsn', 'P2ry6', 'Myo1g')
macrophages_wt_infected <- AddModuleScore(macrophages_wt_infected, features = list(phagocytosis_score), name = 'phagocytosis_score')
FeaturePlot(macrophages_wt_infected, features = 'phagocytosis_score1', reduction = 'wt.infected.mac.umap')

#Inflammatory response
inflammatory_score <- c('Adam8', 'Aif1', 'Alox5', 'Alox5ap', 'Fabp4', 'Apoe', 'App', 'C1qa', 'Ciita', 'C3ar1', 'C5ar1', 'Cd14', 'Cd44', 'Cd68', 'Ccr2', 
'Csf1r', 'Ctsc', 'Cyba', 'Cybb', 'Ednrb', 'Fcer1g', 'Fcgr1', 'Fcgr3', 'Fn1', 'B4galt1', 'Grn', 'Cxcl1', 'Hp', 'Igf1', 'Il6', 'Itgam', 'Itgb2', 'Jun', 'Lpl',
'Lrp1', 'Mif', 'Nfkb1', 'Slc11a1', 'Ptger4', 'Rel', 'Ccl6', 'Cxcl2', 'Thbs1', 'Tnfrsf1b', 'Tyrobp', 
'Pik3cg', 'Cxcl13', 'Ccl24', 'Pf4', 'Pycard', 'Trem2', 'Cd163', 'Gpsm3', 'Tlr7', 'Stab1')
macrophages_wt_infected <- AddModuleScore(macrophages_wt_infected, features = list(inflammatory_score), name = 'inflammatory_score')
FeaturePlot(macrophages_wt_infected, features = 'inflammatory_score1', reduction = 'wt.infected.mac.umap')

#Oxidative stress
oxidative_stress_score <- c('Prdx5', 'Txn1', 'Gsr', 'Ptgs2', 'Ccs', 'Prdx6', 'Gpx4', 'Sesn1', 'Sod3', 'Ltc4s')
macrophages_wt_infected <- AddModuleScore(macrophages_wt_infected, features = list(oxidative_stress_score), name = 'oxidative_stress_score')
FeaturePlot(macrophages_wt_infected, features = 'oxidative_stress_score1', reduction = 'wt.infected.mac.umap')

#ECM Organization
ecm_score <- c('Col1a1', 'Nid1', 'Dpt', 'B4galt1', 'Lum', 'Col3a1', 'Ccdc80', 'Ramp2', 'Serpinh1', 'Ddr2')
macrophages_wt_infected <- AddModuleScore(macrophages_wt_infected, features = list(ecm_score), name = 'ecm_score')
FeaturePlot(macrophages_wt_infected, features = 'ecm_score1', reduction = 'wt.infected.mac.umap') +
  ggtitle('ECM Organization score')

#Cytokin cytokine receptor interacitons
cytokine_score <- c('Ccl6', 'Ccl24', 'Cxcl1', 'Cxcl2', 'Pf4', 'Cxcl13', 'Cxcl12', 'Cxcl16', 'Il6', 'Ccr2', 'Il2rg', 'Csf1r', 'Tnfrsf1b')
macrophages_wt_infected <- AddModuleScore(macrophages_wt_infected, features = list(cytokine_score), name = 'cytokine_score')
FeaturePlot(macrophages_wt_infected, features = 'cytokine_score1', reduction = 'wt.infected.mac.umap') 

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

#Prepare UMAP for macrophages
wt_cerebrum_macrophages_day5 <- prepSeuratObj(wt_cerebrum_macrophages_day5)
ElbowPlot(wt_cerebrum_macrophages_day5, ndims = 40)

#Higher num neighbors for fewer clusters
wt_cerebrum_macrophages_day5 <- prepUmapSeuratObj(wt_cerebrum_macrophages_day5, nDims = 20, reductionName = 'day5_macs_wt',
                                             resolution_value = 0.8)

DimPlot(wt_cerebrum_macrophages_day5, reduction = 'day5_macs_wt')
DimPlot(wt_cerebrum_macrophages_day5, reduction = 'day5_macs_wt', group.by = 'Treatment')

#M1
FeaturePlot(wt_cerebrum_macrophages_day5, features = c('Tnf', 'Il1b', 'Il12a', 'Il23a', 'Nos2', 'Cd68', 'Cd80',
                                                  'Cd86', 'Fcgr1', 'Fcgr2b', 'Fcgr3') , reduction = 'day5_macs_wt')

#M2a
FeaturePlot(wt_cerebrum_macrophages_day5, features = c('Mrc1', 'Cd68', 'Cd163', 'Arg1', 'Il10', 'Tgfb1', 'Il1b') , reduction = 'day5_macs_wt')

#M2b
FeaturePlot(wt_cerebrum_macrophages_day5, features = c('Cd86', 'Cd68', 'Mrc1', 'Il10',
                                                       'Il1b', 'Il6', 'Tnf') , reduction = 'day5_macs_wt')

#M2c
FeaturePlot(wt_cerebrum_macrophages_day5, features = c('Cd163', 'Cd68', 'Mrc1', 'Arg1',
                                                       'Il10', 'Tgfb1') , reduction = 'day5_macs_wt')

#M2d
FeaturePlot(wt_cerebrum_macrophages_day5, features = c('Cd68', 'Mrc1', 'Il10') , reduction = 'day5_macs_wt')

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



