#Load libraries
library(dplyr)
library(Seurat)
library(ggplot2)
library(SingleR)
library(celldex)
library(scran)
library(GSEABase)
library(AUCell)
library(RColorBrewer)

### SUBSET TO JUST WILDTYPE ###
#Source function
source('~/Documents/ÖverbyLab/scripts/langatFunctions.R')

#This is where the 10x data is
setwd('~/Documents/ÖverbyLab/single_nuclei_proj/')

#Read in processed data
sn_integrated_dat <- LoadSeuratRds('~/Documents/ÖverbyLab/single_nuclei_proj/LGTVscCombined.rds')

#Custom pathways
#Gene lists
necroptosis <- c('Tnf', 'Tnfrsf1a', 'Ripk2', 'Mlkl', 'Ripk1', 'Ripk3')
pyroptosis <- c('Gsdmc', 'Nlrp3', 'Aim2', 'Gsdmd', 'Il18', 'Il1b', 'Casp9', 'Casp8', 'Casp6', 'Casp3', 'Casp4', 'Casp1')

#Check data
newCols <-  c(brewer.pal(12, 'Paired'), '#99FFE6', '#CE99FF', '#18662E','#737272',  '#FF8AEF')
DimPlot(sn_integrated_dat, group.by =   'manualAnnotation', cols = newCols, reduction = 'umap.integrated')

#Look at infection and see if labels make sense. Why not more infection in ifnar knockout?
FeaturePlot(sn_integrated_dat, features = 'LGTV', reduction = 'umap.integrated', split.by = 'new_inf')
FeaturePlot(sn_integrated_dat, features = 'rna_LGTV', reduction = 'umap.integrated',  split.by = 'new_genotype')
table(sn_integrated_dat$new_genotype, sn_integrated_dat$new_inf)

#subset to wt
sn_integrated_dat_wt <- subset(sn_integrated_dat, new_genotype %in% c('wt', 'wt (same)'))

sn_integrated_dat_wt <- prepSeuratObj(sn_integrated_dat_wt)
ElbowPlot(sn_integrated_dat_wt) 
sn_integrated_dat_wt <- prepUmapSeuratObj(sn_integrated_dat_wt, nDims = 20,num_neighbors = 30L, 
                                          reductionName = 'wt.umap.integrated',
                                          resolution_value = 0.8)

DimPlot(sn_integrated_dat_wt, group.by =  'manualAnnotation', cols = newCols, reduction = 'wt.umap.integrated')
DimPlot(sn_integrated_dat_wt, group.by =  'infected', cols = newCols, reduction = 'wt.umap.integrated')
FeaturePlot(sn_integrated_dat_wt, features = 'rna_LGTV', reduction = 'wt.umap.integrated')

#Pathway analysis

#Use MAST to test for differences
treatment_markers <- FindAllMarkers(sn_integrated_dat_wt, group.by = 'infected', test.use = 'MAST', 
                                    only.pos = TRUE)
treatment_markers[necroptosis,]
treatment_markers[pyroptosis,]

############Resident pathway enrichment analysis############
############################################################

#Look at differentially expressed pathways from MAST results
#Upregulated in infection

upregulated_infection <- subset(treatment_markers, avg_log2FC > 1 & p_val_adj < 0.01 & cluster == 'rLGTV')

#Set high p value threshold to see all pathways
upregulated_infection_paths <- gprofiler2::gost(query = rownames(upregulated_infection), organism = 'mmusculus', evcodes = TRUE,
                                                user_threshold = 1)

upregulated_infection_paths$result[upregulated_infection_paths$result$source == 'KEGG',]
upregulated_infection_paths$result[upregulated_infection_paths$result$source == 'GO:MF',]
upregulated_infection_paths$result[upregulated_infection_paths$result$source == 'GO:BP',]







