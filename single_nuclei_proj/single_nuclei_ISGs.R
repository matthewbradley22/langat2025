#Load libraries 
library(dplyr)
library(Seurat)
library(ggplot2)
library(SingleR)
library(celldex)
library(scran)
library(RColorBrewer)

### SUBSET TO JUST WILDTYPE ###
#Source function
source('~/Documents/ÖverbyLab/scripts/langatFunctions.R')

#This is where the 10x data is
setwd('~/Documents/ÖverbyLab/single_nuclei_proj/')

#Read in processed data
sn_integrated_dat <- LoadSeuratRds('~/Documents/ÖverbyLab/single_nuclei_proj/LGTVscCombined.rds')

#Check data
newCols <-  c(brewer.pal(12, 'Paired'), '#99FFE6', '#CE99FF', '#18662E','#737272',  '#FF8AEF')
DimPlot(sn_integrated_dat, group.by =   'manualAnnotation', cols = newCols, reduction = 'umap.integrated')

#subset to wt
sn_integrated_dat_wt <- subset(sn_integrated_dat, new_genotype %in% c('wt', 'wt (same)'))

sn_integrated_dat_wt <- prepSeuratObj(sn_integrated_dat_wt)
ElbowPlot(sn_integrated_dat_wt) 
sn_integrated_dat_wt <- prepUmapSeuratObj(sn_integrated_dat_wt, nDims = 20,num_neighbors = 30L, 
                                          reductionName = 'wt.umap.integrated',
                                          resolution_value = 0.8)

DimPlot(sn_integrated_dat_wt, group.by =  'manualAnnotation', cols = newCols, reduction = 'wt.umap.integrated')+
  xlab('Umap1')+
  ylab('Umap1')+
  theme(axis.ticks = element_blank(),
        axis.text=element_blank(),
        legend.text=element_text(size=17))+
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank())+
  ggtitle('single nuclei WT')

#Load in ISGs
#Molecular signatures database
#Should compare using db_species vs not using it
all_gene_sets <- msigdbr(species = "Mus musculus", db_species = "MM")

mouse_gene_sets <- all_gene_sets %>%
  split(x = .$gene_symbol, f = .$gs_name)

#Get total number of positive DEGs as above, but subset to ISGs
ifnA_response <- mouse_gene_sets$HALLMARK_INTERFERON_ALPHA_RESPONSE
ifnA_GOBP_response <- mouse_gene_sets$GOBP_RESPONSE_TO_INTERFERON_ALPHA
type1_response <- mouse_gene_sets$GOBP_RESPONSE_TO_TYPE_I_INTERFERON
#reactome_ifn_antiviral <- mouse_gene_sets$REACTOME_ANTIVIRAL_MECHANISM_BY_IFN_STIMULATED_GENES
all_ISGs_type1 = unique(c(ifnA_response, ifnA_GOBP_response, type1_response))

#ISG module score
sn_integrated_dat_wt <- AddModuleScore(sn_integrated_dat_wt, features = list(all_ISGs_type1), name = 'ISG_score')

#Split by infection for plotting
sn_integrated_dat_wt_inf <- subset(sn_integrated_dat_wt, new_inf == 'LGTV')
sn_integrated_dat_wt_uninf <- subset(sn_integrated_dat_wt, new_inf != 'LGTV')

#Look at some isg dotplots with cell type and infection to see if astrocytes are more active than in single cell data
DotPlot(sn_integrated_dat_wt_inf, features = c(all_ISGs_type1[30:50], 'Rsad2'), group.by = 'manualAnnotation', scale = FALSE)+
  theme(axis.text.x = element_text(angle = 90))+
  ggtitle('Single Nuclei infected ISGs sample')

DotPlot(sn_integrated_dat_wt_uninf, features = c(all_ISGs_type1[30:50], 'Rsad2'), group.by = 'manualAnnotation')+
  theme(axis.text.x = element_text(angle = 90))+
  ggtitle('Single Nuclei uninfected ISGs samples')

sn_isg_treatment_dat <- DotPlot(sn_integrated_dat_wt, features = 'ISG_score1', group.by = 'treatment_celltype')$data
infection_celltype <- str_split_fixed(sn_isg_treatment_dat$id, " ", 2)
colnames(infection_celltype) <- c('infection', 'celltype')
sn_isg_treatment_dat_final <- cbind(sn_isg_treatment_dat, infection_celltype)
ggplot(sn_isg_treatment_dat_final, aes(x = infection, y = celltype,  fill = avg.exp))+
  geom_tile()+
  scale_fill_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white", "lightblue"), 
                        values = c(1.0,0.7,0.4,0.2,0),
                        limits = c(-0.15, 0.55))+
  geom_text(aes(label=round(avg.exp, digits = 2)))+
  ggtitle('Single nuclei ISG scores')

#Arrest and transmigration genes from https://academic.oup.com/view-large/3854487
arrest_only_markers <- c('Vcam1', 'Ccr1')
transmigration_only_markers <- c('Mcam', 'F11r', 'Jam2', 'Pecam1', 'Pvr', 'Cd99', 'Cd99l2', 'Cdh5', 'Epha1', 'Ephb1')
all_arrest_markers <- c('Itgal', 'Itgb2', 'Itgam', 'Vcam1', 'Amica1', 'Selp', 'Ccr1')

#Plot arrest and transmigration markers
DotPlot(sn_integrated_dat_wt, features = arrest_only_markers, scale = FALSE, group.by = 'manualAnnotation')+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))+
  scale_size_continuous(range = c(0.5,6), limits = c(0,100))+
  scale_color_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white"), 
                        values = c(1.0,0.7,0.4,0),
                        limits = c(0,3))+
  ggtitle('Single nuclei arrest markers')

DotPlot(sn_integrated_dat_wt, features = transmigration_only_markers, scale = FALSE, group.by = 'manualAnnotation')+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))+
  scale_size_continuous(range = c(0.5,6), limits = c(0,100))+
  scale_color_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white"), 
                        values = c(1.0,0.7,0.4,0),
                        limits = c(0,3))+
  ggtitle('Single nuclei transmigration markers')
