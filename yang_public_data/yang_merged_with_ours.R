#Merging our data with yang's
#Could subset celltypes to make this faster if it's really slow

library(Seurat)
library(dplyr)
library(ggplot2)
source('~/Documents/ÖverbyLab//scripts/langatFunctions.R')
#Load our data
#Load data
ParseSeuratObj_int <- LoadSeuratRds("~/Documents/ÖverbyLab/data/FilteredRpcaIntegratedDatNoDoublets.rds") 

#Only using cells from gal3 project
wt_cerebrum <-  subset(ParseSeuratObj_int, Treatment %in% c('PBS', 'rLGTV') & Organ == 'Cerebrum' & Genotype == 'WT')

#Plot data
newCols <-  c(brewer.pal(12, 'Paired'), '#99FFE6', '#CE99FF', '#18662E','#737272',  '#FF8AEF')
newCols[11] =  '#FF8AEF'
DimPlot(ParseSeuratObj_int, label = FALSE, group.by = 'manualAnnotation', reduction = 'umap.integrated',
        cols = newCols)

#Load in Yang data
yang_data <- LoadSeuratRds("~/Documents/ÖverbyLab/yang_public_data/yang_rpca_integrated_obj.RDS")
newCols <-  c(brewer.pal(12, 'Paired'), '#99FFE6', '#CE99FF', '#18662E','#737272',  '#FF8AEF')
newCols[11] =  '#FF8AEF'

#Join layers for yang to match parse
yang_data <- JoinLayers(yang_data)

#Plot data
DimPlot(yang_data,reduction = "umap.rpca", label = FALSE, group.by = 'manualAnnotation', cols = newCols)

#Make treatment column names same case before merging
yang_data$Treatment = yang_data$treatment
yang_data$treatment = NULL

#load in single-nuclei data
#Read in processed data
sn_integrated_dat <- LoadSeuratRds('~/Documents/ÖverbyLab/single_nuclei_proj/LGTVscCombined.rds')

#Make treatment column names same case before merging
#Combine infected column for later
sn_integrated_dat$Treatment = sn_integrated_dat$treatment
sn_integrated_dat$treatment = NULL

sn_integrated_dat$infected_grouped = ifelse(sn_integrated_dat$new_inf %in% c('mock', 'none'), 'uninfected', 'infected')
#subset to only wt
sn_integrated_dat_wt <- subset(sn_integrated_dat, new_genotype %in% c('wt', 'wt (same)'))

sn_integrated_dat_wt <- prepSeuratObj(sn_integrated_dat_wt)
ElbowPlot(sn_integrated_dat_wt) 
sn_integrated_dat_wt <- prepUmapSeuratObj(sn_integrated_dat_wt, nDims = 20,num_neighbors = 30L, 
                                          reductionName = 'wt.umap.integrated',
                                          resolution_value = 0.8)

DimPlot(sn_integrated_dat_wt, group.by =  'manualAnnotation', cols = newCols, reduction = 'wt.umap.integrated')

#Create celltype subsets to integrate
parse_endo <- subset(wt_cerebrum, manualAnnotation == 'Endothelial')
parse_macro <- subset(wt_cerebrum, manualAnnotation == 'Macrophage/Monocytes')

yang_endo <- subset(yang_data, manualAnnotation == 'Endothelial')
yang_macro <- subset(yang_data, manualAnnotation == 'Macro/Mono')

sn_endo <- subset(sn_integrated_dat_wt, manualAnnotation == 'Endothelial')
sn_macro <- subset(sn_integrated_dat_wt, manualAnnotation == 'Macrophages')

########### Monocyte / macropahge analyses ########### 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 


#monocyte cells combined
macro.combined <- merge(parse_macro, y = c(yang_macro,sn_macro), add.cell.ids = c("parse", "yang", "sn"), project = "merged_macros")

macro.combined <- NormalizeData(macro.combined)
macro.combined <- FindVariableFeatures(macro.combined)
macro.combined <- ScaleData(macro.combined)
macro.combined <- RunPCA(macro.combined)

#Need to allow greater ram usage to run pca integration
options(future.globals.maxSize = 10000 * 1024^2)

#Integrate data
combined_macros <- IntegrateLayers(
  object = macro.combined, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca",
  verbose = FALSE
)

#Reset max ram to not accidentally use a bunch later without realizing
options(future.globals.maxSize = 500 * 1024^2)

ElbowPlot(combined_macros, ndims = 40)

#Prepare data for plotting
combined_macros <- FindNeighbors(combined_macros, reduction = "integrated.rpca", dims = 1:25)
combined_macros <- FindClusters(combined_macros, resolution = 0.8, cluster.name = "rpca_clusters")
combined_macros <- RunUMAP(combined_macros, reduction = "integrated.rpca", dims = 1:25, reduction.name = "umap.rpca")
combined_macros$manualAnnotation = 'macro/mono'
combined_macros$dataset <- gsub("_.+", '', colnames(combined_macros))
combined_macros$timepoint_or_treatment <- factor(dplyr::case_when(combined_macros$dataset == 'parse' & combined_macros$Treatment != 'PBS' ~ combined_macros$Timepoint,
                                                           combined_macros$dataset == 'parse' & combined_macros$Treatment == 'PBS' ~ combined_macros$Treatment,
                                                           combined_macros$dataset == 'yang' ~ combined_macros$Treatment,
                                                           combined_macros$dataset == 'sn' ~ combined_macros$infected_grouped,),
                                                 levels = c('PBS', 'Day 3', 'Day 4', 'Day 5', 'healthy', 'mild', 'moderate', 'severe', 'uninfected', 'infected'))
combined_macros$dataset_group <- factor(paste(combined_macros$dataset, combined_macros$timepoint_or_treatment, sep = '_'),
                                        levels = c('parse_PBS', 'parse_Day 3', 'parse_Day 4', 'parse_Day 5',
                                                   'yang_healthy', 'yang_mild','yang_moderate', 'yang_severe',
                                                   'sn_uninfected', 'sn_infected'))

#Plot datasets to compare umap locations
pdf('/Users/matthewbradley/Documents/ÖverbyLab/yang_public_data/combined_data_plots/combined_mac_plot.pdf', width = 6, height = 5)
DimPlot(combined_macros, label = FALSE, group.by = 'dataset_group', reduction = 'umap.rpca',
        cols = newCols)+
  ggtitle('Macrophage integrated UMAP')+
  xlab('')+
  ylab('')
dev.off()

#Remove sn_uninfected as these seem to be microglia (not grouped with most macrophages)
combined_macros <- subset(combined_macros, dataset_group != 'sn_uninfected')
DimPlot(combined_macros, label = FALSE, group.by = 'dataset_group', reduction = 'umap.rpca',
        cols = newCols)+
  ggtitle('Macrophage integrated UMAP')+
  xlab('')+
  ylab('')

#Look at seurat clusters, later look into what is separating groups
DimPlot(combined_macros, reduction = 'umap.rpca')
FeaturePlot(combined_macros, reduction = 'umap.rpca', features = 'Nos2')

#Function to make dotplots
comparative_dotplot <- function(data, genes, title, sc_timepoints = TRUE){
  comp_dot_dat <- DotPlot(data, features = genes, 
                          group.by = 'dataset_group', scale = FALSE)$data
  
  comp_meta <- str_split_fixed(comp_dot_dat$id, "_", 2)
  colnames(comp_meta) <- c('dataset', 'group')
  comp_dot_dat <- cbind(comp_dot_dat, comp_meta)
  
  if(sc_timepoints) {
    comp_dot_dat$group <- factor(comp_dot_dat$group, levels = c('PBS', 'Day 3', 'Day 4', 'Day 5',
                                                                'healthy', 'mild', 'moderate', 'severe',
                                                                'uninfected', 'infected'))
  }
  
  legend_max = max(comp_dot_dat$avg.exp.scaled)
  
  #Make plot
  dot_plot <- ggplot(comp_dot_dat, aes(x = group, y = features.plot))+
    facet_grid(~dataset, scales = 'free')+
    geom_point(aes(size = pct.exp, fill = avg.exp.scaled), pch = 21)+
    #coord_flip()+
    scale_size_continuous(range = c(0.5,6), limits = c(0,100))+
    scale_fill_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white"), 
                         values = c(1.0,0.7,0.4,0),
                         limits = c(0,legend_max))+
    ggtitle(title)+
    theme_bw()+
    theme(panel.grid = element_blank(),
          #panel.border = element_rect(colour = "white", fill = NA),
          panel.spacing = unit(0, "line"),
          strip.background = element_blank(),
          axis.text.x = element_text(angle = 90))+
    scale_size(range = c(0,9), limits = c(0, 100))+
    ylab('')+
    xlab('')
  
  plot(dot_plot)
}

pdf('/Users/matthewbradley/Documents/ÖverbyLab/yang_public_data/combined_data_plots/mono_mhc2.pdf', width = 6, height = 5)
comparative_dotplot(data = combined_macros, genes = c('H2-Aa', 'H2-Ab1', 'H2-DMa', 'H2-DMb1', 'H2-DMb2', 'H2-Eb1'), title = 'Monocyte MHC2 genes')
dev.off()

pdf('/Users/matthewbradley/Documents/ÖverbyLab/yang_public_data/combined_data_plots/mono_trafficking_genes.pdf', width = 6, height = 5)
comparative_dotplot(data = combined_macros, genes = c('Pecam1', 'Pvr', 'Cd99l2', 'Epha1', 'Ephb1', 'Itga4', 'Itgb1', 'Ccr1', 'Ccr2', 'Ccr5'), title = 'Monocyte trafficking genes')
dev.off()

#Glance at other genes
comparative_dotplot(data = combined_macros, genes = c('Nos2'), title = 'Nos2')

########### Endothelial analyses ########### 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

#endothelial cells combined
endo.combined <- merge(parse_endo, y = c(yang_endo, sn_endo), add.cell.ids = c("parse", "yang", "sn"), project = "merged_endos")

endo.combined <- NormalizeData(endo.combined)
endo.combined <- FindVariableFeatures(endo.combined)
endo.combined <- ScaleData(endo.combined)
endo.combined <- RunPCA(endo.combined)

#Need to allow greater ram usage to run pca integration
options(future.globals.maxSize = 10000 * 1024^2)

#Integrate data
combined_endos <- IntegrateLayers(
  object = endo.combined, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca",
  verbose = FALSE
)

#Reset max ram to not accidentally use a bunch later without realizing
options(future.globals.maxSize = 500 * 1024^2)

ElbowPlot(combined_endos, ndims = 40)

#Prepare data for plotting
combined_endos <- FindNeighbors(combined_endos, reduction = "integrated.rpca", dims = 1:25)
combined_endos <- FindClusters(combined_endos, resolution = 2, cluster.name = "rpca_clusters")
combined_endos <- RunUMAP(combined_endos, reduction = "integrated.rpca", dims = 1:25, reduction.name = "umap.rpca")
combined_endos$manualAnnotation = 'endothelial'
combined_endos$dataset <- gsub("_.+", '', colnames(combined_endos))
combined_endos$timepoint_or_treatment <- factor(dplyr::case_when(combined_endos$dataset == 'parse' & combined_endos$Treatment != 'PBS' ~ combined_endos$Timepoint,
                                                                  combined_endos$dataset == 'parse' & combined_endos$Treatment == 'PBS' ~ combined_endos$Treatment,
                                                                  combined_endos$dataset == 'yang' ~ combined_endos$Treatment,
                                                                 combined_endos$dataset == 'sn' ~ combined_endos$infected_grouped),
                                                 levels = c('PBS', 'Day 3', 'Day 4', 'Day 5', 'healthy', 'mild', 'moderate', 'severe',
                                                            'uninfected', 'infected'))
combined_endos$dataset_group <- factor(paste(combined_endos$dataset, combined_endos$timepoint_or_treatment, sep = '_'),
                                        levels = c('parse_PBS', 'parse_Day 3', 'parse_Day 4', 'parse_Day 5',
                                                   'yang_healthy', 'yang_mild','yang_moderate', 'yang_severe',
                                                   'sn_uninfected', 'sn_infected'))
#Plot data together
pdf('/Users/matthewbradley/Documents/ÖverbyLab/yang_public_data/combined_data_plots/combined_endothelial_plot.pdf', width = 6, height = 5)
DimPlot(combined_endos, label = FALSE, group.by = 'dataset_group', reduction = 'umap.rpca',
        cols = newCols)+
  ggtitle('Endothelial integrated UMAP')+
  xlab('')+
  ylab('')
dev.off()

#Endothelial trafficking genes
pdf('/Users/matthewbradley/Documents/ÖverbyLab/yang_public_data/combined_data_plots/endo_trafficking_genes.pdf', width = 6, height = 5)
comparative_dotplot(data = combined_endos, genes = c('Vcam1', 'Icam1', 'Icam2', 'Sele', 'Selp', 'Mcam', 'F11r', 'Jam2', 'Pecam1', 'Pvr', 'Cd99l2',
                                                     'Cdh5', 'Cxcl12'), title = 'Endothelial trafficking genes')
dev.off()


########### Single nuclei comparison #############
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##


#Read in single nuclei data for ccl comparison. Only usjng healthy and severe Yang
#Read in processed data
sn_integrated_dat <- LoadSeuratRds('~/Documents/ÖverbyLab/single_nuclei_proj/LGTVscCombined.rds')

#Check data
newCols <-  c(brewer.pal(12, 'Paired'), '#99FFE6', '#CE99FF', '#18662E','#737272',  '#FF8AEF')
DimPlot(sn_integrated_dat, group.by =   'manualAnnotation', cols = newCols, reduction = 'umap.integrated')

#subset to wt for sn data. healthy and severe for yang data
sn_integrated_dat_wt <- subset(sn_integrated_dat, new_genotype %in% c('wt', 'wt (same)'))
yang_data_no_mild <- subset(yang_data, Treatment %in% c('healthy', 'severe'))

#merge wt sn with yang
sn.combined <- merge(sn_integrated_dat_wt, y = yang_data_no_mild, add.cell.ids = c("singleNuclei", "yang"), project = "sn_merged")

sn.combined <- NormalizeData(sn.combined)
sn.combined <- FindVariableFeatures(sn.combined)
sn.combined <- ScaleData(sn.combined)
sn.combined <- RunPCA(sn.combined)

#Need to allow greater ram usage to run pca integration
options(future.globals.maxSize = 20000 * 1024^2)

#Integrate data
sn_merged <- IntegrateLayers(
  object = sn.combined, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca",
  verbose = FALSE
)

#Reset max ram to not accidentally use a bunchlater without realizing
options(future.globals.maxSize = 500 * 1024^2)

ElbowPlot(sn_merged, ndims = 40)

#Prepare data for plotting
sn_merged <- FindNeighbors(sn_merged, reduction = "integrated.rpca", dims = 1:25)
sn_merged <- FindClusters(sn_merged, resolution = 2, cluster.name = "rpca_clusters")
sn_merged <- RunUMAP(sn_merged, reduction = "integrated.rpca", dims = 1:25, reduction.name = "umap.rpca")

#plot merged data by dataset
sn_merged$dataset <- gsub("_.+", '', colnames(sn_merged))

DimPlot(sn_merged, label = FALSE, split.by = 'dataset', group.by = 'manualAnnotation', reduction = 'umap.rpca',
        cols = newCols)+
  ggtitle('single-nuclei and yang')+
  xlab('')+
  ylab('')

sn_merged$timepoint_or_treatment <- (dplyr::case_when(sn_merged$dataset == 'singleNuclei' ~ sn_merged$new_inf,
                                                            sn_merged$dataset == 'yang' ~ sn_merged$Treatment))
sn_merged[[]][sn_merged$new_inf %in% c('none', 'mock'),]$timepoint_or_treatment = 'healthy'

sn_merged$dataset_group <- factor(paste(sn_merged$dataset, sn_merged$timepoint_or_treatment, sep = '_'))

pdf('/Users/matthewbradley/Documents/ÖverbyLab/yang_public_data/combined_data_plots/single_nuclei_and_yang_umap.pdf', width = 9, height = 6)
DimPlot(sn_merged, label = FALSE, group.by = 'dataset_group', split.by = 'dataset', reduction = 'umap.rpca',
        cols = newCols)+
  ggtitle('single-nuclei and yang')+
  xlab('')+
  ylab('')
dev.off()

#plot relevant ccl genes for macrophage recruitment
ccl_genes <- c('Ccl2', 'Ccl5', 'Ccl7', 'Ccl12')

comparative_dotplot(data = sn_merged, genes = ccl_genes, title = '', sc_timepoints = FALSE)

#Want to look at celltypes from each dataset split
sn_merged$dataset_group_celltype <- paste(sn_merged$dataset_group, sn_merged$manualAnnotation, sep = '_')
ccl_celltype_levels <- DotPlot(sn_merged, features = ccl_genes, 
        group.by = 'dataset_group_celltype', scale = FALSE)$data

ccl_meta <- as.data.frame(str_split_fixed(ccl_celltype_levels$id, "_", 3))
colnames(ccl_meta) <- c('dataset', 'treatment', 'celltype')
ccl_meta$dataset_treatment = paste(ccl_meta$dataset, ccl_meta$treatment, sep = '_')
ccl_celltype_levels <- cbind(ccl_celltype_levels, ccl_meta)

ccl_celltype_levels_sn <- ccl_celltype_levels[ccl_celltype_levels$dataset == 'singleNuclei',]
ccl_celltype_levels_yang <- ccl_celltype_levels[ccl_celltype_levels$dataset == 'yang',]

pdf('/Users/matthewbradley/Documents/ÖverbyLab/yang_public_data/combined_data_plots/combined_chemokine_dot.pdf', width = 7, height = 6)
ggplot(ccl_celltype_levels, aes(x = celltype, y = features.plot))+
  facet_wrap(~dataset_treatment, nrow = 1, ncol = 4)+
  geom_point(aes(size = pct.exp, fill = avg.exp.scaled), pch = 21)+
  coord_flip()+
  scale_size_continuous(range = c(0.5,6), limits = c(0,100))+
  scale_fill_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white"), 
                       values = c(1.0,0.7,0.4,0),
                       limits = c(0,5.5))+
  ggtitle('Single-nuclei chemokines')+
  theme_bw()+
  theme(panel.grid = element_blank(),
        #panel.border = element_rect(colour = "white", fill = NA),
        panel.spacing = unit(0, "line"),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90))+
  scale_size(range = c(0,9), limits = c(0, 100))+
  ylab('')+
  xlab('')
dev.off()

pdf('/Users/matthewbradley/Documents/ÖverbyLab/yang_public_data/combined_data_plots/sn_chemokine_dot.pdf', width = 7, height = 6)
ggplot(ccl_celltype_levels_sn, aes(x = celltype, y = features.plot))+
  facet_wrap(~treatment)+
  geom_point(aes(size = pct.exp, fill = avg.exp.scaled), pch = 21)+
  coord_flip()+
  scale_size_continuous(range = c(0.5,6), limits = c(0,100))+
  scale_fill_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white"), 
                       values = c(1.0,0.7,0.4,0),
                       limits = c(0,5.5))+
  ggtitle('Single-nuclei chemokines')+
  theme_bw()+
  theme(panel.grid = element_blank(),
        #panel.border = element_rect(colour = "white", fill = NA),
        panel.spacing = unit(0, "line"),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90))+
  scale_size(range = c(0,9), limits = c(0, 100))+
  ylab('')+
  xlab('')
dev.off()

pdf('/Users/matthewbradley/Documents/ÖverbyLab/yang_public_data/combined_data_plots/yang_chemokine_dot.pdf', width = 7, height = 6)
ggplot(ccl_celltype_levels_yang, aes(x = celltype, y = features.plot))+
  facet_wrap(~treatment)+
  geom_point(aes(size = pct.exp, fill = avg.exp.scaled), pch = 21)+
  coord_flip()+
  scale_size_continuous(range = c(0.5,6), limits = c(0,100))+
  scale_fill_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white"), 
                       values = c(1.0,0.7,0.4,0),
                       limits = c(0,5.5))+
  ggtitle('Yang chemokines')+
  theme_bw()+
  theme(panel.grid = element_blank(),
        #panel.border = element_rect(colour = "white", fill = NA),
        panel.spacing = unit(0, "line"),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90))+
  scale_size(range = c(0,9), limits = c(0, 100))+
  ylab('')+
  xlab('')
dev.off()
