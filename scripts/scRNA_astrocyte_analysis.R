library(Seurat)
library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(tidyr)
source('~/Documents/ÖverbyLab//scripts/langatFunctions.R')

#Load data
ParseSeuratObj_int <- LoadSeuratRds("~/Documents/ÖverbyLab/data/FilteredRpcaIntegratedDatNoDoublets.rds") 
ParseSeuratObj_int$hasVirus = ifelse(ParseSeuratObj_int$virusCountPAdj >= 10, 1, 0)
ParseSeuratObj_int$manualAnnotation[ParseSeuratObj_int$manualAnnotation == 'Macrophage/Monocytes'] = 'Macro/Mono'


#Check data
newCols <-  c(brewer.pal(12, 'Paired'), '#99FFE6', '#CE99FF', '#18662E','#737272',  '#FF8AEF')
newCols[11] =  '#FF8AEF'

DimPlot(ParseSeuratObj_int, label = FALSE, group.by = 'manualAnnotation', reduction = 'umap.integrated',
        cols = newCols)+
  theme(axis.ticks = element_blank(),
        axis.text=element_blank(),
        legend.text=element_text(size=17))+
  xlab('Umap1')+
  ylab('Umap2')+
  guides(color=guide_legend(override.aes=list(size=8)))+
  ggtitle('')

DimPlot(ParseSeuratObj_int, label = FALSE, group.by = 'Treatment', reduction = 'umap.integrated',
        cols = newCols)

#Start looking at overall astrocyte degs, then split by group
astrocytes <- subset(ParseSeuratObj_int, manualAnnotation == 'Astrocytes')

#Recluster astrocytes
astrocytes <- prepSeuratObj(astrocytes)
ElbowPlot(astrocytes, ndims = 40)
astrocytes <- prepUmapSeuratObj(astrocytes, nDims = 20, reductionName = 'astrocytes_umap')

#Add extra column for umap plot
astrocytes$geno_timepoint_treatment = paste(astrocytes$Genotype, astrocytes$Timepoint, astrocytes$Treatment)
DimPlot(astrocytes, reduction = 'astrocytes_umap')
astrocytes$infection_group <- ifelse(astrocytes$Treatment %in% c('rChLGTV', 'rLGTV'), 'infected', 'uninfected')
DimPlot(astrocytes, reduction = 'astrocytes_umap', group.by = 'infection_group')
DimPlot(astrocytes, reduction = 'astrocytes_umap', group.by = 'Genotype')
DimPlot(astrocytes, reduction = 'astrocytes_umap', group.by = 'Timepoint')
DimPlot(astrocytes, reduction = 'astrocytes_umap', group.by = 'geno_timepoint_treatment', cols = newCols)
table(astrocytes$infection_group, astrocytes$Genotype)

#Look at markers

astro_infected_markers <- FindMarkers(astrocytes, ident.1 = 'infected', group.by = 'infection_group',
                             test.use = 'MAST', only.pos = TRUE)
astro_sig_markers <- astro_infected_markers %>% dplyr::filter(avg_log2FC > 1 & p_val_adj < 0.01)
astro_marker_paths <- gprofiler2::gost(query = rownames(astro_sig_markers), organism = 'mmusculus', evcodes = TRUE)
astro_marker_paths$result

#Look into path 558 - seems to be evidence of signalling
astro_marker_paths$result[astro_marker_paths$result$source == 'KEGG',]

#Plot markers - overall, a lot of clear ISGs
#Irf1/8 paper https://pmc.ncbi.nlm.nih.gov/articles/PMC4821649/
FeaturePlot(astrocytes, features = 'Irf1', reduction = 'astrocytes_umap')
DotPlot(astrocytes, features = 'Irf1', scale = FALSE, group.by = 'Treatment')

#Rig 1 type 1 interferon response https://www.sciencedirect.com/science/article/pii/S2211124715004647
FeaturePlot(astrocytes, features = 'Ddx60', reduction = 'astrocytes_umap')
DotPlot(astrocytes, features = 'Ddx60', scale = FALSE, group.by = 'Treatment')

#MHC class 1 https://www.nature.com/articles/nri3339
FeaturePlot(astrocytes, features = 'Nlrc5', reduction = 'astrocytes_umap')
DotPlot(astrocytes, features = 'Nlrc5', scale = FALSE, group.by = 'Treatment')

DotPlot(astrocytes, features = c('Stat1', 'Stat2'), scale = FALSE, group.by = 'Treatment')

FeaturePlot(astrocytes, features = 'Ccl2', reduction = 'astrocytes_umap')

#'Cxcl10' should also be plotted but drowns out other features
DotPlot(astrocytes, features = c('Ccl2', 'Ccl7', 'Ccl12', 'Ccl4',
                                 'Ccl3', 'Ccl5', 'Cxcl9',
                                 'Cxcr3', 'Ccr1', scale = FALSE), 
        scale = FALSE, group.by = 'geno_timepoint_treatment')

table(astrocytes$geno_timepoint_treatment)
FeaturePlot(astrocytes, features = 'Ccl3', reduction = 'astrocytes_umap')
