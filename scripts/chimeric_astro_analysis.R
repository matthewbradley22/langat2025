library(Seurat)
library(RColorBrewer)
library(dplyr)
library(UpSetR)
library(ggplot2)
library(stringr)
library(tidyr)
library(gprofiler2)
source('~/Documents/ÖverbyLab//scripts/langatFunctions.R')

#Load data
ParseSeuratObj_int <- LoadSeuratRds("~/Documents/ÖverbyLab/data/FilteredRpcaIntegratedDatNoDoublets.rds") 

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

#Chimeric day 3 astrocytes
astrocytes <- subset(ParseSeuratObj_int, manualAnnotation == 'Astrocytes' & Treatment != 'rLGTV' & (Timepoint == 'Day 3' | Treatment == 'PBS'))
ips_astro <- subset(astrocytes, Genotype == 'IPS1')
wt_astro <- subset(astrocytes, Genotype == 'WT')

#Plot astros
wt_astro <- prepSeuratObj(wt_astro, use_all_genes = FALSE)
ElbowPlot(wt_astro, ndims = 40)
wt_astro <- prepUmapSeuratObj(wt_astro, nDims = 20, reductionName = 'astrocytes_umap', resolution_value = 0.8)
wt_astro$treatment_organ = paste(wt_astro$Treatment, wt_astro$Organ, sep = '_')
DimPlot(wt_astro, reduction = 'astrocytes_umap', group.by = 'treatment_organ')+
  ggtitle('WT Day 3 Astrocytes')+
  scale_color_manual(values = c('#39B6E3', '#3986E3', '#E35539', '#E36F39'))+
  ylab('Umap2')+
  xlab('Umap1')

ips_astro <- prepSeuratObj(ips_astro, use_all_genes = FALSE)
ElbowPlot(ips_astro, ndims = 40)
ips_astro <- prepUmapSeuratObj(ips_astro, nDims = 20, reductionName = 'astrocytes_umap_ips', resolution_value = 0.8)
ips_astro$treatment_organ = paste(ips_astro$Treatment, ips_astro$Organ, sep = '_')
DimPlot(ips_astro, reduction = 'astrocytes_umap_ips', group.by = 'treatment_organ')+
  scale_color_manual(values = c('#39B6E3', '#3986E3', '#E35539', '#E36F39'))+
  ggtitle('IPS1 Day 3 Astrocytes')+
  ylab('Umap2')+
  xlab('Umap1')

#Find DEGs
wt_astro_degs <- FindMarkers(wt_astro, group.by = 'Treatment', ident.1 = 'rChLGTV', test.use = 'MAST')
wt_astro_degs_sig <- wt_astro_degs[wt_astro_degs$p_val_adj<0.01 & wt_astro_degs$avg_log2FC > 1,]

ips_astro_degs <- FindMarkers(ips_astro, group.by = 'Treatment', ident.1 = 'rChLGTV', test.use = 'MAST')
ips_astro_degs_sig <- ips_astro_degs[ips_astro_degs$p_val_adj<0.01 & ips_astro_degs$avg_log2FC > 1,]

#DO pathyway on wt only, ips only, and intersection
wt_only <- setdiff(rownames(wt_astro_degs_sig), rownames(ips_astro_degs_sig))
wt_paths <- gprofiler2::gost(query = wt_only, organism = 'mmusculus', evcodes = TRUE, sources = c('GO:BP', 'KEGG'))
wt_paths$result$term_name

ips_only <- setdiff(rownames(ips_astro_degs_sig), rownames(wt_astro_degs_sig))
ips_paths <- gprofiler2::gost(query = ips_only, organism = 'mmusculus', evcodes = TRUE, sources = c('GO:BP', 'KEGG'))
ips_paths$result$term_name

both <- intersect(rownames(wt_astro_degs_sig), rownames(ips_astro_degs_sig))
both_paths <- gprofiler2::gost(query = both, organism = 'mmusculus', evcodes = TRUE, sources = c('GO:BP', 'KEGG'))
both_paths$result$term_name

#Look at all cell types to see if trend is the same
ips <- subset(ParseSeuratObj_int, Treatment != 'rLGTV' & (Timepoint == 'Day 3' | Treatment == 'PBS') & Genotype == 'IPS1')
wt <- subset(ParseSeuratObj_int, Treatment != 'rLGTV' & (Timepoint == 'Day 3' | Treatment == 'PBS') & Genotype == 'WT')

#Plot wt
wt <- prepSeuratObj(wt, use_all_genes = FALSE)
ElbowPlot(wt, ndims = 40)
wt <- prepUmapSeuratObj(wt, nDims = 20, reductionName = 'wt_umap', resolution_value = 0.8)
wt$treatment_organ = paste(wt$Treatment, wt$Organ, sep = '_')
wt_umap <- DimPlot(wt, reduction = 'wt_umap', group.by = 'treatment_organ')+
  ggtitle('Day 3 wt')+
  scale_color_manual(values = c('#39B6E3', '#3986E3', '#E35539', '#E36F39'))+
  ylab('Umap2')+
  xlab('Umap1')

#Plot ips
ips <- prepSeuratObj(ips, use_all_genes = FALSE)
ElbowPlot(ips, ndims = 40)
ips <- prepUmapSeuratObj(ips, nDims = 20, reductionName = 'ips_umap', resolution_value = 0.8)
ips$treatment_organ = paste(ips$Treatment, ips$Organ, sep = '_')
ips_umap <- DimPlot(ips, reduction = 'ips_umap', group.by = 'treatment_organ')+
  ggtitle('Day 3 IPS')+
  scale_color_manual(values = c('#39B6E3', '#3986E3', '#E35539', '#E36F39'))+
  ylab('Umap2')+
  xlab('Umap1')

grid.arrange(wt_umap, ips_umap, ncol=2)
