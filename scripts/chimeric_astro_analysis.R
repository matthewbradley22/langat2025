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

#Split astros by day and to just chimeric
astrocytes <- subset(ParseSeuratObj_int, manualAnnotation == 'Astrocytes' & Treatment != 'rLGTV' & (Timepoint == 'Day 3' | Treatment == 'PBS'))
astrocytes_5 <- subset(ParseSeuratObj_int, manualAnnotation == 'Astrocytes' & Treatment != 'rLGTV' & (Timepoint == 'Day 5' | Treatment == 'PBS'))

#Split astros by geno
ips_astro <- subset(astrocytes, Genotype == 'IPS1')
wt_astro <- subset(astrocytes, Genotype == 'WT')
ips_astro_5 <- subset(astrocytes_5, Genotype == 'IPS1')
wt_astro_5 <- subset(astrocytes_5, Genotype == 'WT')

#Plot astros
umap_astro_group <- function(astro_dat, num_dims = 20, returnElbow = FALSE, main = NULL){
  astro_dat <- prepSeuratObj(astro_dat, use_all_genes = FALSE)
  if(returnElbow){
    print(ElbowPlot(astro_dat, ndims = 40))
    return()
  }
  astro_dat <- prepUmapSeuratObj(astro_dat, nDims = num_dims, reductionName = 'astrocytes_umap', resolution_value = 0.8)
  astro_dat$treatment_organ = paste(astro_dat$Treatment, astro_dat$Organ, sep = '_')
  p1 <- DimPlot(astro_dat, reduction = 'astrocytes_umap', group.by = 'treatment_organ')+
    ggtitle(main)+
    scale_color_manual(values = c('#39B6E3', '#3986E3', '#E35539', '#E36F39'))+
    ylab('Umap2')+
    xlab('Umap1')+
    theme(text = element_text(size = 24))
  print(p1)
}

pdf('~/Documents/ÖverbyLab/single_cell_ISG_figures/day3_fig/day3_plots/wt_astro_umap.pdf', width = 6, height = 6)
umap_astro_group(wt_astro)
dev.off()

pdf('~/Documents/ÖverbyLab/single_cell_ISG_figures/day3_fig/day3_plots/ips_astro_umap.pdf', width = 6, height = 6)
umap_astro_group(ips_astro)
dev.off()

png('~/Documents/ÖverbyLab/single_cell_ISG_figures/day4wt_day5ips/day5_astro_wt_umap.png', width = 600, height = 600)
umap_astro_group(wt_astro_5, main = 'WT astros')
dev.off()

png('~/Documents/ÖverbyLab/single_cell_ISG_figures/day4wt_day5ips/day5_astro_ips_umap.png', width = 600, height = 600)
umap_astro_group(ips_astro_5, main = 'IPS astros')
dev.off()

#Find DEGs
wt_astro_degs <- FindMarkers(wt_astro, group.by = 'Treatment', ident.1 = 'rChLGTV', test.use = 'MAST')
wt_astro_degs_sig <- wt_astro_degs[wt_astro_degs$p_val_adj<0.01 & wt_astro_degs$avg_log2FC > 1,]

ips_astro_degs <- FindMarkers(ips_astro, group.by = 'Treatment', ident.1 = 'rChLGTV', test.use = 'MAST')
ips_astro_degs_sig <- ips_astro_degs[ips_astro_degs$p_val_adj<0.01 & ips_astro_degs$avg_log2FC > 1,]

#GO pathyway on wt only, ips only, and intersection
wt_paths_all <- gprofiler2::gost(query = rownames(wt_astro_degs_sig), organism = 'mmusculus', evcodes = TRUE, sources = c('GO:BP', 'KEGG'))
wt_paths_all$result$term_name

ips_paths_all <- gprofiler2::gost(query = rownames(ips_astro_degs_sig), organism = 'mmusculus', evcodes = TRUE, sources = c('GO:BP', 'KEGG'))
ips_paths_all$result$term_name

wt_only <- setdiff(rownames(wt_astro_degs_sig), rownames(ips_astro_degs_sig))
wt_paths <- gprofiler2::gost(query = wt_only, organism = 'mmusculus', evcodes = TRUE, sources = c('GO:BP', 'KEGG'))
wt_paths$result$term_name

ips_only <- setdiff(rownames(ips_astro_degs_sig), rownames(wt_astro_degs_sig))
ips_paths <- gprofiler2::gost(query = ips_only, organism = 'mmusculus', evcodes = TRUE, sources = c('GO:BP', 'KEGG'))
ips_paths$result$term_name

both <- intersect(rownames(wt_astro_degs_sig), rownames(ips_astro_degs_sig))
both_paths <- gprofiler2::gost(query = both, organism = 'mmusculus', evcodes = TRUE, sources = c('GO:BP', 'KEGG'))
both_paths$result$term_name

#Do not produce log file for venn diagram
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

astro_day3 <- VennDiagram::venn.diagram(list(wt = rownames(wt_astro_degs_sig), ips = rownames(ips_astro_degs_sig)), filename = NULL,
                                     fill = brewer.pal(3, "Pastel2")[1:2], cex = 1.5, cat.cex = 1.5)
grid::grid.draw(astro_day3)
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
