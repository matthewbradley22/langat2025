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
astrocytes_4 <- subset(ParseSeuratObj_int, manualAnnotation == 'Astrocytes' & Treatment != 'rLGTV' & (Timepoint == 'Day 4' | Treatment == 'PBS'))
astrocytes_5 <- subset(ParseSeuratObj_int, manualAnnotation == 'Astrocytes' & Treatment != 'rLGTV' & (Timepoint == 'Day 5' | Treatment == 'PBS'))

#Split astros by geno
ips_astro <- subset(astrocytes, Genotype == 'IPS1' | Treatment == 'PBS')
wt_astro <- subset(astrocytes, Genotype == 'WT' | Treatment == 'PBS')
ips_astro_4 <- subset(astrocytes_4, Genotype == 'IPS1' | Treatment == 'PBS')
wt_astro_4 <- subset(astrocytes_4, Genotype == 'WT' | Treatment == 'PBS')
ips_astro_5 <- subset(astrocytes_5, Genotype == 'IPS1' | Treatment == 'PBS')
wt_astro_5 <- subset(astrocytes_5, Genotype == 'WT' | Treatment == 'PBS')

#Plot astros
umap_astro_group <- function(astro_dat, num_dims = 20, returnElbow = FALSE, main = NULL, grouping = 'treatment_organ',
                             group_colors = c('#39B6E3', '#3986E3', '#FFAA70', '#C25610'),
                             ystart = -8, xstart = -14, xlength = 6){
  astro_dat <- prepSeuratObj(astro_dat, use_all_genes = FALSE)
  if(returnElbow){
    print(ElbowPlot(astro_dat, ndims = 40))
    return()
  }
  astro_dat <- prepUmapSeuratObj(astro_dat, nDims = num_dims, reductionName = 'astrocytes_umap', resolution_value = 0.8)
  astro_dat$treatment_organ = paste(astro_dat$Treatment, astro_dat$Organ, sep = '_')
  p1 <- DimPlot(astro_dat, reduction = 'astrocytes_umap', group.by = grouping)+
    ggtitle(main)+
    scale_color_manual(values = group_colors)+
    theme_void()+
    annotate("segment", 
             x =xstart, xend = xstart + xlength, 
             y = ystart, yend = ystart, 
             linewidth = 1.3,
             arrow = arrow(type = "closed", length = unit(10, 'pt'))) +
    annotate("segment", 
             x = xstart, xend = xstart, 
             y = ystart, yend = ystart + 4, 
             linewidth = 1.3,
             arrow = arrow(type = "closed", length = unit(10, 'pt'))) +
    theme(text = element_text(size = 24))+
    guides(color=guide_legend(override.aes=list(size=6)))
  print(p1)
}

astrocytes$treatment_genotype <- paste(astrocytes$Treatment, astrocytes$Genotype, sep = '_')
png('~/Documents/ÖverbyLab/single_cell_ISG_figures/astrocytes_fig/astro_day3_umap.png', width = 3100, height = 2500, res = 400)
umap_astro_group(astrocytes, grouping = 'treatment_genotype', group_colors = c('#FFAA70', '#39B6E3', '#C25610', '#3986E3'))
dev.off()

astrocytes_4$treatment_genotype <- paste(astrocytes_4$Treatment, astrocytes_4$Genotype, sep = '_')
png('~/Documents/ÖverbyLab/single_cell_ISG_figures/astrocytes_fig/astro_day4_umap.png', width = 3100, height = 2500, res = 400)
umap_astro_group(astrocytes_4, grouping = 'treatment_genotype', group_colors = c('#FFAA70', '#39B6E3', '#C25610', '#3986E3'),
                 ystart = -9, xstart = -8, xlength = 5)
dev.off()

astrocytes_5$treatment_genotype <- paste(astrocytes_5$Treatment, astrocytes_5$Genotype, sep = '_')
png('~/Documents/ÖverbyLab/single_cell_ISG_figures/astrocytes_fig/astro_day5_umap.png', width = 3100, height = 2500, res = 400)
umap_astro_group(astrocytes_5, grouping = 'treatment_genotype', group_colors = c('#FFAA70', '#39B6E3', '#C25610', '#3986E3'),
                 ystart = -9, xstart = -8, xlength = 5)
dev.off()

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

#Day 3
wt_astro_degs <- FindMarkers(wt_astro, group.by = 'Treatment', ident.1 = 'rChLGTV', test.use = 'MAST')
wt_astro_degs_sig <- wt_astro_degs[wt_astro_degs$p_val_adj<0.01 & wt_astro_degs$avg_log2FC > 1,]

ips_astro_degs <- FindMarkers(ips_astro, group.by = 'Treatment', ident.1 = 'rChLGTV', test.use = 'MAST')
ips_astro_degs_sig <- ips_astro_degs[ips_astro_degs$p_val_adj<0.01 & ips_astro_degs$avg_log2FC > 1,]

#Day 4
wt_astro_degs_4 <- FindMarkers(wt_astro_4, group.by = 'Treatment', ident.1 = 'rChLGTV', test.use = 'MAST')
wt_astro_degs_4_sig <- wt_astro_degs_4[wt_astro_degs_4$p_val_adj<0.01 & wt_astro_degs_4$avg_log2FC > 1,]

ips_astro_degs_4 <- FindMarkers(ips_astro_4, group.by = 'Treatment', ident.1 = 'rChLGTV', test.use = 'MAST')
ips_astro_degs_4_sig <- ips_astro_degs_4[ips_astro_degs_4$p_val_adj<0.01 & ips_astro_degs_4$avg_log2FC > 1,]

#Day 5
wt_astro_degs_5 <- FindMarkers(wt_astro_5, group.by = 'Treatment', ident.1 = 'rChLGTV', test.use = 'MAST')
wt_astro_degs_5_sig <- wt_astro_degs_5[wt_astro_degs_5$p_val_adj<0.01 & wt_astro_degs_5$avg_log2FC > 1,]

ips_astro_degs_5 <- FindMarkers(ips_astro_5, group.by = 'Treatment', ident.1 = 'rChLGTV', test.use = 'MAST')
ips_astro_degs_5_sig <- ips_astro_degs_5[ips_astro_degs_5$p_val_adj<0.01 & ips_astro_degs_5$avg_log2FC > 1,]

#GO pathyway on wt only, ips only, and intersection
wt_paths_all <- gprofiler2::gost(query = rownames(wt_astro_degs_sig), organism = 'mmusculus', evcodes = TRUE, sources = c('GO:BP', 'KEGG'))
wt_paths_all$result$term_name

pdf('~/Documents/ÖverbyLab/single_cell_ISG_figures/astrocytes_fig/day3_wt_pathways.pdf', width = 6, height = 5)
ggplot(head(wt_paths_all$result), aes(x = -log10(p_value), y = reorder(term_name, p_value)))+
  geom_bar(stat = 'identity', fill = "#B3BFE2")+
  theme_classic()+
  ylab('')+
  theme(text = element_text(size = 22))+
  geom_text(aes(label = intersection_size), hjust = 1)+
  xlim(c(0,68))
dev.off()

ips_paths_all <- gprofiler2::gost(query = rownames(ips_astro_degs_sig), organism = 'mmusculus', evcodes = TRUE, sources = c('GO:BP', 'KEGG'))
ips_paths_all$result$term_name

ips_paths_all$result$term_name[ips_paths_all$result$term_name == 'biological process involved in interspecies interaction between organisms'] = 'interspecies interaction'
pdf('~/Documents/ÖverbyLab/single_cell_ISG_figures/astrocytes_fig/day3_ips_pathways.pdf', width = 6, height = 5)
ggplot(head(ips_paths_all$result), aes(x = -log10(p_value), y = reorder(term_name, p_value)))+
  geom_bar(stat = 'identity', fill = "#FDC0AC")+
  theme_classic()+
  ylab('')+
  theme(text = element_text(size = 22))+
  geom_text(aes(label = intersection_size), hjust = 1)+
  xlim(c(0,68))
dev.off()

wt_only <- setdiff(rownames(wt_astro_degs_sig), rownames(ips_astro_degs_sig))
wt_paths <- gprofiler2::gost(query = wt_only, organism = 'mmusculus', evcodes = TRUE, sources = c('GO:BP', 'KEGG'))
wt_paths$result$term_name

ips_only <- setdiff(rownames(ips_astro_degs_sig), rownames(wt_astro_degs_sig))
ips_paths <- gprofiler2::gost(query = ips_only, organism = 'mmusculus', evcodes = TRUE, sources = c('GO:BP', 'KEGG'))
ips_paths$result$term_name

both <- intersect(rownames(wt_astro_degs_sig), rownames(ips_astro_degs_sig))
both_paths <- gprofiler2::gost(query = both, organism = 'mmusculus', evcodes = TRUE, sources = c('GO:BP', 'KEGG'))
both_paths$result$term_name

#Make venn diagrams of degs at each timepoint
#Day 3
day3_wt_ips_eul <- euler(c('wt' = sum(!rownames(wt_astro_degs_sig) %in% rownames(ips_astro_degs_sig)),
                             'ips' = sum(!rownames(ips_astro_degs_sig) %in% rownames(wt_astro_degs_sig)),
                             'wt&ips' = sum(rownames(wt_astro_degs_sig) %in% rownames(ips_astro_degs_sig))))

day4_wt_ips_eul <- euler(c('wt' = sum(!rownames(wt_astro_degs_4_sig) %in% rownames(ips_astro_degs_4_sig)),
                           'ips' = sum(!rownames(ips_astro_degs_4_sig) %in% rownames(wt_astro_degs_4_sig)),
                           'wt&ips' = sum(rownames(wt_astro_degs_4_sig) %in% rownames(ips_astro_degs_4_sig))))

day5_wt_ips_eul <- euler(c('wt' = sum(!rownames(wt_astro_degs_5_sig) %in% rownames(ips_astro_degs_5_sig)),
                           'ips' = sum(!rownames(ips_astro_degs_5_sig) %in% rownames(wt_astro_degs_5_sig)),
                           'wt&ips' = sum(rownames(wt_astro_degs_5_sig) %in% rownames(ips_astro_degs_5_sig))))

pdf('~/Documents/ÖverbyLab/single_cell_ISG_figures/astrocytes_fig/day3_wt_ips_astro_venn.pdf', width = 6, height = 5)
plot(day3_wt_ips_eul, quantities = list(cex = 1.8), fills = c("#B3BFE2", "#FDC0AC"))       
dev.off()

pdf('~/Documents/ÖverbyLab/single_cell_ISG_figures/astrocytes_fig/day4_wt_ips_astro_venn.pdf', width = 6, height = 5)
plot(day4_wt_ips_eul, quantities = list(cex = 1.8), fills = c("#B3BFE2", "#FDC0AC"))       
dev.off()

pdf('~/Documents/ÖverbyLab/single_cell_ISG_figures/astrocytes_fig/day5_wt_ips_astro_venn.pdf', width = 6, height = 5)
plot(day5_wt_ips_eul, quantities = list(cex = 1.8), fills = c("#B3BFE2", "#FDC0AC"))       
dev.off()

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

#dotplot of relevant isgs across all groups
top_isgs <- head(rownames(wt_astro_degs_sig), n = 8)

astrocytes_all <- subset(ParseSeuratObj_int, manualAnnotation == 'Astrocytes' & Treatment != 'rLGTV')
astrocytes_all$treatment_genotype_time <- paste(astrocytes_all$Treatment, astrocytes_all$Genotype, astrocytes_all$Timepoint, sep = '_')

isg_dot_dat <- DotPlot(astrocytes_all, features = top_isgs, scale = FALSE, group.by = 'treatment_genotype_time')$data

plot_dot_dat <- function(dat){
  dat %>% tidyr::separate(id, into = c('treatment', 'genotype', 'timepoint'), sep = '_') %>% 
    dplyr::mutate(geno_time = paste(genotype, timepoint, sep = ' ')) %>% 
    ggplot(aes(x = geno_time, y = features.plot, size = pct.exp, fill = avg.exp.scaled))+
    geom_point(pch = 21)+
    facet_wrap(~treatment, scales = 'free')+
    theme_classic()+
    scale_fill_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white"), 
                         values = c(1.0,0.7,0.4,0))+
    theme(text = element_text(size = 17),
          axis.text.x = element_text(angle = 45, hjust = 1))
}

pdf('~/Documents/ÖverbyLab/single_cell_ISG_figures/astrocytes_fig/top_isg_dotplot.pdf', width = 7, height = 6)
plot_dot_dat(isg_dot_dat)
dev.off()

#Look at cytokines and mhc II genes
mhc_dot_dat <- DotPlot(astrocytes_all, features = c('H2-Aa', 'H2-Ab1', 'H2-Ea', 'H2-Eb1', 'H2-DMb1'), scale = FALSE, group.by = 'treatment_genotype_time')$data
plot_dot_dat(mhc_dot_dat)

cyto_dot_dat <- DotPlot(astrocytes_all, features = c('Tnf', 'Il6', 'Ccl4', 'Ccl2', 'Ccl5',
                                                     'Il12a', 'Cxcl1', 'Il1a', 'Il2', 'Il15', 'Cxcl10'), scale = FALSE, group.by = 'treatment_genotype_time')$data
plot_dot_dat(cyto_dot_dat)

cyto_recept_dot <- DotPlot(astrocytes_all, features = c('Cxcr4', 'Ccr1', 'Ccr5', 'Cx3cr1'), scale = FALSE, group.by = 'treatment_genotype_time')$data
plot_dot_dat(cyto_recept_dot)

#Look at brain region associated markers from https://www.nature.com/articles/s41593-021-00905-6
region_speific_markers <- DotPlot(astrocytes_all, features = c('Mfge8', 'Igfbp2', 'Nnat', 'Dbi',
                                                               'Agt', 'Nupr1', 'Thrsp', 'Vim', 'Gfap',
                                                               'Slc43a3', 'Prss56', 'Crym'), scale = FALSE, group.by = 'treatment_genotype_time')$data
plot_dot_dat(region_speific_markers)

#Common lps induced gene
serp <- DotPlot(astrocytes_all, features = c('Serpina3n'), scale = FALSE, group.by = 'treatment_genotype_time')$data
plot_dot_dat(serp)

#Unique inflammation genes from same paper as above
inf_genes <- c('Igtp', 'Tap1', 'Stat1', 'Agt', 'Cd34', 'Fbln5',
               'Pgf', 'Cd109', 'Igfbp7', 'Timp1', 'Gpx1', 'Hspb1', 'Gap43')
inf_dot <- DotPlot(astrocytes_all, features = inf_genes, scale = FALSE, group.by = 'treatment_genotype_time')$data
plot_dot_dat(inf_dot)

#Proliferation genes, just general genes not from specific paper - very low percent expressed
prolif_genes <- c('Mki67', 'Top2a')
prolif_dot <- DotPlot(astrocytes_all, features = prolif_genes, scale = FALSE, group.by = 'treatment_genotype_time')$data
plot_dot_dat(prolif_dot)

#Look at clusters within infected astrocytes to explore heterogeneity a bit
infected_astrocytes <- subset(ParseSeuratObj_int, manualAnnotation == 'Astrocytes' & Treatment == 'rChLGTV')

infected_astrocytes <- prepSeuratObj(infected_astrocytes)

ElbowPlot(infected_astrocytes, ndims = 40)

infected_astrocytes <- prepUmapSeuratObj(infected_astrocytes, nDims = 15, reductionName = 'inf_astro_umap', resolution_value = 0.8)

DimPlot(infected_astrocytes, reduction = 'inf_astro_umap', group.by = 'Organ') +
  ggtitle('Infected astrocytes')+
  xlab('')+
  ylab('')


