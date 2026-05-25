library(Seurat)
library(RColorBrewer)
library(dplyr)
library(UpSetR)
library(stringr)
library(ggplot2)
library(tidyr)
library(VennDiagram)
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

#Compare infected ips vs infected wt over all celltypes
chimeric_mock <- subset(ParseSeuratObj_int, Treatment != 'rLGTV')
chimeric_mock$time_geno_treatment <- paste(chimeric_mock$Timepoint, chimeric_mock$Genotype, chimeric_mock$Treatment, sep = '_')

#Alternative pathway genes
alternative_path_genes <- c('Tlr3', 'Tlr7', 'Tlr8', 'Tlr9', 'Ticam1', 'Traf3', 'Myd88', 'Traf6', 'Tram1')
alternative_path_celltype_levels <- list()
for(i in 1:length(unique(chimeric_mock$manualAnnotation))){
  cur_celltype <- unique(chimeric_mock$manualAnnotation)[i]
  cur_dat <- subset(chimeric_mock, manualAnnotation == cur_celltype)
  alternative_dot_dat <- DotPlot(cur_dat, features = alternative_path_genes, group.by = 'time_geno_treatment', scale = FALSE)$data 
  alternative_dot_dat <- alternative_dot_dat %>% tidyr::separate(id, into = c('time', 'geno', 'treatment'), sep = '_')
  alternative_dot_dat$geno <- factor(alternative_dot_dat$geno, levels = c('WT', 'IPS1'))
  alternative_dot_dat$geno_treatment <- factor(paste(alternative_dot_dat$geno, alternative_dot_dat$treatment, sep = '_'),
                                               levels = c('WT_PBS', 'WT_rChLGTV', 'IPS1_PBS', 'IPS1_rChLGTV'))
  alternative_path_celltype_levels[[length(alternative_path_celltype_levels) + 1]] = alternative_dot_dat
  names(alternative_path_celltype_levels)[length(alternative_path_celltype_levels)] = cur_celltype
}

pdf('~/Documents/ÖverbyLab/single_cell_ISG_figures/astrocytes_fig/alternative_path_astro_genes.pdf', width = 9, height = 7)
ggplot(alternative_path_celltype_levels$Astrocytes, aes(x = time, y = features.plot, size = pct.exp, color = avg.exp.scaled))+
  facet_wrap(~geno_treatment, scales = 'free_x')+
  geom_point()+
  scale_color_gradientn(colours = c("#F03C0C","#F0A451","white"), 
                        values = c(1.0,0.6,0))+
  theme_classic()+
  theme(text = element_text(size = 24))+
  scale_size_continuous(range = c(1,9))
dev.off()

#Which genes are differentially expressed 
day4_infected <- subset(chimeric_mock, Treatment == 'rChLGTV' & Timepoint == 'Day 4')
day5_infected <- subset(chimeric_mock, Treatment == 'rChLGTV' & Timepoint == 'Day 5')

#Findmarkers
day4_infected_markers <- FindMarkers(day4_infected, test.use = 'MAST', group.by = 'Genotype', ident.1 = 'WT')
dplyr::filter(day4_infected_markers[alternative_path_genes,], p_val_adj < 0.01 & abs(avg_log2FC) > 1)

day5_infected_markers <- FindMarkers(day5_infected, test.use = 'MAST', group.by = 'Genotype', ident.1 = 'WT')
dplyr::filter(day5_infected_markers[alternative_path_genes,], p_val_adj < 0.01 & abs(avg_log2FC) > 1)

#Are Myd88 or Ticam1 upregulated in infection vs mock?
day5_ips_pbs_astros <- subset(chimeric_mock, ((Timepoint == 'Day 5' & Genotype == 'IPS1') | Treatment == 'PBS') & manualAnnotation == 'Astrocytes')
table(day5_ips_pbs_astros$Timepoint, day5_wt_pbs_astros$Treatment, day5_wt_pbs_astros$Genotype)
day5_wt_pbs_astros <- subset(chimeric_mock, ((Timepoint == 'Day 5' & Genotype == 'WT') | Treatment == 'PBS') & manualAnnotation == 'Astrocytes')
table(day5_wt_pbs_astros$Timepoint, day5_wt_pbs_astros$Treatment, day5_wt_pbs_astros$Genotype)

day5_ips_astro_markers <- FindMarkers(day5_ips_pbs_astros, test.use = 'MAST', group.by = 'Treatment', ident.1 = 'rChLGTV')
day5_wt_astro_markers <- FindMarkers(day5_wt_pbs_astros, test.use = 'MAST', group.by = 'Treatment', ident.1 = 'rChLGTV')

#Myd88 sig up in ips day 5 astros
day5_ips_astro_markers[c('Myd88', 'Ticam1'),]
day5_wt_astro_markers[c('Myd88', 'Ticam1'),]

#Stat genes by celltype
celltypes <- unique(chimeric_mock$manualAnnotation)[unique(chimeric_mock$manualAnnotation) != 'unknown']
stat1_by_celltype_list <- list()
stat2_by_celltype_list <- list()

for(i in 1:length(celltypes)){
  cur_celltype <- subset(chimeric_mock, manualAnnotation == celltypes[[i]])
  cur_celltype_stat1 <- DotPlot(cur_celltype, features = 'Stat1', group.by = 'time_geno_treatment', scale = FALSE)$data %>% 
    tidyr::separate(col = 'id', into = c('time', 'genotype', 'treatment'), sep = '_')
  cur_celltype_stat1$genotype <- factor(cur_celltype_stat1$genotype, levels = c('WT', 'IPS1'))
  p1 <- ggplot(cur_celltype_stat1, aes(x = time, y = treatment, fill = avg.exp.scaled, size = pct.exp))+
    geom_point(pch = 21)+
    facet_wrap(~genotype)+
    scale_fill_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white"), 
                         values = c(1.0,0.7,0.4,0))+
    ggtitle(celltypes[[i]])+
    theme_classic()+
    theme(text = element_text(size = 18))
  stat1_by_celltype_list[[length(stat1_by_celltype_list) + 1]] <- p1
  
  cur_celltype_stat2 <- DotPlot(cur_celltype, features = 'Stat2', group.by = 'time_geno_treatment', scale = FALSE)$data %>% 
    tidyr::separate(col = 'id', into = c('time', 'genotype', 'treatment'), sep = '_')
  cur_celltype_stat2$genotype <- factor(cur_celltype_stat2$genotype, levels = c('WT', 'IPS1'))
  p2 <- ggplot(cur_celltype_stat2, aes(x = time, y = treatment, fill = avg.exp.scaled, size = pct.exp))+
    geom_point(pch = 21)+
    facet_wrap(~genotype)+
    scale_fill_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white"), 
                         values = c(1.0,0.7,0.4,0))+
    ggtitle(celltypes[[i]])+
    theme_classic()+
    theme(text = element_text(size = 18))
  stat2_by_celltype_list[[length(stat2_by_celltype_list) + 1]] <- p2
  
}

names(stat1_by_celltype_list) <- celltypes
names(stat2_by_celltype_list) <- celltypes

pdf('~/Documents/ÖverbyLab/single_cell_ISG_figures/astrocytes_fig/astro_stat1_exp.pdf', width = 8, height = 5)
print(stat1_by_celltype_list$Astrocytes + ggtitle('Stat1'))
dev.off()

pdf('~/Documents/ÖverbyLab/single_cell_ISG_figures/astrocytes_fig/astro_stat2_exp.pdf', width = 8, height = 5)
print(stat2_by_celltype_list$Astrocytes + ggtitle('Stat2'))
dev.off()

#Plot stats in combined celltypes
chimeric_mock_infected <- subset(chimeric_mock, Treatment == 'rChLGTV')
chimeric_mock_infected$time_geno_celltype <- paste(chimeric_mock_infected$Timepoint, chimeric_mock_infected$Genotype, chimeric_mock_infected$manualAnnotation,
                                                   sep = '_')
#Only keep cells with enough sampled
celltypes_with_enough <- chimeric_mock_infected[[]] %>% dplyr::group_by(Timepoint, Genotype, manualAnnotation) %>% 
  dplyr::summarise(cell_count = n()) %>% 
  dplyr::rename(time = Timepoint, genotype = Genotype, celltype = manualAnnotation)

stat1_dat <- DotPlot(chimeric_mock_infected, features = 'Stat1', group.by = 'time_geno_celltype', scale = FALSE)$data %>% 
  tidyr::separate(col = 'id', into = c('time', 'genotype', 'celltype'), sep = '_') %>% 
  dplyr::filter(celltype != 'unknown')

stat1_dat_filtered <- stat1_dat %>% left_join(celltypes_with_enough, by = c('time' = 'time', 'genotype' = 'genotype', 'celltype' = 'celltype')) %>% 
  dplyr::filter(cell_count > 50)
stat1_dat_filtered$genotype = factor(stat1_dat_filtered$genotype, levels = c('WT', 'IPS1'))
stat1_dat_filtered$celltype <- factor(stat1_dat_filtered$celltype, levels = rev(c('Astrocytes', 'Choroid Plexus', 'Endothelial', 'Ependymal', 'Immature Neurons', 
                                                                'Microglia', 'Muscle cells', 'Neurons', 'Oligodendrocytes' ,'Pericytes',
                                                                'B Cells', 'Granulocytes', 'Macrophage/Monocytes', 'Nk cells', 'T cells')))

pdf('~/Documents/ÖverbyLab/single_cell_ISG_figures/timepont_fig/stat1_by_celltype.pdf', width = 8, height = 5)
ggplot(stat1_dat_filtered, aes(x = time, y = celltype, fill = avg.exp.scaled, size = pct.exp))+
  geom_point(pch = 21)+
  facet_wrap(~genotype)+
  scale_fill_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white"), 
                       values = c(1.0,0.7,0.5,0),
                       limits = c(0,3.2))+
  theme_classic()+
  theme(text = element_text(size = 18))+
  scale_size(limits = c(20, 100), range = c(1, 8))+
  ggtitle('Stat1')
dev.off()

stat2_dat <- DotPlot(chimeric_mock_infected, features = 'Stat2', group.by = 'time_geno_celltype', scale = FALSE)$data %>% 
  tidyr::separate(col = 'id', into = c('time', 'genotype', 'celltype'), sep = '_') %>% 
  dplyr::filter(celltype != 'unknown')

stat2_dat_filtered <- stat2_dat %>% left_join(celltypes_with_enough, by = c('time' = 'time', 'genotype' = 'genotype', 'celltype' = 'celltype')) %>% 
  dplyr::filter(cell_count > 50)
stat2_dat_filtered$genotype = factor(stat2_dat_filtered$genotype, levels = c('WT', 'IPS1'))
stat2_dat_filtered$celltype <- factor(stat2_dat_filtered$celltype, levels = rev(c('Astrocytes', 'Choroid Plexus', 'Endothelial', 'Ependymal', 'Immature Neurons', 
                                                                                  'Microglia', 'Muscle cells', 'Neurons', 'Oligodendrocytes' ,'Pericytes',
                                                                                  'B Cells', 'Granulocytes', 'Macrophage/Monocytes', 'Nk cells', 'T cells')))
pdf('~/Documents/ÖverbyLab/single_cell_ISG_figures/timepont_fig/stat2_by_celltype.pdf', width = 8, height = 5)
ggplot(stat2_dat_filtered, aes(x = time, y = celltype, fill = avg.exp.scaled, size = pct.exp))+
  geom_point(pch = 21)+
  facet_wrap(~genotype)+
  scale_fill_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white"), 
                       values = c(1.0,0.7,0.4,0),
                       limits = c(0,3.2))+
  theme_classic()+
  theme(text = element_text(size = 18))+
  scale_size(limits = c(20, 100), range = c(1, 8))+
  ggtitle('Stat2')
dev.off()

pdf('~/Documents/ÖverbyLab/single_cell_ISG_figures/timepont_fig/stats_over_time.pdf', width = 8, height = 5)
ggarrange(stat1_levels, stat2_levels, nrow = 2, common.legend = TRUE, legend = 'right')
dev.off()
#Stat1 and 2, day 3 between genotypes
astrocytes_day3 <- subset(chimeric_mock, manualAnnotation == 'Astrocytes' & Timepoint == 'Day 3' & Treatment == 'rChLGTV')
wt_ips_3_markers <- FindMarkers(astrocytes_day3, group.by = 'Genotype', ident.1 = 'WT', test.use = 'MAST')
wt_ips_3_markers[c('Stat1', 'Stat2'),]

pdf('~/Documents/ÖverbyLab/single_cell_ISG_figures/astrocytes_fig/day_3_inf_astro_stat_vln.pdf', width = 8, height = 5)
VlnPlot(astrocytes_day3, features = c('Stat1', 'Stat2'), group.by = 'Genotype', pt.size = 0)
dev.off()

#Split astros by day and to just chimeric
astrocytes_3 <- subset(ParseSeuratObj_int, manualAnnotation == 'Astrocytes' & Treatment != 'rLGTV' & (Timepoint == 'Day 3' | Treatment == 'PBS'))

#Split astros by geno
ips_astro_3 <- subset(astrocytes, Genotype == 'IPS1' | Treatment == 'PBS')
wt_astro_3 <- subset(astrocytes, Genotype == 'WT' | Treatment == 'PBS')

table(ips_astro_3$Treatment, ips_astro_3$Genotype, ips_astro_3$Timepoint, ips_astro_3$manualAnnotation)
table(wt_astro_3$Treatment, wt_astro_3$Genotype, wt_astro_3$Timepoint, wt_astro_3$manualAnnotation)

#Get upregulated genes 
day3_ips_astro_markers <- FindMarkers(ips_astro_3, group.by = 'Treatment', ident.1 = 'rChLGTV', test.use = 'MAST')
day3_wt_astro_markers <- FindMarkers(wt_astro_3, group.by = 'Treatment', ident.1 = 'rChLGTV', test.use = 'MAST')

day3_ips_astro_markers_sig <- day3_ips_astro_markers %>% as.data.frame() %>% 
  dplyr::filter(p_val_adj < 0.01 & avg_log2FC > 1)
day3_wt_astro_markers_sig <- day3_wt_astro_markers %>% as.data.frame() %>% 
  dplyr::filter(p_val_adj < 0.01 & avg_log2FC > 1)

ips_3_paths <- gprofiler2::gost(query = rownames(day3_ips_astro_markers_sig), organism = 'mmusculus', evcodes = TRUE, sources = c('GO:BP', 'KEGG'))
wt_3_paths <- gprofiler2::gost(query = rownames(day3_wt_astro_markers_sig), organism = 'mmusculus', evcodes = TRUE, sources = c('GO:BP', 'KEGG'))

ips_3_paths$result$term_name %>% head()
wt_3_paths$result$term_name %>% head()

ips_3_paths$result$term_name[ips_3_paths$result$term_name == 'biological process involved in interspecies interaction between organisms'] = 'interaction between organisms'

pdf('~/Documents/ÖverbyLab/single_cell_ISG_figures/astrocytes_fig/ips_day3_paths.pdf', width = 6, height = 5)
ggplot(head(ips_3_paths$result), aes(x = term_name, y = -log10(p_value)))+
  geom_bar(stat = 'identity', fill = '#FDC0AC')+
  coord_flip()+
  ylim(c(0, 67))+
  theme_classic()+
  theme(text = element_text(size = 18))+
  ylab('')
dev.off()

pdf('~/Documents/ÖverbyLab/single_cell_ISG_figures/astrocytes_fig/wt_day3_paths.pdf', width = 6, height = 5)
ggplot(head(wt_3_paths$result), aes(x = term_name, y = -log10(p_value)))+
  geom_bar(stat = 'identity', fill = '#B3BFE2')+
  coord_flip()+
  ylim(c(0, 67))+
  theme_classic()+
  theme(text = element_text(size = 18))+
  ylab('')
dev.off()
