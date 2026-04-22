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

ggplot(alternative_path_celltype_levels$Astrocytes, aes(x = time, y = features.plot, size = pct.exp, color = avg.exp.scaled))+
  facet_wrap(~geno_treatment, scales = 'free_x')+
  geom_point()+
  scale_color_gradientn(colours = c("#F03C0C","#F0A451","white"), 
                        values = c(1.0,0.6,0))+
  theme_classic()

#Which genes are differentially expressed 
day4_infected <- subset(chimeric_mock, Treatment == 'rChLGTV' & Timepoint == 'Day 4')
day5_infected <- subset(chimeric_mock, Treatment == 'rChLGTV' & Timepoint == 'Day 5')

#Findmarkers
day4_infected_markers <- FindMarkers(day4_infected, test.use = 'MAST', group.by = 'Genotype', ident.1 = 'WT')
dplyr::filter(day4_infected_markers[alternative_path_genes,], p_val_adj < 0.01 & abs(avg_log2FC) > 1)

day5_infected_markers <- FindMarkers(day5_infected, test.use = 'MAST', group.by = 'Genotype', ident.1 = 'WT')
dplyr::filter(day5_infected_markers[alternative_path_genes,], p_val_adj < 0.01 & abs(avg_log2FC) > 1)

