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
umap_color_list <- c( "#7047A1", "#B370AE","#292270",  "#166DF0","#6D92F8",  "#6DC3F8", "#8a0000","#F76363", "#FF96A2", "#D6644B", 
                      "#F08C3A", "#fdc087","#074F00", "#208d1f","#7bcd79", 
                      "gray")

DimPlot(ParseSeuratObj_int, label = FALSE, group.by = 'manualAnnotation', reduction = 'umap.integrated',
        cols = umap_color_list)+
  theme(axis.ticks = element_blank(),
        axis.text=element_blank(),
        legend.text=element_text(size=17))+
  xlab('Umap1')+
  ylab('Umap2')+
  guides(color=guide_legend(override.aes=list(size=8)))+
  ggtitle('')

chimeric_mock <- subset(ParseSeuratObj_int, Treatment != 'rLGTV' & manualAnnotation != 'unknown')

chimeric_mock$manualAnnotation <- factor(chimeric_mock$manualAnnotation, 
                                         levels = rev(c( 'unknown',  'T cells',   'Nk cells', 'Macrophage/Monocytes', 
                                                         'Granulocytes', 'B Cells',  'Pericytes', 'Oligodendrocytes','Neurons',
                                                         'Muscle cells', 'Microglia', 'Immature Neurons', 'Ependymal','Endothelial', 
                                                         'Choroid Plexus', 'Astrocytes')))
#Celltype proportions between organs
#PBS WT counts
chimeric_mock[[]] %>% 
  group_by(Genotype, Treatment, Timepoint, Organ, manualAnnotation) %>% 
  dplyr::summarise(cell_count = n()) %>% 
  dplyr::filter(Treatment == 'PBS' & Genotype == 'WT') %>% 
  dplyr::mutate(freq = cell_count / sum(cell_count)) %>% 
  ggplot(aes(x = Timepoint, y = freq, fill = manualAnnotation))+
  geom_bar(stat = 'identity')+
  facet_wrap(~Organ)+
  scale_fill_manual(values = umap_color_list)+
  theme_classic()+
  ggtitle('PBS WT')

#PBS IPS1
chimeric_mock[[]] %>% 
  group_by(Genotype, Treatment, Timepoint, Organ, manualAnnotation) %>% 
  dplyr::summarise(cell_count = n()) %>% 
  dplyr::filter(Treatment == 'PBS' & Genotype == 'IPS1') %>% 
  dplyr::mutate(freq = cell_count / sum(cell_count)) %>% 
  ggplot(aes(x = Timepoint, y = freq, fill = manualAnnotation))+
  geom_bar(stat = 'identity')+
  facet_wrap(~Organ)+
  scale_fill_manual(values = umap_color_list)+
  theme_classic()+
  ggtitle('PBS IPS1')

#Infected wt
chimeric_mock[[]] %>% 
  group_by(Genotype, Treatment, Timepoint, Organ, manualAnnotation) %>% 
  dplyr::summarise(cell_count = n()) %>% 
  dplyr::filter(Treatment == 'rChLGTV' & Genotype == 'WT') %>% 
  dplyr::mutate(freq = cell_count / sum(cell_count)) %>% 
  ggplot(aes(x = Timepoint, y = freq, fill = manualAnnotation))+
  geom_bar(stat = 'identity')+
  facet_wrap(~Organ)+
  scale_fill_manual(values = umap_color_list)+
  theme_classic()+
  ggtitle('rChLGTV WT')

#Infected IPS1
chimeric_mock[[]] %>% 
  group_by(Genotype, Treatment, Timepoint, Organ, manualAnnotation) %>% 
  dplyr::summarise(cell_count = n()) %>% 
  dplyr::filter(Treatment == 'rChLGTV' & Genotype == 'IPS1') %>% 
  dplyr::mutate(freq = cell_count / sum(cell_count)) %>% 
  ggplot(aes(x = Timepoint, y = freq, fill = manualAnnotation))+
  geom_bar(stat = 'identity')+
  facet_wrap(~Organ)+
  scale_fill_manual(values = umap_color_list)+
  theme_classic()+
  ggtitle('rChLGTV WT')


#Upset plots by genotype and time - how similar are degs between organs?
cerebellum <- subset(chimeric_mock, Organ == 'Cerebellum')
cerebrum <- subset(chimeric_mock, Organ == 'Cerebrum')

#Cerebellum deg counts
unique_sample_combos <- unique(paste(cerebellum$Genotype, cerebellum$Treatment, cerebellum$Timepoint))
unique_infected <- unique_sample_combos[!grepl('PBS', x = unique_sample_combos)]
unique_infected_list <- stringr::str_split(unique_infected, pattern = ' ', n = 3)

cerebellum_marker_list <- lapply(unique_infected_list, FUN = function(x){
  cur_genotype = x[1]
  cur_treatment = x[2]
  cur_time = x[3]
  cur_cerebellum <- subset(cerebellum, (Genotype == cur_genotype & Treatment == cur_treatment &
                             Timepoint == cur_time) | Treatment == 'PBS')
  #table(cur_cerebellum$Timepoint, cur_cerebellum$Treatment,cur_cerebellum$Genotype)
  cur_cerebellum_markers <- FindAllMarkers(object = cur_cerebellum, group.by = 'Treatment', test.use = 'MAST', only.pos = TRUE)
})

names(cerebellum_marker_list) = lapply(unique_infected_list, FUN = paste, collapse = '_')

sig_genes_cerebellum <- lapply(cerebellum_marker_list, FUN = function(x){
  sig_genes <- x %>% as.data.frame() %>% 
    dplyr::filter(cluster == 'rChLGTV' & p_val_adj < 0.01 & avg_log2FC > 1)
})

sig_gene_counts_cerebellum <- lapply(cerebellum_marker_list, FUN = function(x){
  sig_genes <- x %>% as.data.frame() %>% 
    dplyr::filter(cluster == 'rChLGTV' & p_val_adj < 0.01 & avg_log2FC > 1) %>% 
    nrow()
})

data.frame(comp = names(sig_gene_counts_cerebellum), counts = unlist(sig_gene_counts_cerebellum)) %>% 
  tidyr::separate(col = comp, into = c('genotype', 'treatment', 'time'), sep = '_') %>% 
  ggplot(aes(x = time, y = counts, fill = genotype))+
  facet_wrap(~genotype)+
  geom_bar(stat = 'identity')

#Cerebrum deg counts
cerebrum_marker_list <- lapply(unique_infected_list, FUN = function(x){
  cur_genotype = x[1]
  cur_treatment = x[2]
  cur_time = x[3]
  cur_cerebrum <- subset(cerebrum, (Genotype == cur_genotype & Treatment == cur_treatment &
                                          Timepoint == cur_time) | Treatment == 'PBS')
  #table(cur_cerebrum$Timepoint, cur_cerebrum$Treatment,cur_cerebrum$Genotype)
  cur_cerebrum_markers <- FindAllMarkers(object = cur_cerebrum, group.by = 'Treatment', test.use = 'MAST', only.pos = TRUE)
  cur_cerebrum_markers
})

names(cerebrum_marker_list) = lapply(unique_infected_list, FUN = paste, collapse = '_')

saveRDS(cerebrum_marker_list, file = '~/Documents/ÖverbyLab/cerebrum_cerebellum_degs/wt_cerebrum_degs.rds')

sig_genes_cerebrum <- lapply(cerebrum_marker_list, FUN = function(x){
  sig_genes <- x %>% as.data.frame() %>% 
    dplyr::filter(cluster == 'rChLGTV' & p_val_adj < 0.01 & avg_log2FC > 1)
})

sig_gene_counts_cerebrum <- lapply(cerebrum_marker_list, FUN = function(x){
  sig_genes <- x %>% as.data.frame() %>% 
    dplyr::filter(cluster == 'rChLGTV' & p_val_adj < 0.01 & avg_log2FC > 1) %>% 
    nrow()
})

data.frame(comp = names(sig_gene_counts_cerebrum), counts = unlist(sig_gene_counts_cerebrum)) %>% 
  tidyr::separate(col = comp, into = c('genotype', 'treatment', 'time'), sep = '_') %>% 
  ggplot(aes(x = time, y = counts, fill = genotype))+
  facet_wrap(~genotype)+
  geom_bar(stat = 'identity')+
  ggtitle('cerebrum')

#Compare degs between both datasets
for(i in 1:length(sig_gene_counts_cerebrum)){
  cur_cerebellum = rownames(sig_genes_cerebellum[[i]])
  cur_cerebrum = rownames(sig_genes_cerebrum[[i]])
  
  cerebllum_only = sum(!cur_cerebellum %in% cur_cerebrum)
  cerebrum_only <- sum(!cur_cerebrum %in% cur_cerebellum)
  both_organs <- sum(cur_cerebellum %in% cur_cerebrum)
  
  print(names(sig_gene_counts_cerebrum)[i])
  print(paste('cerebellum only count:', cerebllum_only))
  print(paste('cerebrum only count:', cerebrum_only))
  print(paste('found in both:', both_organs))
}


