library(Seurat)
library(RColorBrewer)
library(dplyr)
library(UpSetR)
library(ggplot2)
library(tidyr)
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
chimeric_infected <- subset(ParseSeuratObj_int, Treatment == 'rChLGTV')

chimeric_infected <- prepSeuratObj(chimeric_infected)
ElbowPlot(chimeric_infected, ndims = 30)
chimeric_infected <- prepUmapSeuratObj(chimeric_infected, nDims = 20, reductionName = 'ch_umap')
DimPlot(chimeric_infected, group.by = 'Genotype', reduction = 'ch_umap')

#Markers upregulated in wt
wt_up <- FindMarkers(chimeric_infected, group.by = 'Genotype', ident.1 = 'WT', only.pos = TRUE, test.use = 'MAST')
wt_up_sig <- wt_up %>% as.data.frame() %>% dplyr::filter(p_val_adj < 0.01 & avg_log2FC > 1)

#Compare mavs by orig.ident to confirm it lines up with barcoding plate csv we have, 09 - 24 is ips
mavs_dat <- DotPlot(ParseSeuratObj_int, features = 'Mavs', group.by = 'orig.ident', scale = FALSE)$data
mavs_dat$id <- mavs_dat$id %>% as.character() %>% as.numeric()
mavs_dat <- mavs_dat %>% dplyr::mutate(genotype = ifelse(id <= 24, yes = 'IPS1', no = 'WT'))

ggplot(mavs_dat, aes(x = features.plot, y = id, size = pct.exp, color = avg.exp.scaled))+
  facet_wrap(~genotype)+
  geom_point()+ 
  theme(panel.background = element_rect(fill = 'white', colour = 'white'))+
  scale_color_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white"), 
                       values = c(1.0,0.7,0.4,0))+
  ylab('plate id')+
  xlab('')+
  ggtitle('MAVS expression')

#Number of degs wt vs ips by timepoint
deg_counts_by_time <- list()
mock_only <- subset(chimeric_mock, Treatment == 'PBS')

#Number of markers mock vs mock
mock_markers <- FindMarkers(mock_only, test.use = 'MAST', group.by = 'Genotype', ident.1 = 'WT')
mock_markers_up_sig <- mock_markers[mock_markers$p_val_adj < 0.01 & (mock_markers$avg_log2FC) > 1,]
mock_markers_down_sig <- mock_markers[mock_markers$p_val_adj < 0.01 & (mock_markers$avg_log2FC) < -1,]
nrow(mock_markers_up_sig)
nrow(mock_markers_down_sig)

deg_counts_by_time <- c(deg_counts_by_time, mock_wt_up = nrow(mock_markers_up_sig))
deg_counts_by_time <- c(deg_counts_by_time, mock_wt_down = nrow(mock_markers_down_sig))

#Number of markers infected wt vs ips
inf_chimeric <- subset(chimeric_mock, Treatment == 'rChLGTV')
times <- unique(inf_chimeric$Timepoint)

for(i in 1:length(times)){
  cur_timepoint <- subset(inf_chimeric, Timepoint == times[[i]])
  inf_markers <- FindMarkers(cur_timepoint, test.use = 'MAST', group.by = 'Genotype', ident.1 = 'WT')
  
  inf_markers_up_sig <- inf_markers[inf_markers$p_val_adj < 0.01 & (inf_markers$avg_log2FC) > 1,]
  inf_markers_down_sig <- inf_markers[inf_markers$p_val_adj < 0.01 & (inf_markers$avg_log2FC)  < -1,]
  
  up_name = paste(times[[i]], 'up', sep = '_')
  down_name= paste(times[[i]], 'down', sep = '_')
  
  deg_counts_by_time <- c(deg_counts_by_time, up_name = nrow(inf_markers_up_sig))
  deg_counts_by_time <- c(deg_counts_by_time, down_name =  nrow(inf_markers_down_sig))
}

names(deg_counts_by_time) = c('mock_up', 'mock_down', 'day3_up', 'day3_down',
                              'day4_up', 'day4_down', 'day5_up', 'day5_down')

ips_wt_deg_counts <- data.frame(comp = names(deg_counts_by_time), count = unlist(deg_counts_by_time),
                                comp_type = factor(c('mock', 'mock', 'day3', 'day3', 'day4', 'day4', 'day5', 'day5'),
                                                   levels = c('mock', 'day3', 'day4', 'day5')),
                                comp_direction = rep(c('wt_up', 'ips_up'), 4))

pdf('~/Documents/ÖverbyLab/single_cell_ISG_figures/sc_celltype_fig_plots/deg_count_time.pdf', height = 8, width = 10)
ggplot(ips_wt_deg_counts, aes(x = comp_type, y = count, fill = comp_direction)) +
  geom_bar(stat = 'identity', position = 'dodge')+
  ylab('Number of DEGs')+
  xlab('')+
  ggtitle('Number of DEGs')+
  scale_fill_manual(values = c(newCols[4], newCols[7]))+
  theme(axis.text = element_text(size = 24),
        plot.title = element_text(size = 26),
        legend.text = element_text(size = 24),
        axis.title = element_text(size = 24),
        axis.line = element_line(color = 'black'),
        panel.background = element_rect(fill = 'white', colour = 'white'),
        legend.title = element_text(size = 0))
dev.off()

#Number of markers mock vs infected by timepoint split by genotype
#Compare both against all mocks
chimeric_mock_wt <- subset(chimeric_mock, Genotype == 'WT' | Treatment == 'PBS')
chimeric_mock_ips <- subset(chimeric_mock, Genotype == 'IPS1' | Treatment == 'PBS')

deg_counts_by_time_m_vs_i <- list()

times_with_both_treatments = c('Day 3', 'Day 5')

for(i in 1:length(times_with_both_treatments)){
  cur_timepoint_wt <- subset(chimeric_mock_wt, Timepoint == times_with_both_treatments[[i]])
  cur_timepoint_ips <- subset(chimeric_mock_ips, Timepoint == times_with_both_treatments[[i]])
  
  wt_markers <- FindMarkers(cur_timepoint_wt, test.use = 'MAST', group.by = 'Treatment', ident.1 = 'rChLGTV')
  wt_markers_up_sig <- wt_markers[wt_markers$p_val_adj < 0.01 & (wt_markers$avg_log2FC) > 1,]
  wt_markers_down_sig <- wt_markers[wt_markers$p_val_adj < 0.01 & (wt_markers$avg_log2FC) < -1,]
  
  ips_markers <- FindMarkers(cur_timepoint_ips, test.use = 'MAST', group.by = 'Treatment', ident.1 = 'rChLGTV')
  ips_markers_up_sig <- ips_markers[ips_markers$p_val_adj < 0.01 & (ips_markers$avg_log2FC) > 1,]
  ips_markers_down_sig <- ips_markers[ips_markers$p_val_adj < 0.01 & (ips_markers$avg_log2FC) < -1,]
  
  deg_counts_by_time_m_vs_i <- c(deg_counts_by_time_m_vs_i, wt_markers_up_sig, wt_markers_down_sig,
                                 ips_markers_up_sig, ips_markers_down_sig)
}

#saveRDS(deg_counts_by_time_m_vs_i, file = '~/Documents/ÖverbyLab/single_cell_ISG_figures/sc_celltype_fig_plots/mock_infected_deg_lists.rds')

names(deg_counts_by_time_m_vs_i) = c('WT_up_3', 'WT_down_3', 'IPS_up_3', 'IPS_down_3',
                                     'WT_up_5', 'WT_down_5', 'IPS_up_5', 'IPS_down_5')

m_vs_i_deg_df <- data.frame(comp = c('WT_up_3', 'WT_down_3', 'IPS_up_3', 'IPS_down_3',
                                     'WT_up_5', 'WT_down_5', 'IPS_up_5', 'IPS_down_5'),
                            deg_count = unlist(unname(deg_counts_by_time_m_vs_i)),
                            comp_geno = factor(c('WT', 'WT', 'IPS', 'IPS', 'WT', 'WT', 'IPS', 'IPS'), levels = c('WT', 'IPS')),
                            comp_time = c(rep('3', 4), rep('5', 4)),
                            comp_direction = factor(rep(c('up', 'down'), 4), levels = c('up', 'down')))

ggplot(m_vs_i_deg_df, aes(x = comp_time, y = deg_count, fill = comp_direction))+
  facet_wrap(~comp_geno)+
  geom_bar(stat = 'identity', position = 'dodge')


