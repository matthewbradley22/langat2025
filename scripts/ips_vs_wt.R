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
chimeric_infected <- subset(ParseSeuratObj_int, Treatment == 'rChLGTV')

chimeric_infected <- prepSeuratObj(chimeric_infected)
ElbowPlot(chimeric_infected, ndims = 30)
chimeric_infected <- prepUmapSeuratObj(chimeric_infected, nDims = 20, reductionName = 'ch_umap')
DimPlot(chimeric_infected, group.by = 'manualAnnotation', reduction = 'ch_umap', cols = newCols)

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

deg_counts_by_time[[1]] <- mock_markers_up_sig
deg_counts_by_time[[2]] <- mock_markers_down_sig

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
  
  deg_counts_by_time[[length(deg_counts_by_time) + 1]] <- inf_markers_up_sig
  deg_counts_by_time[[length(deg_counts_by_time) + 1]] <- inf_markers_down_sig
}

names(deg_counts_by_time) = c('mock_up', 'mock_down', 'day3_up', 'day3_down',
                              'day4_up', 'day4_down', 'day5_up', 'day5_down')

deg_counts_by_time_counts = lapply(deg_counts_by_time, nrow)
ips_wt_deg_counts <- data.frame(comp = names(deg_counts_by_time), count = unlist(deg_counts_by_time_counts),
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

#Pathway analysis of up in wt and ips at each timepoint
wt_day3_up_paths <- gprofiler2::gost(query = rownames(deg_counts_by_time$day3_up), organism = 'mmusculus', evcodes = TRUE, sources = c('GO:BP', 'KEGG'))
ips_day3_up_paths <- gprofiler2::gost(query = rownames(deg_counts_by_time$day3_down), organism = 'mmusculus', evcodes = TRUE, sources = c('GO:BP', 'KEGG'))
wt_day5_up_paths <- gprofiler2::gost(query = rownames(deg_counts_by_time$day5_up), organism = 'mmusculus', evcodes = TRUE, sources = c('GO:BP', 'KEGG'))
ips_day5_up_paths <- gprofiler2::gost(query = rownames(deg_counts_by_time$day5_down), organism = 'mmusculus', evcodes = TRUE, sources = c('GO:BP', 'KEGG'))

#Number of markers mock vs infected by timepoint split by genotype
#Compare both against all mocks
chimeric_mock_wt <- subset(chimeric_mock, Genotype == 'WT' | Treatment == 'PBS')
chimeric_mock_ips <- subset(chimeric_mock, Genotype == 'IPS1' | Treatment == 'PBS')

deg_counts_by_time_m_vs_i <- list()

#Skip for loop with
deg_counts_by_time_m_vs_i <- readRDS('~/Documents/ÖverbyLab/single_cell_ISG_figures/sc_celltype_fig_plots/mock_infected_deg_lists.rds')

#Can just include day 4 too since we're comparing to all pbs, doesn't need day 4 pbs to be present
times_with_both_treatments = c('Day 3', 'Day 5')

#Can skip
for(i in 1:length(times_with_both_treatments)){
  cur_timepoint_wt <- subset(chimeric_mock_wt, Timepoint == times_with_both_treatments[[i]])
  cur_timepoint_ips <- subset(chimeric_mock_ips, Timepoint == times_with_both_treatments[[i]])
  
  wt_markers <- FindMarkers(cur_timepoint_wt, test.use = 'MAST', group.by = 'Treatment', ident.1 = 'rChLGTV')
  wt_markers_up_sig <- wt_markers[wt_markers$p_val_adj < 0.01 & (wt_markers$avg_log2FC) > 1,]
  wt_markers_down_sig <- wt_markers[wt_markers$p_val_adj < 0.01 & (wt_markers$avg_log2FC) < -1,]
  
  ips_markers <- FindMarkers(cur_timepoint_ips, test.use = 'MAST', group.by = 'Treatment', ident.1 = 'rChLGTV')
  ips_markers_up_sig <- ips_markers[ips_markers$p_val_adj < 0.01 & (ips_markers$avg_log2FC) > 1,]
  ips_markers_down_sig <- ips_markers[ips_markers$p_val_adj < 0.01 & (ips_markers$avg_log2FC) < -1,]
  
  deg_counts_by_time_m_vs_i[[length(deg_counts_by_time_m_vs_i) + 1]] <- wt_markers_up_sig
  deg_counts_by_time_m_vs_i[[length(deg_counts_by_time_m_vs_i) + 1]] <- wt_markers_down_sig
  deg_counts_by_time_m_vs_i[[length(deg_counts_by_time_m_vs_i) + 1]] <- ips_markers_up_sig
  deg_counts_by_time_m_vs_i[[length(deg_counts_by_time_m_vs_i) + 1]] <- ips_markers_down_sig
}


names(deg_counts_by_time_m_vs_i) = c('WT_up_3', 'WT_down_3', 'IPS_up_3', 'IPS_down_3',
                                     'WT_up_5', 'WT_down_5', 'IPS_up_5', 'IPS_down_5')

#saveRDS(deg_counts_by_time_m_vs_i, file = '~/Documents/ÖverbyLab/single_cell_ISG_figures/sc_celltype_fig_plots/mock_infected_deg_lists.rds')

deg_counts <- unname(unlist(lapply(deg_counts_by_time_m_vs_i, nrow)))

m_vs_i_deg_df <- data.frame(comp = names(deg_counts_by_time_m_vs_i),
                            deg_count = deg_counts,
                            comp_geno = factor(c('WT', 'WT', 'IPS', 'IPS', 'WT', 'WT', 'IPS', 'IPS'), levels = c('WT', 'IPS')),
                            comp_time = c(rep('3', 4), rep('5', 4)),
                            comp_direction = factor(rep(c('up', 'down'), 4), levels = c('up', 'down')))

ggplot(m_vs_i_deg_df, aes(x = comp_time, y = deg_count, fill = comp_direction))+
  facet_wrap(~comp_geno)+
  geom_bar(stat = 'identity', position = 'dodge')+
  xlab('')+
  ylab('')+
  ggtitle('DEGs infected vs mock')+
  scale_fill_manual(values = c(newCols[8], newCols[10]))+
  theme(axis.text = element_text(size = 24),
        plot.title = element_text(size = 26),
        legend.text = element_text(size = 24),
        axis.title = element_text(size = 24),
        axis.line = element_line(color = 'black'),
        panel.background = element_rect(fill = 'white', colour = 'white'),
        legend.title = element_text(size = 0))

#Pathway differences between ips and wt across timepoints
#Day 3 up

day_3_up_wt <- deg_counts_by_time_m_vs_i$WT_up_3
day_3_up_ips <- deg_counts_by_time_m_vs_i$IPS_up_3
day_5_up_wt <- deg_counts_by_time_m_vs_i$WT_up_5
day_5_up_ips <- deg_counts_by_time_m_vs_i$IPS_up_5
day_5_down_ips <- deg_counts_by_time_m_vs_i$IPS_down_5

#Do not produce log file
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
vd_day3 <- VennDiagram::venn.diagram(list(wt = rownames(day_3_up_wt), ips = rownames(day_3_up_ips)), filename = NULL,
                                     fill = brewer.pal(3, "Pastel2")[1:2], cex = 1.5, cat.cex = 1.5)
grid::grid.draw(vd_day3)

vd_day5 <- VennDiagram::venn.diagram(list(wt = rownames(day_5_up_wt), ips = rownames(day_5_up_ips)), filename = NULL,
                                     fill = brewer.pal(3, "Pastel2")[1:2], cex = 1.5, cat.cex = 1.5)
grid::grid.draw(vd_day5)

vd_day3_wt_5_ips <- VennDiagram::venn.diagram(list(wt = rownames(day_3_up_wt), ips = rownames(day_5_up_ips)), filename = NULL,
                                     fill = brewer.pal(3, "Pastel2")[1:2], cex = 1.5, cat.cex = 1.5)
grid::grid.draw(vd_day3_wt_5_ips)

#Up day 3 wt but not ips
day_3_only_wt_up <- rownames(day_3_up_wt[!rownames(day_3_up_wt) %in% rownames(day_3_up_ips),])
day_3_only_wt_up_paths <- gprofiler2::gost(query = day_3_only_wt_up, organism = 'mmusculus', evcodes = TRUE, sources = c('GO:BP', 'KEGG'))

plot_go_results <- function(go_dat){
  ggplot(head(go_dat$result, n = 10), aes(x = -log10(p_value), y = reorder(term_name, p_value)))+
    geom_bar(stat = 'identity')+
    ylab('')+
    xlab('')+
    theme(axis.text = element_text(size = 20))
}

plot_go_results(day_3_only_wt_up_paths)

#Day 5 what is ips upregulating that wt is not
day_5_only_ips_up <- rownames(day_5_up_ips[!rownames(day_5_up_ips) %in% rownames(day_5_up_wt),])
day_5_only_ips_up_paths <- gprofiler2::gost(query = day_5_only_ips_up, organism = 'mmusculus', evcodes = TRUE, sources = c('GO:BP', 'KEGG'))
ggplot(head(day_5_only_ips_up_paths$result, n = 10), aes(x = -log10(p_value), y = reorder(term_name, p_value)))+
  geom_bar(stat = 'identity')+
  ylab('')+
  xlab('')+
  theme(axis.text = element_text(size = 20))

#What is wt upregulating that ips is not
day_5_only_wt_up <- rownames(day_5_up_wt[!rownames(day_5_up_wt) %in% rownames(day_5_up_ips),])
day_5_only_wt_up_paths <- gprofiler2::gost(query = day_5_only_wt_up, organism = 'mmusculus', evcodes = TRUE, sources = c('GO:BP', 'KEGG'))

ggplot(head(day_5_only_wt_up_paths$result, n = 10), aes(x = -log10(p_value), y = reorder(term_name, p_value)))+
  geom_bar(stat = 'identity')+
  ylab('')+
  xlab('')+
  theme(axis.text = element_text(size = 20))

#Downregulated differences
day_5_down_wt <- deg_counts_by_time_m_vs_i$WT_down_5
day_5_down_ips <- deg_counts_by_time_m_vs_i$IPS_down_5

vd_day5_down <- VennDiagram::venn.diagram(list(wt = rownames(day_5_down_wt), ips = rownames(day_5_down_ips)), filename = NULL,
                                     fill = brewer.pal(3, "Pastel2")[1:2], cex = 1.5, cat.cex = 1.5)
grid::grid.draw(vd_day5_down)

common_down_5 <- rownames(day_5_down_wt)[rownames(day_5_down_wt) %in% rownames(day_5_down_ips)]

day_5_down_paths <- gprofiler2::gost(query = (common_down_5), organism = 'mmusculus', evcodes = TRUE, sources = c('GO:BP', 'KEGG'))
day_5_down_paths$result

#Downregulated at day 5 ips
day_5_ips_down_paths <- gprofiler2::gost(query = rownames(day_5_down_ips), organism = 'mmusculus', evcodes = TRUE, sources = c('GO:BP', 'KEGG'))
head(day_5_ips_down_paths$result, n = 10)

#look at specific genes across groups
chimeric_mock$time_geno_treatment <- paste(chimeric_mock$Timepoint, chimeric_mock$Genotype, chimeric_mock$Treatment, sep = '_')

plot_gene_across_groups <- function(dat, gene){
  select_gene_data <- DotPlot(dat, features = gene, group.by = 'time_geno_treatment', scale = FALSE)$data 
  select_gene_data_meta <- str_split_fixed(select_gene_data$id, pattern = '_', n = 3)
  colnames(select_gene_data_meta) <- c('timepoint', 'genotype', 'treatment')
  select_gene_data <- cbind(select_gene_data, select_gene_data_meta)
  dot_plot_gene <- ggplot(select_gene_data, aes(x = timepoint, y = treatment, size = pct.exp, color = avg.exp.scaled))+
    facet_wrap(~genotype)+
    geom_point()+
    ggtitle(gene)+
    scale_color_gradientn(colours = c("#F03C0C","#FFF0A3","white"), 
                          values = c(1.0,0.5,0))
  print(dot_plot_gene)
}

plot_gene_across_groups(chimeric_mock, 'Isg20')


#Do above analysis by celltype, seems like trends may differ alot across celltypes
deg_counts_time_m_vs_i_celltype <- list()

#Skip slow for loop
deg_counts_time_m_vs_i_celltype <- readRDS("~/Documents/ÖverbyLab/celltype_venn_data.rds")

celltypes_of_interest = c('Microglia', 'Astrocytes', 'Choroid Plexus', 'Endothelial', 'Ependymal')
times = c('Day 3', 'Day 4', 'Day 5')

#Can skip
for(i in 1:length(celltypes_of_interest)){
  print(paste("starting cell", celltypes_of_interest[i]))
  current_celltype <- celltypes_of_interest[[i]]
  current_cell_data_wt <- subset(chimeric_mock_wt, manualAnnotation == current_celltype)
  current_cell_data_ips <- subset(chimeric_mock_ips, manualAnnotation == current_celltype)
  for(j in 1:length(times)){
    print(paste("starting time", times[j]))
    cur_timepoint_wt <- subset(current_cell_data_wt, Timepoint == times[[j]] | Treatment == 'PBS')
    cur_timepoint_ips <- subset(current_cell_data_ips, Timepoint == times[[j]] | Treatment == 'PBS')
    
    wt_markers <- FindMarkers(cur_timepoint_wt, test.use = 'MAST', group.by = 'Treatment', ident.1 = 'rChLGTV')
    wt_markers_up_sig <- wt_markers[wt_markers$p_val_adj < 0.01 & (wt_markers$avg_log2FC) > 1,]
    wt_markers_down_sig <- wt_markers[wt_markers$p_val_adj < 0.01 & (wt_markers$avg_log2FC) < -1,]
    
    ips_markers <- FindMarkers(cur_timepoint_ips, test.use = 'MAST', group.by = 'Treatment', ident.1 = 'rChLGTV')
    ips_markers_up_sig <- ips_markers[ips_markers$p_val_adj < 0.01 & (ips_markers$avg_log2FC) > 1,]
    ips_markers_down_sig <- ips_markers[ips_markers$p_val_adj < 0.01 & (ips_markers$avg_log2FC) < -1,]
    
    deg_counts_time_m_vs_i_celltype[[length(deg_counts_time_m_vs_i_celltype) + 1]] <- wt_markers_up_sig
    deg_counts_time_m_vs_i_celltype[[length(deg_counts_time_m_vs_i_celltype) + 1]] <- wt_markers_down_sig
    deg_counts_time_m_vs_i_celltype[[length(deg_counts_time_m_vs_i_celltype) + 1]] <- ips_markers_up_sig
    deg_counts_time_m_vs_i_celltype[[length(deg_counts_time_m_vs_i_celltype) + 1]] <- ips_markers_down_sig
    
  }
}

#Create names for new list
name_of_each_comp = c('wt_up', 'wt_down', 'ips_up', 'ips_down')
names_without_comps <- paste(rep(celltypes_of_interest, each = length(times) * length(name_of_each_comp)), rep(times, each = length(name_of_each_comp)), sep = '_')
names_with_comps <- paste(names_without_comps, name_of_each_comp, sep = '_')
names_with_comps <- gsub(' ', '', names_with_comps)
names(deg_counts_time_m_vs_i_celltype) = names_with_comps

#saveRDS(deg_counts_time_m_vs_i_celltype, file = "~/Documents/ÖverbyLab/celltype_venn_data.rds")

#Need to remove space from choroid plexus as this wasn't how list was named
celltypes_of_interest_nospace = c('Microglia', 'Astrocytes', 'ChoroidPlexus', 'Endothelial', 'Ependymal')

#Do not let venndiagram write log files for every comp
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

#Generate ven diagrams for each celltype
for(i in 1:length(celltypes_of_interest_nospace)){
  current_celltype <- celltypes_of_interest_nospace[i]
  for(j in 1:length(times)){
    cur_time <- gsub(' ', '', times[[j]])
    cur_celltype_wt_up_name <- paste(current_celltype, cur_time, 'wt_up', sep = '_')
    cur_celltype_ips_up_name <- paste(current_celltype, cur_time, 'ips_up', sep = '_')
    cur_wt_up_genes <- rownames(deg_counts_time_m_vs_i_celltype[[cur_celltype_wt_up_name]])
    cur_ips_up_genes <- rownames(deg_counts_time_m_vs_i_celltype[[cur_celltype_ips_up_name]])
    
    path_to_save <-  '~/Documents/ÖverbyLab/celltype_venns/'
    current_comp <- paste(current_celltype, cur_time, 'wt_vs_ips_up', sep = '_')
    VennDiagram::venn.diagram(list(wt = cur_wt_up_genes, ips = cur_ips_up_genes), filename = paste0(path_to_save, current_comp, '.png'),
                              fill = brewer.pal(3, "Pastel2")[1:2], cex = 3, cat.cex = 1.5)
  }
}

#Path analysis for select differences
#Microglia day 5
micro_5_wt_only <- rownames(deg_counts_time_m_vs_i_celltype$Microglia_Day5_wt_up)[!rownames(deg_counts_time_m_vs_i_celltype$Microglia_Day5_wt_up) %in% 
                                                       rownames(deg_counts_time_m_vs_i_celltype$Microglia_Day5_ips_up)]
micro_5_wt_only_paths <- gprofiler2::gost(query = micro_5_wt_only, organism = 'mmusculus', evcodes = TRUE, sources = c('GO:BP', 'KEGG'))
plot_go_results(micro_5_wt_only_paths)

#Endothelial day 5
endo_5_wt_only <- rownames(deg_counts_time_m_vs_i_celltype$Endothelial_Day5_wt_up)[!rownames(deg_counts_time_m_vs_i_celltype$Endothelial_Day5_wt_up) %in% 
                                                                                    rownames(deg_counts_time_m_vs_i_celltype$Endothelial_Day5_ips_up)]
endo_5_wt_only_paths <- gprofiler2::gost(query = endo_5_wt_only, organism = 'mmusculus', evcodes = TRUE, sources = c('GO:BP', 'KEGG'))
plot_go_results(endo_5_wt_only_paths)
plot_gene_across_groups(subset(chimeric_mock, manualAnnotation == 'Endothelial'), 'Cx3cl1')

#Ependymal day 5
epen_5_wt_only <- rownames(deg_counts_time_m_vs_i_celltype$Ependymal_Day5_wt_up)[!rownames(deg_counts_time_m_vs_i_celltype$Ependymal_Day5_wt_up) %in% 
                                                                                     rownames(deg_counts_time_m_vs_i_celltype$Ependymal_Day5_ips_up)]
epen_5_wt_only_paths <- gprofiler2::gost(query = epen_5_wt_only, organism = 'mmusculus', evcodes = TRUE, sources = c('GO:BP', 'KEGG'))
plot_go_results(epen_5_wt_only_paths)
plot_gene_across_groups((chimeric_mock), 'Cx3cl1')
