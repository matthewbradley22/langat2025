#Packages and functions
library(gridExtra)
library(ggpubr)
library(Seurat)
library(gprofiler2)
library(UCell)
library(RColorBrewer)
library(rstatix)
library(irGSEA)
library(UCell)
library(msigdbr)
source('~/Documents/ÖverbyLab//scripts/langatFunctions.R')

#Load data
ParseSeuratObj_int <- LoadSeuratRds("~/Documents/ÖverbyLab/data/FilteredRpcaIntegratedDatNoDoublets.rds") 
ParseSeuratObj_int$hasVirus = ifelse(ParseSeuratObj_int$virusCountPAdj >= 10, 1, 0)

#Check data
umap_color_list <- c( "#7047A1", "#B370AE","#292270",  "#166DF0","#6D92F8",  "#6DC3F8", "#8a0000","#F76363", "#FF96A2", 
                      "#D6644B", "#F08C3A", "#fdc087","#074F00", "#208d1f","#7bcd79", 
                                         "gray")

ParseSeuratObj_int$manualAnnotation <- factor(ParseSeuratObj_int$manualAnnotation, 
                                              levels = c('Astrocytes', 'Choroid Plexus', 'Endothelial',
                                                         'Ependymal', 'Immature Neurons', 'Microglia', 'Muscle cells',
                                                         'Neurons', 'Oligodendrocytes', 'Pericytes', 'B Cells',
                                                         'Granulocytes', 'Macrophage/Monocytes', 'Nk cells', 
                                                                      'T cells', 'unknown'))

DimPlot(ParseSeuratObj_int, label = FALSE, group.by = 'manualAnnotation', reduction = 'umap.integrated',
        cols = umap_color_list)

#Wt cerebrum celltypes across times
wt_cerebrum <-  subset(ParseSeuratObj_int, Treatment %in% c('PBS', 'rLGTV') & Organ == 'Cerebrum' & Genotype == 'WT')
table(wt_cerebrum$manualAnnotation, wt_cerebrum$Genotype, wt_cerebrum$Treatment, wt_cerebrum$Timepoint)

#Plot resident cell ratios
wt_cerebrum_resident_mock <- subset(wt_cerebrum, manualAnnotation %in% c('Astrocytes', 'Oligodendrocytes', 'Microglia', 'Endothelial',
                                                                    'Choroid Plexus', 'Immature Neurons', 'Ependymal', 'Pericytes',
                                                                    'Muscle cells', 'Neurons') & Treatment == 'PBS')

wt_cerebrum_resident_inf <- subset(wt_cerebrum, manualAnnotation %in% c('Astrocytes', 'Oligodendrocytes', 'Microglia', 'Endothelial',
                                                                         'Choroid Plexus', 'Immature Neurons', 'Ependymal', 'Pericytes',
                                                                         'Muscle cells', 'Neurons') & Treatment == 'rLGTV')

pdf('~/Documents/ÖverbyLab/scPlots/galectin3_proj/microglia/mock_resident_cell_bars.pdf', width = 5, height = 5)
table(wt_cerebrum_resident_mock$manualAnnotation) %>%  # can include wt_cerebrum_resident_mock$Timepoint to split by time
  as.data.frame() %>% dplyr::mutate(freq_props = Freq/sum(Freq))%>% 
  ggplot(aes(x = 'mock', y = freq_props, fill = Var1))+
  geom_bar(stat = 'identity', position = 'stack', width = 0.4)+
  scale_fill_manual(values = umap_color_list)+
  theme_classic()+
  theme(axis.text.x = element_text(size = 22),
        axis.text.y = element_text(size = 22),
        axis.title=element_text(size=22),
        plot.title = element_text(size = 22))+
  xlab('')+
  ggtitle('Mock resident')+
  ylab('Proportion of cells')
dev.off()

pdf('~/Documents/ÖverbyLab/scPlots/galectin3_proj/microglia/infected_resident_cell_bars.pdf', width = 6, height = 5)
table(wt_cerebrum_resident_inf$Timepoint, wt_cerebrum_resident_inf$manualAnnotation) %>% 
  as.data.frame() %>% dplyr::group_by(Var1) %>% dplyr::mutate(freq_props = Freq/sum(Freq))%>% 
  ggplot(aes(x = Var1, y = freq_props, fill = Var2))+
  geom_bar(stat = 'identity', position = 'stack', width = 0.6)+
  scale_fill_manual(values = umap_color_list)+
  theme_classic()+
  theme(axis.text.x = element_text(size = 22),
        axis.text.y = element_text(size = 22),
        axis.title=element_text(size=22),
        plot.title = element_text(size = 22))+
  xlab('')+
  ggtitle('Infected resident')+
  ylab('Proportion of cells')
dev.off()

#Subset to just microglia for further analysis
wt_cerebrum_microglia <- subset(wt_cerebrum, manualAnnotation == 'Microglia')
table(wt_cerebrum_microglia$Genotype, wt_cerebrum_microglia$Treatment, 
      wt_cerebrum_microglia$Timepoint, wt_cerebrum_microglia$manualAnnotation)

#Microglia UMAP
wt_cerebrum_microglia <- prepSeuratObj(wt_cerebrum_microglia)
ElbowPlot(wt_cerebrum_microglia, ndims = 40)
wt_cerebrum_microglia <- prepUmapSeuratObj(wt_cerebrum_microglia, nDims = 15, reductionName = 'micro.umap',
                                             resolution_value = 0.5)

pdf('~/Documents/ÖverbyLab/scPlots/galectin3_proj/microglia/microglia_treatment_umap.pdf', width = 5, height = 5)
DimPlot(wt_cerebrum_microglia, reduction = 'micro.umap', label = FALSE, group.by = 'Treatment',
        label.size = 6, cols = c('#6DC3F8', '#D6644B', '#208d1f'))+
  ggtitle('Microglia')+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank())+
  ylab('UMAP2')+
  xlab('UMAP1')
dev.off()

wt_cerebrum_microglia$time_by_treatment <- paste(wt_cerebrum_microglia$Treatment, wt_cerebrum_microglia$Timepoint)
pdf('~/Documents/ÖverbyLab/scPlots/galectin3_proj/microglia/microglia_time_umap.pdf', width = 5, height = 5)
DimPlot(wt_cerebrum_microglia, reduction = 'micro.umap', label = FALSE, group.by = 'time_by_treatment',
        label.size = 6, cols = c("#D6644B", "#8a0000",  "#6DC3F8", "#166DF0","#292270"))+
  ggtitle('Microglia')+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank())+
  ylab('UMAP2')+
  xlab('UMAP1')
dev.off()


wt_cerebrum_microglia_markers <- FindAllMarkers(wt_cerebrum_microglia, group.by = 'Treatment', test.use = 'MAST')

top_markers <- wt_cerebrum_microglia_markers %>% dplyr::filter(avg_log2FC > 1 & p_val_adj < 0.01) %>% 
  dplyr::group_by(cluster) %>% dplyr::slice_head(n = 6)

pdf('~/Documents/ÖverbyLab/scPlots/galectin3_proj/microglia/microglia_treatment_genes.pdf', width = 7, height = 5)
DotPlot(wt_cerebrum_microglia, features = top_markers$gene, group.by = 'Treatment', scale = FALSE)$data %>% 
  ggplot(aes(x = features.plot, y = id, size = pct.exp, fill = avg.exp.scaled))+
  geom_point(pch = 21)+
  scale_fill_gradientn(colours = c('white', '#FFD991', '#FF7530', '#FF4024'), 
                       values = c(0, 0.3, 0.6, 1))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        text = element_text(size = 20))+
  ylab('')+
  xlab('')
dev.off()

#Why do mock not create a single cluster?
#Switch cluster numbers so mock clusters are 0, 1, and 2.
wt_cerebrum_microglia$custom_clusters <- dplyr::case_when(wt_cerebrum_microglia$seurat_clusters == "5" ~ "2",
                                                          wt_cerebrum_microglia$seurat_clusters == "2" ~ "5",
                                                          .default = wt_cerebrum_microglia$seurat_clusters)

pdf('~/Documents/ÖverbyLab/scPlots/galectin3_proj/microglia/microglia_UMAP.pdf', width = 7, height = 5)
DimPlot(wt_cerebrum_microglia, reduction = 'micro.umap', label = FALSE, group.by = 'custom_clusters', 
        cols = c('#6B9973', '#A4CE05', '#E7A515', '#EDCF53', '#F78E93', '#CE3B42', '#8A0100'))+
      ggtitle('Microglia')+
      theme(axis.text = element_blank(), axis.ticks = element_blank())+
      ylab('UMAP2')+
  xlab('UMAP1')
dev.off()

wt_cerebrum_microglia[[]] %>% dplyr::group_by(Treatment, Timepoint) %>% 
  dplyr::summarise(mean_count = mean(nCount_RNA), mean_features = mean(nFeature_RNA), mean_mt = mean(percent.mt))

table(wt_cerebrum_microglia$Treatment, wt_cerebrum_microglia$seurat_clusters)

# - - - - - - - - - - - - - - - - - - 
#### Across time infected samples ####
# - - - - - - - - - - - - - - - - - - 

#Look across time at microglia
microglia_infected <- subset(wt_cerebrum_microglia, Treatment == 'rLGTV')

#Infected microglia UMAP
microglia_infected <- prepSeuratObj(microglia_infected)
ElbowPlot(microglia_infected, ndims = 40)
microglia_infected <- prepUmapSeuratObj(microglia_infected, nDims = 15, reductionName = 'micro.inf.umap',
                                           resolution_value = 0.5)

pdf('~/Documents/ÖverbyLab/scPlots/galectin3_proj/microglia/infected_micro_umap.pdf', width = 5, height = 5)
DimPlot(microglia_infected, reduction = 'micro.inf.umap', label = FALSE, group.by = 'Timepoint',
        label.size = 6, cols = c('#6DC3F8', '#D6644B', '#208d1f'))+
  ggtitle('Infected microglia')+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank())+
  ylab('UMAP2')+
  xlab('UMAP1')
dev.off()

pdf('~/Documents/ÖverbyLab/scPlots/galectin3_proj/microglia/infected_micro_cluster_umap.pdf', width = 5, height = 5)
DimPlot(microglia_infected, reduction = 'micro.inf.umap', label = FALSE)+
        #cols = c("#7047A1", "#B370AE","#292270",  "#166DF0","#6DC3F8"))+
  ggtitle('Infected microglia')+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank())+
  ylab('UMAP2')+
  xlab('UMAP1')
dev.off()

#Infected timepoint markers
microglia_infected$time_comb <- ifelse(microglia_infected$Timepoint == 'Day 5', 'Day 5', 'Day 3_4')
infected_time_markers <- FindAllMarkers(microglia_infected, group.by = 'time_comb', test.use = 'MAST')
top_infected_markers <- infected_time_markers %>% dplyr::filter(avg_log2FC > 1 & p_val_adj < 0.01) %>% 
  dplyr::group_by(cluster) %>% 
  dplyr::slice_head(n = 5)

pdf('~/Documents/ÖverbyLab/scPlots/galectin3_proj/microglia/microglia_inf_timepoint_genes.pdf', width = 8, height = 5)
DotPlot(microglia_infected, features = top_infected_markers$gene, group.by = 'Timepoint', scale = FALSE)$data %>% 
  ggplot(aes(x = features.plot, y = id, size = pct.exp, fill = avg.exp.scaled))+
  geom_point(pch = 21)+
  scale_fill_gradientn(colours = c('white', '#FFD991', '#FF7530', '#FF4024'), 
                       values = c(0, 0.3, 0.6, 1))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        text = element_text(size = 20))+
  ylab('')+
  xlab('')+
  ggtitle('Infected Microglia')
dev.off()


#Pathway analysis
early_infect_markers <- infected_time_markers %>% 
  dplyr::filter(cluster == 'Day 3_4' & p_val_adj < 0.01 & avg_log2FC > 1) 
early_infect_markers$gene %>% write.csv(file = '~/Documents/ÖverbyLab/scPlots/galectin3_proj/microglia/early_markers.csv',
                                        row.names = FALSE, quote = FALSE)

early_paths <- gprofiler2::gost(query = early_infect_markers$gene, organism = 'mmusculus', evcodes = TRUE,
                                sources = c('GO:BP', 'KEGG', 'GO:CC', 'GO:MP'))

pdf('~/Documents/ÖverbyLab/scPlots/galectin3_proj/microglia/early_time_paths.pdf', width = 6, height = 5)
ggplot(head(early_paths$result, n = 10), aes(x = -log10(p_value), y = reorder(term_name, desc(p_value))))+
  geom_bar(stat = 'identity', fill = '#6DC3F8', color = 'black')+
  theme_classic()+
  ylab('')+
  xlab('-Log10 p-val')+
  theme(text = element_text(size = 18))+
  ggtitle('Early')
dev.off()

late_infect_markers <- infected_time_markers %>% 
  dplyr::filter(cluster == 'Day 5' & p_val_adj < 0.01 & avg_log2FC > 1)

late_infect_markers$gene %>% write.csv(file = '~/Documents/ÖverbyLab/scPlots/galectin3_proj/microglia/late_markers.csv',
                                       row.names = FALSE, quote = FALSE)

late_paths <- gprofiler2::gost(query = late_infect_markers$gene, organism = 'mmusculus', evcodes = TRUE,
                               sources = c('GO:BP', 'KEGG', 'GO:CC', 'GO:MP'))
pdf('~/Documents/ÖverbyLab/scPlots/galectin3_proj/microglia/late_time_paths.pdf', width = 6, height = 5)
ggplot(head(late_paths$result, n = 10), aes(x = -log10(p_value), y = reorder(term_name, desc(p_value))))+
  geom_bar(stat = 'identity', fill = '#6DC3F8', color = 'black')+
  theme_classic()+
  ylab('')+
  xlab('-Log10 p-val')+
  theme(text = element_text(size = 18))+
  ggtitle('Late')
dev.off()

# - - - - - - - - - - - - - - - - 
#### Across time mock samples ####
# - - - - - - - - - - - - - - - - 

microglia_mock <- subset(wt_cerebrum_microglia, Treatment == 'PBS')

#Mock microglia UMAP
microglia_mock <- prepSeuratObj(microglia_mock)
ElbowPlot(microglia_mock, ndims = 40)
microglia_mock <- prepUmapSeuratObj(microglia_mock, nDims = 15, reductionName = 'micro.mock.umap',
                                        resolution_value = 0.5)

pdf('~/Documents/ÖverbyLab/scPlots/galectin3_proj/microglia/mock_micro_umap.pdf', width = 5, height = 5)
DimPlot(microglia_mock, reduction = 'micro.mock.umap', label = FALSE, group.by = 'Timepoint',
        label.size = 6, cols = c('#6DC3F8', '#D6644B', '#208d1f'))+
  ggtitle('Mock microglia')+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank())+
  ylab('UMAP2')+
  xlab('UMAP1')
dev.off()

#Look at markers separating day 3 and day 5 mock microglia
mock_markers <- FindAllMarkers(microglia_mock, group.by = 'Timepoint', test.use = 'MAST')
mock_markers_3 <- dplyr::filter(mock_markers, p_val_adj < 0.01 & avg_log2FC > 1 & cluster == 'Day 3')
mock_markers_5 <- dplyr::filter(mock_markers, p_val_adj < 0.01 & avg_log2FC > 1 & cluster == 'Day 5')

#No sig pathways for 3 dpi
mock_comp_paths <- gprofiler2::gost(query = mock_markers_5$gene, organism = 'mmusculus', evcodes = TRUE,
                                    sources = c('GO:BP', 'KEGG', 'GO:CC', 'GO:MP'))

mock_comp_paths$result

ggplot(head(mock_comp_paths$result, n = 8), aes(x = -log10(p_value), y = reorder(term_name, -p_value)))+
  geom_bar(stat = 'identity', color = 'black', fill = 'lightgrey')+
  theme_classic()+
  ylab('')+
  theme(text = element_text(size = 18))

FeaturePlot(wt_cerebrum_microglia, features = 'mt-Cytb', reduction = 'micro.umap')
FeaturePlot(microglia_mock, features = 'mt-Cytb', reduction = 'micro.mock.umap')

mock_dot_dat <- DotPlot(microglia_mock, features = head(mock_markers_5$gene, n = 10), group.by = 'Timepoint', scale = FALSE)$data

mock_dot_dat %>% dplyr::mutate(features.plot = factor(features.plot, levels = rev(head(mock_markers_5$gene, n = 10)))) %>% 
  ggplot(aes(x = id, y = features.plot, size = pct.exp, fill = avg.exp.scaled))+
  geom_point(pch = 21)+
  theme_classic()+
  scale_fill_gradientn(colours = c('white', '#FFD991', '#FF7530', '#FF4024'), 
                       values = c(0, 0.3, 0.6, 1))+
  ggtitle('Mock microglia markers')

#Look at mock markers by cluster
DimPlot(microglia_mock, reduction = 'micro.mock.umap', label = FALSE, group.by = 'custom_clusters',
        cols = c('#6B9973', '#A4CE05', '#E7A515', '#EDCF53', '#F78E93', '#CE3B42', '#8A0100'))+
  ggtitle('Mock microglia')+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank())+
  ylab('UMAP2')+
  xlab('UMAP1')

mock_cluster_markers <- FindAllMarkers(microglia_mock, group.by = 'custom_clusters', test.use = 'MAST')
sig_cluster_markers <- mock_cluster_markers %>% dplyr::filter(avg_log2FC > 1 & p_val_adj < 0.01) %>% 
  dplyr::filter(cluster %in% c(0, 1, 2))

top_mock_cluster_genes <- sig_cluster_markers %>% dplyr::group_by(cluster) %>%  dplyr::slice_head(n = 6) %>% 
  dplyr::pull(gene) %>% unique()
  
DotPlot(microglia_mock, features = rev(top_mock_cluster_genes), group.by = 'custom_clusters', scale = FALSE)$data %>% 
  dplyr::filter(id %in% c(0, 1, 2)) %>% 
  ggplot(aes(x = id, y = features.plot, size = pct.exp, fill = avg.exp.scaled))+
  geom_point(pch = 21)+
  theme_classic()+
  scale_fill_gradientn(colours = c('white', '#FFD991', '#FF7530', '#FF4024'), 
                       values = c(0, 0.3, 0.6, 1))+
  ggtitle('')

sig_cluster_markers
sig_cluster_paths <- lapply(c(0, 1, 2), FUN = function(x){
  cur_markers <- dplyr::filter(sig_cluster_markers, cluster == x)
  cur_pathways <- gprofiler2::gost(query = cur_markers$gene, organism = 'mmusculus', evcodes = TRUE,
                                   sources = c('GO:BP', 'KEGG', 'GO:CC', 'GO:MP'))
  cur_pathways
})

ggplot(head(sig_cluster_paths[[2]]$result, n = 8), aes(x = -log10(p_value), y = reorder(term_name, -p_value)))+
  geom_bar(stat = 'identity')+
  theme_classic()+
  ylab('')+
  ggtitle('Cluster 1')+
  theme(text = element_text(size = 16))

ggplot(head(sig_cluster_paths[[3]]$result, n = 8), aes(x = -log10(p_value), y = reorder(term_name, -p_value)))+
  geom_bar(stat = 'identity')+
  theme_classic()+
  ggtitle('Cluster 2')+
  ylab('')+
  theme(text = element_text(size = 16))

FeaturePlot(microglia_mock, reduction = 'micro.mock.umap', features = 'Egr1')

#Look at negative regulators of inflammation in infected samples across time
neg_reg_genes <- c('Socs3', 'Il1rn', 'Il27', 'Lilrb4a', 'Tgfb1', 'Il10', 'Il4', 'Arg1', 'Mrc1', 'Igf1', 'Csf1')

pdf('~/Documents/ÖverbyLab/scPlots/galectin3_proj/microglia/micro_neg_ref_inflammation.pdf', width = 6, height = 5)
DotPlot(microglia_infected, features =neg_reg_genes, group.by = 'Timepoint', scale = FALSE)$data %>% 
  ggplot(aes(x = features.plot, y = id, size = pct.exp, fill = avg.exp.scaled))+
  geom_point(pch = 21)+
  coord_flip()+
  scale_fill_gradientn(colours = c('white', '#FFD991', '#FF7530', '#FF4024'), 
                       values = c(0, 0.3, 0.6, 1))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        text = element_text(size = 20))+
  ylab('')+
  xlab('')+
  ggtitle('Anti-inflammatory')
dev.off()


neg_regulation_dot <- DotPlot(wt_cerebrum_microglia, 
                              features = neg_reg_genes, group.by = 'infected_clusters', scale = FALSE)$data
neg_regulation_dot$id = factor(neg_regulation_dot$id, levels = c('mock', '0', '1', '2', '3', '4'))

ggplot(neg_regulation_dot, aes(x = id, y = features.plot, fill = avg.exp.scaled, size = pct.exp))+
  geom_point(pch = 21)+
  theme_classic()+
  scale_fill_gradientn(colours = c('white', '#FFD991', '#FF7530', '#FF4024'), 
                       values = c(0, 0.3, 0.6, 1))+
  ggtitle('Anti inflammatory')


#- - - - - - - - - - - - - - - - - - - - - 
#### DEGs between infected clusters #####
#- - - - - - - - - - - - - - - - - - - - - 

#Infected microglia cluster markers and paths
infected_clust_markers <- FindAllMarkers(microglia_infected, group.by = 'seurat_clusters', 
                                         test.use = 'MAST', only.pos = TRUE)

inf_micro_path_list <- list()

for(i in 1:length(unique(infected_clust_markers$cluster))-1){
  cur_clust = dplyr::filter(infected_clust_markers, cluster == i & p_val < 0.01 & avg_log2FC > 1)
  cur_genes <- cur_clust$gene
  cur_paths <- gprofiler2::gost(query = cur_genes, organism = 'mmusculus', evcodes = TRUE,
                                sources = c('GO:BP', 'KEGG', 'GO:CC', 'GO:MP'))
  if(length(cur_paths) > 0){
    inf_micro_path_list[[length(inf_micro_path_list) + 1]] = cur_paths$result
    names(inf_micro_path_list)[length(inf_micro_path_list)] = paste0('cluster_', i)
  }
}

#Lots of cell death pathways
pdf('~/Documents/ÖverbyLab/scPlots/galectin3_proj/microglia/inf_mic_clust0_barplot.pdf', width = 10, height = 8)
ggplot(head(inf_micro_path_list$cluster_0, n = 10), aes(x = -log10(p_value), y = reorder(term_name, -(p_value))))+
  geom_bar(stat = 'identity', color = 'black', alpha = 0.5)+
  theme_classic()+
  ylab('')+
  theme(text = element_text(size = 22))+
  ggtitle('Cluster 0')
dev.off()

#Cell cycle pathways
pdf('~/Documents/ÖverbyLab/scPlots/galectin3_proj/microglia/inf_mic_clust2_barplot.pdf', width = 10, height = 8)
ggplot(head(inf_micro_path_list$cluster_2, n = 10), aes(x = -log10(p_value), y = reorder(term_name, -(p_value))))+
  geom_bar(stat = 'identity',  color = 'black', alpha = 0.5)+
  theme_classic()+
  ylab('')+
  theme(text = element_text(size = 18))+
  ggtitle('Cluster 2')
dev.off()

#Immune system response
pdf('~/Documents/ÖverbyLab/scPlots/galectin3_proj/microglia/inf_mic_clust3_barplot.pdf', width = 10, height = 8)
ggplot(head(inf_micro_path_list$cluster_3, n = 10), aes(x = -log10(p_value), y = reorder(term_name, -(p_value))))+
  geom_bar(stat = 'identity', color = 'black', alpha = 0.5)+
  theme_classic()+
  ylab('')+
  theme(text = element_text(size = 18))+
  ggtitle('Cluster 3')
dev.off()

#Also immune system but some negative regulators?
pdf('~/Documents/ÖverbyLab/scPlots/galectin3_proj/microglia/inf_mic_clust4_barplot.pdf', width = 10, height = 8)
ggplot(head(inf_micro_path_list$cluster_4, n = 10), aes(x = -log10(p_value), y = reorder(term_name, -(p_value))))+
  geom_bar(stat = 'identity', color = 'black', alpha = 0.5)+
  theme_classic()+
  ylab('')+
  theme(text = element_text(size = 18))+
  ggtitle('Cluster 4')
dev.off()


FeaturePlot(microglia_infected, features = 'Vcam1', reduction = 'micro.inf.umap')


#Compare some of the specific groups to narrow down roles
cluster_0_3 <- FindMarkers(microglia_infected, group.by = 'seurat_clusters', ident.1 = '0', ident.2 = '3',
                                         test.use = 'MAST')
cluster_0_markers <- dplyr::filter(cluster_0_3, avg_log2FC > 1 & p_val_adj < 0.01)
cluster_3_markers <- dplyr::filter(cluster_0_3, avg_log2FC < -1 & p_val_adj < 0.01)

FeaturePlot(microglia_infected, features = 'Lilrb4a', reduction = 'micro.inf.umap')

clust0_paths <- gprofiler2::gost(query = rownames(cluster_0_markers),organism = 'mmusculus', evcodes = TRUE,
                 sources = c('GO:BP', 'KEGG', 'GO:CC', 'GO:MP'))
clust3_paths <- gprofiler2::gost(query = rownames(cluster_3_markers),organism = 'mmusculus', evcodes = TRUE,
                                 sources = c('GO:BP', 'KEGG', 'GO:CC', 'GO:MP'))
ggplot(head(clust0_paths$result, n = 8), aes(x = -log10(p_value), y = reorder(term_name, -p_value)))+
  geom_bar(stat = 'identity')+
  theme_classic()+
  theme(text = element_text(size = 15))+
  ylab('')

ggplot(head(clust3_paths$result, n = 8), aes(x = -log10(p_value), y = reorder(term_name, -p_value)))+
  geom_bar(stat = 'identity')+
  theme_classic()+
  theme(text = element_text(size = 15))+
  ylab('')

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#### DEGs between each infected cluster and mock #####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

infected_clusters <- microglia_infected[[]]['seurat_clusters']
infected_clusters$cell_id <- rownames(infected_clusters)
colnames(infected_clusters)[1] = 'infected_clusters'

wt_cerebrum_microglia$cell_id <- colnames(wt_cerebrum_microglia)
wt_cerebrum_microglia$infected_clusters = NULL #In case running several times
wt_cerebrum_microglia[[]] <- left_join(wt_cerebrum_microglia[[]], infected_clusters, by = 'cell_id')
wt_cerebrum_microglia$infected_clusters <- factor(wt_cerebrum_microglia$infected_clusters, 
                                                  levels = c(levels(wt_cerebrum_microglia$infected_clusters), 'mock'))
wt_cerebrum_microglia[[]][is.na(wt_cerebrum_microglia$infected_clusters),]$infected_clusters = 'mock'

#Should have 0 through 4 and mock only (if clusters labelled 0-6 check that microglia_infected was run through seuart umap)
DimPlot(wt_cerebrum_microglia, reduction = 'micro.umap', label = FALSE, group.by = 'infected_clusters')+
  ggtitle('Microglia')+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank())+
  ylab('UMAP2')+
  xlab('UMAP1')

#Show proportions of timepoints by each cluster
table(wt_cerebrum_microglia$seurat_clusters, wt_cerebrum_microglia$Timepoint) %>% 
  as.data.frame() %>% dplyr::group_by(Var1) %>% dplyr::mutate(freq_props = Freq/sum(Freq))%>% 
  dplyr::mutate(Var1 = factor(Var1, levels = c('mock', '0', '1', '2', '3', '4'))) %>% 
  ggplot(aes(x = Var1, y = freq_props, fill = Var2))+
  geom_bar(stat = 'identity', position = 'stack', width = 0.6)+
  theme_classic()+
  theme(axis.text.x = element_text(size = 22),
        axis.text.y = element_text(size = 22),
        axis.title=element_text(size=22),
        plot.title = element_text(size = 22))+
  xlab('')+
  ggtitle('')+
  ylab('Proportion of cells')

table(wt_cerebrum_microglia$seurat_clusters, wt_cerebrum_microglia$Timepoint) %>% 
  as.data.frame() %>% dplyr::group_by(Var1) %>% 
  dplyr::mutate(Var1 = factor(Var1, levels = c('mock', '0', '1', '2', '3', '4'))) %>% 
  ggplot(aes(x = Var1, y = Freq, fill = Var2))+
  geom_bar(stat = 'identity', position = 'stack', width = 0.6)+
  theme_classic()+
  theme(axis.text.x = element_text(size = 22),
        axis.text.y = element_text(size = 22),
        axis.title=element_text(size=22),
        plot.title = element_text(size = 22))+
  xlab('')+
  ylab('Proportion of cells')

#Run after subsetting infected microglia but before finding any new clusters
#Show proportions of clusters across each timepoint

pdf('~/Documents/ÖverbyLab/scPlots/galectin3_proj/microglia/infected_micro_time_clusters.pdf', width = 6, height = 6)
table(microglia_infected$Timepoint, microglia_infected$seurat_clusters) %>% 
  as.data.frame() %>% dplyr::group_by(Var1) %>% dplyr::mutate(freq_props = Freq/sum(Freq))%>% 
  ggplot(aes(x = Var1, y = freq_props, fill = Var2))+
  geom_bar(stat = 'identity', position = 'stack', width = 0.6)+
  theme_classic()+
  theme(axis.text.x = element_text(size = 22),
        axis.text.y = element_text(size = 22),
        axis.title=element_text(size=22),
        plot.title = element_text(size = 22))+
  xlab('')+
  ggtitle('')+
  ylab('Proportion of cells')
dev.off()

pdf('~/Documents/ÖverbyLab/scPlots/galectin3_proj/microglia/mock_micro_time_clusters.pdf', width = 5, height = 6)
table(microglia_mock$Timepoint, microglia_mock$seurat_clusters) %>% 
  as.data.frame() %>% dplyr::group_by(Var1) %>% dplyr::mutate(freq_props = Freq/sum(Freq))%>% 
  ggplot(aes(x = Var1, y = freq_props, fill = Var2))+
  geom_bar(stat = 'identity', position = 'stack', width = 0.6)+
  theme_classic()+
  theme(axis.text.x = element_text(size = 22),
        axis.text.y = element_text(size = 22),
        axis.title=element_text(size=22),
        plot.title = element_text(size = 22))+
  xlab('')+
  ggtitle('')+
  ylab('Proportion of cells')
dev.off()

table(microglia_infected$Timepoint, microglia_infected$seurat_clusters) %>% 
  as.data.frame() %>% dplyr::group_by(Var1) %>% 
  ggplot(aes(x = Var1, y = Freq, fill = Var2))+
  geom_bar(stat = 'identity', position = 'stack', width = 0.6)+
  theme_classic()+
  theme(axis.text.x = element_text(size = 22),
        axis.text.y = element_text(size = 22),
        axis.title=element_text(size=22),
        plot.title = element_text(size = 22))+
  xlab('')+
  ylab('Number of cells')

#Find infected markers vs mock
deg_vs_mock_list = list()

for(i in levels(wt_cerebrum_microglia$infected_clusters)){
  if(i != 'mock'){
    print(paste('starting', i))
    clust_markers <- FindMarkers(wt_cerebrum_microglia, group.by = 'infected_clusters', ident.1 = i , ident.2 = 'mock', test.use = 'MAST',
                only.pos = TRUE)
    deg_vs_mock_list[[length(deg_vs_mock_list) + 1]] = clust_markers
  }
}

#Name list elements 
names(deg_vs_mock_list) = paste('cluster',  levels(wt_cerebrum_microglia$infected_clusters)[1:5])

num_sig_genes <- lapply(deg_vs_mock_list, FUN = function(x){
  sig_x <- dplyr::filter(x, p_val_adj < 0.01 & avg_log2FC > 1)
  nrow(sig_x)
})

as.data.frame(num_sig_genes) %>% t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'group') %>% 
  ggplot(aes(x = group, y = V1))+
  geom_bar(stat = 'identity', , color = 'black', fill ='lightgrey')+
  ylab('Num Upregulated DEGs')+
  xlab('Cluster')+
  theme_classic()

#Look at pathways
path_vs_mock_list = list()
for(i in 1:length(deg_vs_mock_list)){
  cur_paths <- deg_vs_mock_list[[i]] %>% as.data.frame() %>% 
    dplyr::filter(p_val_adj < 0.01 & avg_log2FC > 1) %>% 
    rownames() %>% gprofiler2::gost(organism = 'mmusculus', evcodes = TRUE,
                                    sources = c('GO:BP', 'KEGG', 'GO:CC', 'GO:MP'))
  path_vs_mock_list[[length(path_vs_mock_list) + 1]] = cur_paths$result
}
names(path_vs_mock_list) = names(deg_vs_mock_list)

pdf('~/Documents/ÖverbyLab/scPlots/galectin3_proj/microglia/inf_vs_mock_clust_0.pdf', width = 10, height = 8)
ggplot(head(path_vs_mock_list$`cluster 0`, n = 10), aes(x = -log10(p_value), y = reorder(term_name, -(p_value))))+
  geom_bar(stat = 'identity', color = 'black')+
  theme_classic()+
  ylab('')+
  xlab('-Log10 p-val')+
  theme(text = element_text(size = 18))+
  ggtitle('Cluster 0')
dev.off()

for(i in 1:length(path_vs_mock_list)){
  path_bar <- ggplot(head(path_vs_mock_list[[i]], n = 10), aes(x = -log10(p_value), y = reorder(term_name, -(p_value))))+
    geom_bar(stat = 'identity', color = 'black')+
    theme_classic()+
    ylab('')+
    xlab('-Log10 p-val')+
    theme(text = element_text(size = 22))+
    ggtitle(names(path_vs_mock_list)[i])
  
  clust = as.character(i-1)
  file_path = paste0('~/Documents/ÖverbyLab/scPlots/galectin3_proj/microglia/inf_vs_mock_clust_', clust, '.pdf')
  pdf(file_path, width = 11, height = 8)
  print(path_bar)
  dev.off()
}

#Bar plots with grouped clusters
path_vs_mock_list

paths_to_plot <- lapply(path_vs_mock_list, FUN = function(x){
  relevant_paths <- x[x$term_name %in% c('response to stress', 'innate immune response', 
                                         'defense response to other organism', 'immune response',
                                         'response to biotic stimulus', 'defense_response'),]
  relevant_paths[,c('term_name', 'p_value')]
})

do.call(rbind, paths_to_plot) %>% 
  tibble::rownames_to_column(var = 'id') %>% 
  dplyr::mutate(id = substr(id, 1, 9)) %>% 
  dplyr::mutate(id = factor(id, levels = c(paste0('cluster ', 4:0)))) %>% 
  ggplot(aes(x = -log10(p_value), y = term_name, fill = id))+
  geom_bar(stat = 'identity', position = 'dodge', color = 'black')+
  theme_classic()+
  scale_fill_manual(values=c("#7047A1",  "#6DC3F8", "#FF96A2", 
                             "#D6644B", "#F08C3A", "#fdc087","#074F00"))+
  theme(axis.text = element_text(size = 16))+
  guides(fill = guide_legend(reverse = TRUE))

# Look at genes differntially expressed between infected clusters,
#as well as between the given cluster and mock.

sig_vs_all <- list()
for(i in 0:4){
  clust_vs_infected <- dplyr::filter(infected_clust_markers, cluster == i & avg_log2FC > 1 & p_val_adj < 0.01)
  clust_vs_mock <- dplyr::filter(deg_vs_mock_list[[i+1]], avg_log2FC > 1 & p_val_adj < 0.01)
  clust_sig <- clust_vs_infected[rownames(clust_vs_infected) %in% rownames(clust_vs_mock),]
  sig_vs_all[[length(sig_vs_all) + 1]] = clust_sig
}

top_sig_genes <- unlist(lapply(sig_vs_all, FUN = function(x){
  rownames(head(x, n = 3))
}))

top_sig_genes_dot <- DotPlot(wt_cerebrum_microglia, features = top_sig_genes, group.by = 'infected_clusters', scale = FALSE)$data
top_sig_genes_dot$id = factor(top_sig_genes_dot$id, levels = rev(c('mock', '0', '1', '2', '3', '4')))

ggplot(top_sig_genes_dot, aes(x = features.plot, y = id, fill = avg.exp.scaled, size = pct.exp))+
  geom_point(pch = 21)+
  scale_fill_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white"), 
                       values = c(1.0,0.7,0.4,0))+
  theme_classic()

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#### DEGs between each infected cluster and all cells including mock #####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

infected_vs_all_markers <- FindAllMarkers(wt_cerebrum_microglia, group.by = 'seurat_clusters', test.use = 'MAST',
            only.pos = TRUE)

infected_all_sig <- dplyr::filter(infected_vs_all_markers, p_val_adj < 0.01 & avg_log2FC > 1)

#Get upregulated pathways in all clusters
infected_vs_all_path_list <- list()

for(i in 1:length(unique(infected_all_sig$cluster))){
  print(paste('cluster', unique(infected_all_sig$cluster)[i]))
  cur_cluster = infected_all_sig[infected_all_sig$cluster == unique(infected_all_sig$cluster)[i],]
  cur_cluster_paths <- gprofiler2::gost(query = rownames(cur_cluster),organism = 'mmusculus', 
                                        evcodes = TRUE, sources = c('GO:BP', 'KEGG', 'GO:CC', 'GO:MP'))
  infected_vs_all_path_list[[length(infected_vs_all_path_list) + 1]] = cur_cluster_paths$result
}

#No significant pathways in cluster 1
names(infected_vs_all_path_list) = paste0('cluster_', c('0','1', '2', '3', '5', '6'))

for(i in 1:length(infected_vs_all_path_list)){
  
  cur_dat <- infected_vs_all_path_list[[i]]
  
  path_bar <- ggplot(head(cur_dat, n = 10), aes(x = -log10(p_value), y = reorder(term_name, -p_value)))+
    geom_bar(stat = 'identity')+
    theme_classic()+
    theme(text = element_text(size = 18))+
    ylab('')+
    ggtitle(names(infected_vs_all_path_list)[i])
  
  print(path_bar)
}

infected_vs_all_path_list$cluster_0
infected_vs_all_path_list$cluster_1
infected_vs_all_path_list$cluster_2
infected_vs_all_path_list$cluster_3
infected_vs_all_path_list$cluster_5
infected_vs_all_path_list$cluster_6

# - - - - - - - - - - - - 
#### Top deg heatmap #### 
# - - - - - - - - - - - - 

wt_cerebrum_microglia$clusters_with_mock = dplyr::case_when(wt_cerebrum_microglia$custom_clusters %in% c(0, 1, 2) ~ 'mock',
                                                            .default = as.character(wt_cerebrum_microglia$custom_clusters))
  

all_clust_markers <-  FindAllMarkers(wt_cerebrum_microglia, group.by = 'clusters_with_mock', 
                                     test.use = 'MAST', only.pos = TRUE)

# unique_cluster_markers <- all_clust_markers %>% dplyr::group_by(gene) %>% 
#   dplyr::summarise(total = n()) %>% 
#   dplyr::filter(total == 1)

#Making cluster a factor a ordering by that so mock can be next to cluster 0 on heatmap
top_inf_clust_genes <- all_clust_markers %>% 
  dplyr::group_by(cluster) %>% 
  dplyr::slice_head(n = 10) %>% 
  dplyr::mutate(cluster = factor(cluster, levels = c('mock', '3', '4', '5', '6'))) %>% 
  dplyr::arrange(cluster) %>% 
  dplyr::pull(gene)

top_gene_dot_dat <- DotPlot(wt_cerebrum_microglia, features = unique(top_inf_clust_genes), group.by = 'clusters_with_mock', scale = FALSE)$data
unique(top_inf_clust_genes)

gene_order <- c("Fcrls", "Abca9", "P3h2", "Snx29",  "Cx3cr1", "Nav2", "St6gal1", "Selenop", "Csmd3", "Plxdc2", "Diaph3", "Top2a", "Cep128", "Smc2", "Mki67",
                "Rad51b", "Ezh2", "Kif15", "Cit", "Smc4", "Pik3ap1", "Rnf213", "Ifit2", "Trim30a", "Ifi204", "Hsph1", "Herc6", "Ifi207", "Sp100", "Ddx60", "Junb",
                "Nf1", "C1qc", "Wsb1", "Thoc2", "Evl", "Mycbp2", "Atp6v0a2", "Tra2a", "Lrch3","Fmnl2", "Fth1", "Lilrb4a", "Ccl5", "Gm49339", "Il12b",  "Aoah", "Il1rn",
                "Cd40", "Acod1")

pdf('~/Documents/ÖverbyLab/scPlots/galectin3_proj/microglia/cluster_markers.pdf', height = 10, width = 6)
top_gene_dot_dat %>% 
  dplyr::mutate(features.plot = factor(features.plot, levels = rev(gene_order))) %>% 
  dplyr::mutate(id = factor(id, levels = c('mock', '3', '4', '5', '6'))) %>% 
  ggplot(aes(x = id, y = features.plot, fill = avg.exp.scaled)) +
  geom_tile()+
  theme_classic()+
  scale_fill_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white"), 
                       values = c(1.0,0.7,0.4,0))+
  ylab('')+
  xlab('cluster')+
  theme(axis.text = element_text(size = 16))
dev.off()


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#### Microglia groups from https://www.nature.com/articles/s41590-026-02472-z ####
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

surveilance <- c('Hexb', 'Olfml3', 'P2ry12', 'Siglech', 'Sparc', 'Tgfbr1', 'Tmem119')
neuro_protect <- c('Bdnf', 'Gdnf', 'Igf1', 'Lif', 'Tgfb1')
#No Clec7a in data so removed
phagocytosis <- c('Apoe', 'Axl', 'C3ar1', 'Cd9', 'Cd63', 'Cd68', 'Lgals3', 'Myo1e', 'Spp1', 'Trem2', 'Tyrobp')
inflammation <- c('Apod', 'Aoah', 'C5ar1', 'Ccl12', 'Fgr', 'Msr1', 'Rnf169')
cyto_production <- c('Ccl2', 'Cxcl10', 'Il10', 'Il1b', 'Nfkbia', 'Tnf')
antigen_pres <- c('B2m', 'H2-D1', 'H2-Aa', 'H2-Ab1', 'H2-DMa', 'H2-Eb1', 'H2-K1', 'Nlrc5')
ifn_sig <- c('Cd69', 'Cxcl10', 'Ifi204', 'Ifi213', 'Ifit2', 'Ifit3', 'Ifitm3', 'Irf7', 'Isg15', 'Mx1', 'Oasl2', 'Usp18')
proliferation <- c('Birc5', 'Brip1', 'Mcm5', 'Mki67', 'Rad51b', 'Top2a')

#Surveilance
plotList_surv <- lapply(surveilance, featurePlotLight, data = wt_cerebrum_microglia, reduction_choice = 'micro.umap', maxLim = 5.5)
do.call(ggarrange, c(plotList_surv, common.legend = TRUE, legend = 'right'))


wt_cerebrum_microglia <- AddModuleScore_UCell(wt_cerebrum_microglia, 
                                              features=list('surveilance' = surveilance), maxRank = 1200, name = NULL)

pdf('~/Documents/ÖverbyLab/scPlots/galectin3_proj/microglia/surveillance_violin.pdf', height = 5, width = 7)
VlnPlot(object = wt_cerebrum_microglia, features = 'surveilance', group.by = 'custom_clusters', pt.size = 0)
dev.off()

pdf('~/Documents/ÖverbyLab/scPlots/galectin3_proj/microglia/surveillance_dotplot.pdf', height = 5, width = 7)
DotPlot(wt_cerebrum_microglia, features = surveilance, group.by = 'custom_clusters', scale = FALSE)$data %>% 
  ggplot(aes(x = id, y = features.plot, size = pct.exp, fill = avg.exp.scaled))+
  geom_point(pch = 21)+
  theme_classic()+
  scale_fill_gradientn(colours = c('white', '#FFD991', '#FF7530', '#FF4024'), 
                       values = c(0, 0.3, 0.6, 1))+
  ggtitle('Surveillance genes')
dev.off()

#Neuro protection
plotList_np <- lapply(neuro_protect, featurePlotLight, data = wt_cerebrum_microglia, reduction_choice = 'micro.umap', maxLim = 4.5)
do.call(ggarrange, c(plotList_np, common.legend = TRUE, legend = 'right'))

wt_cerebrum_microglia <- AddModuleScore_UCell(wt_cerebrum_microglia, 
                                              features=list('neuro_protection' = neuro_protect), maxRank = 1200, name = NULL)

pdf('~/Documents/ÖverbyLab/scPlots/galectin3_proj/microglia/neuro_protection_violin.pdf', height = 5, width = 7)
VlnPlot(object = wt_cerebrum_microglia, features = 'neuro_protection', group.by = 'custom_clusters', pt.size = 0)
dev.off()

#Phagotcytosis
plotList_phago <- lapply(phagocytosis, featurePlotLight, data = wt_cerebrum_microglia, reduction_choice = 'micro.umap', maxLim = 6)
do.call(ggarrange, c(plotList_phago, common.legend = TRUE, legend = 'right'))

wt_cerebrum_microglia <- AddModuleScore_UCell(wt_cerebrum_microglia, 
                                              features=list('phagocytosis' = phagocytosis), maxRank = 1200, name = NULL)

pdf('~/Documents/ÖverbyLab/scPlots/galectin3_proj/microglia/phagocytosis_violin.pdf', height = 5, width = 7)
VlnPlot(object = wt_cerebrum_microglia, features = 'phagocytosis', group.by = 'custom_clusters', pt.size = 0)
dev.off()

#Inflammation
plotList_inf <- lapply(inflammation, featurePlotLight, data = wt_cerebrum_microglia, reduction_choice = 'micro.umap', maxLim = 4.7)
do.call(ggarrange, c(plotList_inf, common.legend = TRUE, legend = 'right'))

#cytokine production
plotList_cyto <- lapply(cyto_production, featurePlotLight, data = wt_cerebrum_microglia, reduction_choice = 'micro.umap', maxLim = 4.7)
do.call(ggarrange, c(plotList_cyto, common.legend = TRUE, legend = 'right'))

wt_cerebrum_microglia <- AddModuleScore_UCell(wt_cerebrum_microglia, 
                                              features=list('cytokine_production' = cyto_production), maxRank = 1200, name = NULL)

pdf('~/Documents/ÖverbyLab/scPlots/galectin3_proj/microglia/cytokine_production_violin.pdf', height = 5, width = 7)
VlnPlot(object = wt_cerebrum_microglia, features = 'cytokine_production', group.by = 'custom_clusters', pt.size = 0)
dev.off()

#antigen presentation
plotList_antigen <- lapply(antigen_pres, featurePlotLight, data = wt_cerebrum_microglia, reduction_choice = 'micro.umap', maxLim = 5.3)
do.call(ggarrange, c(plotList_antigen, common.legend = TRUE, legend = 'right'))

wt_cerebrum_microglia <- AddModuleScore_UCell(wt_cerebrum_microglia, 
                                              features=list('antigen_presentation' = antigen_pres), maxRank = 1200, name = NULL)

pdf('~/Documents/ÖverbyLab/scPlots/galectin3_proj/microglia/antigen_presentation_violin.pdf', height = 5, width = 7)
VlnPlot(object = wt_cerebrum_microglia, features = 'antigen_presentation', group.by = 'custom_clusters', pt.size = 0)
dev.off()

pdf('~/Documents/ÖverbyLab/scPlots/galectin3_proj/microglia/antigen_presentation_dotplot.pdf', height = 5, width = 7)
DotPlot(wt_cerebrum_microglia, features = antigen_pres, group.by = 'custom_clusters', scale = FALSE)$data %>% 
  ggplot(aes(x = id, y = features.plot, size = pct.exp, fill = avg.exp.scaled))+
  geom_point(pch = 21)+
  theme_classic()+
  scale_fill_gradientn(colours = c('white', '#FFD991', '#FF7530', '#FF4024'), 
                       values = c(0, 0.3, 0.6, 1))+
  ggtitle('Antigen presentation genes')
dev.off()

#ifn sig
plotList_ifn <- lapply(ifn_sig, featurePlotLight, data = wt_cerebrum_microglia, reduction_choice = 'micro.umap', maxLim = 5.5)
do.call(ggarrange, c(plotList_ifn, common.legend = TRUE, legend = 'right'))

#proliferation
plotList_prolif <- lapply(proliferation, featurePlotLight, data = microglia_infected, reduction_choice = 'micro.inf.umap', maxLim = 4.5)
pdf('~/Documents/ÖverbyLab/scPlots/galectin3_proj/microglia/proliferation_genes.pdf', width = 6, height = 5)
do.call(ggarrange, c(plotList_prolif, common.legend = TRUE, legend = 'right'))
dev.off()

wt_cerebrum_microglia <- AddModuleScore_UCell(wt_cerebrum_microglia, 
                                              features=list('proliferation' = proliferation), maxRank = 1200, name = NULL)

pdf('~/Documents/ÖverbyLab/scPlots/galectin3_proj/microglia/proliferation_production_violin.pdf', height = 5, width = 7)
VlnPlot(object = wt_cerebrum_microglia, features = 'proliferation', group.by = 'custom_clusters', pt.size = 0)
dev.off()

#Function for plotting multiple featureplots together on one scale
featurePlotLight <- function(gene, data, reduction_choice, scale = FALSE, minLim = 0, maxLim = 5){
  dat = FeaturePlot(data, gene, reduction = reduction_choice)$data
  colnames(dat) = c('umap1', 'umap2', 'ident', 'expression')
  ggplot(dat, aes(x = umap1, y = umap2, color = expression))+
    geom_point(size = 0.1)+  
    theme(line = element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_blank(),
          panel.background = element_rect(fill = '#F2F2F2', color = '#F2F2F2'))+
    scale_color_gradient(low = 'lightgrey', high = 'blue', limits = c(minLim,maxLim))+
    ggtitle(gene)
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#### DAMs https://www.cell.com/cell/fulltext/S0092-8674(17)30578-0#mmc1 ####
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#Downregulated in dams
featurePlotLight('P2ry12', data = wt_cerebrum_microglia, reduction_choice = 'micro.umap', maxLim = NA)
featurePlotLight('Tmem119', data = wt_cerebrum_microglia, reduction_choice = 'micro.umap', maxLim = NA)
featurePlotLight('Selplg', data = wt_cerebrum_microglia, reduction_choice = 'micro.umap', maxLim = NA)

#Upregulated in dams
featurePlotLight('Cst7', data = wt_cerebrum_microglia, reduction_choice = 'micro.umap', maxLim = NA)
featurePlotLight('Lpl', data = wt_cerebrum_microglia, reduction_choice = 'micro.umap', maxLim = NA)
featurePlotLight('Apoe', data = wt_cerebrum_microglia, reduction_choice = 'micro.umap', maxLim = NA)
featurePlotLight('Spp1', data = wt_cerebrum_microglia, reduction_choice = 'micro.umap', maxLim = NA)

#Just infected
featurePlotLight('P2ry12', data = microglia_infected, reduction_choice = 'micro.inf.umap', maxLim = NA)
featurePlotLight('Tmem119', data = microglia_infected, reduction_choice = 'micro.inf.umap', maxLim = NA)
featurePlotLight('Cst7', data = microglia_infected, reduction_choice = 'micro.inf.umap', maxLim = NA)
featurePlotLight('Itgax', data = microglia_infected, reduction_choice = 'micro.inf.umap', maxLim = NA)
featurePlotLight('Hexb', data = microglia_infected, reduction_choice = 'micro.inf.umap', maxLim = NA)

#Top dam genes (from supp table 2 in paper, may change to supp table 3)
dam_signatures <- list(dam_down = c('P2ry12','Tmem119','Cx3cr1','Selplg','Serinc3', 'Marcks',
                                     'Malat1','Glul','Lgmn','4632428N05Rik','Txnip','Rhob',
                                     'Sparc','Slco2b1','Csf1r'),
                        dam_up = c('Cst7','Lpl', 'Clec7a', 'Itgax', 'Spp1', 'Igf1',
                                   'Apoe','Axl','Ank','Ch25h', 'Ctsd', 'Ccl3','Baiap2l2','Csf1'))

dam_signatures_supp_3 <- list(dam_up = c('Cst7', 'Lyz2', 'Lpl', 'Ctsb', 'Ctsd', 'Apoe', 'B2m', 'Gnas',
                                         'Cd9', 'Ank', 'H2-D1', 'Tyrobp', 'Axl', 'Ctsz','Itgax'),
                              dom_down = c('P2ry12', 'Tmem119', 'Cx3cr1', 'Selplg', 'Marcks', 'Serinc3', 'Sparc',
                                           'Cd164', 'Maf', 'P2ry13','Ccr5','Glul','Rhob', 'Siglech','Olfml3'))
wt_cerebrum_microglia <- AddModuleScore_UCell(wt_cerebrum_microglia, 
                                      features=dam_signatures, maxRank = 1200, name = NULL)

plotList_dam <- lapply(c('dam_down', 'dam_up'), featurePlotLight, 
                       data = wt_cerebrum_microglia, reduction_choice = 'micro.umap', maxLim = 1)
do.call(ggarrange, c(plotList_dam, common.legend = TRUE, legend = 'right'))

VlnPlot(object = wt_cerebrum_microglia, features = 'dam_up', group.by = 'seurat_clusters', pt.size = 0)
VlnPlot(object = wt_cerebrum_microglia, features = 'dam_down', group.by = 'seurat_clusters', pt.size = 0)

#heatmap instead
dam_dot_up <- DotPlot(wt_cerebrum_microglia, features = dam_signatures_supp_3$dam_up, group.by = 'infected_clusters', scale = FALSE)$data
dam_dot_up$id = factor(dam_dot_up$id, levels = c('mock', '0', '1', '2', '3', '4'))

#Set levels based on which data from the paper is being used
feature_levels_supp2 <- rev(c('Itgax',  'Ch25h','Ccl3', 'Cst7','Axl','Csf1', 'Ctsd', 'Lpl', 
                              'Spp1', 'Igf1', 'Apoe', 'Ank', 'Baiap2l2'))
feature_levels_supp3 <- rev(c('Itgax', 'H2-D1', 'Ctsz','Axl', 'Cst7', 'Ctsb', 'Ctsd', 'B2m', 'Gnas',
                              'Cd9', 'Ank', 'Tyrobp', 'Lyz2', 'Lpl', 'Apoe'))

dam_dot_up$features.plot = factor(dam_dot_up$features.plot, 
                                  levels = feature_levels_supp3)
ggplot(dam_dot_up, aes(x = id, y = features.plot, fill = avg.exp.scaled, size = pct.exp))+
  geom_point(pch = 21)+
  theme_classic()+
  scale_fill_gradientn(colours = c('white', '#FFD991', '#FF7530', '#FF4024'), 
                       values = c(0, 0.3, 0.6, 1))+
  ggtitle('DAM Up genes')

#Dam down
dam_dot_down <- DotPlot(wt_cerebrum_microglia, features = dam_signatures_supp_3$dom_down, group.by = 'infected_clusters', scale = FALSE)$data
dam_dot_down$id = factor(dam_dot_down$id, levels = c('mock', '0', '1', '2', '3', '4'))
dam_dot_down$features.plot = factor(dam_dot_down$features.plot, 
                                  levels = rev(c('P2ry12', 'Tmem119', 'Cx3cr1','Glul', 'Siglech','Selplg', 'Marcks', 'Serinc3', 'Sparc',
                                                 'Cd164', 'Maf', 'P2ry13','Ccr5','Rhob', 'Olfml3')))
ggplot(dam_dot_down, aes(x = id, y = features.plot, fill = avg.exp.scaled, size = pct.exp))+
  geom_point(pch = 21)+
  theme_classic()+
  scale_fill_gradientn(colours = c('white', '#FFD991', '#FF7530', '#FF4024'), 
                       values = c(0, 0.3, 0.6, 1))+
  ggtitle('DAM down genes')

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#### Microglia groups from https://link.springer.com/article/10.15252/embr.201846171#Sec2 ####
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

#Homeostatic genes
featurePlotLight('Olfml3', data = wt_cerebrum_microglia, reduction_choice = 'micro.umap', maxLim = NA)
featurePlotLight('Fcrls', data = wt_cerebrum_microglia, reduction_choice = 'micro.umap', maxLim = NA)
featurePlotLight('Tmem119', data = wt_cerebrum_microglia, reduction_choice = 'micro.umap', maxLim = NA)
featurePlotLight('Gpr34', data = wt_cerebrum_microglia, reduction_choice = 'micro.umap', maxLim = NA)
featurePlotLight('Mef2c', data = wt_cerebrum_microglia, reduction_choice = 'micro.umap', maxLim = NA)

homeo_markers <- DotPlot(wt_cerebrum_microglia, features = c('Olfml3', 'Fcrls', 'Tmem119', 'Gpr34', 'Mef2c'), group.by = 'seurat_clusters', scale = FALSE)$data

ggplot(homeo_markers, aes(x = id, y = features.plot, fill = avg.exp.scaled, size = pct.exp))+
  geom_point(pch = 21)+
  theme_classic()+
  scale_fill_gradientn(colours = c('white', '#FFD991', '#FF7530', '#FF4024'), 
                       values = c(0, 0.3, 0.6, 1))+
  ggtitle('Homeostatic genes')

#Phagocytotic genes
featurePlotLight('Tyrobp', data = wt_cerebrum_microglia, reduction_choice = 'micro.umap', maxLim = NA)
featurePlotLight('Trem2', data = wt_cerebrum_microglia, reduction_choice = 'micro.umap', maxLim = NA)

#Anti inflammatory
featurePlotLight('Mrc1', data = wt_cerebrum_microglia, reduction_choice = 'micro.umap', maxLim = NA)
featurePlotLight('Arg1', data = wt_cerebrum_microglia, reduction_choice = 'micro.umap', maxLim = NA)

#Pro inflammatory
featurePlotLight('Il1b', data = wt_cerebrum_microglia, reduction_choice = 'micro.umap', maxLim = NA)
featurePlotLight('Tnf', data = wt_cerebrum_microglia, reduction_choice = 'micro.umap', maxLim = NA)
featurePlotLight('Ccl2', data = wt_cerebrum_microglia, reduction_choice = 'micro.umap', maxLim = NA)

DotPlot(wt_cerebrum_microglia, features =c('Il1b', 'Tnf', 'Ccl2'), group.by = 'infected_clusters', scale = FALSE)$data %>% 
  dplyr::mutate(id = factor(id, levels = c('mock', '0', '1', '2', '3', '4'))) %>% 
  ggplot(aes(x = id, y = features.plot, size = pct.exp, fill = avg.exp.scaled))+
  geom_point(pch = 21)+
  scale_fill_gradientn(colours = c('white', '#FFD991', '#FF7530', '#FF4024'), 
                       values = c(0, 0.3, 0.6, 1))+
  theme_classic()+
  ylab('')+
  xlab('')+
  ggtitle('pro-inflammatory')

#LPS inflammatory markers
lps_markers <- c('Ms4a6c', 'Msr1', 'Igsf6', 'Ms4a6b', 'Ms4a6d','Srgn', 'Ccl12',
                 'Slfn2', 'D17H6S56E-5','Gpr84','Saa3', 'Ctsc','Pilra', 'Ifi204',
                 'Nfkbia', 'mt-Rnr1','Ifitm3', 'Rps2','AW112010', 'C3ar1', 'Cpd', 'Cd52')
lps_dot <- DotPlot(wt_cerebrum_microglia, features = lps_markers, group.by = 'seurat_clusters', scale = FALSE)$data
#lps_dot$id = factor(lps_dot$id, levels = c())

ggplot(lps_dot, aes(x = id, y = features.plot, fill = avg.exp.scaled, size = pct.exp))+
  geom_point(pch = 21)+
  theme_classic()+
  scale_fill_gradientn(colours = c('white', '#FFD991', '#FF7530', '#FF4024'), 
                       values = c(0, 0.3, 0.6, 1))+
  ggtitle('IAM (inflammatory markers)')

wt_cerebrum_microglia <- AddModuleScore_UCell(wt_cerebrum_microglia, 
                                              features=list('inflammatory_score' = lps_markers), maxRank = 1200, name = NULL)

pdf('~/Documents/ÖverbyLab/scPlots/galectin3_proj/microglia/inflammatory_score_violin.pdf', height = 5, width = 7)
VlnPlot(object = wt_cerebrum_microglia, features = 'inflammatory_score', group.by = 'custom_clusters', pt.size = 0)
dev.off()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#### Microglia groups from https://www.nature.com/articles/s41586-019-0924-x ####
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#Gpnmb also mentioned here as phagocytotic https://www.nature.com/articles/s41467-026-73003-5
DotPlot(wt_cerebrum_microglia, features = 'Gpnmb', group.by = 'infected_clusters', scale = FALSE)$data %>% 
  dplyr::mutate(id = factor(id, levels = c('mock', '0', '1', '2', '3', '4'))) %>% 
  ggplot(aes(x = id, y = features.plot, fill = avg.exp.scaled, size = pct.exp))+
  geom_point(pch = 21)+
  theme_classic()+
  scale_fill_gradientn(colours = c('white', '#FFD991', '#FF7530', '#FF4024'), 
                       values = c(0, 0.3, 0.6, 1))+
  ggtitle('Gpnmb')

DotPlot(wt_cerebrum_microglia, features = 'Cybb', group.by = 'infected_clusters', scale = FALSE)$data %>% 
  dplyr::mutate(id = factor(id, levels = c('mock', '0', '1', '2', '3', '4'))) %>% 
  ggplot(aes(x = id, y = features.plot, fill = avg.exp.scaled, size = pct.exp))+
  geom_point(pch = 21)+
  theme_classic()+
  scale_fill_gradientn(colours = c('white', '#FFD991', '#FF7530', '#FF4024'), 
                       values = c(0, 0.3, 0.6, 1))+
  ggtitle('Cybb')

# - - - -- - - - - - - - - - - - - -
#### Hallmark apoptosis genes ####
# - - - - - - - - - - - - - - - - - - 
gene_sets <- msigdbr(species = "mouse", db_species = "MM")
apoptosis_set <- gene_sets[gene_sets$gs_name == 'HALLMARK_APOPTOSIS',]$gene_symbol

wt_cerebrum_microglia <- AddModuleScore_UCell(wt_cerebrum_microglia, 
                                              features= list(apoptosis = apoptosis_set), maxRank = 1200, name = NULL)

FeaturePlot(wt_cerebrum_microglia, features = 'apoptosis', reduction = 'micro.umap')+
  ggtitle('Hallmark apoptosis')

VlnPlot(wt_cerebrum_microglia, features = 'apoptosis', pt.size = 0, group.by = 'infected_clusters')+
  ggtitle('Hallmark apoptosis')

# - - - - - - - - - - - - - - - - - - - - - - 
#### Random genes of interest to infection ####
# - - - - - - - - - - - - - - - - - - - - - - 
#chemokines
chemokine_markers <- c('Ccl2', 'Ccl5', 'Ccl12', 'Cxcl10')
chemo_dot <- DotPlot(wt_cerebrum_microglia, features = chemokine_markers, group.by = 'infected_clusters', scale = FALSE)$data
chemo_dot$id = factor(chemo_dot$id, levels = c('mock', '0', '1', '2', '3', '4'))

ggplot(chemo_dot, aes(x = id, y = features.plot, fill = avg.exp.scaled, size = pct.exp))+
  geom_point(pch = 21)+
  theme_classic()+
  scale_fill_gradientn(colours = c('white', '#FFD991', '#FF7530', '#FF4024'), 
                       values = c(0, 0.3, 0.6, 1))+
  ggtitle('Chemokines')

FeaturePlot(wt_cerebrum_microglia, features = 'Ccl12', reduction = 'micro.umap')+
  ggtitle('Ccl12')

#mhc 2 proteins
mhc2 = c("H2-Aa",  "H2-Ab1",  "H2-DMa",  "H2-DMb1", "H2-DMb2" ,"H2-Eb1" )

mhc2_dot <- DotPlot(wt_cerebrum_microglia, features = mhc2, group.by = 'treatment_time', scale = FALSE)$data
mhc2_dot$id = factor(mhc2_dot$id, levels = c('mock', '0', '1', '2', '3', '4'))

ggplot(mhc2_dot, aes(x = id, y = features.plot, fill = avg.exp.scaled, size = pct.exp))+
  geom_point(pch = 21)+
  theme_classic()+
  scale_fill_gradientn(colours = c('white', '#FFD991', '#FF7530', '#FF4024'), 
                       values = c(0, 0.3, 0.6, 1))+
  ggtitle('mhc2')

FeaturePlot(wt_cerebrum_microglia, features = 'H2-D1', reduction = 'micro.umap')+
  ggtitle('H2-D1')

#ISGs (some from https://link.springer.com/article/10.1186/s12974-025-03471-x#Sec27)
isg_list = c('Tnf', 'Il1b', 'Isg15')

isg_dot <- DotPlot(wt_cerebrum_microglia, features = isg_list, group.by = 'infected_clusters', scale = FALSE)$data
isg_dot$id = factor(isg_dot$id, levels = c('mock', '0', '1', '2', '3', '4'))

ggplot(isg_dot, aes(x = id, y = features.plot, fill = avg.exp.scaled, size = pct.exp))+
  geom_point(pch = 21)+
  theme_classic()+
  scale_fill_gradientn(colours = c('white', '#FFD991', '#FF7530', '#FF4024'), 
                       values = c(0, 0.3, 0.6, 1))+
  ggtitle('isgs')

FeaturePlot(wt_cerebrum_microglia, features = 'Tnf', reduction = 'micro.umap')+
  ggtitle('Tnf')

#F480 plot
pdf('~/Documents/ÖverbyLab/scPlots/galectin3_proj/microglia/f480_feature_plot.pdf', width = 6, height = 5)
FeaturePlot(wt_cerebrum_microglia, features = 'Adgre1', reduction = 'micro.umap')+
  ggtitle('F4/80')+
  ylab('')+
  xlab('')
dev.off()

pdf('~/Documents/ÖverbyLab/scPlots/galectin3_proj/microglia/f480_dotplot.pdf', width = 6, height = 5)
DotPlot(wt_cerebrum_microglia, features = 'Adgre1', group.by = 'clusters_with_mock', scale = FALSE)$data %>% 
  dplyr::mutate(id = factor(id, levels = c('mock', '3', '4', '5', '6'))) %>% 
  ggplot(aes(x = id, y = features.plot, fill = avg.exp.scaled, size = pct.exp))+
  geom_point(pch = 21)+
  theme_classic()+
  scale_fill_gradientn(colours = c('white', '#FFD991', '#FF7530', '#FF4024'), 
                       values = c(0, 0.3, 0.6, 1))+
  ggtitle('F4/80')
dev.off()
#Microglia cluster 14 score
microglia_clust_14 <- dplyr::filter(hsv_gene_list, cluster == 14 & celltype == 'Microglia') %>% dplyr::pull(gene)

wt_cerebrum_microglia <- AddModuleScore_UCell(wt_cerebrum_microglia, features = list(clust14 = microglia_clust_14), name = 'hsv_clust_14_score')
wt_cerebrum_microglia$treatment_time <- paste(wt_cerebrum_microglia$Treatment, wt_cerebrum_microglia$Timepoint, sep = '_')

DotPlot(wt_cerebrum_microglia, features = 'clust14hsv_clust_14_score', group.by = 'treatment_time', scale = FALSE)$data %>% 
  ggplot(aes(x = id, y = features.plot, fill = avg.exp.scaled))+
  geom_tile()+
  theme_classic()+
  scale_fill_gradientn(colours = c('white', '#FFD991', '#FF7530', '#FF4024'), 
                       values = c(0, 0.3, 0.6, 1))+
  ggtitle('Clust 14 scores')

VlnPlot(wt_cerebrum_microglia,  features = 'clust14hsv_clust_14_score', group.by = 'treatment_time', pt.size = 0)
FeaturePlot(wt_cerebrum_microglia, reduction = 'micro.umap', features = 'clust14hsv_clust_14_score')
  ggtitle('Microglia')+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank())+
  ylab('UMAP2')+
  xlab('UMAP1')

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#### Pathways from https://www.sciencedirect.com/science/article/pii/S0896627325008037#sec2 ####
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

#Homeostatic
plotList_homeo <- lapply(c('Tmem119', 'Cx3cr1'), featurePlotLight, data = wt_cerebrum_microglia, reduction_choice = 'micro.umap', maxLim = 5.5)
do.call(ggarrange, c(plotList_homeo, common.legend = TRUE, legend = 'right'))

#inflammatory
plotList_inf <- lapply(c('Nfkb1', 'Il1b'), featurePlotLight, data = wt_cerebrum_microglia, reduction_choice = 'micro.umap', maxLim = 5.5)
do.call(ggarrange, c(plotList_inf, common.legend = TRUE, legend = 'right'))

#stress
plotList_stress <- lapply(c('Hsph1', 'Hspd1'), featurePlotLight, data = wt_cerebrum_microglia, reduction_choice = 'micro.umap', maxLim = 5.5)
do.call(ggarrange, c(plotList_stress, common.legend = TRUE, legend = 'right'))

FeaturePlot(wt_cerebrum_microglia, features = 'C5ar1', reduction = 'micro.umap')

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#### Virus associated microglia (VAM) from https://link.springer.com/article/10.1186/s12974-025-03471-x#Sec27 ####
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

#Import gene lists from sup table 3 (last excel sheet in supp)
subcluster_gene_list <- read.csv("~/Documents/ÖverbyLab/scPlots/galectin3_proj/microglia/vam_gene_list.csv")
vam_genes <- subcluster_gene_list[subcluster_gene_list$celltype == 'Microglia' & subcluster_gene_list$cluster == 14,]$gene
top_vam_genes <- head(vam_genes, n = 12)

vam_dot <- DotPlot(wt_cerebrum_microglia, features = top_vam_genes, group.by = 'custom_clusters', scale = FALSE)$data

#Dotplot of vam genes
pdf('~/Documents/ÖverbyLab/scPlots/galectin3_proj/microglia/Vam_dotplot.pdf', width = 6, height = 5)
ggplot(vam_dot, aes(x = id, y = features.plot, fill = avg.exp.scaled, size = pct.exp))+
  geom_point(pch = 21)+
  theme_classic()+
  scale_fill_gradientn(colours = c('white', '#FFD991', '#FF7530', '#FF4024'), 
                       values = c(0, 0.3, 0.6, 1))+
  ggtitle('VAM Markers')
dev.off()

#Score and violin plot of vam genes
wt_cerebrum_microglia <- AddModuleScore_UCell(wt_cerebrum_microglia, 
                                              features=list('vam' = vam_genes), maxRank = 1200, name = NULL)

pdf('~/Documents/ÖverbyLab/scPlots/galectin3_proj/microglia/vam_score_violin.pdf', height = 5, width = 7)
VlnPlot(object = wt_cerebrum_microglia, features = 'vam', group.by = 'custom_clusters', pt.size = 0)
dev.off()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#### Homeostasis markers from fig 3 https://www.nature.com/articles/s41380-026-03529-z#MOESM5 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
hm_markers <- c('Selplg', 'P2ry12', 'Arsb', 'Gatm', 'Tspo', 'Uba52', 'Tmem119', 'Elmo1',
                'Adap2', 'Sall1', 'Dleu2', 'Gm15564', 'C1qc', 'C1qa', 'Hexb', 'Aif1', 'Cx3cr1', 'Rps10')

hm_dot <- DotPlot(wt_cerebrum_microglia, features = hm_markers, group.by = 'custom_clusters', scale = FALSE)$data

pdf('~/Documents/ÖverbyLab/scPlots/galectin3_proj/microglia/homeostasis_dotplot.pdf', width = 6, height = 5)
ggplot(hm_dot, aes(x = id, y = features.plot, fill = avg.exp.scaled, size = pct.exp))+
  geom_point(pch = 21)+
  theme_classic()+
  scale_fill_gradientn(colours = c('white', '#FFD991', '#FF7530', '#FF4024'), 
                       values = c(0, 0.3, 0.6, 1))+
  ggtitle('Homeostasis Markers')
dev.off()

#Score and violin plot of hm genes
wt_cerebrum_microglia <- AddModuleScore_UCell(wt_cerebrum_microglia, 
                                              features=list('hm_score' = hm_markers), maxRank = 1200, name = NULL)

pdf('~/Documents/ÖverbyLab/scPlots/galectin3_proj/microglia/homeostasis_score_violin.pdf', height = 5, width = 7)
VlnPlot(object = wt_cerebrum_microglia, features = 'hm_score', group.by = 'custom_clusters', pt.size = 0)
dev.off()
# - - - - - - - - - - - - 
#### GSEA line plots ####
# - - - - - - - - - - - - 

#Load in gene lists
response_to_stress <- read.table("~/Documents/ÖverbyLab/gene_lists_for_gsea/respones_to_stress.csv", quote="\"", comment.char="")
response_to_stress <- unique(response_to_stress$V1)

#Create clusterprofiler gsea objects. each infected cluster vs mock
create_gsea_object <- function(dat, min.p.val = 0.05){
  dat_markers <- FindAllMarkers(dat, group.by = 'infected_clusters', logfc.threshold = 0, 
                             min.pct = 0)
  dat_markers <- dat_markers %>% as.data.frame() %>% dplyr::arrange(desc(avg_log2FC)) 
  dat_markers_list <- as.list(dat_markers$avg_log2FC)
  names(dat_markers_list) <- rownames(dat_markers)
  #Probably a better way to do this but the above creates a list of lists, so get rid of that now
  dat_markers_list <- unlist(dat_markers_list)
  
  dat_gsea<- gseGO(geneList = dat_markers_list,
                   OrgDb = org.Mm.eg.db,
                   ont   = "BP",
                   minGSSize = 100,
                   maxGSSize = 5000,
                   pvalueCutoff = min.p.val,
                   verbose = FALSE, 
                   keyType = 'SYMBOL')
}

#Need to make this a lapply statement where we run each group vs mock through
#function
microglia_gsea <- create_gsea_object(wt_cerebrum_microglia)

