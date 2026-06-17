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

wt_cerebrum_microglia <- subset(wt_cerebrum, manualAnnotation == 'Microglia')
table(wt_cerebrum_microglia$Genotype, wt_cerebrum_microglia$Treatment, 
      wt_cerebrum_microglia$Timepoint, wt_cerebrum_microglia$manualAnnotation)

#Microglia UMAP
wt_cerebrum_microglia <- prepSeuratObj(wt_cerebrum_microglia)
ElbowPlot(wt_cerebrum_microglia, ndims = 40)
wt_cerebrum_microglia <- prepUmapSeuratObj(wt_cerebrum_microglia, nDims = 15, reductionName = 'micro.umap',
                                             resolution_value = 0.8)

pdf('~/Documents/ÖverbyLab/scPlots/galectin3_proj/microglia/microglia_treatment_umap.pdf', width = 5, height = 5)
DimPlot(wt_cerebrum_microglia, reduction = 'micro.umap', label = FALSE, group.by = 'Treatment',
        label.size = 6, cols = c('#6DC3F8', '#D6644B', '#208d1f'))+
  ggtitle('Microglia')+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank())+
  ylab('UMAP2')+
  xlab('UMAP1')
dev.off()

pdf('~/Documents/ÖverbyLab/scPlots/galectin3_proj/microglia/microglia_time_umap.pdf', width = 5, height = 5)
DimPlot(wt_cerebrum_microglia, reduction = 'micro.umap', label = FALSE, group.by = 'Timepoint',
        label.size = 6, cols = c('#6DC3F8', '#D6644B', '#208d1f'))+
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
DimPlot(wt_cerebrum_microglia, reduction = 'micro.umap', label = FALSE, group.by = 'seurat_clusters')+
      ggtitle('Microglia')+
      theme(axis.text = element_blank(), axis.ticks = element_blank())+
      ylab('UMAP2')+
  xlab('UMAP1')

mock_markers <- FindMarkers(wt_cerebrum_microglia, ident.1 = 5, ident.2 = 0, group.by = 'seurat_clusters', test.use = 'MAST')
mock_markers <- dplyr::filter(mock_markers, p_val_adj < 0.01 & avg_log2FC > 1)
FeaturePlot(wt_cerebrum_microglia, features = 'Egr1', reduction = 'micro.umap')

mock_comp_paths <- gprofiler2::gost(query = rownames(mock_markers), organism = 'mmusculus', evcodes = TRUE,
                 sources = c('GO:BP', 'KEGG', 'GO:CC', 'GO:MP'))

# - - - - - - - - - - - - - - - - - - 
#### Across time infected samples ####
# - - - - - - - - - - - - - - - - - - 

#Look across time at microglia
microglia_infected <- subset(wt_cerebrum_microglia, Treatment == 'rLGTV')
microglia_mock <- subset(wt_cerebrum_microglia, Treatment == 'PBS')

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

#Mock microglia UMAP
microglia_mock <- prepSeuratObj(microglia_mock)
ElbowPlot(microglia_mock, ndims = 40)
microglia_mock <- prepUmapSeuratObj(microglia_mock, nDims = 15, reductionName = 'micro.mock.umap',
                                        resolution_value = 0.8)

pdf('~/Documents/ÖverbyLab/scPlots/galectin3_proj/microglia/mock_micro_umap.pdf', width = 5, height = 5)
DimPlot(microglia_mock, reduction = 'micro.mock.umap', label = FALSE, group.by = 'Timepoint',
        label.size = 6, cols = c('#6DC3F8', '#D6644B', '#208d1f'))+
  ggtitle('Mock microglia')+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank())+
  ylab('UMAP2')+
  xlab('UMAP1')
dev.off()

#Timepoint markers microglia
#Infected

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

#PBS
mock_time_markers <- FindAllMarkers(microglia_mock, group.by = 'Timepoint', test.use = 'MAST')
top_mock_markers <- mock_time_markers %>% dplyr::filter(avg_log2FC > 1 & p_val_adj < 0.01) %>% 
  dplyr::group_by(cluster) %>% 
  dplyr::slice_head(n = 5)

DotPlot(microglia_mock, features = top_mock_markers$gene, group.by = 'Timepoint', scale = FALSE)$data %>% 
  ggplot(aes(x = features.plot, y = id, size = pct.exp, fill = avg.exp.scaled))+
  geom_point(pch = 21)+
  scale_fill_gradientn(colours = c('white', '#FFD991', '#FF7530', '#FF4024'), 
                       values = c(0, 0.3, 0.6, 1))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        text = element_text(size = 20))+
  ylab('')+
  xlab('')+
  ggtitle('Mock Microglia')

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

#Look at negative regulators of inflammation in infected samples across time
neg_reg_genes <- c('Socs3', 'Il1rn', 'Il27', 'Lilrb4a')

pdf('~/Documents/ÖverbyLab/scPlots/galectin3_proj/microglia/micro_neg_ref_inflammation.pdf', width = 6, height = 5)
DotPlot(microglia_infected, features =neg_reg_genes, group.by = 'Timepoint', scale = FALSE)$data %>% 
  ggplot(aes(x = features.plot, y = id, size = pct.exp, fill = avg.exp.scaled))+
  geom_point(pch = 21)+
  scale_fill_gradientn(colours = c('white', '#FFD991', '#FF7530', '#FF4024'), 
                       values = c(0, 0.3, 0.6, 1))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        text = element_text(size = 20))+
  ylab('')+
  xlab('')+
  ggtitle('Negative reg inflammation')
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


#- - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#### DEGs between infected clusters #####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - 

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

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#### DEGs between each infected cluster and all cells including mock #####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

infected_vs_all_markers <- FindAllMarkers(wt_cerebrum_microglia, group.by = 'infected_clusters', test.use = 'MAST',
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

#No DEGs in cluster 1
names(infected_vs_all_path_list) = paste0('cluster_', c('0', '2', '3', '4', 'mock'))

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
infected_vs_all_path_list$cluster_2
infected_vs_all_path_list$cluster_3
infected_vs_all_path_list$cluster_4
infected_vs_all_path_list$cluster_mock



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

#Neuro protection
plotList_np <- lapply(neuro_protect, featurePlotLight, data = wt_cerebrum_microglia, reduction_choice = 'micro.umap', maxLim = 4.5)
do.call(ggarrange, c(plotList_np, common.legend = TRUE, legend = 'right'))

#Phagotcytosis
plotList_phago <- lapply(phagocytosis, featurePlotLight, data = wt_cerebrum_microglia, reduction_choice = 'micro.umap', maxLim = 6)
do.call(ggarrange, c(plotList_phago, common.legend = TRUE, legend = 'right'))

#Inflammation
plotList_inf <- lapply(inflammation, featurePlotLight, data = wt_cerebrum_microglia, reduction_choice = 'micro.umap', maxLim = 4.7)
do.call(ggarrange, c(plotList_inf, common.legend = TRUE, legend = 'right'))

#cytokine production
plotList_cyto <- lapply(cyto_production, featurePlotLight, data = wt_cerebrum_microglia, reduction_choice = 'micro.umap', maxLim = 4.7)
do.call(ggarrange, c(plotList_cyto, common.legend = TRUE, legend = 'right'))

#antigen presentation
plotList_antigen <- lapply(antigen_pres, featurePlotLight, data = wt_cerebrum_microglia, reduction_choice = 'micro.umap', maxLim = 5.3)
do.call(ggarrange, c(plotList_antigen, common.legend = TRUE, legend = 'right'))

#ifn sig
plotList_ifn <- lapply(ifn_sig, featurePlotLight, data = wt_cerebrum_microglia, reduction_choice = 'micro.umap', maxLim = 5.5)
do.call(ggarrange, c(plotList_ifn, common.legend = TRUE, legend = 'right'))

#proliferation
plotList_prolif <- lapply(proliferation, featurePlotLight, data = microglia_infected, reduction_choice = 'micro.inf.umap', maxLim = 4.5)
pdf('~/Documents/ÖverbyLab/scPlots/galectin3_proj/microglia/proliferation_genes.pdf', width = 6, height = 5)
do.call(ggarrange, c(plotList_prolif, common.legend = TRUE, legend = 'right'))
dev.off()

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
lps_dot <- DotPlot(wt_cerebrum_microglia, features = lps_markers, group.by = 'infected_clusters', scale = FALSE)$data
lps_dot$id = factor(lps_dot$id, levels = c('mock', '0', '1', '2', '3', '4'))

ggplot(lps_dot, aes(x = id, y = features.plot, fill = avg.exp.scaled, size = pct.exp))+
  geom_point(pch = 21)+
  theme_classic()+
  scale_fill_gradientn(colours = c('white', '#FFD991', '#FF7530', '#FF4024'), 
                       values = c(0, 0.3, 0.6, 1))+
  ggtitle('IAM (inflammatory markers)')

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

# - - - - - - - - - - - - - - - - - -
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
FeaturePlot(wt_cerebrum_microglia, features = 'Adgre1', reduction = 'micro.umap')+
  ggtitle('F4/80')

DotPlot(wt_cerebrum_microglia, features = 'Adgre1', group.by = 'infected_clusters', scale = FALSE)$data %>% 
  dplyr::mutate(id = factor(id, levels = c('mock', '0', '1', '2', '3', '4'))) %>% 
  ggplot(aes(x = id, y = features.plot, fill = avg.exp.scaled, size = pct.exp))+
  geom_point(pch = 21)+
  theme_classic()+
  scale_fill_gradientn(colours = c('white', '#FFD991', '#FF7530', '#FF4024'), 
                       values = c(0, 0.3, 0.6, 1))+
  ggtitle('F4/80')

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