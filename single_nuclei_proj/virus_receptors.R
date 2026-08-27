
#Source functions
source('~/Documents/ÖverbyLab/scripts/langatFunctions.R')

#Read in processed data
sn_integrated_dat <- LoadSeuratRds('~/Documents/ÖverbyLab/single_nuclei_proj/LGTVscCombined.rds')

newCols <-  c(brewer.pal(12, 'Paired'), '#99FFE6', '#CE99FF', '#18662E','#737272',  '#FF8AEF')

#Check data
DimPlot(sn_integrated_dat, group.by =  'manualAnnotation', cols = newCols, reduction = 'umap.integrated')+
  xlab('Umap1')+
  ylab('Umap1')+
  theme(axis.ticks = element_blank(),
        axis.text=element_blank(),
        legend.text=element_text(size=17))+
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank())+
  ggtitle('single nuclei WT')

#Clean up genotype column
sn_integrated_dat$geno_simplified <- case_when(sn_integrated_dat$new_genotype %in% c('KO', 'KO (same)') ~ 'KO',
                                               sn_integrated_dat$new_genotype %in% c('wt', 'wt (same)') ~ 'wt')

#Look at various flavivirus receptors, and compare to infection levels + lgtv to compare infection level with receptor
receptors <- c('Axl', 'Havcr1', 'Havcr2', 'Timd4', 'Tyro3', 'Lrp8', 'Lrp1',
               'Lrp4', 'Cd209a', 'Cd209b', 'Cd209c', 'Itgb4', 'Hspa5', 'Ncam1', 'Hspa1a', 'Vim', 'Itgav',
               'Itgb3', 'Cldn1', 'Clec5a', 'Mrc1', 'Mer', 'Scarb1')

#Variables of interest
table(sn_integrated_dat$manualAnnotation, sn_integrated_dat$geno_simplified, sn_integrated_dat$infected)

#Group variables for plotting later
sn_integrated_dat$geno_celltype <- paste(sn_integrated_dat$geno_simplified, sn_integrated_dat$manualAnnotation, sep = '_')
sn_integrated_dat$treatment_celltype <- paste(sn_integrated_dat$infected, sn_integrated_dat$manualAnnotation, sep = '_')

# - - - - - - - - - - - - - - - - - - - - - -
#### Infection and receptors by celltype ####
# - - - - - - - - - - - - - - - - - - - - - -

pdf('~/Documents/ÖverbyLab/single_nuclei_proj/langat_receptors/receptor_dotplot.pdf', width = 10, height = 5)
DotPlot(sn_integrated_dat, features = receptors, group.by = 'treatment_celltype', scale = FALSE)$data %>% 
  tidyr::separate(col = 'id', into = c('treatment', 'celltype'), sep = '_') %>% 
  #create 0 data for havcr1 since we don't have it
  add_row(avg.exp = 0, pct.exp = NA, features.plot = 'Havcr1', treatment = 'FALSE', celltype = 'Astrocytes', avg.exp.scaled = NA) %>% 
  dplyr::mutate(infected = factor(ifelse(treatment == FALSE, yes = 'Mock', no = 'LGTV'), levels = c('Mock', 'LGTV'))) %>% 
  ggplot(aes(x = features.plot, y = celltype, fill = avg.exp.scaled, size = pct.exp))+
  geom_point(pch = 21)+
  facet_wrap(~infected)+
  theme_classic()+
  scale_fill_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white"),
                       values = c(1.0,0.7,0.4,0))+
  xlab('')+
  ylab('')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(fill = NA, color = "black", linetype = "dashed"))
dev.off()

#Split by infection for plots
infected_samples <- subset(sn_integrated_dat, infected == TRUE)

pdf('~/Documents/ÖverbyLab/single_nuclei_proj/langat_receptors/lgtv_dotplot.pdf', width = 10, height = 5)
DotPlot(infected_samples, features = 'rna_LGTV', group.by = 'geno_celltype', scale = FALSE)$data %>% 
  tidyr::separate(col = 'id', into = c('genotype', 'celltype'), sep = '_') %>% 
  dplyr::mutate(genotype = factor(genotype, levels = c('wt', 'KO'))) %>% 
  ggplot(aes(x = features.plot, y = celltype, fill = avg.exp.scaled, size = pct.exp))+
  geom_point(pch = 21)+
  facet_wrap(~genotype)+
  theme_classic()+
  scale_fill_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white"),
                       values = c(1.0,0.7,0.4,0))+
  ggtitle('Infected samples')+
  xlab('')+
  ylab('')+
  scale_size(range = c(0, 8))

dev.off()  

#Correlations with virus
#Look at correlations between receptors and viruses by celltype
present_receptors <- receptors[receptors %in% rownames(sn_integrated_dat)]

infected_split_celltype <- SplitObject(infected_samples, split.by = "manualAnnotation")

lapply(infected_split_celltype, FUN = function(x){
  receptor_expression <- x[['RNA']]$data[c(present_receptors, 'LGTV'),] 
  
  t(as.matrix(receptor_expression)) %>% 
    cor() %>% 
    as.data.frame() %>% 
    dplyr::select(LGTV) %>% 
    dplyr::arrange(desc(LGTV))
})



# - - - - - - - - - - - - - - - - - - - - - - - -
#### Infection across neuronal subclusters ####
# - - - - - - - - - - - - - - - - - - - - - - - -
neurons <- subset(sn_integrated_dat, manualAnnotation %in% c('Ex Neurons', 'In Neurons'))
infected_neurons <- subset(infected_samples, manualAnnotation %in% c('Ex Neurons', 'In Neurons'))

neurons <- prepSeuratObj(neurons)
ElbowPlot(neurons, ndims = 40)
neurons <- prepUmapSeuratObj(neurons, nDims = 20, reductionName = 'neuron.umap', resolution_value = 0.4)

DimPlot(neurons, reduction = 'neuron.umap', label = FALSE, group.by = 'manualAnnotation',
        label.size = 6)

DimPlot(neurons, reduction = 'neuron.umap', label = FALSE, group.by = 'new_genotype',
        label.size = 6)+
  ggtitle('Neurons')

DimPlot(neurons, reduction = 'neuron.umap', label = FALSE, group.by = 'infected',
        label.size = 6)+
  ggtitle('Neurons')

DimPlot(neurons, reduction = 'neuron.umap', label = TRUE,
        label.size = 6)

#Infected neurons only
infected_neurons <- prepSeuratObj(infected_neurons)
ElbowPlot(infected_neurons, ndims = 40)
infected_neurons <- prepUmapSeuratObj(infected_neurons, nDims = 20, reductionName = 'inf_neuron.umap', resolution_value = 0.8)

DimPlot(infected_neurons, reduction = 'inf_neuron.umap', label = FALSE, group.by = 'manualAnnotation',
        label.size = 6)

DimPlot(infected_neurons, reduction = 'inf_neuron.umap', label = FALSE, group.by = 'new_genotype',
        label.size = 6)+
  ggtitle('Neurons')

pdf('~/Documents/ÖverbyLab/single_nuclei_proj/langat_receptors/neuron_clusters.pdf', width = 6, height = 6)
DimPlot(infected_neurons, reduction = 'inf_neuron.umap', label = TRUE,
        label.size = 6)
dev.off()

FeaturePlot(infected_neurons, reduction = 'inf_neuron.umap', features = 'rna_LGTV')

#Need to split by genotype
pdf('~/Documents/ÖverbyLab/single_nuclei_proj/langat_receptors/neuron_clusters_lgtv.pdf', width = 8, height = 7)
DotPlot(infected_neurons, features = 'rna_LGTV', group.by = 'seurat_clusters', scale = FALSE)$data %>% 
  ggplot(aes(x = features.plot, y = id, fill = avg.exp.scaled, size = pct.exp))+
  geom_point(pch = 21)+
  theme_classic()+
  scale_fill_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white"),
                       values = c(1.0,0.7,0.4,0))+
  ggtitle('Infected neuron lgtv by cluster')+
  xlab('')+
  ylab('')+
  scale_size(range = c(0, 8))
dev.off()

#Facet by genotype but first subset so points only shown up in respective genotype if there are enough cells in that group/geno
cluster_geno_counts = infected_neurons[[]] %>% group_by(seurat_clusters, geno_simplified) %>% 
  dplyr::summarise(num = n()) %>% 
  dplyr::mutate(enough_cells = ifelse(num > 25, yes = TRUE, no = FALSE))

num_clusters <- max(as.numeric((infected_neurons$seurat_clusters)))

infected_neurons$genotype_clusters <- paste(infected_neurons$geno_simplified, infected_neurons$seurat_clusters, sep = '_')

pdf('~/Documents/ÖverbyLab/single_nuclei_proj/langat_receptors/neuron_clusters_receptors.pdf', width = 9, height = 8)
DotPlot(infected_neurons, features = receptors, group.by = 'genotype_clusters', scale = FALSE)$data %>% 
  tidyr::separate(col = 'id', into = c('geno_simplified', 'seurat_clusters'), sep = '_') %>% 
  dplyr::left_join(cluster_geno_counts, by = c('geno_simplified', 'seurat_clusters')) %>% 
  dplyr::filter(enough_cells == TRUE) %>% 
  #create 0 data for havcr1 since we don't have it
  add_row(avg.exp = 0, pct.exp = NA, features.plot = 'Havcr1', geno_simplified = 'wt', seurat_clusters = factor(0), avg.exp.scaled = NA) %>% 
  dplyr::mutate(seurat_clusters = factor(seurat_clusters, levels = rev(as.character(seq(0,num_clusters))))) %>% 
  dplyr::mutate(geno_simplified = factor(geno_simplified, levels = c('wt', 'KO'))) %>% 
  ggplot(aes(x = features.plot, y = seurat_clusters, fill = avg.exp.scaled, size = pct.exp))+
  geom_point(pch = 21)+
  theme_classic()+
  facet_wrap(~geno_simplified)+
  scale_fill_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white"),
                       values = c(1.0,0.7,0.4,0))+
  xlab('')+
  ylab('')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  geom_hline(yintercept =  20.5, linetype = 'dashed', col = 'blue')+
  geom_hline(yintercept =  4.5, linetype = 'dashed', col = 'blue')+
  geom_hline(yintercept =  16.5, linetype = 'dashed', col = 'blue')+
  ggtitle('Infected neuron clusters')
dev.off()

#Look at wt neurons split by infected and uninfected, and see which receptors change
infected_neurons$lgtv_present <- infected_neurons[['RNA']]$data['LGTV',] > 0

pdf('~/Documents/ÖverbyLab/single_nuclei_proj/langat_receptors/neurons_bystander_infected_dot.pdf', width = 7, height = 6)
DotPlot(infected_neurons, features = receptors, group.by = 'lgtv_present', scale = FALSE)$data %>% 
  ggplot(aes(x = id, y = features.plot, fill = avg.exp.scaled, size = pct.exp))+
  geom_point(pch = 21)+
  theme_classic()+
  scale_fill_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white"),
                       values = c(1.0,0.7,0.4,0))+
  ggtitle('Neurons bystander vs infected')+
  xlab('')+
  ylab('')+
  scale_size(range = c(0, 8))
dev.off()

# - - - - - - - - - - - - - - - - 
#### Macropahge infection ####
# - - - - - - - - - - - - - - - - 
infected_macro <- subset(infected_samples, manualAnnotation %in% c('Macrophages'))
infected_macro$lgtv_present <- infected_macro[['RNA']]$data['LGTV',] > 0

pdf('~/Documents/ÖverbyLab/single_nuclei_proj/langat_receptors/macrophage_bystander_infected_dot.pdf', width = 7, height = 6)
DotPlot(infected_macro, features = receptors, group.by = 'lgtv_present', scale = FALSE)$data %>% 
  ggplot(aes(x = id, y = features.plot, fill = avg.exp.scaled, size = pct.exp))+
  geom_point(pch = 21)+
  theme_classic()+
  scale_fill_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white"),
                       values = c(1.0,0.7,0.4,0))+
  ggtitle('Macrophage bystander vs infected')+
  xlab('')+
  ylab('')+
  scale_size(range = c(0, 8))
dev.off()
