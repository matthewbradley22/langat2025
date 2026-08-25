
#Source function
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
receptors <- c('Axl', 'Havcr1', 'Timd4', 'Tyro3', 'Lrp8', 'Lrp1',
               'Lrp4', 'Cd209a')

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
  ylab('')
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

mock_samples <- subset(sn_integrated_dat, infected == FALSE)

DotPlot(mock_samples, features = receptors, group.by = 'geno_celltype', scale = FALSE)$data %>% 
  tidyr::separate(col = 'id', into = c('genotype', 'celltype'), sep = '_') %>% 
  ggplot(aes(x = features.plot, y = celltype, fill = avg.exp.scaled, size = pct.exp))+
  geom_point(pch = 21)+
  facet_wrap(~genotype)+
  theme_classic()+
  scale_fill_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white"),
                       values = c(1.0,0.7,0.4,0))



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
infected_neurons <- prepUmapSeuratObj(infected_neurons, nDims = 20, reductionName = 'inf_neuron.umap', resolution_value = 0.4)

DimPlot(infected_neurons, reduction = 'neuron.umap', label = FALSE, group.by = 'manualAnnotation',
        label.size = 6)

DimPlot(infected_neurons, reduction = 'neuron.umap', label = FALSE, group.by = 'new_genotype',
        label.size = 6)+
  ggtitle('Neurons')

DimPlot(infected_neurons, reduction = 'neuron.umap', label = TRUE,
        label.size = 6)

FeaturePlot(infected_neurons, reduction = 'neuron.umap', features = 'rna_LGTV')

#Need to split by genotype
DotPlot(infected_neurons, features = 'rna_LGTV', group.by = 'seurat_clusters', scale = FALSE)$data %>% 
  ggplot(aes(x = features.plot, y = id, fill = avg.exp.scaled, size = pct.exp))+
  geom_point(pch = 21)+
  theme_classic()+
  scale_fill_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white"),
                       values = c(1.0,0.7,0.4,0))+
  ggtitle('Infected samples')+
  xlab('')+
  ylab('')+
  scale_size(range = c(0, 8))

#Facet by genotype but first subset so points only shown up in respective genotype if there are enough cells in that group/geno
cluster_geno_counts = infected_neurons[[]] %>% group_by(seurat_clusters, geno_simplified) %>% 
  dplyr::summarise(num = n()) %>% 
  dplyr::mutate(enough_cells = ifelse(num > 25, yes = TRUE, no = FALSE))

infected_neurons$genotype_clusters <- paste(infected_neurons$geno_simplified, infected_neurons$seurat_clusters, sep = '_')
DotPlot(infected_neurons, features = receptors, group.by = 'genotype_clusters', scale = FALSE)$data %>% 
  tidyr::separate(col = 'id', into = c('geno_simplified', 'seurat_clusters'), sep = '_') %>% 
  dplyr::left_join(cluster_geno_counts, by = c('geno_simplified', 'seurat_clusters')) %>% 
  dplyr::filter(enough_cells == TRUE) %>% 
  #create 0 data for havcr1 since we don't have it
  add_row(avg.exp = 0, pct.exp = NA, features.plot = 'Havcr1', geno_simplified = 'wt', seurat_clusters = factor(0), avg.exp.scaled = NA) %>% 
  dplyr::mutate(seurat_clusters = factor(seurat_clusters, levels = rev(as.character(seq(0,17))))) %>% 
  dplyr::mutate(geno_simplified = factor(geno_simplified, levels = c('wt', 'KO'))) %>% 
  ggplot(aes(x = features.plot, y = seurat_clusters, fill = avg.exp.scaled, size = pct.exp))+
  geom_point(pch = 21)+
  theme_classic()+
  facet_wrap(~geno_simplified)+
  scale_fill_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white"),
                       values = c(1.0,0.7,0.4,0))+
  xlab('')+
  ylab('')


