library(SingleR)
library(celldex)
library(readr)
library(scDblFinder)
library(stringr)
source('~/Documents/ÖverbyLab//scripts/langatFunctions.R')

#Annotate yang japanese encephelitis cells
yang_data <- LoadSeuratRds("~/Documents/ÖverbyLab/yang_public_data/yang_rpca_integrated_obj.RDS")
DimPlot(yang_data,reduction = "umap.rpca", label = TRUE)+
  theme(legend.position = 'none')

#Look at cell type markers

#Some neuron markers
FeaturePlot(yang_data, 'Snap25', reduction = 'umap.rpca')
FeaturePlot(yang_data, 'Rbfox3', reduction = 'umap.rpca')
FeaturePlot(yang_data, 'Dpp10', reduction = 'umap.rpca')
FeaturePlot(yang_data, 'Syt1', reduction = 'umap.rpca')

#Ex neurons
FeaturePlot(yang_data, 'Slc17a7', reduction = 'umap.rpca')
FeaturePlot(yang_data, 'Sv2b', reduction = 'umap.rpca')
FeaturePlot(yang_data, 'Arpp21', reduction = 'umap.rpca')

#in neurons
FeaturePlot(yang_data, 'Gad1', reduction = 'umap.rpca')
FeaturePlot(yang_data, 'Gad2', reduction = 'umap.rpca')
FeaturePlot(yang_data, 'Dlx6os1', reduction = 'umap.rpca')


#Endothelial cells
#brain endo marker Pglyrp1 from https://www.sciencedirect.com/science/article/pii/S0092867420300623
FeaturePlot(yang_data, 'Flt1', reduction = 'umap.rpca')
FeaturePlot(yang_data, 'Cd93', reduction = 'umap.rpca')
FeaturePlot(yang_data, 'Vwf', reduction = 'umap.rpca')
FeaturePlot(yang_data, 'Cldn5', reduction = 'umap.rpca')
FeaturePlot(yang_data, 'Pglyrp1', reduction = 'umap.rpca')
FeaturePlot(yang_data, 'Emcn', reduction = 'umap.rpca')
FeaturePlot(yang_data, 'Pecam1', reduction = 'umap.rpca')
FeaturePlot(yang_data, 'Adgrl4', reduction = 'umap.rpca')
FeaturePlot(yang_data, 'Slco1a4', reduction = 'umap.rpca')


#Astrocytes
FeaturePlot(yang_data, 'Gfap', reduction = 'umap.rpca')
FeaturePlot(yang_data, 'Aqp4', reduction = 'umap.rpca')
FeaturePlot(yang_data, 'Fgfr3', reduction = 'umap.rpca')
FeaturePlot(yang_data, 'Aldh1l1', reduction = 'umap.rpca')
FeaturePlot(yang_data, 'Slc1a3', reduction = 'umap.rpca')
FeaturePlot(yang_data, 'Gli3', reduction = 'umap.rpca')
FeaturePlot(yang_data, 'Slc39a12', reduction = 'umap.rpca')

#Microglia
FeaturePlot(yang_data, 'Ctss', reduction = 'umap.rpca')
FeaturePlot(yang_data, 'Csf1r', reduction = 'umap.rpca')
FeaturePlot(yang_data, 'Tmem119', reduction = 'umap.rpca')
FeaturePlot(yang_data, 'P2ry12', reduction = 'umap.rpca')
FeaturePlot(yang_data, 'Cx3cr1', reduction = 'umap.rpca')
FeaturePlot(yang_data, 'Itgam', reduction = 'umap.rpca')

#Macrophage/monocyte markers
FeaturePlot(yang_data, 'Ptprc', reduction = 'umap.rpca')
FeaturePlot(yang_data, 'Ccr2', reduction = 'umap.rpca')
FeaturePlot(yang_data, 'Lyz2', reduction = 'umap.rpca')
#Ly6c2 used in yng paper
FeaturePlot(yang_data, 'Ly6c2', reduction = 'umap.rpca')

#Oligo
FeaturePlot(yang_data, 'Mag',reduction = 'umap.rpca')
FeaturePlot(yang_data, 'Mog', reduction = 'umap.rpca')
FeaturePlot(yang_data, 'Plp1', reduction = 'umap.rpca')
FeaturePlot(yang_data, 'Gjc3',reduction = 'umap.rpca') #Also high in OPCs
FeaturePlot(yang_data, 'Mbp', reduction = 'umap.rpca')
FeaturePlot(yang_data, 'Olig1', reduction = 'umap.rpca')

#OPC
FeaturePlot(yang_data, 'Pdgfra', reduction = 'umap.rpca')
FeaturePlot(yang_data, 'Stk32a', reduction = 'umap.rpca')
FeaturePlot(yang_data, 'Cspg4', reduction = 'umap.rpca')

#T cells
#Panglao and this paper https://www.nature.com/articles/s41467-022-32627-z/figures/1
FeaturePlot(yang_data, 'Trbc2', reduction = 'umap.rpca')
FeaturePlot(yang_data, 'Cd3g', reduction = 'umap.rpca')
FeaturePlot(yang_data, 'Cd3d', reduction = 'umap.rpca')
FeaturePlot(yang_data, 'Cd3e', reduction = 'umap.rpca')

#NK cells same sources as above
#This paper supp table 2 also has some https://academic.oup.com/bioinformatics/article/38/3/785/6390798
FeaturePlot(yang_data, 'Nkg7', reduction = 'umap.rpca')
FeaturePlot(yang_data, 'Klrd1', reduction = 'umap.rpca')
FeaturePlot(yang_data, 'Ncr1', reduction = 'umap.rpca')
FeaturePlot(yang_data, 'Klrb1b', reduction = 'umap.rpca')
FeaturePlot(yang_data, 'Clnk', reduction = 'umap.rpca')


#Pericytes
#Along with PangloaDB, this paper has pericyte markers: 
#https://www.sciencedirect.com/science/article/pii/S1537189124001605
FeaturePlot(yang_data, 'Vtn', reduction = 'umap.rpca')
FeaturePlot(yang_data, 'Abcc9', reduction = 'umap.rpca')
FeaturePlot(yang_data, 'Kcnj8', reduction = 'umap.rpca')
FeaturePlot(yang_data, 'Atp13a5', reduction = 'umap.rpca')
FeaturePlot(yang_data, 'Anpep', reduction = 'umap.rpca')
FeaturePlot(yang_data, 'Ace2', reduction = 'umap.rpca')


#Label celltypes
#Saved some cell labels in yang_mic_mac_analysis.R script, so running this again will actually remove those
#updated macrophage microglia labels. Updating script below to include that analysis, but just something to be cognizant of
yang_data[[]] <- yang_data[[]] %>% dplyr::mutate(manualAnnotation = case_when(seurat_clusters %in% c(19, 10, 11, 20, 5, 14, 4, 41)~'Neurons',
                                                                              seurat_clusters %in% c(9, 39, 27)~'Endothelial',
                                                                              seurat_clusters %in% c(18, 46)~'Astrocytes',
                                                                              seurat_clusters %in% c(2)~'Microglia',
                                                                              seurat_clusters %in% c(0, 1, 17, 24, 40)~'Macrophage/Monocytes',
                                                                              seurat_clusters %in% c(21, 22, 13, 7, 6, 43, 23, 28)~'Oligo',
                                                                              seurat_clusters %in% c(3, 44, 37, 44)~'OPC',
                                                                              seurat_clusters %in% c(8, 33)~'T cells',
                                                                              seurat_clusters %in% c(26)~'Nk cells',
                                                                              seurat_clusters %in% c(47, 16, 35)~'Peri',
                                                                              seurat_clusters %in% c(38, 25)~'Immature Neurons', #Based on genes below
                                                                              .default = 'unknown'))


#Plot assignments
newCols <-  c(brewer.pal(12, 'Paired'), '#99FFE6', '#CE99FF', '#18662E','#737272',  '#FF8AEF')
newCols[11] =  '#FF8AEF'

DimPlot(yang_data, group.by = 'manualAnnotation', reduction = "umap.rpca",
        label = FALSE, cols = newCols)+
  xlab('')+
  ylab('')

#Check my annotations vs singleR
#Mouse rna seq reference
data_for_labelling <- c(paste0("yang_data[['RNA']]@layers$data.", c(1:6)))
ref <- fetchReference("mouse_rnaseq", "2024-02-26")

#Run separately on each layer, not sure better way for this
#Creating this vector is unnecessary and creates pretty large object, fix later
data_layers <- c(yang_data[['RNA']]$data.1, yang_data[['RNA']]$data.2, yang_data[['RNA']]$data.3, yang_data[['RNA']]$data.4, yang_data[['RNA']]$data.5, yang_data[['RNA']]$data.6)

prediction_list <- list()
for(i in 1:length(data_layers)){
  pred <- SingleR(test=data_layers[[i]], ref=ref, labels=ref$label.main)
  prediction_list[[i]] = pred
}

cell_labels <- as.data.frame(do.call(rbind, prediction_list))
cell_labels$id <- rownames(cell_labels)
yang_data$id <- colnames(yang_data)
yang_data[[]] <- left_join(yang_data[[]], cell_labels[c('id', 'pruned.labels')], by = 'id')
yang_data$single_r_labels <- yang_data$pruned.labels
yang_data$pruned.labels <- NULL

plotScoreHeatmap(prediction_list[[2]])

#UMAP with cell assignments
DimPlot(yang_data, group.by = 'single_r_labels', reduction = "umap.rpca",
        label = TRUE)+
  theme(legend.position = 'none')
yang_dat_merged_layers <- prepseu

########### Allen brain atlas cell assignments ###############
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
yang_dat_merged_layers <- JoinLayers(yang_data)
clust16_47_35 <- subset(yang_dat_merged_layers, seurat_clusters %in% c(16, 47, 35))
clust16_47_35[['RNA']]$counts %>% t() %>% write.csv(file = '~/Documents/ÖverbyLab/yang_public_data/gene_counts_for_allen_atlas/16_47_35.csv', row.names = TRUE)

cluster16_47_35Map <- read_csv("~/Documents/ÖverbyLab/yang_public_data/gene_counts_for_allen_atlas/16_47_35_allen/16_47_35csv_10xWholeMouseBrain(CCN20230722)_HierarchicalMapping_UTC_1771585059472.csv", 
                         skip = 4)
table(cluster16_47_35Map$subclass_name) %>% sort()

#Look at cluster 12
clust12 <- subset(yang_dat_merged_layers, seurat_clusters %in% c(12))
clust12[['RNA']]$counts %>% t() %>% write.csv(file = '~/Documents/ÖverbyLab/yang_public_data/gene_counts_for_allen_atlas/12.csv', row.names = TRUE)

cluster12Map <- read_csv("~/Documents/ÖverbyLab/yang_public_data/gene_counts_for_allen_atlas/12csv_10xWholeMouseBrain(CCN20230722)_HierarchicalMapping_UTC_1771599174687/12csv_10xWholeMouseBrain(CCN20230722)_HierarchicalMapping_UTC_1771599174687.csv", 
                               skip = 4)
table(cluster12Map$subclass_name) %>% sort()
cluster12Map_labels <- cluster12Map[c('cell_id', 'subclass_name')]
yang_data[[]] = left_join(yang_data[[]], cluster12Map_labels, by = c('id' = 'cell_id'))
DimPlot(yang_data, reduction = "umap.rpca", group.by = 'subclass_name')
yang_data$manualAnnotation[yang_data$subclass_name == '334 Microglia NN'] = 'Microglia' 

#Check for doublets
#dbr.sd set to 1 since i'm not sure expected doublet rate with parse
sce_dbl <- scDblFinder(yang_dat_merged_layers[['RNA']]$counts, dbr.sd=1)
hist(sce_dbl$scDblFinder.score)

yang_data$sc_dbl_labels <- sce_dbl$scDblFinder.class
DimPlot(yang_data, group.by = 'sc_dbl_labels', reduction = "umap.rpca")

#Look at all neurons
#Based on markers, clusters 10 and 19 appear excitatory, and the rest a bit inconclusive
neurons <- subset(yang_dat_merged_layers, manualAnnotation == 'Neuron')
neurons[['RNA']]$counts %>% t() %>% write.csv(file = '~/Documents/ÖverbyLab/yang_public_data/gene_counts_for_allen_atlas/yang_neurons.csv', row.names = TRUE)

neuron_map <- read_csv("~/Documents/ÖverbyLab/yang_public_data/gene_counts_for_allen_atlas/yang_neuron_map/yang_neuronscsv_10xWholeMouseBrain(CCN20230722)_HierarchicalMapping_UTC_1773742365225.csv", 
                         skip = 4)
table(neuron_map$class_name) %>% sort()
neuron_map <- neuron_map[c('cell_id', 'class_name')]
colnames(neuron_map) = c('id', 'neuron_labelling_mapMy')
neuron_map['neuron_labelling_mapMy'] <-  stringr::word(neuron_map$neuron_labelling_mapMy, -1)
yang_data[[]] <- left_join(yang_data[[]], neuron_map)
DimPlot(yang_data,reduction = "umap.rpca", group.by = 'neuron_labelling_mapMy', cols = newCols)

#Update neuron labels and remove group of likely doublets
yang_data <- subset(yang_data, seurat_clusters != 4)
yang_data[[]] <- yang_data[[]] %>% mutate(manualAnnotation = case_when(neuron_labelling_mapMy == 'GABA' ~ 'In Neurons',
                                                                                                 neuron_labelling_mapMy == 'Glut' ~ 'Ex Neurons',
                                                                                                 .default = manualAnnotation))
DimPlot(yang_data, group.by = 'manualAnnotation', reduction = 'umap.rpca', cols = newCols)

#Look at cluster 4 markers to see what celltype it could be, or possibly doublets
#Compare to oligo and then to neurons. comparing to all celltypes not super helpful
#Higher expression of Mag and Mog than neurons, but don't appear to be oligo... Either weird infected cell type
#or doublets, maybe remove or label as unknown for now
table(yang_data$seurat_clusters, yang_dat_merged_layers$neuron_labelling_mapMy)
clust4_markers <- FindMarkers(yang_dat_merged_layers, group.by = 'seurat_clusters', ident.1 = 4, ident.2 = 7,
                              only.pos = TRUE,
                              test.use = 'MAST')
clust4_markers_vs_neurons <- FindMarkers(yang_dat_merged_layers, group.by = 'seurat_clusters', ident.1 = 4, ident.2 = 10,
                              only.pos = TRUE,
                              test.use = 'MAST')
clust4_markers %>% dplyr::mutate(pct.diff = pct.1 - pct.2) %>% dplyr::arrange(desc(pct.diff))
clust4_markers_vs_neurons %>% dplyr::mutate(pct.diff = pct.1 - pct.2) %>% dplyr::arrange(desc(pct.diff))

#Look at cluster 25
clust25_markers <- FindMarkers(yang_dat_merged_layers, group.by = 'seurat_clusters', ident.1 = 25,
                              only.pos = TRUE,
                              test.use = 'MAST')

#Could be immature neurons based on marker mapping to allen atlas
head(clust25_markers)
FeaturePlot(yang_data, reduction = 'umap.rpca', features = 'Igfbpl1') 

#Immature neuron markers
FeaturePlot(yang_data, 'Mex3a', reduction = 'umap.rpca')
FeaturePlot(yang_data, 'Draxin', reduction = 'umap.rpca')
FeaturePlot(yang_data, 'Dcx', reduction = 'umap.rpca')

############# Macrophage labelling ############# 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 


#Subset to just macrophage and microglia
mic_mac <- subset(yang_data, seurat_clusters %in% c(0, 1, 40, 24, 17, 15, 12, 2, 29))
#Clustering will change but want to be able to refer back to original clusters
mic_mac$previous_clusters = mic_mac$seurat_clusters

#Prep for umap
mic_mac <- prepSeuratObj(mic_mac)
ElbowPlot(mic_mac, ndims = 30)
mic_mac <- prepUmapSeuratObj(mic_mac, nDims = 20, reductionName = 'mac_umap', num_neighbors = 30L)

DimPlot(mic_mac, reduction = "mac_umap", label = FALSE, group.by = 'manualAnnotation', cols = newCols)
#Remove columns that match name of later left_join with allan data, shouldn't be here
mic_mac$subclass_name <- NULL
mic_mac$neuron_labelling_mapMy <- NULL
#Look at important gene markers
#Microglia
FeaturePlot(mic_mac, features = 'Tmem119', reduction = 'mac_umap')
FeaturePlot(mic_mac, features = 'Cx3cr1', reduction = 'mac_umap')
FeaturePlot(mic_mac, 'Csf1r', reduction = 'mac_umap')

#Macrophages
FeaturePlot(mic_mac, features = 'Ly6c2', reduction = 'mac_umap')
FeaturePlot(mic_mac, 'Ccr2', reduction = 'mac_umap')
FeaturePlot(mic_mac, 'Flt1', reduction = 'mac_umap') #From yang paper

#Look at clusters in mapMyCells
mic_mac <- JoinLayers(mic_mac)
mic_mac[['RNA']]$counts %>% t() %>% write.csv(file = '~/Documents/ÖverbyLab/yang_public_data/gene_counts_for_allen_atlas/mic_mac_cells.csv', row.names = TRUE)

#Check mapMyCells results
mic_mac_countMy <- read_csv("~/Documents/ÖverbyLab/yang_public_data/gene_counts_for_allen_atlas/mic_mac_mapMyCells/mic_mac_countMyCells.csv", 
                            skip = 4)
mic_mac_countMy <- dplyr::filter(mic_mac_countMy, class_name %in% c('34 Immune'))
mic_mac_countMy_labels <- mic_mac_countMy[c('cell_id', 'subclass_name')]
mic_mac[[]] = left_join(mic_mac[[]], mic_mac_countMy_labels, by = c('id' = 'cell_id'))

#These labels look correct to me based on gene markers
DimPlot(mic_mac, reduction = "mac_umap", group.by = 'subclass_name', cols = newCols)

#Match manual annotation to subclass names
mic_mac$manualAnnotation = 'unknown'
mic_mac[[]][mic_mac$subclass_name == '334 Microglia NN' & !is.na(mic_mac$subclass_name),]$manualAnnotation = 'Microglia'
#I find mapmy to generally label macrophages as dendritic cells due to its low count of actual macrophges in the data
mic_mac[[]][mic_mac$subclass_name %in% c('336 Monocytes NN', '337 DC NN', '335 BAM NN') & !is.na(mic_mac$subclass_name),]$manualAnnotation = 'Macrophage/Monocytes'

#Check results back in main UMAP
yang_data[[]][yang_data$id %in% mic_mac$id,]$manualAnnotation = mic_mac$manualAnnotation
DimPlot(yang_data,reduction = "umap.rpca", label = FALSE, group.by = 'manualAnnotation', cols = newCols)

#Save data labels
SaveSeuratRds(yang_data, file = "~/Documents/ÖverbyLab/yang_public_data/yang_rpca_integrated_obj.RDS")

