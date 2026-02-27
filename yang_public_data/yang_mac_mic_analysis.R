library(Seurat)
library(dplyr)
source('~/Documents/ÖverbyLab//scripts/langatFunctions.R')

#Load in data
yang_data <- LoadSeuratRds("~/Documents/ÖverbyLab/yang_public_data/yang_rpca_integrated_obj.RDS")
newCols <-  c(brewer.pal(12, 'Paired'), '#99FFE6', '#CE99FF', '#18662E','#737272',  '#FF8AEF')
newCols[11] =  '#FF8AEF'
DimPlot(yang_data,reduction = "umap.rpca", label = FALSE, group.by = 'manualAnnotation', cols = newCols)

#Subset to just macrophage and microglia
mic_mac <- subset(yang_data, seurat_clusters %in% c(0, 1, 40, 24, 17, 15, 12, 2, 29))

#Prep for umap
mic_mac <- prepSeuratObj(mic_mac)
ElbowPlot(mic_mac, ndims = 30)
mic_mac <- prepUmapSeuratObj(mic_mac, nDims = 20, reductionName = 'mac_umap', num_neighbors = 30L)

DimPlot(mic_mac, reduction = "mac_umap", label = FALSE, group.by = 'manualAnnotation', cols = newCols)

#Look at important gene markers
#Microglia
FeaturePlot(mic_mac, features = 'Tmem119', reduction = 'mac_umap')
FeaturePlot(mic_mac, features = 'Cx3cr1', reduction = 'mac_umap')
FeaturePlot(mic_mac, 'Csf1r', reduction = 'mac_umap')

#Macrophages
FeaturePlot(mic_mac, features = 'Ly6c2', reduction = 'mac_umap')
FeaturePlot(mic_mac, 'Ccr2', reduction = 'mac_umap')

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
mic_mac[[]][mic_mac$subclass_name %in% c('336 Monocytes NN', '337 DC NN', '335 BAM NN') & !is.na(mic_mac$subclass_name),]$manualAnnotation = 'Macro/Mono'

#Check results back in main UMAP
yang_data[[]][yang_data$id %in% mic_mac$id,]$manualAnnotation = mic_mac$manualAnnotation
DimPlot(yang_data,reduction = "umap.rpca", label = FALSE, group.by = 'manualAnnotation', cols = newCols)

