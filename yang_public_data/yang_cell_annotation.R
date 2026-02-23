library(SingleR)
library(celldex)
library(readr)

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
FeaturePlot(yang_data, 'Ccr2', reduction = 'umap.rpca')


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
yang_data[[]] <- yang_data[[]] %>% dplyr::mutate(manualAnnotation = case_when(seurat_clusters %in% c(19, 10, 11, 20, 5, 14, 4, 41)~'Neuron',
                                                                              seurat_clusters %in% c(9, 39, 27)~'Endothelial',
                                                                              seurat_clusters %in% c(18, 46)~'Astrocytes',
                                                                              seurat_clusters %in% c(2)~'Microglia',
                                                                              seurat_clusters %in% c(0, 1, 17, 24, 40)~'Macro/Mono',
                                                                              seurat_clusters %in% c(21, 22, 13, 7, 6)~'Oligo',
                                                                              seurat_clusters %in% c(3, 44, 37, 44)~'OPC',
                                                                              seurat_clusters %in% c(8)~'T cells',
                                                                              seurat_clusters %in% c(26)~'Nk cells',
                                                                              seurat_clusters %in% c(47, 16, 35)~'Pericytes',
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
