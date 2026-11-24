library(Seurat)
library(ggplot2)
library(SingleR)
library(celldex)
library(pheatmap)
library(RColorBrewer)
library(dplyr)

ParseSeuratObj_int <- LoadSeuratRds("~/Documents/ÖverbyLab/data/FilteredRpcaIntegratedDatNoDoublets.rds")

# - - - - - - - - - - - - - - 
#### Automatic annotation ####
# - - - - - - - - - - - - - - 

#Initial labelling with singleR in order to get an idea of what we have

#Mouse rna seq reference
ref <- fetchReference("mouse_rnaseq", "2024-02-26")

#Make singleR label predictions
pred <- SingleR(test=ParseSeuratObj_int[['RNA']]$data, ref=ref, labels=ref$label.main)
ParseSeuratObj_int$singleR_labels <- pred$labels

#Appear to have high confidence groups of astrocytes, endothelial cells, microglia, and clear infiltrating immune clusters
plotScoreHeatmap(pred)

#UMAP with cell assignments
DimPlot(ParseSeuratObj_int, group.by = 'singleR_labels', reduction = "umap.integrated",
        label = TRUE)


# - - - - - - - - - - - - - - - - - 
#### Manual curation of labels ####
# - - - - - - - - - - - - - - - - - 

#Look at known marker genes and use to make final cell labels

#neuron markers, markers are consistent
FeaturePlot(ParseSeuratObj_int, 'Snap25', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Rbfox3', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Syt1', reduction = 'umap.integrated')

#Endothelial cells, markers are consistent
FeaturePlot(ParseSeuratObj_int, 'Flt1', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Cd93', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Vwf', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Cldn5', reduction = 'umap.integrated')

#Astrocytes, markers less consistent but match singleR labels nicely
FeaturePlot(ParseSeuratObj_int, 'Gfap', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Aqp4', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Fgfr3', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Slc1a3', reduction = 'umap.integrated')

#Microglia, markers agree and match with singleR
FeaturePlot(ParseSeuratObj_int, 'Csf1r', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Tmem119', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'P2ry12', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Cx3cr1', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Itgam', reduction = 'umap.integrated')

#Macrophage/monocyte markers
FeaturePlot(ParseSeuratObj_int, 'Ccr2', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Lyz2', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Ptprc', reduction = 'umap.integrated')

#Oligodendocytes, very consistent markers
FeaturePlot(ParseSeuratObj_int, 'Mag',reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Mog', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Plp1', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Mbp', reduction = 'umap.integrated')

#OPC, quite inconsistent and inconclusive markers
FeaturePlot(ParseSeuratObj_int, 'Pdgfra', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Stk32a', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Cspg4', reduction = 'umap.integrated')

#T cells, consistent markers
FeaturePlot(ParseSeuratObj_int, 'Trbc2', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Cd3g', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Cd3d', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Cd3e', reduction = 'umap.integrated')


#NK cells, consistent
FeaturePlot(ParseSeuratObj_int, 'Nkg7', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Klrd1', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Ncr1', reduction = 'umap.integrated')

#Pericytes, also cosnsitent markers
FeaturePlot(ParseSeuratObj_int, 'Vtn', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Abcc9', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Kcnj8', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Atp13a5', reduction = 'umap.integrated')

#Smooth muscle cells, consistent next to pericytes
FeaturePlot(ParseSeuratObj_int, 'Acta2', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Tagln', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Cnn1', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Myh11', reduction = 'umap.integrated')

#Epithelial cells - not helpful
FeaturePlot(ParseSeuratObj_int, 'Krt14', reduction = 'umap.integrated') #Not present
FeaturePlot(ParseSeuratObj_int, 'Epcam', reduction = 'umap.integrated')

#Ependymal cells, consistent markers
FeaturePlot(ParseSeuratObj_int, 'Foxj1', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Dynlrb2', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Rabl2', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Cfap54', reduction = 'umap.integrated')

#Fibroblast markers, not consistent and lineup somewhat with pericytes
FeaturePlot(ParseSeuratObj_int, 'Vim', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Pdgfrb', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Smoc2', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Pdgfra', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Fbln1', reduction = 'umap.integrated')

#Choroid plexus epithilial cells, consistent markers
FeaturePlot(ParseSeuratObj_int, 'Ttr', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Aqp1', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Kl', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Clic6', reduction = 'umap.integrated')

#Neutrophils and other granluocytes, pretty consistent markers
FeaturePlot(ParseSeuratObj_int, 'S100a8', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Il1r2', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'S100a9', reduction = 'umap.integrated')

#B cells, appears to be a very small population
FeaturePlot(ParseSeuratObj_int, 'Ms4a1', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Cd19', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Ms4a1', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Cd79a', reduction = 'umap.integrated')

# - - - - - - - - - - 
#### Map my cells ####
# - - - - - - - - - - 

#Further confirmation of some clusters was checked at 
# the allen institute's 'map my cells' tool

#Write data to csv for mapMyCells
clust17 <- subset(ParseSeuratObj_int, seurat_clusters == 17)
clust17[['RNA']]$counts %>% t() %>% write.csv(file = './data/cluster17Counts.csv', row.names = TRUE)

#Read in data from allen mapMyCells
cluster17MapMy <- read_csv("data/cluster17Countscsv_10xWholeMouseBrain(CCN20230722)_HierarchicalMapping_UTC_1757599685232/cluster17Countscsv_10xWholeMouseBrain(CCN20230722)_HierarchicalMapping_UTC_1757599685232.csv", 
                           skip = 4)
#Dominated by Bergmann glial cells
table(cluster17MapMy$supertype_name)

#Write data to csv for mapMyCells
clust35 <- subset(ParseSeuratObj_int, seurat_clusters == 35)
clust35[['RNA']]$counts %>% t() %>% write.csv(file = './data/cluster35Counts.csv', row.names = TRUE)

#Read in data from allen mapMyCells
cluster35Map <- read_csv("data/cluster35Countscsv_10xWholeMouseBrain(CCN20230722)_HierarchicalMapping_UTC_1757928986340/cluster35Countscsv_10xWholeMouseBrain(CCN20230722)_HierarchicalMapping_UTC_1757928986340.csv", 
                         skip = 4)

#Write data to csv for mapMyCells
clust25 <- subset(ParseSeuratObj_int, seurat_clusters == 25)
clust25[['RNA']]$counts %>% t() %>% write.csv(file = './data/cluster25Counts.csv', row.names = TRUE)

#Read in data from allen mapMyCells
cluster25Map <- read_csv("data/cluster25Countscsv_10xWholeMouseBrain(CCN20230722)_HierarchicalMapping_UTC_1757663135149/cluster25Countscsv_10xWholeMouseBrain(CCN20230722)_HierarchicalMapping_UTC_1757663135149.csv", 
                         skip = 4)
cluster25Map

# - - - - - - - - - - - - 
#### Final annotation ####
# - - - - - - - - - - - - 

#Using all methods above, finalize manual labels for cell type of each cluster

ParseSeuratObj_int$manualAnnotation <- 
  case_when(ParseSeuratObj_int$seurat_clusters %in% c(26, 40, 44) &
              ParseSeuratObj_int$singleR_labels == 'Neurons' ~ 'Neurons',
            ParseSeuratObj_int$seurat_clusters %in% c(25) &
              ParseSeuratObj_int$singleR_labels != 'Astrocytes' ~ 'Immature Neurons',
            ParseSeuratObj_int$seurat_clusters %in% c(3, 4, 45, 15, 17) & 
              ParseSeuratObj_int$singleR_labels == 'Astrocytes' ~ 'Astrocytes',
            ParseSeuratObj_int$seurat_clusters == 31 ~ 'Pericytes',
            ParseSeuratObj_int$seurat_clusters == 30 ~ 'Muscle cells',
            ParseSeuratObj_int$seurat_clusters %in% c(18, 24, 20)~ 'Choroid Plexus',
            ParseSeuratObj_int$seurat_clusters %in% c(2, 1, 8, 6, 21, 11) &
              ParseSeuratObj_int$singleR_labels == 'Microglia' ~ 'Microglia',
            ParseSeuratObj_int$seurat_clusters %in% c(0, 12, 7, 37)~ 'Macrophage/Monocytes',
            ParseSeuratObj_int$seurat_clusters %in% c(5, 14, 10)~ 'Endothelial',
            ParseSeuratObj_int$seurat_clusters %in% c(13, 23)~ 'Oligodendrocytes',
            ParseSeuratObj_int$seurat_clusters %in% c(9)~ 'Ependymal',
            ParseSeuratObj_int$seurat_clusters %in% c(22) ~ 'T cells',
            ParseSeuratObj_int$seurat_clusters %in% c(29) ~ 'Nk cells',
            ParseSeuratObj_int$seurat_clusters %in% c(28) ~ 'Granulocytes',
            ParseSeuratObj_int$seurat_clusters %in% c(42) ~ 'B Cells',
            ParseSeuratObj_int$seurat_clusters %in% c() ~ 'Fibroblasts',
            ParseSeuratObj_int$seurat_clusters %in% c() ~ 'OPCs (?)',
            ParseSeuratObj_int$seurat_clusters %in% c(16, 32, 19, 27, 38) ~ 'Doublets',
            
            .default = 'unknown') 

#SaveSeuratRds(ParseSeuratObj_int, "./data/seuratSingletsAnnotated.rds")

#Plot relevant gene markers for all celltypes
no_unknown <- subset(ParseSeuratObj_int,  manualAnnotation != 'unknown')

no_unknown$manualAnnotation <- factor(no_unknown$manualAnnotation, 
                                              levels = c( 'T cells',   'Nk cells', 'Macrophage/Monocytes', 
                                                          'Granulocytes', 'B Cells',  'Pericytes', 'Oligodendrocytes','Neurons',
                                                          'Muscle cells', 'Microglia', 'Immature Neurons', 'Ependymal','Endothelial', 
                                                          'Choroid Plexus', 'Astrocytes'))

DotPlot(no_unknown, features = c('Aqp4', 'Fgfr3','Gfap', 'Ttr','Kl', 
                                       'Flt1', 'Pecam1','Vwf',
                                       'Cfap54' ,'Nnat', 'Mia',
                                       'Sox11', 'Celf4','Csf1r', 'Cx3cr1', 
                                       'Tmem119', 'Acta2', 'Tagln',
                                       'Snap25', 'Syt1',  
                                       'Mag', 'Mog', 'Abcc9', 'Vtn', 
                                       'Cd19', 'Ms4a1', 'S100a9', 'Il1r2', 
                                       'Klra2', 'Ccr2', 'Lyz2', 
                                       'Clnk','Nkg7','Cd3e', 'Cd3d'), 
        group.by = 'manualAnnotation', assay = 'RNA')+
  theme(axis.text.x = element_text(angle = 75, vjust = 0.5),
        axis.text = element_text(size = 16))+
  scale_color_gradient2(low = 'blue', mid = 'white', high = 'red', labels = c(-1, 0, 1,2),
                        breaks = c(-1, 0, 1,2), , limits = c(-2,2.5))+
  # scale_color_gradientn(colours = c('lightblue','white', '#FFD991', '#FF4024'), 
  #              values = c(0, 0.35, 0.5, 1),
  #             name = 'Average Expression')+
  scale_size_continuous(range = c(0,8), limits = c(0,100))+
  theme(legend.position = "bottom",
        legend.justification = "center",
        legend.direction = "horizontal",
        legend.title = element_text(hjust = 0.3),
        legend.spacing.x = unit(2, "cm"))+
  guides(size = guide_legend(title.position = "top", title = 'Percent Expressed', theme = theme(legend.key.width  = unit(1, "cm"))),
         color = guide_colorbar(title.position = "top", title = 'Average Scaled Expression', theme = theme(legend.key.width  = unit(4, "cm"))))



