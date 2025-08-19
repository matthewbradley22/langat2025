#CELL ANNOTATION 
#Load in libraries
library(Seurat)
library(SingleR)
library(celldex)
library(pheatmap)
library(dplyr)

#Load in data
#ParseSeuratObj_int <- LoadSeuratRds('./ccaIntegratedDat.rds')
ParseSeuratObj_int <- LoadSeuratRds("./data/FilteredRpcaIntegratedDat.rds")

#Load parse well label map
wellMap <- data.frame(well = c(paste0('A', seq(1,12)), paste0('B', seq(1,12)),
                               paste0('C', seq(1,12)),   paste0('D', seq(1,12))))

####Annotate cell types ####
#Will try automatic annotation
#Mouse rna seq reference
ref <- fetchReference("mouse_rnaseq", "2024-02-26")
pred <- SingleR(test=ParseSeuratObj_int[['RNA']]$data, ref=ref, labels=ref$label.main)
ParseSeuratObj_int$singleR_labels <- pred$labels
plotScoreHeatmap(pred)

#UMAP with cell assignments
DimPlot(ParseSeuratObj_int, group.by = 'singleR_labels', reduction = "umap.integrated",
        label = TRUE)
clusterAssignments <- table(Assigned=pred$pruned.labels, Cluster= ParseSeuratObj_int$seurat_clusters)
pheatmap(log2(clusterAssignments+10), color=colorRampPalette(c("white", "blue"))(101))

####Data is loaded in with doublets, so can start here ####
#Look at doublets in integrated data
DimPlot(ParseSeuratObj_int, group.by = 'scDblFinderLabel', reduction = "umap.integrated")
ParseSeuratObj_int <- subset(ParseSeuratObj_int, scDblFinderLabel == 'singlet')

#Plot viral counts vs well/treatment. 
#removing doublets doesn't get rid of viral contamination in PBS samples
ggplot(ParseSeuratObj_int[[]], aes(x = orig.ident, y = virusCount, col = Treatment))+
  geom_point() +
  scale_x_discrete(labels= wellMap$well)+
  theme(axis.text.x = element_text(angle = 90))

#Count of virus presence
ParseSeuratObj_int[[]] %>% mutate(virusPresent = ifelse(virusCount>0, 'yes', 'no')) %>% 
  group_by(Treatment, virusPresent) %>% dplyr::summarise(count = n()) %>% 
  ggplot(aes(x = Treatment, y = count, fill = virusPresent))+
  geom_bar(stat = 'identity', position = 'dodge')

#Plot of virus presence over small amount
ParseSeuratObj_int[[]] %>% mutate(virusHight = ifelse(virusCount>3, 'yes', 'no')) %>% 
  group_by(Treatment, virusHight) %>% dplyr::summarise(count = n()) %>% 
  ggplot(aes(x = Treatment, y = count, fill = virusHight))+
  geom_bar(stat = 'identity', position = 'dodge')

#Look at virus presence across celltypes
ParseSeuratObj_int[[]] <- ParseSeuratObj_int[[]] %>% mutate(virusPresence = ifelse(virusCount > 0, 'yes', 'no')) 
DimPlot(ParseSeuratObj_int, group.by = 'virusPresence', reduction = 'umap.integrated')

ParseSeuratObj_int[[]] %>% group_by(singleR_labels, virusPresence) %>% dplyr::summarise(count = n()) %>% 
  ggplot(aes(x = singleR_labels, y = count, fill = virusPresence))+
  geom_bar(stat = 'identity', position = 'dodge')+
  theme(axis.text.x = element_text(angle = 90))

#Look at cell counts per well
table(ParseSeuratObj_int$Well) %>% barplot()

#Most pbs cells with virus have very low amount, could just filter by amount
pbsWithVirus <- subset(ParseSeuratObj_int, Treatment == 'PBS' & virusCountPAdj > 0)
summary(pbsWithVirus$virusCountPAdj)
hist(pbsWithVirus$virusCountPAdj)

#### Manual annotiation#### 

#Paper has some canonical markers for cell types
#https://umu.diva-portal.org/smash/get/diva2:1897514/FULLTEXT01.pdf

#This paper https://link.springer.com/article/10.1007/s10571-021-01159-3
#Heterogeneity and Molecular Markers for CNS Glial Cells Revealed by Single-Cell Transcriptomics also 
#has lists of markers

#Some neuron markers
FeaturePlot(ParseSeuratObj_int, 'Snap25', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Rbfox3', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Dpp10', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Syt1', reduction = 'umap.integrated')


#Endothelial cells
#brain endo marker Pglyrp1 from https://www.sciencedirect.com/science/article/pii/S0092867420300623
#Pretty much all markers showing up in 3, 7, 22, 23.
FeaturePlot(ParseSeuratObj_int, 'Flt1', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Cd93', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Vwf', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Cldn5', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Pglyrp1', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Emcn', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Pecam1', reduction = 'umap.integrated')

#Astrocytes
FeaturePlot(ParseSeuratObj_int, 'Gfap', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Aqp4', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Fgfr3', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Aldh1l1', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Slc1a3', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Gli3', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Slc39a12', reduction = 'umap.integrated')

#Microglia
FeaturePlot(ParseSeuratObj_int, 'Ctss', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Csf1r', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Tmem119', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'P2ry12', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Cx3cr1', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Itgam', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'C1qa', reduction = 'umap.integrated')


#Oligo
FeaturePlot(ParseSeuratObj_int, 'Mag',reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Mog', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Mbp', reduction = 'umap.integrated')

#OPC
FeaturePlot(ParseSeuratObj_int, 'Pdgfra', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Cspg4', reduction = 'umap.integrated')

#Macrophage markers
FeaturePlot(ParseSeuratObj_int, 'Ptprc', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Ccr2', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Lyz2', reduction = 'umap.integrated')

#T cells
#Panglao and this paper https://www.nature.com/articles/s41467-022-32627-z/figures/1
FeaturePlot(ParseSeuratObj_int, 'Trbc2', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Cd3g', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Cd3d', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Cd3e', reduction = 'umap.integrated')

#NK cells same sources as above
#This paper supp table 2 also has some https://academic.oup.com/bioinformatics/article/38/3/785/6390798
FeaturePlot(ParseSeuratObj_int, 'Nkg7', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Klrd1', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Ncr1', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Klrb1b', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Clnk', reduction = 'umap.integrated')

#Pericytes
#Along with PangloaDB, this paper has pericyte markers: 
#https://www.sciencedirect.com/science/article/pii/S1537189124001605
FeaturePlot(ParseSeuratObj_int, 'Vtn', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Abcc9', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Kcnj8', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Atp13a5', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Anpep', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Ace2', reduction = 'umap.integrated')

#Smooth muscle cells paper here for markers: https://www.sciencedirect.com/science/article/pii/S1534580722006852
#Should express these, and not fibroblast markers (Pdgfra), or endothelial markers (Pecam1)
FeaturePlot(ParseSeuratObj_int, 'Acta2', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Tagln', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Cnn1', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Myh11', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Myl9', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Rgs5', reduction = 'umap.integrated')

#Epithelial cells - not helpful
FeaturePlot(ParseSeuratObj_int, 'Krt14', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Epcam', reduction = 'umap.integrated')

#Ependymal cells
#Markers from here https://www.frontiersin.org/journals/cellular-neuroscience/articles/10.3389/fncel.2021.703951/full
#And panglao and from here https://ars.els-cdn.com/content/image/1-s2.0-S1534580723000035-gr1.jpg
FeaturePlot(ParseSeuratObj_int, 'Foxj1', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Dynlrb2', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Rabl2', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Cfap54', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Nnat', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Mia', reduction = 'umap.integrated')

#Cardiomocytes markers
FeaturePlot(ParseSeuratObj_int, 'Fgf2', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Tcap', reduction = 'umap.integrated')
#Hopx marker from here https://pmc.ncbi.nlm.nih.gov/articles/PMC6220122/
FeaturePlot(ParseSeuratObj_int, 'Hopx', reduction = 'umap.integrated')

#Cardiomocytes markers from https://pmc.ncbi.nlm.nih.gov/articles/PMC10905351/
FeaturePlot(ParseSeuratObj_int, 'Tnnt1', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Actn1', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Tpm3', reduction = 'umap.integrated')

#Fibroblast markers
#Good fribroblast paper https://pmc.ncbi.nlm.nih.gov/articles/PMC10318398/
#Says some of these markers cannot differentiate between fibro and muscle cells
FeaturePlot(ParseSeuratObj_int, 'Vim', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Pdgfrb', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Pdgfra', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Fbln1', reduction = 'umap.integrated')

#Choroid plexus cells
#TTR only really expressed in choroid plexus looks like https://www.proteinatlas.org/ENSG00000118271-TTR/brain
#More markers here https://www.nature.com/articles/s41380-021-01416-3 and from Pangao
FeaturePlot(ParseSeuratObj_int, 'Ttr', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Aqp1', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Kl', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Clic6', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Folr1', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Prlr', reduction = 'umap.integrated')

#Neutrophils and other granluocytes
#Some markers from a 10x list linked here https://www.10xgenomics.com/support/software/cell-ranger/latest/tutorials/cr-tutorial-neutrophils
FeaturePlot(ParseSeuratObj_int, 'Csf3r', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'S100a8', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Il1r2', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'S100a9', reduction = 'umap.integrated')

#B cells. some markers in supplement of https://www.sciencedirect.com/science/article/pii/S0304383524000582#appsec1
FeaturePlot(ParseSeuratObj_int, 'Ms4a1', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Cd19', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Ms4a1', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Cd79a', reduction = 'umap.integrated')


#Dimplots for convenience
DimPlot(ParseSeuratObj_int, label = TRUE, reduction = 'umap.integrated')
DimPlot(ParseSeuratObj_int, label = TRUE, group.by = 'singleR_labels', reduction = 'umap.integrated')

#Unsure about cluster 23, has both endothelial and pericyte/smooth muscle markers expressed

#### Look at DEGs of select groups ####
#think 31 is pericytes based on degs and top markers matching up w pericytes in panglaodb
markers31 <- FindMarkers(ParseSeuratObj_int, group.by = 'seurat_clusters', ident.1 = 31,
                       only.pos = TRUE)

#Looks like could be smooth muscle cells based on plugging markers in pangaodb and marker expression above
markers32 <- FindMarkers(ParseSeuratObj_int, group.by = 'seurat_clusters', ident.1 = 32,
                       only.pos = TRUE)

#16, as well as 26 and 17 may be Choroid plexus cells based on Ttr, Kl, and Folr1 expression
markers16 <- FindMarkers(ParseSeuratObj_int, group.by = 'seurat_clusters', ident.1 = 16,
                         only.pos = TRUE)
markers26 <- FindMarkers(ParseSeuratObj_int, group.by = 'seurat_clusters', ident.1 = 26,
                         only.pos = TRUE)
markers17 <- FindMarkers(ParseSeuratObj_int, group.by = 'seurat_clusters', ident.1 = 26,
                         only.pos = TRUE)

#Check that cluster 6 expresses microglial degs
#Cluster 6 seems to express a lot of macrophage markers compared to other microglia clusters, not sure what to label
markers6 <- FindMarkers(ParseSeuratObj_int, group.by = 'seurat_clusters', ident.1 = 6,
                         only.pos = TRUE)
markers2 <- FindMarkers(ParseSeuratObj_int, group.by = 'seurat_clusters', ident.1 = 2,
                        only.pos = TRUE)

#Check macrophage/monocyte clusters, 1, 9, 13  look like macrophages/monocytes
markers13 <- FindMarkers(ParseSeuratObj_int, group.by = 'seurat_clusters', ident.1 = 13,
                        only.pos = TRUE)

#Check 14 and 33 clusters
markers14 <- FindMarkers(ParseSeuratObj_int, group.by = 'seurat_clusters', ident.1 = 14,
                         only.pos = TRUE)

#Glance at granulocyte cluster. Looks like neutrophil markers popping up
markers30 <- FindMarkers(ParseSeuratObj_int, group.by = 'seurat_clusters', ident.1 = 30,
                         only.pos = TRUE)

#Sort of looks like macrophage (Lyz2, Adgre1)
markers36 <- FindMarkers(ParseSeuratObj_int, group.by = 'seurat_clusters', ident.1 = 36,
                         only.pos = TRUE)
#34, is it astrocytes? S100b high up, indicates astrocytes
#but also some oligodendrocyte genes (Mobp)
markers34 <- FindMarkers(ParseSeuratObj_int, group.by = 'seurat_clusters', ident.1 = 34,
                         only.pos = TRUE)

#Look at lower half of cluster 28 which is split across macrophages
umapCoords <- ParseSeuratObj_int@reductions$umap.integrated@cell.embeddings %>% as.data.frame()
umapCoords[umapCoords$umapintegrated_1]

#Opened cell selector and chose lower part of cluster 28, below macrophages, which is labelled Microglia
#Only 52 cells which shouldn't matter much, but going to check their markers
intUMAP <- DimPlot(ParseSeuratObj_int, reduction = 'umap.integrated')
cells.located <- CellSelector(plot = intUMAP)
susMicroglia <- ParseSeuratObj_int[[]][cells.located,] %>% filter(singleR_labels == 'Microglia'&
                                                    manualAnnotation == 'Microglia')
ParseSeuratObj_int$susMicroglia <- ifelse(rownames(ParseSeuratObj_int[[]]) %in% rownames(susMicroglia), 'yes', 'no') 
DimPlot(ParseSeuratObj_int, reduction = 'umap.integrated', group.by = 'susMicroglia')


#Custom annotation 
#Created on singlet data
ParseSeuratObj_int$manualAnnotation <- 
  case_when(ParseSeuratObj_int$seurat_clusters %in% c('27', '24',
                                                      '37','44') &
              ParseSeuratObj_int$singleR_labels == 'Neurons' ~ 'Neurons',
            ParseSeuratObj_int$seurat_clusters %in% c('4', '5', '12', '15', '18', '21', '35', '40',
                                                      '25', '34', '38', '39', '33', '14', '43',
                                                      '27', '24')&
              ParseSeuratObj_int$singleR_labels == 'Astrocytes' ~ 'Astrocytes',
            ParseSeuratObj_int$seurat_clusters == '31'~ 'Pericytes',
            ParseSeuratObj_int$seurat_clusters == '32'~ 'Muscle cells',
            ParseSeuratObj_int$seurat_clusters %in% c(16,26,17)~ 'Choroid Plexus',
            ParseSeuratObj_int$seurat_clusters %in% c(2, 11, 0, 10, 6)~ 'Microglia',
            ParseSeuratObj_int$seurat_clusters %in% c(28) &
              ParseSeuratObj_int$singleR_labels == 'Microglia' ~ 'Microglia',
            ParseSeuratObj_int$seurat_clusters %in% c(28) &
              ParseSeuratObj_int$singleR_labels == 'Macrophages' ~ 'Macrophage/Monocytes',
            ParseSeuratObj_int$seurat_clusters %in% c(1, 9, 13, 36)~ 'Macrophage/Monocytes',
            ParseSeuratObj_int$seurat_clusters %in% c(3, 7, 22)~ 'EC',
            ParseSeuratObj_int$seurat_clusters %in% c(8, 19, 25)~ 'Oligodendrocytes',
            ParseSeuratObj_int$seurat_clusters %in% c(14,33)~ 'Ependymal',
            ParseSeuratObj_int$seurat_clusters %in% c(20) ~ 'T cells',
            ParseSeuratObj_int$seurat_clusters %in% c(29) ~ 'Nk cells',
            ParseSeuratObj_int$seurat_clusters %in% c(30) ~ 'Granulocytes',
            ParseSeuratObj_int$seurat_clusters %in% c(42) ~ 'B Cells',
            .default = 'unkown')


newCols <-  c(brewer.pal(12, 'Paired'), '#99FFE6', '#CE99FF', '#737272')
DimPlot(ParseSeuratObj_int, label = FALSE, group.by = 'manualAnnotation', reduction = 'umap.integrated',
        cols = newCols)

#There is an odd group under macrophages that is labelled microglia, need to look more

#SaveSeuratRds(ParseSeuratObj_int, "./data/seuratSingletsAnnotated.rds")

