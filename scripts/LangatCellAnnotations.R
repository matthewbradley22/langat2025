#CELL ANNOTATION 
#Load in libraries
library(Seurat)
library(ggplot2)
library(SingleR)
library(celldex)
library(pheatmap)
library(RColorBrewer)
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
#ParseSeuratObj_int <- subset(ParseSeuratObj_int, scDblFinderLabel == 'singlet')

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
FeaturePlot(ParseSeuratObj_int, 'Flt1', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Cd93', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Vwf', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Cldn5', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Pglyrp1', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Emcn', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Pecam1', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Adgrl4', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Slco1a4', reduction = 'umap.integrated')

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

#Micro/macrophage markers from Allen atlas https://knowledge.brain-map.org/celltypes/CCN202002013
FeaturePlot(ParseSeuratObj_int, 'Hexb', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Inpp5d', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Ms4a4a', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Cd74', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Klra2', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Cd209a', reduction = 'umap.integrated')

#Macrophage/monocyte markers
FeaturePlot(ParseSeuratObj_int, 'Ptprc', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Ccr2', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Lyz2', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Ccr2', reduction = 'umap.integrated')

#Oligo
FeaturePlot(ParseSeuratObj_int, 'Mag',reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Mog', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Plp1', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Gjc3',reduction = 'umap.integrated') #Also high in OPCs
FeaturePlot(ParseSeuratObj_int, 'Mbp', reduction = 'umap.integrated')

#OPC
FeaturePlot(ParseSeuratObj_int, 'Pdgfra', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Stk32a', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Cspg4', reduction = 'umap.integrated')

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
FeaturePlot(ParseSeuratObj_int, 'Smoc2', reduction = 'umap.integrated')
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


#### Look at DEGs of select groups ####

#39, fibroblasts?
markers39 <- FindMarkers(ParseSeuratObj_int, group.by = 'seurat_clusters', ident.1 = 39,
                         only.pos = TRUE)
head(markers39, n = 20)

#27 showing signs of endothelial + pericyte and smooth muscle
#Rgs5 pericyte marker, Adgrf5 endothelial
markers27 <- FindMarkers(ParseSeuratObj_int, group.by = 'seurat_clusters', ident.1 = 27,
                         only.pos = TRUE)
head(markers27, n = 20)

#17 - astrocyte and oligo markers.
#Coexpressing several neuronal genes with neurons in top degs.. Confusing
#Gria4, Kcnj3
#Not sure that A2m is expressed in neurons though
#Pax3 discussed in allen brain paper, expressed in astro and neurons
markers17 <- FindMarkers(ParseSeuratObj_int, group.by = 'seurat_clusters', ident.1 = 17,
                         only.pos = TRUE)
head(markers17, n = 30)
clust17 <- subset(ParseSeuratObj_int, seurat_clusters == 17)
clust17[['RNA']]$counts %>% t() %>% write.csv(file = './data/cluster17Counts.csv', row.names = TRUE)
#35 - doublets? Some oligodendrocyte markers but also some astrocyte etc...
markers35 <- FindMarkers(ParseSeuratObj_int, group.by = 'seurat_clusters', ident.1 = 35,
                         only.pos = TRUE)
head(markers35, n = 20)

#34, macro or micro
markers34 <- FindMarkers(ParseSeuratObj_int, group.by = 'seurat_clusters', ident.1 = 34,
                         only.pos = TRUE)
head(markers34, n = 20)

#25, not sure what these are. Neurons + a few astrocytes?
#Several top genes involved in neural development, immature neurons? Most genes also align w neurons in panglao
#Dscaml1, Meis2, Sox11, Nol4
#Sox11 are critical for neural precursor survival https://www.science.org/doi/full/10.1126/sciadv.abc6093
#Celf4: critical role of CELF in neurodevelopment https://www.sciencedirect.com/science/article/pii/S0969996124001244
#Bclla: https://www.mdpi.com/2079-7737/13/2/126
markers25 <- FindMarkers(ParseSeuratObj_int, group.by = 'seurat_clusters', ident.1 = 25,
                         only.pos = TRUE)
head(markers25, n = 30)
clust25 <- subset(ParseSeuratObj_int, seurat_clusters == 25)
clust25[['RNA']]$counts %>% t() %>% write.csv(file = './data/cluster25Counts.csv', row.names = TRUE)

#Read in data from allen mapMyCells
cluster25Map <- read_csv("data/cluster25Countscsv_10xWholeMouseBrain(CCN20230722)_HierarchicalMapping_UTC_1757663135149/cluster25Countscsv_10xWholeMouseBrain(CCN20230722)_HierarchicalMapping_UTC_1757663135149.csv", 
                                                                                   skip = 4)
cluster25Map
#Paper https://www.nature.com/articles/s41467-019-08453-1 has midbrain neuronal markers
FeaturePlot(ParseSeuratObj_int, 'Sbk1', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Dcx', reduction = 'umap.integrated')

#Allen brain atlas paper also has immature markers https://www.nature.com/articles/s41586-023-06812-z
FeaturePlot(ParseSeuratObj_int, 'Mex3a', reduction = 'umap.integrated')
FeaturePlot(ParseSeuratObj_int, 'Draxin', reduction = 'umap.integrated')



#Look at lower half of cluster 28 which is split across macrophages
umapCoords <- ParseSeuratObj_int@reductions$umap.integrated@cell.embeddings %>% as.data.frame()
umapCoords[umapCoords$umapintegrated_1]

#Opened cell selector and chose lower part of cluster 34, with macrophages, which is labelled Microglia
intUMAP <- DimPlot(ParseSeuratObj_int, reduction = 'umap.integrated')
cells.located <- CellSelector(plot = intUMAP)

#REDOING THIS CURRENTLY
#Custom annotation 
#Created on singlet data
ParseSeuratObj_int$manualAnnotation <- 
  case_when(ParseSeuratObj_int$seurat_clusters %in% c(26, 40, 44) &
              ParseSeuratObj_int$singleR_labels == 'Neurons' ~ 'Neurons',
            ParseSeuratObj_int$seurat_clusters %in% c(25) &
              ParseSeuratObj_int$singleR_labels != 'Astrocytes' ~ 'Immature Neurons',
            ParseSeuratObj_int$seurat_clusters %in% c(3, 4, 45, 15) & 
              ParseSeuratObj_int$singleR_labels == 'Astrocytes' ~ 'Astrocytes',
            ParseSeuratObj_int$seurat_clusters == 31 ~ 'Pericytes',
            ParseSeuratObj_int$seurat_clusters == 30 ~ 'Muscle cells',
            ParseSeuratObj_int$seurat_clusters %in% c(18, 24, 20)~ 'Choroid Plexus',
            ParseSeuratObj_int$seurat_clusters %in% c(2, 1, 8, 6, 21, 11) &
              ParseSeuratObj_int$singleR_labels == 'Microglia' ~ 'Microglia',
            ParseSeuratObj_int$seurat_clusters %in% c(0, 12, 7, 37)~ 'Macrophage/Monocytes', #Particularly convinced about the macrophages by the allan atlas markers above
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
            
            #ParseSeuratObj_int$seurat_clusters == 26 & 
           # ParseSeuratObj_int$scDblFinderLabel == '',
            .default = 'unknown') 

#ParseSeuratObj_int[[]][cells.located,]$manualAnnotation = 'Macrophage/Monocytes'

newCols <-  c(brewer.pal(12, 'Paired'), '#99FFE6', '#CE99FF', '#18662E',  '#FF8AEF','#737272')
DimPlot(ParseSeuratObj_int, label = FALSE, group.by = 'manualAnnotation', reduction = 'umap.integrated',
        cols = newCols)


#Load in doublet removed data to plot
ParseSeuratObj_int <- LoadSeuratRds("./data/FilteredRpcaIntegratedDatNoDoublets.rds")


ParseSeuratObj_int$manualAnnotation = factor(ParseSeuratObj_int$manualAnnotation, levels = 
                                               rev(c('Neurons', 'Immature Neurons', 'Microglia', 'Astrocytes',
                                                 'Macrophage/Monocytes', 'Choroid Plexus',
                                                 'Endothelial', 'Oligodendrocytes',
                                                 'Ependymal', 'B Cells', 'Pericytes',
                                                 'Muscle cells', 'T cells', 'Granulocytes',
                                                 'Nk cells', 'unknown')))
 
DotPlot(ParseSeuratObj_int, features = c('Snap25', 'Syt1', 'Csf1r', 'Cx3cr1', 'Tmem119',
                                         'Aqp4', 'Fgfr3','Gfap', 'Ccr2' , 'Ttr',
                                         'Kl', 'Mag', 'Mog', 'Nnat', 'Cfap54' ,'Mia' ,
                                         'Cd19', 'Ms4a1', 'Vtn', 'Abcc9', 'Acta2', 'Tagln',
                                         'Cd3e', 'Cd3d', 'S100a9', 'Il1r2', 'Clnk', 'Nkg7'
                                         ),
        group.by = 'manualAnnotation', assay = 'RNA')+
  theme(axis.text.x = element_text(angle = 75, vjust = 0.5))

#Look at rna features across clusters
ParseSeuratObj_int[[]] %>% dplyr::group_by(seurat_clusters) %>% 
  dplyr::summarise(meanFeatures = mean(nFeature_RNA)) %>% arrange(desc(meanFeatures)) %>% 
  ggplot(aes(x = seurat_clusters, y = meanFeatures))+
  geom_bar(stat = 'identity')+
  theme(axis.text.x = element_text(angle = 90))


#SaveSeuratRds(ParseSeuratObj_int, "./data/seuratSingletsAnnotated.rds")

