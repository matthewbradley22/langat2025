#Load libraries
library(Seurat)
library(ggplot2)
library(tidyr)
library(dplyr)
source('./langatFunctions.R')

#This is where the 10x data is
setwd('~/Documents/ÖverbyLab/single_nuclei_proj/')

#Read in seurat data, or process from raw data below
sn_integrated_dat <- LoadSeuratRds('./LGTVscCombined.rds')


# - - - - - - - - - - - - - - - - - - - - - - - - - 
#### Process raw data and create Seurat object ####
# - - - - - - - - - - - - - - - - - - - - - - - - - 


#Start here to perform integration
outputDirs <- list.dirs(path = ".", full.names = TRUE, recursive = TRUE)
outputDirs <- outputDirs[grep('filtered_feature_bc_matrix', outputDirs)]
metaData <- read.delim("./samplemeta.csv.correct")

#Begin processing data
scDat <- lapply(outputDirs, FUN = function(x){
  dat <- Read10X(data.dir =x)
  dat
})

scDatSeu <-  lapply(scDat, FUN = function(x){
  CreateSeuratObject(counts = x, project = "hupAnalysis", min.cells = 3, min.features = 300)
})

#Merge for v5 ingtegration
toMerge <- c(scDatSeu[[2]], scDatSeu[[3]], scDatSeu[[4]],
             scDatSeu[[5]], scDatSeu[[6]], scDatSeu[[7]])

#Merge data before integration
scDatObj <- merge(scDatSeu[[1]], y = toMerge, 
                  add.cell.ids = paste0('seuObj', seq(1,7)), project = "single_nuclei")



#Label mitochondrial gene expression
scDatObj[["percent.mt"]] <- PercentageFeatureSet(scDatObj, pattern = "^mt-")

#Look at data features
Idents(scDatObj) <- "all"  #Stops violin plot from grouping by seurat cluster
VlnPlot(scDatObj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3, pt.size = 0) 

get_library_labels <- do.call(rbind, stringr::str_split(string = rownames(scDatObj[[]]), pattern = '_'))
lib_labels <- get_library_labels[,1]
scDatObj$library <- lib_labels
metaData$library = paste0('seuObj', seq(1:7))
scDatObj[[]] <-  left_join(scDatObj[[]], metaData, by = 'library')

#Remove likely contamination/doublets. should have low mito reads since this is nuclei
scDatObj<- subset(scDatObj, subset = nFeature_RNA > 300 & nFeature_RNA < 5000 & percent.mt < 10)

#Run through data processing and visualization before integration
scDatObj <- NormalizeData(scDatObj)
scDatObj <- FindVariableFeatures(scDatObj)
scDatObj <- ScaleData(scDatObj)
scDatObj <- RunPCA(scDatObj)
ElbowPlot(scDatObj, ndims = 50)
scDatObj <- FindNeighbors(scDatObj, dims = 1:30, reduction = "pca")

#Don't think resolution here matters since we rerun findclusters post integration. Tried with multiple
scDatObj <- FindClusters(scDatObj, resolution = 2, cluster.name = "unintegrated_clusters") 
scDatObj <- RunUMAP(scDatObj, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(scDatObj)

#Integrate data
#choosing integration method - rpca fastest, cca takes long
scDatObj_int <- IntegrateLayers(object = scDatObj, method = RPCAIntegration,
                                orig.reduction = "pca",  new.reduction = "integrated.rpca", 
                                verbose = TRUE)

#re-join layers after integration
scDatObj_int[["RNA"]] <- JoinLayers(scDatObj_int[["RNA"]])
ElbowPlot(scDatObj_int, ndims = 30)
scDatObj_int <- FindNeighbors(scDatObj_int, reduction = "integrated.rpca", dims = 1:20)
scDatObj_int <- FindClusters(scDatObj_int, resolution = 1)
scDatObj_int <- RunUMAP(scDatObj_int, dims = 1:20, reduction = "integrated.rpca", 
                        reduction.name = "umap.integrated")

#Visualize post integration
DimPlot(scDatObj_int, reduction = "umap.integrated", group.by = c("seurat_clusters"), label = TRUE)+
  NoLegend()+
  ggtitle('Integrated UMAP')

#remove cells with both sex genes
sexGenes <- scDatObj_int[['RNA']]$data[c('Xist', 'Eif2s3y'),]
sexGenePresence <- colSums(sexGenes > 0)
scDatObj_int$sexGenePresence <- case_when(sexGenePresence == 0 ~ 'None',
                                          sexGenePresence == 1 ~ 'One',
                                          sexGenePresence == 2 ~ 'Two')

table(scDatObj_int$sexGenePresence)
scDatObj_int$sexGenePresence <- subset(scDatObj_int$sexGenePresence, sexGenePresence != 'Two')

#SaveSeuratRds(scDatObj_int, './LGTVscCombined.rds')


# - - - - - - - - - - - - 
#### Cell annotation ####
# - - - - - - - - - - - - 


#Ex neurons
FeaturePlot(sn_integrated_dat, 'Slc17a7', reduction = 'umap.integrated')
FeaturePlot(sn_integrated_dat, 'Sv2b', reduction = 'umap.integrated')
FeaturePlot(sn_integrated_dat, 'Arpp21', reduction = 'umap.integrated')

#in neurons
FeaturePlot(sn_integrated_dat, 'Gad1', reduction = 'umap.integrated')
FeaturePlot(sn_integrated_dat, 'Gad2', reduction = 'umap.integrated')
FeaturePlot(sn_integrated_dat, 'Dlx6os1', reduction = 'umap.integrated')

#Microglia. Seems to align with uninfected section of micro/macro cluster
FeaturePlot(sn_integrated_dat, 'Csf1r', reduction = 'umap.integrated')
FeaturePlot(sn_integrated_dat, 'Tmem119', reduction = 'umap.integrated')
FeaturePlot(sn_integrated_dat, 'P2ry12', reduction = 'umap.integrated')
FeaturePlot(sn_integrated_dat, 'Cx3cr1', reduction = 'umap.integrated')

#Macrophage/monocyte markers
FeaturePlot(sn_integrated_dat, 'Ptprc', reduction = 'umap.integrated')
FeaturePlot(sn_integrated_dat, 'Ccr2', reduction = 'umap.integrated')
FeaturePlot(sn_integrated_dat, 'Lyz2', reduction = 'umap.integrated')
FeaturePlot(sn_integrated_dat, 'Ccr2', reduction = 'umap.integrated')

#Astrocytes
FeaturePlot(sn_integrated_dat, 'Aqp4', reduction = 'umap.integrated')
FeaturePlot(sn_integrated_dat, 'Rorb', reduction = 'umap.integrated')
FeaturePlot(sn_integrated_dat, 'Fgfr3', reduction = 'umap.integrated')
FeaturePlot(sn_integrated_dat, 'Gfap', reduction = 'umap.integrated')

#oligo
FeaturePlot(sn_integrated_dat, 'Mag', reduction = 'umap.integrated')
FeaturePlot(sn_integrated_dat, 'Mog', reduction = 'umap.integrated')
FeaturePlot(sn_integrated_dat, 'Plp1', reduction = 'umap.integrated')

#OPC
FeaturePlot(sn_integrated_dat, 'Pdgfra', reduction = 'umap.integrated')
FeaturePlot(sn_integrated_dat, 'Vcan', reduction = 'umap.integrated')

#VLMCs
FeaturePlot(sn_integrated_dat, 'Dcn', reduction = 'umap.integrated')
FeaturePlot(sn_integrated_dat, 'Col1a1', reduction = 'umap.integrated')

#Pericytes
FeaturePlot(sn_integrated_dat, 'Abcc9', reduction = 'umap.integrated')
FeaturePlot(sn_integrated_dat, 'Pdgfrb', reduction = 'umap.integrated')
FeaturePlot(sn_integrated_dat, 'Vtn', reduction = 'umap.integrated')

#ChP
FeaturePlot(sn_integrated_dat, 'Ttr', reduction = 'umap.integrated')
FeaturePlot(sn_integrated_dat, 'Aqp1', reduction = 'umap.integrated')
FeaturePlot(sn_integrated_dat, 'Kl', reduction = 'umap.integrated')

#Endothelial cells
FeaturePlot(sn_integrated_dat, 'Flt1', reduction = 'umap.integrated')
FeaturePlot(sn_integrated_dat, 'Cd93', reduction = 'umap.integrated')
FeaturePlot(sn_integrated_dat, 'Vwf', reduction = 'umap.integrated')
FeaturePlot(sn_integrated_dat, 'Pecam1', reduction = 'umap.integrated')

#NK cells / CD8 cells
FeaturePlot(sn_integrated_dat, 'Cd3e', reduction = 'umap.integrated')
FeaturePlot(sn_integrated_dat, 'Klrd1', reduction = 'umap.integrated')
FeaturePlot(sn_integrated_dat, 'Gzmb', reduction = 'umap.integrated')
FeaturePlot(sn_integrated_dat, 'Ncr1', reduction = 'umap.integrated')

#Ependymal cells
FeaturePlot(sn_integrated_dat, 'Cfap54', reduction = 'umap.integrated')
FeaturePlot(sn_integrated_dat, 'Nnat', reduction = 'umap.integrated')
FeaturePlot(sn_integrated_dat, 'Mia', reduction = 'umap.integrated')

#Export expression values from several clusters to be uploaded to 
#MapMyCells to double check labels
#Load data into mapmycells to check
clust22 <- subset(sn_integrated_dat, seurat_clusters == 22)
clust22[['RNA']]$counts %>% t() %>% 
  write.csv(file = './mapMyCellsData/clust22.csv', row.names = TRUE)


clust19 <- subset(sn_integrated_dat, seurat_clusters == 19)
clust19[['RNA']]$counts %>% t() %>% 
  write.csv(file = './mapMyCellsData/clust19.csv', row.names = TRUE)

clust11 <- subset(sn_integrated_dat, seurat_clusters == 11)
clust11[['RNA']]$counts %>% t() %>% 
  write.csv(file = './mapMyCellsData/clust11.csv', row.names = TRUE)

clust0_31_29 <- subset(sn_integrated_dat, seurat_clusters %in% c(0, 29, 31))
clust0_31_29[['RNA']]$counts %>% t() %>%
  write.csv(file = './mapMyCellsData/clust0_29_31.csv', row.names = TRUE)

#Read in MapMyCells data and label 
cluster22Map <- readr::read_csv("~/Documents/ÖverbyLab/single_nuclei_proj/mapMyCellsData/clust22_snData_mapMy/clust22csv_10xWholeMouseBrain(CCN20230722)_HierarchicalMapping_UTC_1764931715314.csv", 
                                skip = 4)
table(cluster22Map$class_name)
neuron_labels <- cluster22Map %>% dplyr::select(cell_id, class_name) %>% 
  dplyr::mutate(neu_type = case_when(grepl('GABA', class_name) ~ 'In Neurons',
                                     grepl('Glut', class_name) ~ 'Ex Neurons',
                                     .default = 'unknown'))

cluster19Map <- readr::read_csv("~/Documents/ÖverbyLab/single_nuclei_proj/mapMyCellsData/clust19_snData_mapMy/clust19csv_10xWholeMouseBrain(CCN20230722)_HierarchicalMapping_UTC_1764940170384.csv", 
                                skip = 4)
table(cluster19Map$class_name)
neuron_labels_19 <- cluster19Map %>% dplyr::select(cell_id, class_name) %>% 
  dplyr::mutate(neu_type = case_when(grepl('GABA', class_name) ~ 'In Neurons',
                                     grepl('Glut', class_name) ~ 'Ex Neurons',
                                     .default = 'unknown'))

cluster11Map <- readr::read_csv("~/Documents/ÖverbyLab/single_nuclei_proj/mapMyCellsData/clust11_snData_mapMy/clust11csv_10xWholeMouseBrain(CCN20230722)_HierarchicalMapping_UTC_1765206280627.csv", 
                                skip = 4)
table(cluster11Map$class_name)
neuron_labels_11 <- cluster11Map %>% dplyr::select(cell_id, class_name) %>% 
  dplyr::mutate(neu_type = case_when(grepl('GABA', class_name) ~ 'In Neurons',
                                     grepl('Glut', class_name) ~ 'Ex Neurons',
                                     grepl('Astro', class_name) ~ 'Astrocytes',
                                     .default = 'unknown'))

cluster0_29_31Map <- readr::read_csv("~/Documents/ÖverbyLab/single_nuclei_proj/mapMyCellsData/clust0_29_31csv_10xWholeMouseBrain(CCN20230722)_HierarchicalMapping_UTC_1773656729060/clust0_29_31csv_10xWholeMouseBrain(CCN20230722)_HierarchicalMapping_UTC_1773656729060.csv", 
                                     skip = 4)
table(cluster0_29_31Map$class_name)

#Label cells
#Macrophages and microglia clustering together,
#PBS samples appear to be microglia and infected Macrophages
sn_integrated_dat$manualAnnotation <- 
  case_when(sn_integrated_dat$seurat_clusters %in% c(2, 12, 30) ~ 'Oligo',
            sn_integrated_dat$seurat_clusters %in% c(16) ~ 'OPC',
            sn_integrated_dat$seurat_clusters %in% c(18, 10, 8, 19) ~ 'In Neurons',
            sn_integrated_dat$seurat_clusters %in% c(1, 5, 6, 7, 9, 25, 28, 13, 27, 23, 15, 14, 4,
                                                     17, 26, 11) ~ 'Ex Neurons',
            sn_integrated_dat$seurat_clusters %in% c(3, 21) & sn_integrated_dat$new_inf != 'LGTV' ~ 'Microglia',
            sn_integrated_dat$seurat_clusters %in% c(3, 21) & sn_integrated_dat$new_inf == 'LGTV' ~ 'Macrophages',
            sn_integrated_dat$seurat_clusters %in% c(0) ~ 'Astrocytes',
            sn_integrated_dat$seurat_clusters %in% c() ~ 'ChP',
            sn_integrated_dat$seurat_clusters %in% c(24) ~ 'Endothelial',
            sn_integrated_dat$seurat_clusters %in% c(29) ~ 'Peri',
            sn_integrated_dat$seurat_clusters %in% c(20) ~ 'VLMCs',
            .default = 'unknown')

#Change astrocytes with high ttr to ChP cells
sn_integrated_dat$ttr_exp <-  FetchData(object = sn_integrated_dat, vars = c("Ttr"), layer = "data")
sn_integrated_dat[[]][sn_integrated_dat$manualAnnotation == 'Astrocytes' & sn_integrated_dat$ttr_exp > 1,]$manualAnnotation = 'ChP'

sn_integrated_dat$manualAnnotation[neuron_labels$cell_id] = neuron_labels$neu_type
sn_integrated_dat$manualAnnotation[neuron_labels_19$cell_id] = neuron_labels_19$neu_type
sn_integrated_dat$manualAnnotation[neuron_labels_11$cell_id] = neuron_labels_11$neu_type

sn_integrated_dat <- subset(sn_integrated_dat, manualAnnotation != 'unknown')


# - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#### Analyze relevant chemokines across cell types ####
# - - - - - - - - - - - - - - - - - - - - - - - - - - - 

sn_integrated_dat$treatment_celltype <- paste(sn_integrated_dat$infected, sn_integrated_dat$manualAnnotation, sep = '_')
sn_integrated_dat_wt <- subset(sn_integrated_dat, new_genotype %in% c('wt', 'wt (same)'))

mono_chemokines <- c('Ccl2', 'Ccl5', 'Ccl7')

mono_dot_dat <- DotPlot(sn_integrated_dat_wt, features = genes_of_interest, group.by = 'treatment_celltype', scale = FALSE)$data
mono_dot_dat <- mono_dot_dat %>% 
  tidyr::separate(id, into = c('infected', 'celltype'), sep = '_')

mono_dot_dat$celltype = factor(mono_dot_dat$celltype, levels = rev(c('Astrocytes', 'ChP', 'Endothelial', 'Ex Neurons',
                                                                   'In Neurons', 'Microglia', 'Oligo', 'OPC', 'Peri', 'VLMCs',
                                                                   'Cd8+NK', 'Macrophages')))

ggplot(mono_dot_dat, aes(x = features.plot, y = celltype))+
  geom_point(aes(fill = avg.exp.scaled, size = pct.exp), pch = 21)+
  facet_grid(cols = vars(infected), scales = "free",
             labeller = as_labeller(c("FALSE" = "Uninfected", "TRUE" = "LGTV")))+
  #coord_flip()+
  scale_size_continuous(range = c(0,9), limits = c(0,100))+
  scale_fill_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white"), 
                       values = c(1.0,0.7,0.4,0),
                       limits = c(0,3.5))+
  ggtitle('Single-nuclei chemokines genes')+
  theme_bw()+
  theme(panel.grid = element_blank(),
        #panel.border = element_rect(colour = "white", fill = NA),
        panel.spacing = unit(0, "line"),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90))+
  ylab('')+
  xlab('')


# - - - - - - - - - - - - - - - - - - - -
#### Chemokines across neuron types ####
# - - - - - - - - - - - - - - - - - - - -



