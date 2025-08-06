#Load in libraries
library(Seurat)
library(SingleR)
library(celldex)

#List out our data files, and read them into R
parseOutput <- list.files("./unfilteredParseOutput")

ParseMatricies <- lapply(parseOutput, FUN = function(x){
  ReadParseBio(paste0("unfilteredParseOutput/", x))
})

#Convert data to Seurat objects  
ParseSeuratObj <- lapply(ParseMatricies, FUN = function(x){
  CreateSeuratObject(counts = x, project = "Langat", 
                     min.cells = 1, min.features = 20)
})

#Can also combine this data using split-pipe combined, instead of doing it here
toMerge <- c(ParseSeuratObj[[2]], ParseSeuratObj[[3]], ParseSeuratObj[[4]],
             ParseSeuratObj[[5]], ParseSeuratObj[[6]], ParseSeuratObj[[7]],
             ParseSeuratObj[[8]])

#Merge data before integration
ParseSeuratObj <- merge(ParseSeuratObj[[1]], y = toMerge, 
                  add.cell.ids = paste0('seuObj', seq(1,8)), project = "ParseLangat")

#Use Parse barcoding plate to label cells
#Going to try just left joining and seeing if virus expression matches up
ParseBarcodePlate <- read_csv("ParseBarcodePlate.csv")
ParseBarcodePlate$orig.ident <- c(paste0('0', as.character(seq(1:9))), seq(10,48))
ParseSeuratObj[[]] <- left_join(ParseSeuratObj[[]], ParseBarcodePlate, by = 'orig.ident')

#Assign viral reads from labels
#List files with viral reads
virusCountFiles <- list.files('./viralReadCountsUnadjusted/')
virusCountFiles <- virusCountFiles[grep('P36207', virusCountFiles)]

viralCountsUnadjusted <- lapply(virusCountFiles, FUN = function(x){
  subLibViralCounts <- read_table(paste0("viralReadCountsUnadjusted/", x), 
                                 col_names = FALSE)
  subLibViralCounts
})

#For each sublibrary create a df with counts of virus per cell barcode
#Also adds a sublibrary column for joining purposes afterwards
for(i in 1:length(viralCountsUnadjusted)){
  dat = viralCountsUnadjusted[[i]]
  dat$subLib <- paste0('seuObj', as.character(i))
  dat$cell <- substr(dat$X2, start = 6, stop = 14)
  dat <- dat[c(1,3,4)]
  colnames(dat)[1] = 'virusCount'
  viralCountsUnadjusted[[i]] <- dat
}

viralCountsUnadjusted <- bind_rows(viralCountsUnadjusted)

ParseSeuratObj$cell <- substr(colnames(ParseSeuratObj), start = 9, stop = 16)
ParseSeuratObj$subLib <- substr(colnames(ParseSeuratObj), start = 1, stop = 7)

ParseSeuratObj[[]] <- left_join(ParseSeuratObj[[]], viralCountsUnadjusted, by = c('subLib', 'cell'))
ParseSeuratObj$virusCount[is.na(ParseSeuratObj$virusCount)] = 0

wellMap <- data.frame(well = c(paste0('A', seq(1,12)), paste0('B', seq(1,12)),
                               paste0('C', seq(1,12)),   paste0('D', seq(1,12))))

#Plot viral counts vs well/treatment. Looks right
ggplot(ParseSeuratObj[[]], aes(x = orig.ident, y = virusCount, col = Treatment))+
  geom_point() +
  scale_x_discrete(labels= wellMap$well)+
  theme(axis.text.x = element_text(angle = 90))

#Check percent of cells in each group expressing virus


#Label mitochondrial gene expression
ParseSeuratObj[["percent.mt"]] <- PercentageFeatureSet(ParseSeuratObj, pattern = "^mt-")

#Subset for integration for now but should reexamine cutoff
ParseSeuratObj <- subset(ParseSeuratObj, nFeature_RNA > 200)

#Look at data features
Idents(ParseSeuratObj) <- "all"  #Stops violin plot from grouping by seurat cluster
VlnPlot(ParseSeuratObj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3, pt.size = 0)

ParseSeuratObj <- subset(ParseSeuratObj, nFeature_RNA < 7500 & nCount_RNA < 100000)

#Run through data processing and visualization before integration
ParseSeuratObj <- NormalizeData(ParseSeuratObj)
ParseSeuratObj <- FindVariableFeatures(ParseSeuratObj)
ParseSeuratObj <- ScaleData(ParseSeuratObj)
ParseSeuratObj <- RunPCA(ParseSeuratObj)
ElbowPlot(ParseSeuratObj, ndims = 40)
ParseSeuratObj <- FindNeighbors(ParseSeuratObj, dims = 1:30, reduction = "pca")
ParseSeuratObj <- FindClusters(ParseSeuratObj, resolution = 2, cluster.name = "unintegrated_clusters")  
ParseSeuratObj <- RunUMAP(ParseSeuratObj, dims = 1:30, reduction = "pca")
DimPlot(ParseSeuratObj, group.by = 'orig.ident')

#Integrate data
#Need to choose integration method - rpca fastest
#Can look at running in parellel to increase speed. CCA takes 5+ hours. RPCA much much faster
ParseSeuratObj_int <- IntegrateLayers(object = ParseSeuratObj, method = CCAIntegration,
                                      orig.reduction = "pca",  new.reduction = "integrated.cca", 
                                      verbose = TRUE)

# re-join layers after integration
ParseSeuratObj_int[["RNA"]] <- JoinLayers(ParseSeuratObj_int[["RNA"]])

ParseSeuratObj_int <- FindNeighbors(ParseSeuratObj_int, reduction = "integrated.cca", dims = 1:30)
ParseSeuratObj_int <- FindClusters(ParseSeuratObj_int, resolution = 1)

ParseSeuratObj_int <- RunUMAP(ParseSeuratObj_int, dims = 1:30, reduction = "integrated.cca")
DimPlot(ParseSeuratObj_int, reduction = "umap", group.by = c("seurat_clusters"), label = TRUE)


####Annotate cell types ####
#Will try automatic annotation
#Mouse rna seq reference
ref <- fetchReference("mouse_rnaseq", "2024-02-26")
pred <- SingleR(test=ParseSeuratObj_int[['RNA']]$data, ref=ref, labels=ref$label.main)
ParseSeuratObj_int$singleR_labels <- pred$labels
plotScoreHeatmap(pred)
DimPlot(ParseSeuratObj_int, group.by = 'singleR_labels', label = TRUE)

#### Look at canonical gene markers for annotation ####
#Paper has some canonical markers for cell types
#https://umu.diva-portal.org/smash/get/diva2:1897514/FULLTEXT01.pdf

#This paper https://link.springer.com/article/10.1007/s10571-021-01159-3
#Heterogeneity and Molecular Markers for CNS Glial Cells Revealed by Single-Cell Transcriptomics also 
#has lists of markers

#Some neuron markers
FeaturePlot(ParseSeuratObj_int, 'Snap25')
FeaturePlot(ParseSeuratObj_int, 'Pcp2')
FeaturePlot(ParseSeuratObj_int, 'Rbfox3')
FeaturePlot(ParseSeuratObj_int, 'Dpp10')

FeaturePlot(ParseSeuratObj_int, 'Foxp1')

#Endothelial cells
FeaturePlot(ParseSeuratObj_int, 'Flt1')

#Astrocytes
FeaturePlot(ParseSeuratObj_int, 'Gfap')
FeaturePlot(ParseSeuratObj_int, 'Aqp4')
FeaturePlot(ParseSeuratObj_int, 'Fgfr3')
FeaturePlot(ParseSeuratObj_int, 'Aldh1l1')
FeaturePlot(ParseSeuratObj_int, 'Slc1a3')

#Microglia
FeaturePlot(ParseSeuratObj_int, 'Ctss')
FeaturePlot(ParseSeuratObj_int, 'Csf1r')
FeaturePlot(ParseSeuratObj_int, 'Cx3cr1')
FeaturePlot(ParseSeuratObj_int, 'C1qa')
FeaturePlot(ParseSeuratObj_int, 'Tmem119')
FeaturePlot(ParseSeuratObj_int, 'P2ry12')

#Oligo
FeaturePlot(ParseSeuratObj_int, 'Mag')
FeaturePlot(ParseSeuratObj_int, 'Mog')

#OPC
FeaturePlot(ParseSeuratObj_int, 'Pdgfra')
FeaturePlot(ParseSeuratObj_int, 'Cspg4')

#Macrophage markers
FeaturePlot(ParseSeuratObj_int, 'Ptprc')


#Look at specific clusters top markers to confirm cell types
#Cluster 0 uprefulated w Chil3, Ms4a8a, Saa3, GM15056. Looks like macrophages
possibleMacrophages <- FindMarkers(ParseSeuratObj_int, ident.1 = 0, group.by = 'seurat_clusters', only.pos = TRUE)



