library(Seurat)
library(dplyr)
library(scDblFinder)
library(readr)
library(ggplot2)
library(gridExtra)
library(Matrix)

#Retrieve data and read into R
parseOutput <- list.files("FilteredParseOutput/")

ParseMatricies <- lapply(parseOutput, FUN = function(x){
  ReadParseBio(paste0("FilteredParseOutput/", x))
})

#Convert data to Seurat objects 
ParseSeuratObj <- lapply(ParseMatricies, FUN = function(x){
  CreateSeuratObject(counts = x, project = "Langat", min.cells = 3)
})

#Combine Parse sublibraries into one seurat object before integration
toMerge <- c(ParseSeuratObj[[2]], ParseSeuratObj[[3]], ParseSeuratObj[[4]],
             ParseSeuratObj[[5]], ParseSeuratObj[[6]], ParseSeuratObj[[7]],
             ParseSeuratObj[[8]])


ParseSeuratObj <- merge(ParseSeuratObj[[1]], y = toMerge, 
                        add.cell.ids = paste0('seuObj', seq(1,8)), project = "ParseLangat")

#Quick section to look at raw data ifn values (set min.cells to 0 for this only)
countDats <- list(ParseSeuratObj[['RNA']]$counts.1, ParseSeuratObj[['RNA']]$counts.2,
                  ParseSeuratObj[['RNA']]$counts.3, ParseSeuratObj[['RNA']]$counts.4,
                  ParseSeuratObj[['RNA']]$counts.5, ParseSeuratObj[['RNA']]$counts.6,
                  ParseSeuratObj[['RNA']]$counts.7, ParseSeuratObj[['RNA']]$counts.8)

#Label mitochondrial gene expression
ParseSeuratObj[["percent.mt"]] <- PercentageFeatureSet(ParseSeuratObj, pattern = "^mt-")

#Look at data features
Idents(ParseSeuratObj) <- "all"  #Stops violin plot from grouping by seurat cluster
VlnPlot(ParseSeuratObj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3, pt.size = 0) 

#Filter out cells with high rna features and mitochondrial gene counts
ParseSeuratObj <- subset(ParseSeuratObj, nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 5)

#Use Parse barcoding plate to attach metadata to cells
ParseBarcodePlate <- read_csv("~/Documents/ÖverbyLab/ParseBarcodePlate.csv")
ParseBarcodePlate$orig.ident <- c(paste0('0', as.character(seq(1:9))), seq(10,48))
ParseSeuratObj[[]] <- left_join(ParseSeuratObj[[]], ParseBarcodePlate, by = 'orig.ident')

#Label sublibraries 
ParseSeuratObj$subLib <- substr(colnames(ParseSeuratObj), start = 1, stop = 7)

#Standard seurat processing before integration
#Normalize data
ParseSeuratObj <- NormalizeData(ParseSeuratObj)

#Scale data
ParseSeuratObj <- FindVariableFeatures(ParseSeuratObj)
ParseSeuratObj <- ScaleData(ParseSeuratObj)

#PCA and clustering functions
ParseSeuratObj <- RunPCA(ParseSeuratObj)
ElbowPlot(ParseSeuratObj, ndims = 50)
ParseSeuratObj <- FindNeighbors(ParseSeuratObj, dims = 1:30, reduction = "pca")

#Clusters and UMAP are run post integration so these are not important, just visualizing before integration
ParseSeuratObj <- FindClusters(ParseSeuratObj, resolution = 2, cluster.name = "unintegrated_clusters") 
ParseSeuratObj <- RunUMAP(ParseSeuratObj, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(ParseSeuratObj)

#Create SCE object to view doublets
Parse_sce <- ParseSeuratObj %>% JoinLayers() %>%  as.SingleCellExperiment()

#Use sublibraries as samples for scdblfinder
Parse_sce <- scDblFinder(Parse_sce, dbr.sd=1, samples = 'subLib')
table(Parse_sce$scDblFinder.class)
ParseSeuratObj$scDblFinderLabel <- Parse_sce$scDblFinder.class
DimPlot(ParseSeuratObj, group.by = 'scDblFinderLabel')

#Doublets are filtered out later in 02_doublet_removal.R script

#Integrate data
#choosing integration method - rpca fastest, cca takes 5+ hours
ParseSeuratObj_int <- IntegrateLayers(object = ParseSeuratObj, method = RPCAIntegration,
                                      orig.reduction = "pca",  new.reduction = "integrated.rpca", 
                                      verbose = TRUE)

#re-join layers after integration
ParseSeuratObj_int[["RNA"]] <- JoinLayers(ParseSeuratObj_int[["RNA"]])

#Recluster and view integrated data
ParseSeuratObj_int <- FindNeighbors(ParseSeuratObj_int, reduction = "integrated.rpca", dims = 1:30)
ParseSeuratObj_int <- FindClusters(ParseSeuratObj_int, resolution = 1)
ParseSeuratObj_int <- RunUMAP(ParseSeuratObj_int, dims = 1:30, reduction = "integrated.rpca", 
                              reduction.name = "umap.integrated")

intUMAP <- DimPlot(ParseSeuratObj_int, reduction = "umap.integrated", group.by = c("seurat_clusters"), label = TRUE)+
  NoLegend()+
  ggtitle('Integrated UMAP')
unIntUmap <- DimPlot(ParseSeuratObj_int, reduction = "umap.unintegrated", group.by = c("seurat_clusters"), label = TRUE)+
  NoLegend()+
  ggtitle('Unintegrated UMAP')

grid.arrange(intUMAP, unIntUmap, ncol=2)

#remove cells that have both male and female genes, as this indicates doublets
sexGenes <- ParseSeuratObj_int[['RNA']]$data[c('Xist', 'Eif2s3y'),]
sexGenePresence <- colSums(sexGenes > 0)
ParseSeuratObj_int$sexGenePresence <- case_when(sexGenePresence == 0 ~ 'None',
                                                sexGenePresence == 1 ~ 'One',
                                                sexGenePresence == 2 ~ 'Two')

ParseSeuratObj_int <- subset(ParseSeuratObj_int, sexGenePresence != 'Two')

#SaveSeuratRds(ParseSeuratObj_int, "./data/FilteredRpcaIntegratedDat.rds")
