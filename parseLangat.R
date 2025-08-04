#are these samples pooled and how many samples are there?
#Load in libraries
library(Seurat)

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



#Label mitochondrial gene expression
ParseSeuratObj[["percent.mt"]] <- PercentageFeatureSet(ParseSeuratObj, pattern = "^mt-")

#Subset for integration for now but should reexamine cutoff
ParseSeuratObj <- subset(ParseSeuratObj, nFeature_RNA > 200)

#Look at data features
Idents(ParseSeuratObj) <- "all"  #Stops violin plot from grouping by seurat cluster
VlnPlot(ParseSeuratObj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3, pt.size = 0)

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
#Can look at running in parellel to increase speed. Also at running v4 integration
ParseSeuratObj_int <- IntegrateLayers(object = ParseSeuratObj, method = CCAIntegration,
                                      orig.reduction = "pca",  new.reduction = "integrated.cca", 
                                      verbose = TRUE)

# re-join layers after integration
ParseSeuratObj_int[["RNA"]] <- JoinLayers(ParseSeuratObj_int[["RNA"]])

ParseSeuratObj_int <- FindNeighbors(ParseSeuratObj_int, reduction = "integrated.rpca", dims = 1:30)
ParseSeuratObj_int <- FindClusters(ParseSeuratObj_int, resolution = 1)

ParseSeuratObj_int <- RunUMAP(ParseSeuratObj_int, dims = 1:30, reduction = "integrated.rpca")
DimPlot(ParseSeuratObj_int, reduction = "umap", group.by = c("seurat_clusters"))

#Look at canonical gene markers
#Not sure best cell markers, this thread offers some https://www.biostars.org/p/9571136/

#This paper https://link.springer.com/article/10.1007/s10571-021-01159-3
#Heterogeneity and Molecular Markers for CNS Glial Cells Revealed by Single-Cell Transcriptomics has
#lists of markers

#Some neuron markers
FeaturePlot(ParseSeuratObj_int, 'Snap25')

#Endothelial cells
FeaturePlot(ParseSeuratObj_int, 'Flt1')

#Astrocytes
FeaturePlot(ParseSeuratObj_int, 'Gfap')
FeaturePlot(ParseSeuratObj_int, 'Aldh1l1')
FeaturePlot(ParseSeuratObj_int, 'Slc1a3')

#Microglia
FeaturePlot(ParseSeuratObj_int, 'C1qa')
FeaturePlot(ParseSeuratObj_int, 'Tmem119')
FeaturePlot(ParseSeuratObj_int, 'P2ry12')




