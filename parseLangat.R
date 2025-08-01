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
                     min.cells = 3, min.features = 200)
})

toMerge <- c(ParseSeuratObj[[2]], ParseSeuratObj[[3]], ParseSeuratObj[[4]],
             ParseSeuratObj[[5]], ParseSeuratObj[[6]], ParseSeuratObj[[7]],
             ParseSeuratObj[[8]])
#Merge data before integration
ParseSeuratObj <- merge(ParseSeuratObj[[1]], y = toMerge, 
                  add.cell.ids = seq(1,8), project = "ParseLangat")

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
#Need to choose integration method
ParseSeuratObj_int <- IntegrateLayers(object = ParseSeuratObj, method = CCAIntegration,
                                      orig.reduction = "pca", new.reduction = "integrated.cca", verbose = FALSE)

# re-join layers after integration
ifnb[["RNA"]] <- JoinLayers(ifnb[["RNA"]])

ifnb <- FindNeighbors(ifnb, reduction = "integrated.cca", dims = 1:30)
ifnb <- FindClusters(ifnb, resolution = 1)

ifnb <- RunUMAP(ifnb, dims = 1:30, reduction = "integrated.cca")
DimPlot(ifnb, reduction = "umap", group.by = c("stim", "seurat_annotations"))
