#Script analyzing cerebellum

#Load packages
library(Seurat)

#Load in data
ParseSeuratObj_int <- LoadSeuratRds("./FilteredRpcaIntegratedDat.rds")

#Subset to cerebellum data
cerebellumObj <- subset(ParseSeuratObj_int, Organ == 'Cerebellum')
cerebellumObj$isInfected = ifelse(cerebellumObj$virusCountPAdj>0, 'yes', 'no')
#Rerun through data processing and visualization
cerebellumObj <- NormalizeData(cerebellumObj)
cerebellumObj <- FindVariableFeatures(cerebellumObj)
cerebellumObj <- ScaleData(cerebellumObj)
cerebellumObj <- RunPCA(cerebellumObj)
ElbowPlot(cerebellumObj, ndims = 40)
cerebellumObj <- FindNeighbors(cerebellumObj, dims = 1:30, reduction = "pca")
cerebellumObj <- FindClusters(cerebellumObj, resolution = 2, cluster.name = "unintegrated_clusters")  
cerebellumObj <- RunUMAP(cerebellumObj, dims = 1:30, reduction = "pca", reduction.name = "umap")

DimPlot(cerebellumObj, reduction = 'umap', group.by = 'singleR_labels', label = TRUE)
DimPlot(cerebellumObj, reduction = 'umap', group.by = 'Timepoint', label = TRUE)
DimPlot(cerebellumObj, reduction = 'umap', group.by = 'Treatment', label = TRUE)


#Timepoint differences
#



