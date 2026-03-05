#Merging our data with yang's
#Could subset celltypes to make this faster if it's really slow

#Load our data
#Load data
ParseSeuratObj_int <- LoadSeuratRds("~/Documents/ÖverbyLab/data/FilteredRpcaIntegratedDatNoDoublets.rds") 

#Plot data
newCols <-  c(brewer.pal(12, 'Paired'), '#99FFE6', '#CE99FF', '#18662E','#737272',  '#FF8AEF')
newCols[11] =  '#FF8AEF'
DimPlot(ParseSeuratObj_int, label = FALSE, group.by = 'manualAnnotation', reduction = 'umap.integrated',
        cols = newCols)

#Load in Yang data
yang_data <- LoadSeuratRds("~/Documents/ÖverbyLab/yang_public_data/yang_rpca_integrated_obj.RDS")
newCols <-  c(brewer.pal(12, 'Paired'), '#99FFE6', '#CE99FF', '#18662E','#737272',  '#FF8AEF')
newCols[11] =  '#FF8AEF'

#Join layers for yang to match parse
yang_data <- JoinLayers(yang_data)

#Plot data
DimPlot(yang_data,reduction = "umap.rpca", label = FALSE, group.by = 'manualAnnotation', cols = newCols)

#Create celltype subsets to integrate
parse_endo <- subset(ParseSeuratObj_int, manualAnnotation == 'Endothelial')
parse_macro <- subset(ParseSeuratObj_int, manualAnnotation == 'Macrophage/Monocytes')

yang_endo <- subset(yang_data, manualAnnotation == 'Endothelial')
yang_macro <- subset(yang_data, manualAnnotation == 'Macro/Mono')

#monocyte cells combined
macro.combined <- merge(parse_macro, y = yang_macro, add.cell.ids = c("parse", "yang"), project = "merged_macros")

macro.combined <- NormalizeData(macro.combined)
macro.combined <- FindVariableFeatures(macro.combined)
macro.combined <- ScaleData(macro.combined)
macro.combined <- RunPCA(macro.combined)

#Need to allow greater ram usage to run pca integration
options(future.globals.maxSize = 10000 * 1024^2)

#Integrate data
combined_macros <- IntegrateLayers(
  object = macro.combined, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca",
  verbose = FALSE
)

#Reset max ram to not accidentally use a bunch later without realizing
options(future.globals.maxSize = 500 * 1024^2)

ElbowPlot(combined_macros, ndims = 30)

#Prepare data for plotting
combined_macros <- FindNeighbors(combined_macros, reduction = "integrated.rpca", dims = 1:20)
combined_macros <- FindClusters(combined_macros, resolution = 2, cluster.name = "rpca_clusters")
combined_macros$manualAnnotation = 'macro/mono'
combined_macros$dataset <- gsub("_.+", ' ', colnames(combined_macros))

DimPlot(combined_macros, label = FALSE, group.by = 'dataset', reduction = 'integrated.rpca',
        cols = newCols)

