#Packages and functions
source('./scripts/langatFunctions.R')

#Load data
ParseSeuratObj_int <- LoadSeuratRds("./data/FilteredRpcaIntegratedDatNoDoublets.rds") 
ParseSeuratObj_int$hasVirus = ifelse(ParseSeuratObj_int$virusCountPAdj >= 10, 1, 0)

#Check data
newCols <-  c(brewer.pal(12, 'Paired'), '#99FFE6', '#8A5730', '#18662E',  '#FF8AEF','#737272')
DimPlot(ParseSeuratObj_int, label = FALSE, group.by = 'manualAnnotation', reduction = 'umap.integrated',
        cols = newCols)

wt_cerebrum <- subset(ParseSeuratObj_int, Treatment %in% c('PBS', 'rLGTV') & Organ == 'Cerebrum')
#Immune cell subset
#Need to finish labelling last mo/microglia cluster
immune <- subset(wt_cerebrum, manualAnnotation %in% c('Microglia', 'Macrophage/Monocytes', 'B Cells',
                                                             'T cells', 'Granulocytes', 'Nk cells'))

#ReUMAP
immune <- prepSeuratObj(immune)
ElbowPlot(immune, ndims = 40)

immune <- prepUmapSeuratObj(immune, nDims = 30, reductionName = 'immune.umap')

dimPlot(immune, reduction = 'immune.umap')