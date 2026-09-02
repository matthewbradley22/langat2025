library(Seurat)
library(RColorBrewer)
library(UCell)

source("~/Documents/ÖverbyLab//scripts/langatFunctions.R")

#Load data
ParseSeuratObj_int <- LoadSeuratRds("~/Documents/ÖverbyLab/data/FilteredRpcaIntegratedDatNoDoublets.rds") 
ParseSeuratObj_int$hasVirus = ifelse(ParseSeuratObj_int$virusCountPAdj >= 10, 1, 0)

#Check data
newCols <-  c(brewer.pal(12, 'Paired'), '#99FFE6', '#CE99FF', '#18662E','#737272',  '#FF8AEF')
newCols[11] =  '#FF8AEF'
DimPlot(ParseSeuratObj_int, label = FALSE, group.by = 'manualAnnotation', reduction = 'umap.integrated',
        cols = newCols)

#Subset to just monocytes/microglia since dcs will likely be similar
mono_micro <- subset(ParseSeuratObj_int, manualAnnotation %in% c('Microglia', 'Macrophage/Monocytes'))

mono_micro <- prepSeuratObj(mono_micro)
ElbowPlot(mono_micro, ndims = 40)
mono_micro <- prepUmapSeuratObj(mono_micro, nDims = 20, reductionName = 'mono.micro.umap',
                                resolution_value = 0.5)

DimPlot(mono_micro, reduction = 'mono.micro.umap', label = FALSE, cols = newCols)
DimPlot(mono_micro, reduction = 'mono.micro.umap', label = FALSE, cols = newCols, group.by = 'manualAnnotation')

#Marker list thrown together, no specific source
#DC markers
micro_markers <- c('Tmem119', 'P2ry12', 'Sall1', 'Slc2a5', 'Hexb', 'Siglech', 'Gpr34', 'Olfml3', 'Fcrls', 'Csf1r')
mono_micro <- AddModuleScore_UCell(mono_micro, features = list(micro_markers), name = 'micro_list')

FeaturePlot(mono_micro, features = 'signature_1micro_list', reduction = 'mono.micro.umap')

infil_macro_markers <- c('Ccr2', 'Ly6c2', 'Plac8', 'Vim', 'Chil3', 'Itgal')
mono_micro <- AddModuleScore_UCell(mono_micro, features = list(infil_macro_markers), name = 'infil_macro')

FeaturePlot(mono_micro, features = 'signature_1infil_macro', reduction = 'mono.micro.umap')

#Small population of border macrophages
FeaturePlot(mono_micro, features = 'Mrc1', reduction = 'mono.micro.umap')

#DCs - should lack Tmem119, P2ry12
FeaturePlot(mono_micro, features = 'Tmem119', reduction = 'mono.micro.umap')
FeaturePlot(mono_micro, features = 'P2ry12', reduction = 'mono.micro.umap')

pdc_sig <- c('Siglech', 'Bst2', 'Klk1', 'Ccr9', 'Tcf4', 'Bcl11a', 'Ly6d', 'Cox6a2', 'Spib')
mono_micro <- AddModuleScore_UCell(mono_micro, features = list(pdc_sig), name = 'pdc_sig')

#Weaker signature and generally overlaps with tmem119 positive cells
FeaturePlot(mono_micro, features = 'signature_1pdc_sig', reduction = 'mono.micro.umap')

#Some genes from https://www.nature.com/articles/s41593-019-0393-4
FeaturePlot(mono_micro, features = 'Ccr9', reduction = 'mono.micro.umap')
FeaturePlot(mono_micro, features = 'Pacsin1', reduction = 'mono.micro.umap')
FeaturePlot(mono_micro, features = 'Siglech', reduction = 'mono.micro.umap') #overlaps with tmem119

#Some other dendritic/mac markers. Itgax expression is interesting
FeaturePlot(mono_micro, features = 'Itgax', reduction = 'mono.micro.umap')
FeaturePlot(mono_micro, features = 'Cd209a', reduction = 'mono.micro.umap')


