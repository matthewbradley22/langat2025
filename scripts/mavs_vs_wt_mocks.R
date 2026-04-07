#Load libraries 
library(dplyr)
library(Seurat)
library(ggplot2)
library(RColorBrewer)

### SUBSET TO JUST WILDTYPE ###
#Source function
source('~/Documents/ÖverbyLab/scripts/langatFunctions.R')

#Read in data
ParseSeuratObj_int <- LoadSeuratRds("~/Documents/ÖverbyLab/data/FilteredRpcaIntegratedDatNoDoublets.rds") 

#Check data
newCols <-  c(brewer.pal(12, 'Paired'), '#99FFE6', '#CE99FF', '#18662E','#737272',  '#FF8AEF')
newCols[11] =  '#FF8AEF'
DimPlot(ParseSeuratObj_int, label = FALSE, group.by = 'manualAnnotation', reduction = 'umap.integrated',
        cols = newCols)

#Subset to just mock samples
mock_samples <- subset(ParseSeuratObj_int, Treatment == 'PBS')

#Prep data and plot
mock_samples <- prepSeuratObj(mock_samples)
ElbowPlot(mock_samples, ndims = 40)
mock_samples <- prepUmapSeuratObj(mock_samples, nDims = 30, reductionName = 'mock_umap')

DimPlot(mock_samples, reduction = 'mock_umap', group.by = 'Genotype')

#Compare mock astrocytes in cerebrum
mock_astros <- subset(mock_samples, manualAnnotation == 'Astrocytes' & Organ == 'Cerebrum')

#Prep data and plot
mock_astros <- prepSeuratObj(mock_astros)
ElbowPlot(mock_astros, ndims = 40)
mock_astros <- prepUmapSeuratObj(mock_astros, nDims = 20, reductionName = 'mock_astro_umap')

DimPlot(mock_astros, reduction = 'mock_astro_umap', group.by = 'Genotype')

#Are isgs/ccls upregulated in mavs mock... Doesn't look like it
mock_astros_markers <- FindAllMarkers(mock_astros, test.use = 'MAST', only.pos = TRUE, group.by = 'Genotype')
mavs_markers <- dplyr::filter(mock_astros_markers, avg_log2FC > 0.1 & p_val_adj < 0.01 & cluster == 'IPS1')
mavs_paths <- gprofiler2::gost(query = mavs_markers$gene, organism = 'mmusculus', evcodes = TRUE, sources = c('GO:BP', 'KEGG'))

DotPlot(mock_astros, features = c(paste0('Ccl', 1:16)), group.by = 'Genotype', scale = FALSE)+
  theme(axis.text.x = element_text(angle = 90))

DotPlot(mock_astros, features = c('Rsad2', 'Irf7', 'Ifitm3', 'Oas1a', 'Ifit1', 'Ifit2', 'Ifit3'), group.by = 'Genotype', scale = FALSE)+
  theme(axis.text.x = element_text(angle = 90))

#Look at mock comp in all celltypes
mock_markers <- FindMarkers(mock_samples, test.use = 'MAST', group.by = 'Genotype', ident.1 = 'WT')
mock_markers_up_sig <- mock_markers[mock_markers$p_val_adj < 0.01 & (mock_markers$avg_log2FC) > 1,]
mock_markers_down_sig <- mock_markers[mock_markers$p_val_adj < 0.01 & (mock_markers$avg_log2FC) < -1,]
nrow(mock_markers_up_sig)
nrow(mock_markers_down_sig)

DotPlot(mock_samples, features = 'Fcrls', group.by = 'Genotype', scale = FALSE)
wt_mock_up <- gprofiler2::gost(query = rownames(mock_markers_up_sig), organism = 'mmusculus', evcodes = TRUE, sources = c('GO:BP', 'KEGG'))
wt_mock_down <- gprofiler2::gost(query = rownames(mock_markers_down_sig), organism = 'mmusculus', evcodes = TRUE, sources = c('GO:BP', 'KEGG'))
wt_mock_up$result
wt_mock_down$result
