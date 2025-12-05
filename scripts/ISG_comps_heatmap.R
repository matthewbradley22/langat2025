library(Seurat)
library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(tidyr)

#Load data
ParseSeuratObj_int <- LoadSeuratRds("~/Documents/ÖverbyLab/data/FilteredRpcaIntegratedDatNoDoublets.rds") 
ParseSeuratObj_int$hasVirus = ifelse(ParseSeuratObj_int$virusCountPAdj >= 10, 1, 0)


#Check data
newCols <-  c(brewer.pal(12, 'Paired'), '#99FFE6', '#CE99FF', '#18662E','#737272',  '#FF8AEF')
newCols[11] =  '#FF8AEF'

DimPlot(ParseSeuratObj_int, label = FALSE, group.by = 'manualAnnotation', reduction = 'umap.integrated',
        cols = newCols)+
  theme(axis.ticks = element_blank(),
        axis.text=element_blank(),
        legend.text=element_text(size=17))+
  xlab('Umap1')+
  ylab('Umap2')+
  guides(color=guide_legend(override.aes=list(size=8)))+
  ggtitle('')

#subset data (no LGTV)
chimeric_mock <- subset(ParseSeuratObj_int, Treatment != 'rLGTV')

#Genes of interest for heatmap
chemokines <- c('Ccl7', 'Ccl2', 'Ccl12', 'Ccl4', 'Ccl3', 'Ccl5',
                'Cxcl10', 'Cxcl9', 'Cxcr3', 'Ccr1', 'Ccr2', 'Ccr5', 'Ccr3', 'Cxcl11', 'Cx3cr1')

chemokine_expression <- FetchData(object = chimeric_mock, vars = chemokines, layer = "data") %>% 
  rownames_to_column(var = 'cell_id') %>% 
  tidyr::pivot_longer(!cell_id, names_to = 'gene', values_to = 'expression')

