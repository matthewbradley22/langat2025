library(Seurat)
library(RColorBrewer)
library(msigdbr)
library(dplyr)
library(matrixStats)
library(ggplot2)
library(gt)
library(ggrepel)
library(purrr)
library(forcats)
source('~/Documents/ÖverbyLab//scripts/langatFunctions.R')

#Load data
ParseSeuratObj_int <- LoadSeuratRds("~/Documents/ÖverbyLab/data/FilteredRpcaIntegratedDatNoDoublets.rds") 
ParseSeuratObj_int$hasVirus = ifelse(ParseSeuratObj_int$virusCountPAdj >= 10, 1, 0)


#Check data
newCols <-  c(brewer.pal(12, 'Paired'), '#99FFE6', '#CE99FF', '#18662E','#737272',  '#FF8AEF')
newCols[11] =  '#FF8AEF'

DimPlot(ParseSeuratObj_int, label = FALSE, group.by = 'manualAnnotation', reduction = 'umap.integrated',
        cols = newCols)

#Comparing mock wt and mock ips
mock_for_comparison <- subset(ParseSeuratObj_int, Treatment == 'PBS')

#Control for organ and timepoint when comparing genotypes
mock_comp_bulk <- createPseudoBulk(mock_for_comparison, variables = c('Genotype', 'Organ', 'Timepoint'))
mock_comp_bulk <- DESeq(mock_comp_bulk)
mock_comp_bulk_genotype_res <- results(mock_comp_bulk, name = 'Genotype_WT_vs_IPS1')

#Only 66 DEGs found this way
mock_comp_bulk_genotype_res %>% as.data.frame() %>% dplyr::filter(padj < 0.01 & abs(log2FoldChange) > 1) %>% 
  dplyr::arrange(padj)

plotCounts(mock_comp_bulk, gene = 'Mavs', intgroup = 'Genotype')

#Also look at MAST test 
#Maybe should first subset by day and organ to control for other variables
mock_for_comparison_markers <- FindAllMarkers(mock_for_comparison, only.pos = TRUE, test.use = 'MAST', assay = 'RNA',
                                              group.by = 'Genotype')
mock_for_comparison_markers %>% dplyr::filter(p_val_adj < 0.01 & avg_log2FC > 1)
