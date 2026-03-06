#Look at genes upregulated in astrocytes upon infection across datasets
#Start with single cell data

#Packages and functions
library(Seurat)
library(UCell)
library(RColorBrewer)
source('~/Documents/ÖverbyLab//scripts/langatFunctions.R')

############ Loading in all data (single-cell, bulk, single-nuclei) ############ 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

## Load single-cell data
ParseSeuratObj_int <- LoadSeuratRds("~/Documents/ÖverbyLab/data/FilteredRpcaIntegratedDatNoDoublets.rds") 
ParseSeuratObj_int$hasVirus = ifelse(ParseSeuratObj_int$virusCountPAdj >= 10, 1, 0)

#Check data
newCols <-  c(brewer.pal(12, 'Paired'), '#99FFE6', '#CE99FF', '#18662E','#737272',  '#FF8AEF')
newCols[11] =  '#FF8AEF'
DimPlot(ParseSeuratObj_int, label = FALSE, group.by = 'manualAnnotation', reduction = 'umap.integrated',
        cols = newCols)
sc_astros <- subset(ParseSeuratObj_int, manualAnnotation == 'Astrocytes' & Genotype == 'WT' & Treatment != 'rLGTV')

## Loading bulk data takes quite a few lines, so will just run bulk_astrocyte_heatmaps.R to load both datasets
#Just check that they're loaded
dds_mavs
dds_wt
plotPCA(vsd_mavs, intgroup=c("treatment_time"))
plotPCA(vsd_wt, intgroup=c("treatment_time"))

########## Get lists of DEGs for respective datasets ########## 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

#SC markers
sc_astro_markers <- FindAllMarkers(sc_astros, group.by = 'Treatment',test.use = 'MAST', only.pos = TRUE)
sc_astro_markers_chlgtv <- subset(sc_astro_markers, cluster == 'rChLGTV' & p_val_adj < 0.01 & avg_log2FC > 1)

#Bulk markers

#Compare mock to average of others
#Changing dds formula for this comparison, just compare between treatment groups 
dds_wt <- DESeqDataSetFromTximport(txi_wt, metadata_wt, ~Treatment)
dds_wt <- DESeq(dds_wt)
resultsNames(dds_wt)

dds_wt_res <- results(dds_wt, name = "Treatment_chLGTV_vs_Mock")
dds_wt_res_sig <- dds_wt_res %>% as.data.frame() %>% dplyr::filter(padj < 0.0001 & log2FoldChange > 1) %>% 
  dplyr::arrange(padj)
plotCounts(dds_wt, gene = 'ENSMUSG00000035863', intgroup = 'Treatment')

########## Comapre up and downregulated genes between datasets ########## 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

#Convert bulk deg ensembl names to symbols


