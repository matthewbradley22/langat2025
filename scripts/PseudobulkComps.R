library(Seurat)
library(ggrepel)

#Source functions
source('./scripts/langatFunctions.R')

#Load data
ParseSeuratObj_int <- LoadSeuratRds("./data/FilteredRpcaIntegratedDatNoDoublets.rds") 
ParseSeuratObj_int$hasVirus = ifelse(ParseSeuratObj_int$virusCountPAdj >= 10, 1, 0)

#Check data
newCols <-  c(brewer.pal(12, 'Paired'), '#99FFE6', '#CE99FF', '#18662E',  '#FF8AEF','#737272')
DimPlot(ParseSeuratObj_int, label = FALSE, group.by = 'manualAnnotation', reduction = 'umap.integrated',
        cols = newCols)

#Create pseudobulk object
ParseSeuratObj_int$hasVirus <- factor(ParseSeuratObj_int$hasVirus)
langatPseudoBulk <- createPseudoBulk(ParseSeuratObj_int, c('Genotype', 'Treatment', 'Timepoint', 'Organ', 'hasVirus'))

#DeSeq analysis
rld <- rlog(langatPseudoBulk, blind=TRUE)

#PCA data df
pca_data_condition <- plotPCA(rld, intgroup=c("Organ"), returnData = TRUE) 

ggplot(pca_data_condition, aes(x = PC1, y = PC2, color = Organ)) +
  geom_point() + 
  theme_classic() +
  xlab(paste0("PC1: ", round(attr(pca_data_condition, "percentVar")[1] * 100), "% variance")) +
  ylab(paste0("PC2: ", round(attr(pca_data_condition, "percentVar")[2] * 100), "% variance")) 


