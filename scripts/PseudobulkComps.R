library(Seurat)
library(ggrepel)
library(pheatmap)

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
#Run through again with only infected cells?
ParseSeuratObj_int$hasVirus <- factor(ParseSeuratObj_int$hasVirus)
langatPseudoBulk <- createPseudoBulk(ParseSeuratObj_int, c('Genotype', 'Treatment', 'Timepoint', 'Organ', 'hasVirus'))

infected = subset(ParseSeuratObj_int, hasVirus == 1)
langatPseudoBulk <- createPseudoBulk(infected, c('Genotype', 'Treatment', 'Timepoint', 'Organ'))

#DeSeq analysis
rld <- rlog(langatPseudoBulk, blind=TRUE)

#PCA data df
pca_data_condition <- plotPCA(rld, intgroup=c("Organ"), returnData = TRUE) 

ggplot(pca_data_condition, aes(x = PC1, y = PC2, color = Timepoint)) +
  geom_point() + 
  theme_classic() +
  xlab(paste0("PC1: ", round(attr(pca_data_condition, "percentVar")[1] * 100), "% variance")) +
  ylab(paste0("PC2: ", round(attr(pca_data_condition, "percentVar")[2] * 100), "% variance")) 

# Calculate sample correlation
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)

# Plot heatmap
pheatmap(rld_cor)

#Create deseq2 obj
dds <- DESeq(langatPseudoBulk)
plotDispEsts(dds)
resultsNames(dds)

res <- lfcShrink(dds, 
                 coef = "hasVirus_1_vs_0",
                 type = "apeglm")
res[grep('Ifn',rownames(res)),] %>% as.data.frame()
dplyr::arrange(as.data.frame(res), padj)

res['Rsad2',]
res['Lrp8',]
plotCounts(dds, gene='Rsad2', intgroup="hasVirus")
#Look at mavs levels
plotCounts(dds, gene='Mavs', intgroup="Genotype")

VlnPlot(subset(ParseSeuratObj_int, hasVirus == 1), features = 'Mavs', group.by = 'Genotype')
mean(subset(ParseSeuratObj_int, hasVirus == 1 & Genotype == 'WT')[['RNA']]$data['Mavs',] > 0)
mean(subset(ParseSeuratObj_int, hasVirus == 1 & Genotype == 'IPS1')[['RNA']]$data['Mavs',] > 0)
table(ParseSeuratObj_int$Genotype, ParseSeuratObj_int$hasVirus)
