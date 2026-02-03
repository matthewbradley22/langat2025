library(tximport)
library(EnsDb.Mmusculus.v79)
library(stringr)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(Seurat)
library(RColorBrewer)
library(tidyr)
source('~/Documents/ÖverbyLab//scripts/langatFunctions.R')


#Load bulk data
#Load in viral levels from salmon quant
salmonPath = '~/Documents/ÖverbyLab/bulkSalmonData/filtered_transcripts_quant/'
salmonQuantFiles <- list.files(salmonPath, recursive = TRUE)
salmonQuantFiles <- salmonQuantFiles[grep('customQuant.sf', salmonQuantFiles)]

#Remove sample not passing qc
salmonQuantFiles = salmonQuantFiles[!salmonQuantFiles == 'Sample_YE-4324-WT-LGTV-48-4/_salmonQuant/customQuant.sf']

#Remove sample not passing qc
#salmonQuantFiles <- salmonQuantFiles[!salmonQuantFiles == 'Sample_YE-4324-WT-LGTV-48-4/_salmonQuant/noVersionQuant.sf']
salmonQuantFiles <- paste0(salmonPath, salmonQuantFiles)
names(salmonQuantFiles) <- str_extract(salmonQuantFiles, 'Sample.+(\\d)')

#Mavs ch + mock
mavsQuantFiles <- salmonQuantFiles[grep('MAVS-ch|MAVS-mock', names(salmonQuantFiles))]

#WT ch + mock
wtQuantFiles <- salmonQuantFiles[grep('WT-ch|WT-Mock', names(salmonQuantFiles))]

#turn transcript labels to gene names
#txdb = TxDb.Mmusculus.UCSC.mm39.refGene
ensDb <- EnsDb.Mmusculus.v79
ensDbDf <- transcripts(ensDb, return.type = 'DataFrame')
k <- keys(ensDb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(ensDb, keys = k, "GENEID", keytype = "TXNAME")

#Add virus to tx2gene - doesnt work because then rows don't match in original files
tx2gene[nrow(tx2gene)+1,] =c('EU790644', 'LGTV', 'EU790644.1')
tx2gene[nrow(tx2gene)+1,] =c('DQ401140', 'TBEV', 'DQ401140.3')

txi_mavs <- tximport(mavsQuantFiles, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
txi_wt <- tximport(wtQuantFiles, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)

#Prepare data for DESeq2
sampleNames_mavs <- colnames(txi_mavs$abundance)
sampleNames_wt <- colnames(txi_wt$abundance)

metadata_mavs <- str_extract(sampleNames_mavs, '[MW].+') %>% str_split('-')
metadata_mavs <- as.data.frame(do.call("rbind", metadata_mavs))
colnames(metadata_mavs) <- c('Genotype','Treatment', 'Timepoint', 'Replicate')

metadata_wt <- str_extract(sampleNames_wt, '[MW].+') %>% str_split('-')
metadata_wt <- as.data.frame(do.call("rbind", metadata_wt))
colnames(metadata_wt) <- c('Genotype','Treatment', 'Timepoint', 'Replicate')

#Method of labelling doesn't work for ifn samples with no timepoint info, so need to make changes
metadata_mavs[metadata_mavs$Timepoint %in% c(1,2,3,4),]$Replicate = metadata_mavs[metadata_mavs$Timepoint %in% c(1,2,3,4),]$Timepoint
metadata_wt[metadata_wt$Timepoint %in% c(1,2,3,4),]$Replicate = metadata_wt[metadata_wt$Timepoint %in% c(1,2,3,4),]$Timepoint

metadata_mavs[metadata_mavs$Timepoint %in% c(1,2,3,4),]$Timepoint = 0
metadata_wt[metadata_wt$Timepoint %in% c(1,2,3,4),]$Timepoint = 0

#metadata[metadata$Timepoint %in% c(1,2,3,4),]$Timepoint = NA
rownames(metadata_mavs) = sampleNames_mavs
rownames(metadata_wt) = sampleNames_wt

#Ensure mocks are capitalized the same way
metadata_mavs[metadata_mavs$Treatment == 'mock' |metadata_mavs$Treatment == 'Mock',]$Treatment = 'Mock'
metadata_wt[metadata_wt$Treatment == 'mock' |metadata_wt$Treatment == 'Mock',]$Treatment = 'Mock'

metadata_mavs$Treatment <- factor(metadata_mavs$Treatment)
metadata_mavs$Treatment <- relevel(metadata_mavs$Treatment,"Mock")

metadata_wt$Treatment <- factor(metadata_wt$Treatment)
metadata_wt$Treatment <- relevel(metadata_wt$Treatment,"Mock")

metadata_mavs$treatment_time <- paste(metadata_mavs$Treatment, metadata_mavs$Timepoint, sep = '_')
metadata_wt$treatment_time <- paste(metadata_wt$Treatment, metadata_wt$Timepoint, sep = '_')

dds_mavs <- DESeqDataSetFromTximport(txi_mavs, metadata_mavs, ~0 + treatment_time)
plotCounts(dds_mavs, gene = 'ENSMUSG00000037523', intgroup = 'treatment_time') #Confirm mavs gene is much lower in knockout than wt below
dds_mavs <- DESeq(dds_mavs)
resultsNames(dds_mavs)

dds_wt <- DESeqDataSetFromTximport(txi_wt, metadata_wt, ~0 + treatment_time)
plotCounts(dds_wt, gene = 'ENSMUSG00000037523', intgroup = 'treatment_time')
dds_wt <- DESeq(dds_wt)
resultsNames(dds_wt)

vsd_mavs <- vst(dds_mavs)
vsd_wt <- vst(dds_wt)

#pca
plotPCA(vsd_mavs, intgroup=c("treatment_time"))
plotPCA(vsd_wt, intgroup=c("treatment_time"))

#Load single cell data
#Load data
ParseSeuratObj_int <- LoadSeuratRds("~/Documents/ÖverbyLab/data/FilteredRpcaIntegratedDatNoDoublets.rds") 
ParseSeuratObj_int$hasVirus = ifelse(ParseSeuratObj_int$virusCountPAdj >= 10, 1, 0)
ParseSeuratObj_int$manualAnnotation[ParseSeuratObj_int$manualAnnotation == 'Macrophage/Monocytes'] = 'Macro/Mono'

#Create new column for use in pseudobulk later
ParseSeuratObj_int$Treatment_celltype <- paste(ParseSeuratObj_int$Treatment, ParseSeuratObj_int$manualAnnotation, sep = '_')

#Add extra column for umap plot
ParseSeuratObj_int$geno_timepoint_treatment = paste(ParseSeuratObj_int$Genotype, ParseSeuratObj_int$Timepoint, ParseSeuratObj_int$Treatment)

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

DimPlot(ParseSeuratObj_int, label = FALSE, group.by = 'Treatment', reduction = 'umap.integrated',
        cols = newCols)

#Start looking at overall astrocyte degs, then split by group
astrocytes <- subset(ParseSeuratObj_int, manualAnnotation == 'Astrocytes')

#Recluster astrocytes
#This was true before so can switch to get old umap
astrocytes <- prepSeuratObj(astrocytes, use_all_genes = FALSE)
ElbowPlot(astrocytes, ndims = 40)
astrocytes <- prepUmapSeuratObj(astrocytes, nDims = 20, reductionName = 'astrocytes_umap', resolution_value = 0.8)

#UMAPs
DimPlot(astrocytes, reduction = 'astrocytes_umap')
astrocytes$infection_group <- ifelse(astrocytes$Treatment %in% c('rChLGTV', 'rLGTV'), 'infected', 'uninfected')

#Compare gene lists here with those from single cell astros
#from scRNA_astrocyte_analysis.R script





#PBS
pbs_sc_ips <- ips_astrocytes_cerebrum_markers[ips_astrocytes_cerebrum_markers$avg_log2FC > 1 & 
                                                ips_astrocytes_cerebrum_markers$p_val_adj < 0.01 &
                                                ips_astrocytes_cerebrum_markers$cluster == 'PBS',]
pbs_mavs_bulk = gene_lists_for_pathways_mavs[[4]]
pbs_mavs_bulk = pbs_mavs_bulk[pbs_mavs_bulk$pvalue < 0.01 & pbs_mavs_bulk$log2FoldChange > 1,]
pbs_mavs_bulk_genes <- ensembldb::select(EnsDb.Mmusculus.v79, keys= pbs_mavs_bulk$input, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
pbs_sc_ips[which(pbs_sc_ips$gene %in% pbs_mavs_bulk_genes$SYMBOL),]
