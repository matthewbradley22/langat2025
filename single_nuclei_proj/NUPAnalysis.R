#Load libraries
library(dplyr)
library(Seurat)
library(ggplot2)
library(SingleR)
library(celldex)
library(scran)
library(GSEABase)
library(AUCell)

#Source function
source('~/Documents/ÖverbyLab/scripts/langatFunctions.R')

#This is where the 10x data is
setwd('~/Documents/ÖverbyLab/single_nuclei_proj/')

#Read in processed data
sn_integrated_dat <- LoadSeuratRds('./LGTVscCombined.rds')

#Start here to perform integration
outputDirs <- list.dirs(path = ".", full.names = TRUE, recursive = TRUE)
outputDirs <- outputDirs[grep('filtered_feature_bc_matrix', outputDirs)]
metaData <- read.delim("./samplemeta.csv.correct")
  
#Begin processing data
scDat <- lapply(outputDirs, FUN = function(x){
  dat <- Read10X(data.dir =x)
  dat
})

scDatSeu <-  lapply(scDat, FUN = function(x){
  CreateSeuratObject(counts = x, project = "hupAnalysis", min.cells = 3, min.features = 200)
})

for(i in 1: length(scDatSeu)){
  scDatSeu[[i]]$ind <- metaData[i,]$X
  scDatSeu[[i]]$new_genotype <- metaData[i,]$new_genotype
  scDatSeu[[i]]$treatment <- metaData[i,]$treatment
  scDatSeu[[i]]$infected <- metaData[i,]$infected
}

scDatSeu <-  lapply(scDatSeu, FUN = function(x){
  x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = "^mt-")
  x
})

scDatSeu <-  lapply(scDatSeu, FUN = function(x){
  x <- subset(x, subset = nFeature_RNA > 100 & nFeature_RNA < 3500 & percent.mt < 15)
  x
})

scDatSeu <-  lapply(scDatSeu, FUN = function(x){
  x <-NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(x)
  x <- ScaleData(x, features = all.genes)
  x <- RunPCA(x, features = VariableFeatures(object = x))
})

features <- SelectIntegrationFeatures(object.list = scDatSeu)
anchors <- FindIntegrationAnchors(object.list = scDatSeu, anchor.features = features)

scCombined <- IntegrateData(anchorset = anchors)

# Run the standard workflow for visualization and clustering
scCombined <- ScaleData(scCombined, verbose = FALSE)
scCombined <- RunPCA(scCombined, verbose = FALSE)
ElbowPlot(scCombined, ndims = 30)
scCombined <- FindNeighbors(scCombined, reduction = "pca", dims = 1:30)
scCombined <- FindClusters(scCombined, resolution = 0.5)
scCombined <- RunUMAP(scCombined, reduction = "pca", dims = 1:30)
DimPlot(scCombined)

scCombined <- JoinLayers(scCombined, assay = 'RNA')

#remove cells with both sex genes
sexGenes <- scCombined[['RNA']]$data[c('Xist', 'Eif2s3y'),]
sexGenePresence <- colSums(sexGenes > 0)
scCombined$sexGenePresence <- case_when(sexGenePresence == 0 ~ 'None',
                                                sexGenePresence == 1 ~ 'One',
                                                sexGenePresence == 2 ~ 'Two')

table(scCombined$sexGenePresence)
scCombined <- subset(scCombined, sexGenePresence != 'Two')

#SaveSeuratRds(scCombined, './LGTVscCombined.rds')


##############Begin analysis##############


###########Cell annotation##############
DefaultAssay(sn_integrated_dat) <- "RNA"
DimPlot(sn_integrated_dat, label = TRUE)

#Ex neurons
FeaturePlot(sn_integrated_dat, 'Slc17a7')
FeaturePlot(sn_integrated_dat, 'Sv2b')
FeaturePlot(sn_integrated_dat, 'Arpp21')

#in neurons
FeaturePlot(sn_integrated_dat, 'Gad1')
FeaturePlot(sn_integrated_dat, 'Dlx6os1')

#microglia/monocytes
FeaturePlot(sn_integrated_dat, 'Ctss')
FeaturePlot(sn_integrated_dat, 'Cx3cr1')
FeaturePlot(sn_integrated_dat, 'Csf1r')

#Astrocytes
FeaturePlot(sn_integrated_dat, 'Aqp4')
FeaturePlot(sn_integrated_dat, 'Rorb')
FeaturePlot(sn_integrated_dat, 'Fgfr3')
FeaturePlot(sn_integrated_dat, 'Gfap')

#oligo
FeaturePlot(sn_integrated_dat, 'Mag')
FeaturePlot(sn_integrated_dat, 'Mog')
FeaturePlot(sn_integrated_dat, 'Plp1')

#OPC
FeaturePlot(sn_integrated_dat, 'Pdgfra')
FeaturePlot(sn_integrated_dat, 'Vcan')

#VLMCs
FeaturePlot(sn_integrated_dat, 'Dcn')
FeaturePlot(sn_integrated_dat, 'Col1a1')

#Pericytes
FeaturePlot(sn_integrated_dat, 'Abcc9')
FeaturePlot(sn_integrated_dat, 'Pdgfrb')
FeaturePlot(sn_integrated_dat, 'Vtn')

#ChP
FeaturePlot(sn_integrated_dat, 'Ttr')
FeaturePlot(sn_integrated_dat, 'Aqp1')

#Label cells
sn_integrated_dat$manualAnnotation <- 
  case_when(sn_integrated_dat$seurat_clusters %in% c(1) ~ 'Oligo',
            sn_integrated_dat$seurat_clusters %in% c(13) ~ 'OPC',
            sn_integrated_dat$seurat_clusters %in% c(14, 12, 19) ~ 'In Neurons',
            sn_integrated_dat$seurat_clusters %in% c(6,7, 9, 4, 17, 20, 16, 0, 5, 10, 8, 18, 22, 15) ~ 'Ex Neurons',
            sn_integrated_dat$seurat_clusters %in% c(3) ~ 'Micro/MO',
            sn_integrated_dat$seurat_clusters %in% c(2) ~ 'Astrocytes',
            sn_integrated_dat$seurat_clusters %in% c(26) ~ 'ChP',
            sn_integrated_dat$seurat_clusters %in% c(21) ~ 'VLMCs',
            .default = 'unknown') 


#Look at other variables
newCols <-  c(brewer.pal(12, 'Paired'), '#99FFE6', '#CE99FF', '#18662E','#737272',  '#FF8AEF')
pdf("~/Documents/ÖverbyLab/single_nuclei_proj/sn_plots/celltype_umap.pdf", width = 8, height = 6)
DimPlot(sn_integrated_dat, group.by =   'manualAnnotation', cols = newCols)
dev.off()

DimPlot(sn_integrated_dat, group.by = 'treatment')
DimPlot(sn_integrated_dat, group.by = 'infected')
DimPlot(sn_integrated_dat, group.by = 'new_genotype')


###########Necroptosis and pyroptosis analysis###############
#Gene lists
necroptosis <- c('Tnf', 'Tnfrsf1a', 'Ripk2', 'Mlkl', 'Ripk1', 'Ripk3')
pyoptosis <- c('Gsdmc', 'Nlrp3', 'Aim2', 'Gsdmd', 'Il18', 'Il1b', 'Casp9', 'Casp8', 'Casp6', 'Casp3', 'Casp4', 'Casp1')

#Remove ifnar for this analysis, I am guessing this is anything labelled KO in new_genotype, need to confirm
table(sn_integrated_dat$new_genotype)
cells_for_necroptosis <- subset(sn_integrated_dat, new_genotype == 'wt' | new_genotype == 'wt (same)')
table(cells_for_necroptosis$new_genotype)

#Split by treatment to plot each separately


pdf("~/Documents/ÖverbyLab/single_nuclei_proj/sn_plots/sn_necroptosis_dotplot.pdf", width = 9, height = 5)
DotPlot(sn_integrated_dat, features = c(necroptosis, pyoptosis), group.by = 'manualAnnotation', scale = FALSE)+
  coord_flip()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
  ggtitle("Single Nuclei necroptosis and pyoptosis genes")
dev.off()


###########Pseudobulk Comparison of treatment###############
############################################################
sn_integrated_dat_wt <- subset(sn_integrated_dat, new_genotype %in% c('wt', 'wt (same)'))
table(sn_integrated_dat_wt$new_genotype)

#Bool variable doesn't work well with pseudobulk pipeline, just make factor
sn_integrated_dat_wt$infected <- factor(sn_integrated_dat_wt$infected)

sn_wt_bulk <- createPseudoBulk(sn_integrated_dat_wt, variables = c('infected', 'manualAnnotation'))
sn_wt_bulk <- DESeq(sn_wt_bulk)
resultsNames(sn_wt_bulk)
infected_vs_uninfected <- results(sn_wt_bulk, name="infected_TRUE_vs_FALSE")
infected_vs_uninfected_upregulated <- subset(infected_vs_uninfected, padj < 0.01 & log2FoldChange > 1)
infected_vs_uninfected_upregulated
infected_vs_uninfected_up_paths<- gprofiler2::gost(query = rownames(infected_vs_uninfected_upregulated), organism = 'mmusculus', evcodes = TRUE)

#Look at pathway results by path type
infected_vs_uninfected_up_paths$result[infected_vs_uninfected_up_paths$result$source == 'KEGG',]
infected_vs_uninfected_up_paths$result[infected_vs_uninfected_up_paths$result$source == 'GO:MF',]
infected_vs_uninfected_up_paths$result[infected_vs_uninfected_up_paths$result$source == 'GO:BP',]

#Downregulated genes too

#######NUP Analysis###########
###
FeaturePlot(scCombined, features = 'Nup98')
FeaturePlot(scCombined, features = 'LGTV')

lgtvCounts = scCombined[['RNA']]$data[grep('LGTV', rownames(scCombined[['RNA']]$data), ignore.case = TRUE),] 
sum(lgtvCounts > 0) / length(lgtvCounts)

scCombined[[]] <- scCombined[[]] %>% mutate('lgtvCounts' = lgtvCounts)
VlnPlot(scCombined, features = 'lgtvCounts', group.by = 'infected')

scCombined[[]] <- scCombined[[]] %>% mutate(category = case_when(infected == FALSE ~ 'uninfected',
                                                                 infected == TRUE & lgtvCounts == 0 ~ 'bystander',
                                                                 infected == TRUE & lgtvCounts > 0 ~ 'infected'))

scCombined[[]] <- scCombined[[]] %>% mutate(presence = case_when(lgtvCounts > 0 ~ 'present',
                                                                 lgtvCounts == 0 ~ 'absent'))

scCombined[[]] %>% subset(treatment == 'LGTV') %>% group_by(ind, presence) %>% summarise(totals = n())
table(scCombined$infected, scCombined$category)

VlnPlot(scCombined, features = c('Nup98', 'Nup153'), group.by = 'category', assay = 'RNA')

DotPlot(object = scCombined, features = c('Nup98', 'Nup153'), group.by = 'category', assay = 'RNA')

FeatureScatter(object = scCombined, feature1 = 'LGTV', feature2 = 'Nup153', group.by = 'category', plot.cor = FALSE) +
  xlab('LGTV levels') +
  ylab('Nup98 levels')


table(scCombined$infected, scCombined$presence)
corDat <- scCombined[['RNA']]$counts[c('LGTV', 'Nup98', 'Nup153'),] %>% t() %>% as.data.frame()
presenceCorDat <- corDat %>%
  mutate(across(where(is.numeric), ~ifelse(.x >0, 1, 0)))

table(presenceCorDat$LGTV, presenceCorDat$Nup153)

#Label cell types with singleR
DefaultAssay(object = scCombined) <- "RNA"

ref.mouse <- fetchReference("mouse_rnaseq", "2024-02-26")
sce <- as.SingleCellExperiment(scCombined)
pred.sce <- SingleR(test = sce, ref = ref.mouse, labels = ref.mouse$label.fine)

scCombined[["SingleR.labels"]] <- pred.sce$labels

DotPlot(scCombined, features = c('Nup98', 'Nup153'), group.by = 'SingleR.labels')

table(scCombined$SingleR.labels)
DimPlot(scCombined, group.by = 'SingleR.labels')
DotPlot(scCombined, features = 'LGTV', group.by = 'SingleR.labels') + theme(legend.position = 'None')

neurons <- subset(scCombined, SingleR.labels == 'Neurons')
DotPlot(neurons, features = c('Nup98', 'Nup153'), group.by = 'category')

activeNeurons <- subset(scCombined, SingleR.labels == 'Neurons activated')
microgliaAct <- subset(scCombined, SingleR.labels == 'Microglia activated')
DotPlot(activeNeurons, features = c('Nup98', 'Nup153'), group.by = 'category')
DotPlot(microgliaAct, features = c('Nup98', 'Nup153'), group.by = 'category')
#Attempt manual annotation based on marker expression
#Following along https://bioconductor.org/books/3.13/OSCA.basic/cell-type-annotation.html
brainRef <- ZeiselBrainData()
brainRef <- logNormCounts(brainRef)
wilcox.z <- pairwiseWilcox(brainRef, brainRef$level1class, 
                           lfc=1, direction="up")
markers.z <- getTopMarkers(wilcox.z$statistics, wilcox.z$pairs,
                           pairwise=FALSE, n=50)

all.sets <- lapply(names(markers.z), function(x) {
  GeneSet(markers.z[[x]], setName=x)        
})
all.sets <- GeneSetCollection(all.sets)
sceC <- as.SingleCellExperiment(scCombined, assay = 'counts')
rankings <- AUCell_buildRankings(counts(sce),
                                 plotStats=FALSE, verbose=FALSE)
cell.aucs <- AUCell_calcAUC(all.sets, rankings)
results <- t(assay(cell.aucs))
new.labels <- colnames(results)[max.col(results)]
scCombined[['AUC_labels']] <- new.labels
DimPlot(scCombined, group.by = 'AUC_labels')
table(scCombined$SingleR.labels, scCombined$AUC_labels)

VlnPlot(scCombined, features = 'LGTV', group.by = 'AUC_labels') + theme(legend.position = 'None')
table(scCombined$presence, scCombined$AUC_labels)

#Look at LGTV levels in pyramidal neurons
pyramidal_SS <- subset(scCombined, AUC_labels == 'pyramidal SS')
pyramidal_CA1 <-subset(scCombined, AUC_labels == 'pyramidal CA1')


DotPlot(pyramidal_SS, features = c('Nup98', 'Nup153'), group.by = 'category')
DotPlot(pyramidal_CA1, features = c('Nup98', 'Nup153'), group.by = 'category')

#Excitatory neurons
FeaturePlot(scCombined, features = c('Hpcal1', 'Rorb', 'Fezf2', 'Themis', 'Nxph4'))
exMArkerFeats <- list(c('Hpcal1', 'Rorb', 'Fezf2', 'Themis', 'Nxph4')) # from Identification of visual cortex cell types and species differences using single-cell RNA sequencing
exMarkerFeats2 <- list(c('Slc17a7', 'Grm4')) #from Integrative single-cell analysis of transcriptional and epigenetic states in the human adult brain
scCombined <- AddModuleScore(scCombined, name = 'exMarkers',  features = exMArkerFeats, assay = 'RNA')
scCombined <- AddModuleScore(scCombined, name = 'exMarkersTwo',  features = exMarkerFeats2, assay = 'RNA')
FeaturePlot(scCombined, features =  'exMarkersTwo1')
VlnPlot(scCombined, features = 'exMarkers1', group.by = 'SingleR.labels')
VlnPlot(scCombined, features = 'exMarkersTwo1', group.by = 'SingleR.labels')
summary(scCombined$exMarkers1)
scCombined[[]] <- scCombined[[]] %>% mutate(excitatory = case_when(exMarkers1 > 0.1 ~ 'ex',
                                                                   exMarkers1 < 0.1 ~ 'notEx'))
excitatoryNeurons <- subset(scCombined, SingleR.labels == 'Neurons' & excitatory == 'ex')


DotPlot(excitatoryNeurons, features = c('Nup98', 'Nup153'), group.by = 'category')
VlnPlot(excitatoryNeurons, features = c('Nup98', 'Nup153'), group.by = 'category')



#SaveSeuratRds(scCombined, '~/LGTVscCombined.rds')

