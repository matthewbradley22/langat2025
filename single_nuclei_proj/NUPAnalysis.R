#Load libraries
library(dplyr)
library(Seurat)
library(ggplot2)
library(SingleR)
library(celldex)
library(scran)
library(GSEABase)
library(AUCell)
library(gridExtra)

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
  CreateSeuratObject(counts = x, project = "hupAnalysis", min.cells = 3, min.features = 300)
})

#Merge for v5 ingtegration
toMerge <- c(scDatSeu[[2]], scDatSeu[[3]], scDatSeu[[4]],
             scDatSeu[[5]], scDatSeu[[6]], scDatSeu[[7]])

#Merge data before integration
scDatObj <- merge(scDatSeu[[1]], y = toMerge, 
                        add.cell.ids = paste0('seuObj', seq(1,7)), project = "single_nuclei")



#Label mitochondrial gene expression
scDatObj[["percent.mt"]] <- PercentageFeatureSet(scDatObj, pattern = "^mt-")

#Look at data features
Idents(scDatObj) <- "all"  #Stops violin plot from grouping by seurat cluster
VlnPlot(scDatObj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3, pt.size = 0) 


ggplot(scDatObj[[]], aes(x = 1, y = nFeature_RNA))+
  geom_violin() +
  geom_hline(yintercept = 300)

get_library_labels <- do.call(rbind, stringr::str_split(string = rownames(scDatObj[[]]), pattern = '_'))
lib_labels <- get_library_labels[,1]
scDatObj$library <- lib_labels
metaData$library = paste0('seuObj', seq(1:7))
scDatObj[[]] <-  left_join(scDatObj[[]], metaData, by = 'library')

#Remove likely contamination/doublets. should have low mito reads since this is nuclei
scDatObj<- subset(scDatObj, subset = nFeature_RNA > 100 & nFeature_RNA < 5000 & percent.mt < 10)

#Run through data processing and visualization before integration
scDatObj <- NormalizeData(scDatObj)
scDatObj <- FindVariableFeatures(scDatObj)
scDatObj <- ScaleData(scDatObj)
scDatObj <- RunPCA(scDatObj)
ElbowPlot(scDatObj, ndims = 50)
scDatObj <- FindNeighbors(scDatObj, dims = 1:30, reduction = "pca")

#Don't think resolution here matters since we rerun findclusters post integration. Tried with multiple
scDatObj <- FindClusters(scDatObj, resolution = 2, cluster.name = "unintegrated_clusters") 
scDatObj <- RunUMAP(scDatObj, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(scDatObj)

#Integrate data
#choosing integration method - rpca fastest, cca takes long
scDatObj_int <- IntegrateLayers(object = scDatObj, method = RPCAIntegration,
                                      orig.reduction = "pca",  new.reduction = "integrated.rpca", 
                                      verbose = TRUE)

#re-join layers after integration
scDatObj_int[["RNA"]] <- JoinLayers(scDatObj_int[["RNA"]])
ElbowPlot(scDatObj_int, ndims = 30)
scDatObj_int <- FindNeighbors(scDatObj_int, reduction = "integrated.rpca", dims = 1:20)
scDatObj_int <- FindClusters(scDatObj_int, resolution = 1)
scDatObj_int <- RunUMAP(scDatObj_int, dims = 1:20, reduction = "integrated.rpca", 
                              reduction.name = "umap.integrated")

#Visualize post integration
intUMAP <- DimPlot(scDatObj_int, reduction = "umap.integrated", group.by = c("seurat_clusters"), label = TRUE)+
  NoLegend()+
  ggtitle('Integrated UMAP')
unIntUmap <- DimPlot(scDatObj_int, reduction = "umap.unintegrated", group.by = c("seurat_clusters"), label = TRUE)+
  NoLegend()+
  ggtitle('Unintegrated UMAP')

grid.arrange(intUMAP, unIntUmap, ncol=2)

#remove cells with both sex genes
sexGenes <- scDatObj_int[['RNA']]$data[c('Xist', 'Eif2s3y'),]
sexGenePresence <- colSums(sexGenes > 0)
scDatObj_int$sexGenePresence <- case_when(sexGenePresence == 0 ~ 'None',
                                                sexGenePresence == 1 ~ 'One',
                                                sexGenePresence == 2 ~ 'Two')

table(scDatObj_int$sexGenePresence)
scDatObj_int$sexGenePresence <- subset(scDatObj_int$sexGenePresence, sexGenePresence != 'Two')

#SaveSeuratRds(scDatObj_int, './LGTVscCombined.rds')


##############Begin analysis##############

###########Cell annotation##############

#Integrated assay for finding celltypes
DimPlot(sn_integrated_dat, label = TRUE, reduction = 'umap.integrated')

#Ex neurons
FeaturePlot(sn_integrated_dat, 'Slc17a7', reduction = 'umap.integrated')
FeaturePlot(sn_integrated_dat, 'Sv2b', reduction = 'umap.integrated')
FeaturePlot(sn_integrated_dat, 'Arpp21', reduction = 'umap.integrated')

#in neurons
FeaturePlot(sn_integrated_dat, 'Gad1', reduction = 'umap.integrated')
FeaturePlot(sn_integrated_dat, 'Dlx6os1', reduction = 'umap.integrated')

#Microglia
FeaturePlot(sn_integrated_dat, 'Csf1r', reduction = 'umap.integrated')
FeaturePlot(sn_integrated_dat, 'Tmem119', reduction = 'umap.integrated')
FeaturePlot(sn_integrated_dat, 'P2ry12', reduction = 'umap.integrated')
FeaturePlot(sn_integrated_dat, 'Cx3cr1', reduction = 'umap.integrated')
FeaturePlot(sn_integrated_dat, 'Itgam', reduction = 'umap.integrated')

#microglia/monocytes
FeaturePlot(sn_integrated_dat, 'Ctss', reduction = 'umap.integrated')
FeaturePlot(sn_integrated_dat, 'Cx3cr1', reduction = 'umap.integrated')
FeaturePlot(sn_integrated_dat, 'Csf1r', reduction = 'umap.integrated')

#Macrophage/monocyte markers
FeaturePlot(sn_integrated_dat, 'Ptprc', reduction = 'umap.integrated')
FeaturePlot(sn_integrated_dat, 'Ccr2', reduction = 'umap.integrated')
FeaturePlot(sn_integrated_dat, 'Lyz2', reduction = 'umap.integrated')
FeaturePlot(sn_integrated_dat, 'Ccr2', reduction = 'umap.integrated')

#Astrocytes
FeaturePlot(sn_integrated_dat, 'Aqp4', reduction = 'umap.integrated')
FeaturePlot(sn_integrated_dat, 'Rorb', reduction = 'umap.integrated')
FeaturePlot(sn_integrated_dat, 'Fgfr3', reduction = 'umap.integrated')
FeaturePlot(sn_integrated_dat, 'Gfap', reduction = 'umap.integrated')

#oligo
FeaturePlot(sn_integrated_dat, 'Mag', reduction = 'umap.integrated')
FeaturePlot(sn_integrated_dat, 'Mog', reduction = 'umap.integrated')
FeaturePlot(sn_integrated_dat, 'Plp1', reduction = 'umap.integrated')

#OPC
FeaturePlot(sn_integrated_dat, 'Pdgfra', reduction = 'umap.integrated')
FeaturePlot(sn_integrated_dat, 'Vcan', reduction = 'umap.integrated')

#VLMCs
FeaturePlot(sn_integrated_dat, 'Dcn', reduction = 'umap.integrated')
FeaturePlot(sn_integrated_dat, 'Col1a1', reduction = 'umap.integrated')

#Pericytes
FeaturePlot(sn_integrated_dat, 'Abcc9', reduction = 'umap.integrated')
FeaturePlot(sn_integrated_dat, 'Pdgfrb', reduction = 'umap.integrated')
FeaturePlot(sn_integrated_dat, 'Vtn', reduction = 'umap.integrated')

#ChP
FeaturePlot(sn_integrated_dat, 'Ttr', reduction = 'umap.integrated')
FeaturePlot(sn_integrated_dat, 'Aqp1', reduction = 'umap.integrated')

#Endothelial cells
#brain endo marker Pglyrp1 from https://www.sciencedirect.com/science/article/pii/S0092867420300623
FeaturePlot(sn_integrated_dat, 'Flt1', reduction = 'umap.integrated')
FeaturePlot(sn_integrated_dat, 'Cd93', reduction = 'umap.integrated')
FeaturePlot(sn_integrated_dat, 'Vwf', reduction = 'umap.integrated')
FeaturePlot(sn_integrated_dat, 'Cldn5', reduction = 'umap.integrated')
FeaturePlot(sn_integrated_dat, 'Pglyrp1', reduction = 'umap.integrated')
FeaturePlot(sn_integrated_dat, 'Emcn', reduction = 'umap.integrated')
FeaturePlot(sn_integrated_dat, 'Pecam1', reduction = 'umap.integrated')
FeaturePlot(sn_integrated_dat, 'Adgrl4', reduction = 'umap.integrated')
FeaturePlot(sn_integrated_dat, 'Slco1a4', reduction = 'umap.integrated')

#Label cells
sn_integrated_dat$manualAnnotation <- 
  case_when(sn_integrated_dat$seurat_clusters %in% c(2, 12, 30) ~ 'Oligo',
            sn_integrated_dat$seurat_clusters %in% c(16) ~ 'OPC',
            sn_integrated_dat$seurat_clusters %in% c(18, 10, 8) ~ 'In Neurons',
            sn_integrated_dat$seurat_clusters %in% c(1, 5, 6, 7, 9, 25, 28, 13, 27, 23, 15, 14, 4,
                                                     17, 26) ~ 'Ex Neurons',
            sn_integrated_dat$seurat_clusters %in% c(3, 21) ~ 'Micro/MO',
            sn_integrated_dat$seurat_clusters %in% c(0) ~ 'Astrocytes',
            sn_integrated_dat$seurat_clusters %in% c() ~ 'ChP',
            sn_integrated_dat$seurat_clusters %in% c() ~ 'VLMCs',
            sn_integrated_dat$seurat_clusters %in% c(24) ~ 'Endothelial',
            .default = 'unknown')


#Look at other variables
newCols <-  c(brewer.pal(12, 'Paired'), '#99FFE6', '#CE99FF', '#18662E','#737272',  '#FF8AEF')
pdf("~/Documents/ÖverbyLab/single_nuclei_proj/sn_plots/celltype_umap.pdf", width = 8, height = 6)
DimPlot(sn_integrated_dat, group.by =   'manualAnnotation', cols = newCols, reduction = 'umap.integrated')
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

#Barplot pathways
infected_vs_uninfected_up_paths$result[infected_vs_uninfected_up_paths$result$term_name %in% 
                                   c('Apoptosis', 'TNF signaling pathway', 'cytokine production',
                                     'response to cytokine', 'cell death',
                                     'necroptotic process', 'Regulation of necroptotic cell death'),] %>% 
  ggplot(aes(x = -log10(p_value), y = term_name, fill = source))+
  geom_bar(stat = 'identity', position = 'dodge')+
  ggtitle('Single nuclei data upregulated in infection')

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

