#Very few samples in wells 16, 17, 18
#Load in libraries
library(Seurat)
library(dplyr)
library(scDblFinder)
library(readr)
library(ggplot2)

#List out our data files, and read them into R
parseOutput <- list.files("./FilteredParseOutput/")

ParseMatricies <- lapply(parseOutput, FUN = function(x){
  ReadParseBio(paste0("./FilteredParseOutput/", x))
})

#Convert data to Seurat objects  
ParseSeuratObj <- lapply(ParseMatricies, FUN = function(x){
  CreateSeuratObject(counts = x, project = "Langat", 
                     min.cells = 1, min.features = 20)
})

#Can also combine this data using split-pipe combined, instead of doing it here
toMerge <- c(ParseSeuratObj[[2]], ParseSeuratObj[[3]], ParseSeuratObj[[4]],
             ParseSeuratObj[[5]], ParseSeuratObj[[6]], ParseSeuratObj[[7]],
             ParseSeuratObj[[8]])

#Merge data before integration
ParseSeuratObj <- merge(ParseSeuratObj[[1]], y = toMerge, 
                        add.cell.ids = paste0('seuObj', seq(1,8)), project = "ParseLangat")


#Use Parse barcoding plate to label cells
#Going to try just left joining and seeing if virus expression matches up
ParseBarcodePlate <- read_csv("ParseBarcodePlate.csv")
ParseBarcodePlate$orig.ident <- c(paste0('0', as.character(seq(1:9))), seq(10,48))
ParseSeuratObj[[]] <- left_join(ParseSeuratObj[[]], ParseBarcodePlate, by = 'orig.ident')

#Label sublibraries 
ParseSeuratObj$subLib <- substr(colnames(ParseSeuratObj), start = 1, stop = 7)

#Assign viral reads from labels
#List files with viral reads
virusCountFiles <- list.files('./viralReadCountsUnadjusted/')
virusCountFiles <- virusCountFiles[grep('P36207', virusCountFiles)]

viralCountsUnadjusted <- lapply(virusCountFiles, FUN = function(x){
  subLibViralCounts <- read_table(paste0("viralReadCountsUnadjusted/", x), 
                                  col_names = FALSE)
  subLibViralCounts
})

#For each sublibrary create a df with counts of virus per cell barcode
#Also adds a sublibrary column for joining purposes afterwards
for(i in 1:length(viralCountsUnadjusted)){
  dat = viralCountsUnadjusted[[i]]
  dat$subLib <- paste0('seuObj', as.character(i))
  dat$cell <- substr(dat$X2, start = 6, stop = 14)
  dat <- dat[c(1,3,4)]
  colnames(dat)[1] = 'virusCount'
  viralCountsUnadjusted[[i]] <- dat
}

viralCountsUnadjusted <- bind_rows(viralCountsUnadjusted)

ParseSeuratObj$cell <- substr(colnames(ParseSeuratObj), start = 9, stop = 16)

ParseSeuratObj[[]] <- left_join(ParseSeuratObj[[]], viralCountsUnadjusted, by = c('subLib', 'cell'))
ParseSeuratObj$virusCount[is.na(ParseSeuratObj$virusCount)] = 0

wellMap <- data.frame(well = c(paste0('A', seq(1,12)), paste0('B', seq(1,12)),
                               paste0('C', seq(1,12)),   paste0('D', seq(1,12))))

#Plot viral counts vs well/treatment. Looks right
ggplot(ParseSeuratObj[[]], aes(x = orig.ident, y = virusCount, col = Treatment))+
  geom_point() +
  scale_x_discrete(labels= wellMap$well)+
  theme(axis.text.x = element_text(angle = 90))

#Check percent of cells in each group expressing virus
ParseSeuratObj[[]] %>% mutate(virusPresent = ifelse(virusCount>0, 'yes', 'no')) %>% 
  group_by(Treatment, virusPresent) %>% dplyr::summarise(count = n()) %>% 
  ggplot(aes(x = Treatment, y = count, fill = virusPresent))+
  geom_bar(stat = 'identity', position = 'dodge')

#Label mitochondrial gene expression
ParseSeuratObj[["percent.mt"]] <- PercentageFeatureSet(ParseSeuratObj, pattern = "^mt-")


#Look at data features
Idents(ParseSeuratObj) <- "all"  #Stops violin plot from grouping by seurat cluster
VlnPlot(ParseSeuratObj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3, pt.size = 0)

#Run through data processing and visualization before integration
ParseSeuratObj <- NormalizeData(ParseSeuratObj)
ParseSeuratObj <- FindVariableFeatures(ParseSeuratObj)
ParseSeuratObj <- ScaleData(ParseSeuratObj)
ParseSeuratObj <- RunPCA(ParseSeuratObj)
ElbowPlot(ParseSeuratObj, ndims = 40)
ParseSeuratObj <- FindNeighbors(ParseSeuratObj, dims = 1:30, reduction = "pca")
ParseSeuratObj <- FindClusters(ParseSeuratObj, resolution = 2, cluster.name = "unintegrated_clusters")  
ParseSeuratObj <- RunUMAP(ParseSeuratObj, dims = 1:30, reduction = "pca")
DimPlot(ParseSeuratObj, group.by = 'orig.ident')

#Look at potential doublets

Parse_sce <- ParseSeuratObj %>% JoinLayers() %>%  as.SingleCellExperiment()

#Use sublibraries as samples for scdblfinder
Parse_sce <- scDblFinder(Parse_sce, dbr.sd=1, samples = 'subLib')
table(Parse_sce$scDblFinder.class)
ParseSeuratObj$scDblFinderLabel <- Parse_sce$scDblFinder.class

ParseSeuratObj <- subset(ParseSeuratObj, scDblFinderLabel == 'singlet')


#Integrate data
#Need to choose integration method - rpca fastest
#Can look at running in parellel to increase speed. CCA takes 5+ hours. RPCA much much faster
ParseSeuratObj_int <- IntegrateLayers(object = ParseSeuratObj, method = CCAIntegration,
                                      orig.reduction = "pca",  new.reduction = "integrated.cca", 
                                      verbose = TRUE)

# re-join layers after integration
ParseSeuratObj_int[["RNA"]] <- JoinLayers(ParseSeuratObj_int[["RNA"]])

ParseSeuratObj_int <- FindNeighbors(ParseSeuratObj_int, reduction = "integrated.cca", dims = 1:30)
ParseSeuratObj_int <- FindClusters(ParseSeuratObj_int, resolution = 1)

ParseSeuratObj_int <- RunUMAP(ParseSeuratObj_int, dims = 1:30, reduction = "integrated.cca")
DimPlot(ParseSeuratObj_int, reduction = "umap", group.by = c("seurat_clusters"), label = TRUE)


