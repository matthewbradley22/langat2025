#Very few samples in wells 16, 17, 18
#Load in libraries
library(Seurat)
library(dplyr)
library(scDblFinder)
library(readr)
library(ggplot2)
library(gridExtra)

#List out our data files, and read them into R
parseOutput <- list.files("./data/FilteredParseOutput/")

ParseMatricies <- lapply(parseOutput, FUN = function(x){
  ReadParseBio(paste0("./data/FilteredParseOutput/", x))
})

#Convert data to Seurat objects  
ParseSeuratObj <- lapply(ParseMatricies, FUN = function(x){
  CreateSeuratObject(counts = x, project = "Langat")
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
#importing both undajusted and partially adjusted (after fragmentation) viral counts
virusCountFiles <- list.files('./viralReadCountsUnadjusted/')
virusCountFilesP <- list.files('./viralReadCountsPartiallyAdjusted/')

virusCountFiles <- virusCountFiles[grep('P36207', virusCountFiles)]

viralCountsUnadjusted <- lapply(virusCountFiles, FUN = function(x){
  subLibViralCounts <- read_table(paste0("viralReadCountsUnadjusted/", x), 
                                  col_names = FALSE)
  subLibViralCounts
})

viralCountsPartAdjusted <- lapply(virusCountFilesP, FUN = function(x){
  subLibViralCounts <- read_table(paste0("viralReadCountsPartiallyAdjusted/", x), 
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

for(i in 1:length(viralCountsPartAdjusted)){
  dat = viralCountsPartAdjusted[[i]]
  dat$subLib <- paste0('seuObj', as.character(i))
  dat$cell <- substr(dat$X2, start = 6, stop = 14)
  dat <- dat[c(1,3,4)]
  colnames(dat)[1] = 'virusCountPAdj'
  viralCountsPartAdjusted[[i]] <- dat
}



viralCountsUnadjusted <- bind_rows(viralCountsUnadjusted)
viralCountsPartAdjusted <- bind_rows(viralCountsPartAdjusted)


ParseSeuratObj$cell <- substr(colnames(ParseSeuratObj), start = 9, stop = 16)

ParseSeuratObj[[]] <- left_join(ParseSeuratObj[[]], viralCountsUnadjusted, by = c('subLib', 'cell'))
ParseSeuratObj[[]] <- left_join(ParseSeuratObj[[]], viralCountsPartAdjusted, by = c('subLib', 'cell'))

ParseSeuratObj$virusCount[is.na(ParseSeuratObj$virusCount)] = 0
ParseSeuratObj$virusCountPAdj[is.na(ParseSeuratObj$virusCountPAdj)] = 0

wellMap <- data.frame(well = c(paste0('A', seq(1,12)), paste0('B', seq(1,12)),
                               paste0('C', seq(1,12)),   paste0('D', seq(1,12))))

#Plot viral counts vs well/treatment. Looks right
ggplot(ParseSeuratObj[[]], aes(x = orig.ident, y = virusCountPAdj, col = Treatment))+
  geom_point() +
  scale_x_discrete(labels= wellMap$well)+
  theme(axis.text.x = element_text(angle = 90))

ggplot(ParseSeuratObj[[]], aes(x = virusCount, y = virusCountPAdj))+
  geom_point()+
  geom_abline(slope=1)

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
ParseSeuratObj <- RunUMAP(ParseSeuratObj, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(ParseSeuratObj)

#Look at potential doublets
Parse_sce <- ParseSeuratObj %>% JoinLayers() %>%  as.SingleCellExperiment()

#Use sublibraries as samples for scdblfinder
Parse_sce <- scDblFinder(Parse_sce, dbr.sd=1, samples = 'subLib')
table(Parse_sce$scDblFinder.class)
ParseSeuratObj$scDblFinderLabel <- Parse_sce$scDblFinder.class
DimPlot(ParseSeuratObj, group.by = 'scDblFinderLabel')

#Keep doublets for now, see if they group together
#ParseSeuratObj <- subset(ParseSeuratObj, scDblFinderLabel == 'singlet')

#Plot virus locations based on treatment
table(ParseSeuratObj$Organ, ParseSeuratObj$Treatment) %>% as.data.frame() %>% 
  ggplot(aes(x = Var1, y = Freq, fill = Var2))+
  geom_bar(stat = 'identity', position = 'dodge')+
  xlab('Organ')+
  ylab('Number of cells')

#Integrate data
#choosing integration method - rpca fastest, cca takes 5+ hours
ParseSeuratObj_int <- IntegrateLayers(object = ParseSeuratObj, method = RPCAIntegration,
                                      orig.reduction = "pca",  new.reduction = "integrated.rpca", 
                                      verbose = TRUE)

#re-join layers after integration
ParseSeuratObj_int[["RNA"]] <- JoinLayers(ParseSeuratObj_int[["RNA"]])

ParseSeuratObj_int <- FindNeighbors(ParseSeuratObj_int, reduction = "integrated.rpca", dims = 1:30)
ParseSeuratObj_int <- FindClusters(ParseSeuratObj_int, resolution = 1)

ParseSeuratObj_int <- RunUMAP(ParseSeuratObj_int, dims = 1:30, reduction = "integrated.rpca", 
                              reduction.name = "umap.integrated")

intUMAP <- DimPlot(ParseSeuratObj_int, reduction = "umap.integrated", group.by = c("seurat_clusters"), label = TRUE)+
  NoLegend()+
  ggtitle('Integrated UMAP')
unIntUmap <- DimPlot(ParseSeuratObj_int, reduction = "umap.unintegrated", group.by = c("seurat_clusters"), label = TRUE)+
  NoLegend()+
  ggtitle('Unintegrated UMAP')

grid.arrange(intUMAP, unIntUmap, ncol=2)

#Add neo reads in same way viral reads were added 

#Neo reads
neoFilesP <- list.files('./neoCountsPartiallyAdjusted/')

neoCountsPartAdjusted <- lapply(neoFilesP, FUN = function(x){
  subLibNeoCounts <- read_table(paste0("neoCountsPartiallyAdjusted/", x), 
                                col_names = FALSE)
  subLibNeoCounts
})

# Convert neo reads in same way as viral above
for(i in 1:length(neoCountsPartAdjusted)){
  dat = neoCountsPartAdjusted[[i]]
  dat$subLib <- paste0('seuObj', as.character(i))
  dat$cell <- substr(dat$X2, start = 6, stop = 14)
  dat <- dat[c(1,3,4)]
  colnames(dat)[1] = 'neoCountPAdj'
  neoCountsPartAdjusted[[i]] <- dat
}

neoCountsPartAdjusted <- bind_rows(neoCountsPartAdjusted)

ParseSeuratObj_int[[]] <- left_join(ParseSeuratObj_int[[]], 
                                    neoCountsPartAdjusted, by = c('subLib', 'cell'))

ParseSeuratObj_int$neoCountPAdj[is.na(ParseSeuratObj_int$neoCountPAdj)] = 0

#plot neo counts
ParseSeuratObj_int$neoPresence <- ifelse(ParseSeuratObj_int$neoCountPAdj > 0, 'yes', 'no')
ParseSeuratObj_int[[]] %>% group_by(orig.ident) %>% dplyr::summarise(totalNeo = sum(neoCountPAdj),
                                                                     Genotype = unique(Genotype)) %>% 
  ggplot(aes(x = orig.ident, y = totalNeo, fill = Genotype))+
  geom_bar(stat = 'identity')+
  scale_x_discrete(labels= wellMap$well)+
  theme(axis.text.x = element_text(angle = 90))+
  ylab('Total Neo reads')

#proportion of cells expressing neo
ParseSeuratObj_int[[]] %>% group_by(orig.ident) %>% dplyr::summarise(neoProp = mean(neoPresence == 'yes'),
                                                                     Genotype = unique(Genotype)) %>% 
  ggplot(aes(x = orig.ident, y = neoProp, fill = Genotype))+
  geom_bar(stat = 'identity')+
  ylab('Proportion of cells expressing Neo')+
  scale_x_discrete(labels= wellMap$well)+
  theme(axis.text.x = element_text(angle = 90))

#Plot total cells as well
table(ParseSeuratObj_int$Genotype ,ParseSeuratObj_int$orig.ident) %>% as.data.frame() %>% 
  ggplot(aes(x = Var2, y = Freq, fill = Var1))+
  geom_bar(stat = 'identity')+
  scale_x_discrete(labels= wellMap$well)+
  theme(axis.text.x = element_text(angle = 90))+
  ylab('Total Cells in Sample')+
  xlab('Sample Well')

ggplot(ParseSeuratObj_int[[]], aes(x = orig.ident, y = neoPresence, fill = neoPresence))+
  geom_bar(stat = 'identity', position = 'stack') +
  scale_x_discrete(labels= wellMap$well)+
  theme(axis.text.x = element_text(angle = 90))

#Proportion of infection by genotype
ParseSeuratObj_int[[]] %>% mutate(virusPresence = ifelse(virusCountPAdj > 0, 'yes', 'no')) %>% 
  group_by(Genotype) %>% 
  dplyr::summarise(virusPresenceProp = mean(virusPresence == 'yes')) %>% 
  ggplot(aes(x = Genotype, y = virusPresenceProp, fill = Genotype)) +
  geom_bar(stat = 'identity')+
  ylab('Virus presence (threshold 1 viral read)')

ParseSeuratObj_int[[]] %>% mutate(virusPresence = ifelse(virusCountPAdj > 4, 'yes', 'no')) %>% 
  group_by(Genotype) %>% 
  dplyr::summarise(virusPresenceProp = mean(virusPresence == 'yes')) %>% 
  ggplot(aes(x = Genotype, y = virusPresenceProp, fill = Genotype)) +
  geom_bar(stat = 'identity')+
  ylab('Virus presence (threshold 5 viral reads)')

#Which cells have a lot of virus reads
highVirus <- subset(ParseSeuratObj_int, virusCountPAdj > 400)
table(highVirus$singleR_labels)

#Plot viral counts vs well/treatment. Looks right
ggplot(ParseSeuratObj_int[[]], aes(x = orig.ident, y = virusCountPAdj, col = Treatment))+
  geom_point() +
  scale_x_discrete(labels= wellMap$well)+
  theme(axis.text.x = element_text(angle = 90))

ggplot(highVirus[[]], aes(x = orig.ident, y = virusCountPAdj, col = singleR_labels))+
  geom_point() +
  scale_x_discrete(labels= wellMap[table(highVirus$orig.ident) %>% names() %>% as.numeric(),])+
  theme(axis.text.x = element_text(angle = 90))

#Virus dot plot
DotPlot(ParseSeuratObj_int, features = 'virusCountPAdj', group.by = 'singleR_labels')

ParseSeuratObj_int[[]] %>% mutate(virusPresence = ifelse(virusCountPAdj > 2, 1, 0)) %>% 
  group_by(Organ) %>% summarise(virusProp = mean(virusPresence))
