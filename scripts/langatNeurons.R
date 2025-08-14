#Explore neurons in data
#Load packages
library(Seurat)
library(ggplot2)
library(dplyr)

#Source functions
source('./scripts/langatFunctions.R')

#Load in data
ParseSeuratObj_int <- LoadSeuratRds("./data/FilteredRpcaIntegratedDat.rds")

#filter to singlets. went through analysis and filtering made sense, doublets mostly clustered together 
ParseSeuratObj_int <- subset(ParseSeuratObj_int, scDblFinderLabel == 'singlet')

DimPlot(ParseSeuratObj_int, reduction = 'umap.integrated', label = TRUE) +
  NoLegend()

DimPlot(ParseSeuratObj_int, reduction = 'umap.integrated', group.by = 'singleR_labels',
        cols = c(rep('gray', 14), 'red', rep('gray', 3))) 

#Neurons appear to make up clusters 24, 24, 37, and 44
table(subset(ParseSeuratObj_int, singleR_labels == 'Neurons')$seurat_clusters)

#Can look at proportions of neurons in each cluster, if low proportion then consider looking closer
table(subset(ParseSeuratObj_int, seurat_clusters == '27')$singleR_labels)

#Can start with conservative approach of subsetting by celltype and cluster, removing doublets
neurons <- subset(ParseSeuratObj_int, singleR_labels == 'Neurons' & seurat_clusters %in% c('27', '24',
                                                                                          '37', '44'))


table(neurons$Genotype)
table(neurons$Treatment) %>% barplot(main = 'Neurons across treatments') 
table(neurons$Timepoint)%>% barplot(main = 'Neurons across times') 

#More neurons in cerebellum is interesting
table(neurons$Organ)%>% barplot(main = 'Neurons across organs') 

#Check virus presence. Make this slightly stricter later on (need > 1 virus) for closer analysis
neurons$virusPresent <- ifelse(neurons$virusCountPAdj > 0, 1, 0)
table(neurons$virusPresent) %>% barplot(main = 'Neurons with viral reads') 
table(neurons$virusPresent, neurons$Organ) %>% barplot(legend.text = TRUE)

#Rerun through data processing and visualization
neurons <- NormalizeData(neurons)
neurons <- FindVariableFeatures(neurons)
neurons <- ScaleData(neurons)
neurons <- RunPCA(neurons)
ElbowPlot(neurons, ndims = 40)
neurons <- FindNeighbors(neurons, dims = 1:30, reduction = "pca")
neurons <- FindClusters(neurons, resolution = 2, cluster.name = "neuron_clusters")  
neurons <- RunUMAP(neurons, dims = 1:30, reduction = "pca", reduction.name = "umap")

DimPlot(neurons, reduction = 'umap', label = TRUE)

#How many neurons from each brain region
DimPlot(neurons, reduction = 'umap', group.by = 'Organ')

#Split between wt and ips infection levels
neurons[[]] %>% mutate(virusPresence = ifelse(virusCountPAdj > 4, 'yes', 'no')) %>% 
  group_by(Genotype) %>% 
  dplyr::summarise(virusPresenceProp = mean(virusPresence == 'yes')) %>% 
  ggplot(aes(x = Genotype, y = virusPresenceProp, fill = Genotype)) +
  geom_bar(stat = 'identity')+
  ylab('Virus presence (threshold 1 viral read)')

#Using allen cell atlas: https://knowledge.brain-map.org/celltypes for celltype markers
#and panglao  https://panglaodb.se/markers.html?cell_type=%27Purkinje%20neurons%27

#Inhibitory neuron markers
FeaturePlot(neurons, 'Gad1', reduction = 'umap')
FeaturePlot(neurons, 'Gad2', reduction = 'umap')
FeaturePlot(neurons, 'Dlx6os1', reduction = 'umap')
FeaturePlot(neurons, 'Slc6a1', reduction = 'umap')

#Excitatory markers
FeaturePlot(neurons, 'Sv2b', reduction = 'umap')
FeaturePlot(neurons, 'Arpp21', reduction = 'umap')
FeaturePlot(neurons, 'Slc7a7', reduction = 'umap')
FeaturePlot(neurons, 'Slc7a6', reduction = 'umap')

#Purkinje markers from https://pmc.ncbi.nlm.nih.gov/articles/PMC9497131/
FeaturePlot(neurons, 'Rora', reduction = 'umap')
FeaturePlot(neurons, 'Foxp2', reduction = 'umap')
FeaturePlot(neurons, 'Car8', reduction = 'umap')

#Other markers here https://panglaodb.se/markers.html?cell_type=%27Purkinje%20neurons%27
#and from here, good source https://www.nature.com/articles/s41593-021-00872-y#MOESM1
FeaturePlot(neurons, 'Calb1', reduction = 'umap')
FeaturePlot(neurons, 'Skor2', reduction = 'umap')
FeaturePlot(neurons, 'Itpr1', reduction = 'umap')
FeaturePlot(neurons, 'Pcp2', reduction = 'umap')

#Look at deg markers present.
neuronMarkers <- FindAllMarkers(neurons, assay = 'RNA')
neuronMarkers

#### NUP Expression ####
#split this by brain organ
#Look for NUP
FeaturePlot(neurons, 'Nup98', reduction = 'umap')
FeaturePlot(neurons, 'Nup153', reduction = 'umap')
FeaturePlot(neurons, 'virusPresent', reduction = 'umap')

#Compare between cells w and without virus
#Set virus cutoff greater than 4
neurons$virusPresent <- ifelse(neurons$virusCountPAdj > 2, 1, 0)
table(neurons$Treatment, neurons$virusPresent)

neurons = subset(neurons, Treatment != 'PBS' | virusPresent == 0)
table(neurons$Treatment, neurons$virusPresent)

Nup98Exp <- neurons[['RNA']]$data['Nup98',]
Nup153Exp <- neurons[['RNA']]$data['Nup153',]
nupDat <- data.frame(virus = neurons$virusCountPAdj, virusPresent = neurons$virusPresent,
                     nup98 = Nup98Exp, nup153 = Nup153Exp, organ = neurons$Organ)
nupDat <- nupDat %>% mutate(nup98Present = ifelse(nup98 > 0, 1, 0))
nupDat <- nupDat %>% mutate(nup153Present = ifelse(nup153 > 0, 1, 0))

ggplot(nupDat, aes(x = factor(virusPresent), y = nup98))+
  geom_boxplot() +
  geom_point()

nupDat %>% group_by(virusPresent) %>% dplyr::summarise(nup98Mean = mean(nup98))
table(nupDat$virusPresent)

#Quick function for plotting barplots
nupBarPlot <- function(dat, nupDat, xtitle = '', ytitle = ''){
  dat %>% group_by(virusPresent) %>% dplyr::summarise(nupProportions = mean(!! sym(nupDat))) %>% 
    ggplot(aes(x = factor(virusPresent), y = nupProportions, fill = virusPresent))+
    geom_bar(stat = 'identity')+
    ylab(ytitle) +
    xlab(xtitle) +
    theme(legend.position = 'None')+
    ylim(0,1)
}

nupBarPlot(dat = nupDat, nupDat = 'nup98Present', xtitle = 'Virus Presence', ytitle = 'Proportion of cells expressing Nup98')


ggplot(nupDat, aes(x = factor(virusPresent), y = nup153))+
  geom_boxplot()

nupBarPlot(nupDat, 'nup153Present', xtitle = 'Virus Presence', ytitle = 'Proportion of cells expressing Nup153')


#Split by organ
cerebellum <- nupDat[nupDat$organ == 'Cerebellum',]
cerebrum <- nupDat[nupDat$organ == 'Cerebrum',]

ggplot(cerebellum, aes(x = factor(virusPresent), y = nup98))+
  geom_boxplot() +
  geom_point() +
  xlab('Virus Presence')+
  ylab('Nup98 expression')

nupBarPlot(cerebellum, 'nup98Present', xtitle = 'Virus Presence', ytitle = 'Proportion of cells expressing Nup98')

ggplot(cerebrum, aes(x = factor(virusPresent), y = nup98))+
  geom_boxplot() +
  geom_point() +
  xlab('Virus Presence')+
  ylab('Nup98 expression')


nupBarPlot(cerebrum, 'nup98Present', xtitle = 'Virus Presence', ytitle = 'Proportion of cells expressing Nup98')

#NUP153 barplots

nupBarPlot(cerebellum, 'nup153Present', xtitle = 'Virus Presence', ytitle = 'Proportion of cells expressing Nup153')
nupBarPlot(cerebrum, 'nup153Present', xtitle = 'Virus Presence', ytitle = 'Proportion of cells expressing Nup153')


#Compare excitatory and inhibitory viral expression
neurons$subtype = case_when(neurons$neuron_clusters %in% c(1,0,4,16,9,6,17,7,8) ~ 'excitatory',
                            neurons$neuron_clusters %in%  c(2,3,5,10,11)~'inhibitory')
neurons[[]] %>% mutate(virusPresent = ifelse(neurons$virusCountPAdj > 3, 1, 0)) %>% group_by(subtype) %>% 
  summarise(virusProp = mean(virusPresent)) 

#Subset to just excitatory neuron clusters based on marker genes
exNeurons <- subset(neurons, neuron_clusters %in% c(1,0,4,16,9,6,17,7,8))
FeaturePlot(exNeurons, 'Sv2b', reduction = 'umap')

Nup98ExpEx <- exNeurons[['RNA']]$data['Nup98',]
Nup153ExpEx <- exNeurons[['RNA']]$data['Nup153',]
nupDatEx <- data.frame(virus = exNeurons$virusCountPAdj, virusPresent = exNeurons$virusPresent,
                     nup98 = Nup98ExpEx, nup153 = Nup153ExpEx, organ = exNeurons$Organ)

nupDatEx <- nupDatEx %>% mutate(nup98Present = ifelse(nup98 > 0, 1, 0))
nupDatEx <- nupDatEx %>% mutate(nup153Present = ifelse(nup153 > 0, 1, 0))

cerebellumEx <- nupDatEx[nupDatEx$organ == 'Cerebellum',]
cerebrumEx <- nupDatEx[nupDatEx$organ == 'Cerebrum',]

nupBarPlot(cerebellumEx, 'nup98Present', xtitle = 'Virus Presence', ytitle = 'Proportion of cells expressing Nup98')
nupBarPlot(cerebrumEx, 'nup98Present', xtitle = 'Virus Presence', ytitle = 'Proportion of cells expressing Nup98')

nupBarPlot(cerebellumEx, 'nup153Present', xtitle = 'Virus Presence', ytitle = 'Proportion of cells expressing Nup153')
nupBarPlot(cerebrumEx, 'nup153Present', xtitle = 'Virus Presence', ytitle = 'Proportion of cells expressing Nup153')


# Pretty small group sizes for pseudobulk, might not be able to do too much with these
#Remember pbs neurons that had >2 viral reads were filtered
exNeurons$virusPresent = factor(exNeurons$virusPresent)
neuronPB <- createPseudoBulk(exNeurons, c("Genotype", "Treatment", "Timepoint","Organ", 'virusPresent'))
neuronPB <- DESeq(neuronPB)
resultsNames(neuronPB)


#nothing significant
res <- results(neuronPB, contrast = c('virusPresent', '1', '0'))
res[order(res$padj),] %>% as.data.frame() %>% arrange(padj) %>% head(n = 20)

#Percent neurons infected each time point
neurons[[]] %>% group_by(Timepoint) %>% mutate(virusPresent = ifelse(virusCountPAdj>4, 1, 0)) %>% 
  summarise(virusProp = mean(virusPresent)) %>% ggplot(aes(x = Timepoint, y = virusProp,
                                                           fill = Timepoint))+
  geom_bar(stat = 'identity')+
  theme(legend.position = 'None')

cerebellum <- subset(neurons, Organ == 'Cerebellum')
cerebellum[[]] %>% group_by(Timepoint) %>% mutate(virusPresent = ifelse(virusCountPAdj>2, 1, 0)) %>% 
  summarise(virusProp = mean(virusPresent)) %>% ggplot(aes(x = Timepoint, y = virusProp, fill = Timepoint))+
  geom_bar(stat = 'identity')+
  theme(legend.position = 'None')+
  ylab('Proportion cells with virus')


cerebrum <- subset(neurons, Organ == 'Cerebrum')
cerebrum[[]] %>% group_by(Timepoint) %>% mutate(virusPresent = ifelse(virusCountPAdj>2, 1, 0)) %>% 
  summarise(virusProp = mean(virusPresent)) %>% ggplot(aes(x = Timepoint, y = virusProp,
                                                           fill = Timepoint))+
  geom_bar(stat = 'identity')+
  theme(legend.position = 'None')+
  ylab('Proportion cells with virus')
