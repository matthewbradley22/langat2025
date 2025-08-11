#Explore neurons in data

#Load packages
library(Seurat)
library(ggplot2)
library(dplyr)

#Load in data
ParseSeuratObj_int <- LoadSeuratRds("./data/FilteredRpcaIntegratedDat.rds")

#filter to singlets. went through analysis and filtering made sense, doublets mostly clustered together 
ParseSeuratObj_int <- subset(ParseSeuratObj_int, scDblFinderLabel == 'singlet')

DimPlot(ParseSeuratObj_int, reduction = 'umap.integrated', label = TRUE) +
  NoLegend()

DimPlot(ParseSeuratObj_int, reduction = 'umap.integrated', label = TRUE, group.by = 'singleR_labels') +
  NoLegend()

#Neurons appear to make up clusters 24, 24, 37, and 44
table(subset(ParseSeuratObj_int, singleR_labels == 'Neurons')$seurat_clusters)

#Look at proportions of neurons in each cluster, if low proportion then consider looking closer
table(subset(ParseSeuratObj_int, seurat_clusters == '27')$singleR_labels)
table(subset(ParseSeuratObj_int, seurat_clusters == '24')$singleR_labels)
table(subset(ParseSeuratObj_int, seurat_clusters == '37')$singleR_labels)
table(subset(ParseSeuratObj_int, seurat_clusters == '44')$singleR_labels)

#Look at markers in neurons
neuronMarkers <- FindMarkers(ParseSeuratObj_int, group.by = 'singleR_labels',
                             ident.1 = 'Neurons')

#Can start with conservative approach of subsetting by celltype and cluster, removing doublets
neurons <- subset(ParseSeuratObj_int, singleR_labels == 'Neurons' & seurat_clusters %in% c('27', '24',
                                                                                          '37', '44'))


table(neurons$Genotype)
table(neurons$Treatment)
table(neurons$Timepoint)

#More neurons in cerebellum is interesting
table(neurons$Organ)

#Check virus presence. Make this slightly stricter later on (need > 1 virus) for closer analysis
neurons$virusPresent <- ifelse(neurons$virusCountPAdj > 0, 1, 0)
table(neurons$virusPresent)
table(neurons$virusPresent, neurons$Organ)

#Rerun through data processing and visualization
neurons <- NormalizeData(neurons)
neurons <- FindVariableFeatures(neurons)
neurons <- ScaleData(neurons)
neurons <- RunPCA(neurons)
ElbowPlot(neurons, ndims = 40)
neurons <- FindNeighbors(neurons, dims = 1:30, reduction = "pca")
neurons <- FindClusters(neurons, resolution = 2, cluster.name = "neuron_clusters")  
neurons <- RunUMAP(neurons, dims = 1:30, reduction = "pca", reduction.name = "umap")

DimPlot(neurons, reduction = 'umap')

#How many neurons from each brain region
table(neurons$Organ)

#Using allen cell atlas: https://knowledge.brain-map.org/celltypes for celltype markers
#and panglao  https://panglaodb.se/markers.html?cell_type=%27Purkinje%20neurons%27

#Inhibitory neuron markers - don't really show up
FeaturePlot(neurons, 'Gad1', reduction = 'umap')
FeaturePlot(neurons, 'Gad2', reduction = 'umap')
FeaturePlot(neurons, 'Dlx6os1', reduction = 'umap')
FeaturePlot(neurons, 'Dlx6os1', reduction = 'umap')
FeaturePlot(neurons, 'Slc6a1', reduction = 'umap')

#Excitatory markers. Looks like they're all excitatory
FeaturePlot(neurons, 'Sv2b', reduction = 'umap')
FeaturePlot(neurons, 'Arpp21', reduction = 'umap')
FeaturePlot(neurons, 'Slc7a7', reduction = 'umap')
FeaturePlot(neurons, 'Slc7a6', reduction = 'umap')

#Purkinje markers from https://pmc.ncbi.nlm.nih.gov/articles/PMC9497131/
FeaturePlot(neurons, 'Rora', reduction = 'umap')
FeaturePlot(neurons, 'Foxp2', reduction = 'umap')
FeaturePlot(neurons, 'Car8', reduction = 'umap')

#Other markers here https://panglaodb.se/markers.html?cell_type=%27Purkinje%20neurons%27
FeaturePlot(neurons, 'Calb1', reduction = 'umap')


#Look at deg markers present.
neuronMarkers <- FindAllMarkers(neurons, assay = 'RNA')
neuronMarkers

#### NUP Expression ####
#split this by brain organ
#Look for NUP
FeaturePlot(neurons, 'Nup98', reduction = 'umap')
FeaturePlot(neurons, 'Nup153', reduction = 'umap')

#Compare between cells w and without virus
#Set virus cutoff greater than 1
neurons$virusPresent <- ifelse(neurons$virusCountPAdj > 1, 1, 0)
table(neurons$Treatment, neurons$virusPresent)

neurons = subset(neurons, Treatment != 'PBS' | virusPresent == 0)
table(neurons$Treatment, neurons$virusPresent)

Nup98Exp <- neurons[['RNA']]$data['Nup98',]
Nup153Exp <- neurons[['RNA']]$data['Nup153',]
nupDat <- data.frame(virus = neurons$virusCountPAdj, virusPresent = neurons$virusPresent,
                     nup98 = Nup98Exp, nup153 = Nup153Exp)
nupDat <- nupDat %>% mutate(nup98Present = ifelse(nup98 > 0, 1, 0))
nupDat <- nupDat %>% mutate(nup153Present = ifelse(nup153 > 0, 1, 0))
#Median of nup98 much higher in viral infected cells
ggplot(nupDat, aes(x = factor(virusPresent), y = nup98))+
  geom_boxplot() +
  geom_point()

nupDat %>% group_by(virusPresent) %>% dplyr::summarise(nup98Mean = mean(nup98))

nupDat %>% group_by(virusPresent) %>% dplyr::summarise(nupProportions = mean(nup98Present)) %>% 
  ggplot(aes(x = factor(virusPresent), y = nupProportions, fill = virusPresent))+
  geom_bar(stat = 'identity')+
  ylab('Proportion of cells expressing Nup98') +
  xlab('Virus Presence') +
  theme(legend.position = 'None')+
  ylim(0,1)

table(nupDat$virusPresent)

ggplot(nupDat, aes(x = factor(virusPresent), y = nup153))+
  geom_boxplot()

nupDat %>% group_by(virusPresent) %>% dplyr::summarise(nup153 <- mean(nup153))

nupDat %>% group_by(virusPresent) %>% dplyr::summarise(nupProportions = mean(nup153Present)) %>% 
  ggplot(aes(x = factor(virusPresent), y = nupProportions, fill = virusPresent))+
  geom_bar(stat = 'identity')+
  ylab('Proportion of cells expressing Nup153') +
  xlab('Virus Presence') +
  theme(legend.position = 'None')+
  ylim(0,1)


