#Merging our data with yang's
#Could subset celltypes to make this faster if it's really slow

#Load our data
#Load data
ParseSeuratObj_int <- LoadSeuratRds("~/Documents/ÖverbyLab/data/FilteredRpcaIntegratedDatNoDoublets.rds") 

#Only using cells from gal3 project
wt_cerebrum <-  subset(ParseSeuratObj_int, Treatment %in% c('PBS', 'rLGTV') & Organ == 'Cerebrum' & Genotype == 'WT')

#Plot data
newCols <-  c(brewer.pal(12, 'Paired'), '#99FFE6', '#CE99FF', '#18662E','#737272',  '#FF8AEF')
newCols[11] =  '#FF8AEF'
DimPlot(ParseSeuratObj_int, label = FALSE, group.by = 'manualAnnotation', reduction = 'umap.integrated',
        cols = newCols)

#Load in Yang data
yang_data <- LoadSeuratRds("~/Documents/ÖverbyLab/yang_public_data/yang_rpca_integrated_obj.RDS")
newCols <-  c(brewer.pal(12, 'Paired'), '#99FFE6', '#CE99FF', '#18662E','#737272',  '#FF8AEF')
newCols[11] =  '#FF8AEF'

#Join layers for yang to match parse
yang_data <- JoinLayers(yang_data)

#Plot data
DimPlot(yang_data,reduction = "umap.rpca", label = FALSE, group.by = 'manualAnnotation', cols = newCols)

#Make treatment column names same case before merging
yang_data$Treatment = yang_data$treatment
yang_data$treatment = NULL

#Create celltype subsets to integrate
parse_endo <- subset(wt_cerebrum, manualAnnotation == 'Endothelial')
parse_macro <- subset(wt_cerebrum, manualAnnotation == 'Macrophage/Monocytes')

yang_endo <- subset(yang_data, manualAnnotation == 'Endothelial')
yang_macro <- subset(yang_data, manualAnnotation == 'Macro/Mono')


########### Begin monocyte / macropahge analyses ########### 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 


#monocyte cells combined
macro.combined <- merge(parse_macro, y = yang_macro, add.cell.ids = c("parse", "yang"), project = "merged_macros")

macro.combined <- NormalizeData(macro.combined)
macro.combined <- FindVariableFeatures(macro.combined)
macro.combined <- ScaleData(macro.combined)
macro.combined <- RunPCA(macro.combined)

#Need to allow greater ram usage to run pca integration
options(future.globals.maxSize = 10000 * 1024^2)

#Integrate data
combined_macros <- IntegrateLayers(
  object = macro.combined, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca",
  verbose = FALSE
)

#Reset max ram to not accidentally use a bunch later without realizing
options(future.globals.maxSize = 500 * 1024^2)

ElbowPlot(combined_macros, ndims = 40)

#Prepare data for plotting
combined_macros <- FindNeighbors(combined_macros, reduction = "integrated.rpca", dims = 1:25)
combined_macros <- FindClusters(combined_macros, resolution = 2, cluster.name = "rpca_clusters")
combined_macros$manualAnnotation = 'macro/mono'
combined_macros$dataset <- gsub("_.+", '', colnames(combined_macros))
combined_macros$timepoint_or_treatment <- factor(dplyr::case_when(combined_macros$dataset == 'parse' & combined_macros$Treatment != 'PBS' ~ combined_macros$Timepoint,
                                                           combined_macros$dataset == 'parse' & combined_macros$Treatment == 'PBS' ~ combined_macros$Treatment,
                                                           combined_macros$dataset == 'yang' ~ combined_macros$Treatment),
                                                 levels = c('PBS', 'Day 3', 'Day 4', 'Day 5', 'healthy', 'mild', 'moderate', 'severe'))
combined_macros$dataset_group <- factor(paste(combined_macros$dataset, combined_macros$timepoint_or_treatment, sep = '_'),
                                        levels = c('parse_PBS', 'parse_Day 3', 'parse_Day 4', 'parse_Day 5',
                                                   'yang_healthy', 'yang_mild','yang_moderate', 'yang_severe'))

#Plot datasets to compare umap locations
pdf('/Users/matthewbradley/Documents/ÖverbyLab/yang_public_data/combined_data_plots/combined_mac_plot.pdf', width = 6, height = 5)
DimPlot(combined_macros, label = FALSE, split.by = 'dataset', group.by = 'timepoint_or_treatment', reduction = 'integrated.rpca',
        cols = newCols)+
  ggtitle('Macrophage integrated UMAP')+
  xlab('')+
  ylab('')
dev.off()

pdf('/Users/matthewbradley/Documents/ÖverbyLab/yang_public_data/combined_data_plots/combined_mac_mhc2.pdf', width = 6, height = 5)
DotPlot(combined_macros, features = c('H2-Aa', 'H2-Ab1', 'H2-DMa', 'H2-DMb1', 'H2-DMb2', 'H2-Eb1'), group.by = 'dataset_group', scale = FALSE)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
  ggtitle('mhc-2 mac markers')+
  coord_flip()+
  xlab('')+
  ylab('')+
  scale_color_gradientn(colours = c('white', '#FFD991', '#FF7530', '#FF4024'), 
                        values = c(0, 0.3, 0.6, 1),
                        name = 'Average Expression')
dev.off()

#Monocyte trafficking genes
mono_trafficking_dot_dat <- DotPlot(combined_macros, features = c('Pecam1', 'Pvr', 'Cd99l2', 'Epha1', 'Ephb1', 'Itga4', 'Itgb1', 'Ccr1', 'Ccr2', 'Ccr5'), 
                                    group.by = 'dataset_group', scale = FALSE)$data

pdf('/Users/matthewbradley/Documents/ÖverbyLab/yang_public_data/combined_data_plots/mono_trafficking_genes.pdf', width = 6, height = 5)
ggplot(mono_trafficking_dot_dat, aes(x = id, y = features.plot))+
  geom_point(aes(size = pct.exp, fill = avg.exp.scaled), pch = 21)+
  #coord_flip()+
  scale_size_continuous(range = c(0.5,6), limits = c(0,100))+
  scale_fill_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white"), 
                       values = c(1.0,0.7,0.4,0),
                       limits = c(0,2.5))+
  ggtitle('Monocyte trafficking genes')+
  theme_bw()+
  theme(panel.grid = element_blank(),
        #panel.border = element_rect(colour = "white", fill = NA),
        panel.spacing = unit(0, "line"),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90))+
  scale_size(range = c(0,9), limits = c(0, 100))+
  ylab('')+
  xlab('')
dev.off()

########### Begin endothelial analyses ########### 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

#endothelial cells combined
endo.combined <- merge(parse_endo, y = yang_endo, add.cell.ids = c("parse", "yang"), project = "merged_endos")

endo.combined <- NormalizeData(endo.combined)
endo.combined <- FindVariableFeatures(endo.combined)
endo.combined <- ScaleData(endo.combined)
endo.combined <- RunPCA(endo.combined)

#Need to allow greater ram usage to run pca integration
options(future.globals.maxSize = 10000 * 1024^2)

#Integrate data
combined_endos <- IntegrateLayers(
  object = endo.combined, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca",
  verbose = FALSE
)

#Reset max ram to not accidentally use a bunch later without realizing
options(future.globals.maxSize = 500 * 1024^2)

ElbowPlot(combined_endos, ndims = 40)

#Prepare data for plotting
combined_endos <- FindNeighbors(combined_endos, reduction = "integrated.rpca", dims = 1:25)
combined_endos <- FindClusters(combined_endos, resolution = 2, cluster.name = "rpca_clusters")
combined_endos$manualAnnotation = 'endothelial'
combined_endos$dataset <- gsub("_.+", '', colnames(combined_endos))
combined_endos$timepoint_or_treatment <- factor(dplyr::case_when(combined_endos$dataset == 'parse' & combined_endos$Treatment != 'PBS' ~ combined_endos$Timepoint,
                                                                  combined_endos$dataset == 'parse' & combined_endos$Treatment == 'PBS' ~ combined_endos$Treatment,
                                                                  combined_endos$dataset == 'yang' ~ combined_endos$Treatment),
                                                 levels = c('PBS', 'Day 3', 'Day 4', 'Day 5', 'healthy', 'mild', 'moderate', 'severe'))
combined_endos$dataset_group <- factor(paste(combined_endos$dataset, combined_endos$timepoint_or_treatment, sep = '_'),
                                        levels = c('parse_PBS', 'parse_Day 3', 'parse_Day 4', 'parse_Day 5',
                                                   'yang_healthy', 'yang_mild','yang_moderate', 'yang_severe'))
#Plot data together
pdf('/Users/matthewbradley/Documents/ÖverbyLab/yang_public_data/combined_data_plots/combined_endothelial_plot.pdf', width = 6, height = 5)
DimPlot(combined_endos, label = FALSE, split.by = 'dataset', group.by = 'timepoint_or_treatment', reduction = 'integrated.rpca',
        cols = newCols)+
  ggtitle('Endothelial integrated UMAP')+
  xlab('')+
  ylab('')
dev.off()

