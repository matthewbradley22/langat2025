#Look at necroptosis/pyroptosis genes from email
source('~/Documents/ÖverbyLab//scripts/langatFunctions.R')

#Load data
ParseSeuratObj_int <- LoadSeuratRds("~/Documents/ÖverbyLab/data/FilteredRpcaIntegratedDatNoDoublets.rds") 
ParseSeuratObj_int$hasVirus = ifelse(ParseSeuratObj_int$virusCountPAdj >= 10, 1, 0)

#Check data
newCols <-  c(brewer.pal(12, 'Paired'), '#99FFE6', '#CE99FF', '#18662E','#737272',  '#FF8AEF')
newCols[11] =  '#FF8AEF'
DimPlot(ParseSeuratObj_int, label = FALSE, group.by = 'manualAnnotation', reduction = 'umap.integrated',
        cols = newCols)
ParseSeuratObj_int$time_treatment <- paste(ParseSeuratObj_int$Timepoint, ParseSeuratObj_int$Treatment, sep = '_')
#Gene lists
necroptosis <- c('Tnf', 'Tnfrsf1a', 'Ripk2', 'Mlkl', 'Ripk1', 'Ripk3')
pyoptosis <- c('Gsdmc', 'Nlrp3', 'Aim2', 'Gsdmd', 'Il18', 'Il1b', 'Casp9', 'Casp8', 'Casp6', 'Casp3', 'Casp4', 'Casp1')


DotPlot(ParseSeuratObj_int, features = c(necroptosis, pyoptosis), group.by = 'Timepoint', scale = FALSE)+
  coord_flip()

#Subset to same cells as in gal3 project
wt_cerebrum <- subset(ParseSeuratObj_int, Treatment %in% c('PBS', 'rLGTV') & Organ == 'Cerebrum' & Genotype == 'WT')

wt_cerebrum <- prepSeuratObj(wt_cerebrum)
ElbowPlot(wt_cerebrum, ndims = 40)
wt_cerebrum <- prepUmapSeuratObj(wt_cerebrum, nDims = 20, reductionName = 'wt.cerebrum.umap')

DimPlot(wt_cerebrum, reduction = 'wt.cerebrum.umap', label = TRUE)

#Split by treatment
wt_cerebrum_pbs <- subset(wt_cerebrum, Treatment == 'PBS')
#Don't plot celltypes with less than 50 cells
celltypes_for_plot <- names(table(wt_cerebrum_pbs$manualAnnotation)[table(wt_cerebrum_pbs$manualAnnotation) > 50])
wt_cerebrum_pbs_subset <- subset(wt_cerebrum_pbs, manualAnnotation %in% celltypes_for_plot)
pdf("~/Documents/ÖverbyLab/scPlots/necroptosis_genes_sc/wt_cerebrum_pbs_necroptosis_dot.pdf", width = 9, height = 6)
DotPlot(wt_cerebrum_pbs_subset, features = c(necroptosis, pyoptosis), group.by = 'manualAnnotation', scale = FALSE)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
  ggtitle("WT cerebrum PBS necroptosis and pyoptosis genes")
dev.off()

wt_cerebrum_lgtv <- subset(wt_cerebrum, Treatment == 'rLGTV')
celltypes_for_plot_lgtv <- names(table(wt_cerebrum_lgtv$manualAnnotation)[table(wt_cerebrum_lgtv$manualAnnotation) > 50])
wt_cerebrum_lgtv_subset <- subset(wt_cerebrum_lgtv, manualAnnotation %in% celltypes_for_plot_lgtv)
pdf("~/Documents/ÖverbyLab/scPlots/necroptosis_genes_sc/wt_cerebrum_LGTV_necroptosis_dot.pdf", width = 9, height = 6)
DotPlot(wt_cerebrum_lgtv, features = c(necroptosis, pyoptosis), group.by = 'manualAnnotation', scale = FALSE)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
  ggtitle("WT cerebrum LGTV necroptosis and pyoptosis genes")
dev.off()

mac_mono <- subset(wt_cerebrum, manualAnnotation == 'Macrophage/Monocytes')
DotPlot(mac_mono, features = c(necroptosis, pyoptosis), group.by = 'time_treatment', scale = FALSE)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
  ggtitle("Monocyte macrophage LGTV necroptosis and pyoptosis genes")
table(mac_mono$time_treatment)

#Split up data to look at gene expression across variables
