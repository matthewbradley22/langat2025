#Packages and functions
library(Seurat)
library(gprofiler2)
library(UCell)
library(rstatix)
library(msigdbr)
source('~/Documents/ÖverbyLab//scripts/langatFunctions.R')

#Load data
ParseSeuratObj_int <- LoadSeuratRds("~/Documents/ÖverbyLab/data/FilteredRpcaIntegratedDatNoDoublets.rds") 
ParseSeuratObj_int$hasVirus = ifelse(ParseSeuratObj_int$virusCountPAdj >= 10, 1, 0)

#Check data
umap_color_list <- c( "#7047A1", "#B370AE","#292270",  "#166DF0","#6D92F8",  "#6DC3F8", "#8a0000","#F76363", "#FF96A2", 
                      "#D6644B", "#F08C3A", "#fdc087","#074F00", "#208d1f","#7bcd79", 
                      "gray")

ParseSeuratObj_int$manualAnnotation <- factor(ParseSeuratObj_int$manualAnnotation, 
                                              levels = c('Astrocytes', 'Choroid Plexus', 'Endothelial',
                                                         'Ependymal', 'Immature Neurons', 'Microglia', 'Muscle cells',
                                                         'Neurons', 'Oligodendrocytes', 'Pericytes', 'B Cells',
                                                         'Granulocytes', 'Macrophage/Monocytes', 'Nk cells', 
                                                         'T cells', 'unknown'))

DimPlot(ParseSeuratObj_int, label = FALSE, group.by = 'manualAnnotation', reduction = 'umap.integrated',
        cols = umap_color_list)


#Read in fully corrected viral levels
adjusted_viral_reads <- list.files('~/Documents/ÖverbyLab/viralReadCountsFullyAdjusted/', full.names = TRUE)
adjusted_viral_reads_df <- lapply(adjusted_viral_reads, FUN = function(x){
  dat <- readr::read_table(x, col_names = FALSE)
  colnames(dat) = c('virusCountFinalAdj', 'barcode')
  dat$barcode = substring(dat$barcode, 6)
  dat
})

for(i in 1:length(adjusted_viral_reads_df)){
  adjusted_viral_reads_df[[i]]$sublib = paste0('seuObj', i)
}

adjusted_viral_reads_binded <- do.call(rbind, adjusted_viral_reads_df)

ParseSeuratObj_int[[]] <- left_join(ParseSeuratObj_int[[]], adjusted_viral_reads_binded, by = c('cell' = 'barcode', 'subLib' = 'sublib'))
ParseSeuratObj_int$virusCountFinalAdj <- ifelse(is.na(ParseSeuratObj_int$virusCountFinalAdj), 0, ParseSeuratObj_int$virusCountFinalAdj)
ParseSeuratObj_int[[]][which(ParseSeuratObj_int$virusCountFinalAdj != ParseSeuratObj_int$virusCountPAdj),][c('virusCount', 'virusCountPAdj', 'virusCountFinalAdj')]

#Normalize viral count by cell
parse_counts <- ParseSeuratObj_int[['RNA']]$counts
viral_read_matrix <- matrix(data = ParseSeuratObj_int$virusCountFinalAdj, ncol = length(ParseSeuratObj_int$virusCountFinalAdj))
rownames(viral_read_matrix) = 'lgtv_final'
parse_counts_with_virus <- rbind(parse_counts, viral_read_matrix)

parse_counts_with_virus_norm <- NormalizeData(parse_counts_with_virus)

ParseSeuratObj_int$virus_count_normalized <- parse_counts_with_virus_norm['lgtv_final',]
#value 2 comes from viralCoverages.R script
ParseSeuratObj_int$virusCountFinalAdj_corrected = ParseSeuratObj_int$virus_count_normalized - 2


#Look at various flavivirus receptors, and compare to infection levels + lgtv to compare infection level with receptor
receptors <- c('Axl', 'Havcr1', 'Havcr2', 'Timd4', 'Tyro3', 'Lrp8', 'Lrp1',
               'Lrp4', 'Cd209a', 'Cd209b', 'Cd209c', 'Itgb4', 'Hspa5', 'Ncam1', 'Hspa1a', 'Vim', 'Itgav',
               'Itgb3', 'Cldn1', 'Clec5a', 'Mrc1', 'Mer', 'Scarb1')

#First look at wt immature neurons
wt_im_neurons <- subset(ParseSeuratObj_int, manualAnnotation == 'Immature Neurons' & Genotype == 'WT')

wt_im_neurons[[]] <- wt_im_neurons[[]] %>% dplyr::mutate(virus_present = ifelse(virus_count_normalized > 2, 'yes', 'no'))


DotPlot(wt_im_neurons, features = receptors, group.by = 'virus_present', scale = FALSE) +
  coord_flip()

ParseSeuratObj_int$treatment_celltype <- paste(ParseSeuratObj_int$Treatment, ParseSeuratObj_int$manualAnnotation, sep = '_')

DotPlot(ParseSeuratObj_int, features = receptors, group.by = 'treatment_celltype', scale = FALSE)$data %>% 
  tidyr::separate(col = 'id', into = c('treatment', 'celltype'), sep = '_') %>% 
  dplyr::filter(treatment != 'rChLGTV') %>% 
  #create 0 data for havcr1 since we don't have it
  #add_row(avg.exp = 0, pct.exp = NA, features.plot = 'Havcr1', treatment = 'FALSE', celltype = 'Astrocytes', avg.exp.scaled = NA) %>% 
  dplyr::mutate(infected = factor(ifelse(treatment == 'PBS', yes = 'Mock', no = 'LGTV'), levels = c('Mock', 'LGTV'))) %>% 
  ggplot(aes(x = features.plot, y = celltype, fill = avg.exp.scaled, size = pct.exp))+
  geom_point(pch = 21)+
  facet_wrap(~infected)+
  theme_classic()+
  scale_fill_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white"),
                       values = c(1.0,0.7,0.4,0))+
  xlab('')+
  ylab('')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(fill = NA, color = "black", linetype = "dashed"))


#Look at receptors in all wt celltypes
wt_mock <- subset(ParseSeuratObj_int, Genotype == 'WT' & Treatment == 'PBS')

DotPlot(wt_mock, features = receptors, group.by = 'manualAnnotation', scale = FALSE)$data %>% 
  dplyr::filter(id != 'unknown') %>% 
  ggplot(aes(x = features.plot, y = id, fill = avg.exp.scaled, size = pct.exp))+
  geom_point(pch = 21)+
  theme_classic()+
  scale_fill_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white"),
                       values = c(1.0,0.7,0.4,0))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 13))+
  ylab('')+
  xlab('')

#Microglia subcluster receptors
wt_mock_micro <- subset(wt_mock, manualAnnotation == 'Microglia')

wt_mock_micro <- prepSeuratObj(wt_mock_micro)
ElbowPlot(wt_mock_micro, ndims = 40)
wt_mock_micro <- prepUmapSeuratObj(wt_mock_micro, nDims = 15, reductionName = 'micro.umap', resolution_value = 0.6)

DimPlot(wt_mock_micro, reduction = 'micro.umap', label = FALSE, label.size = 6)

DotPlot(wt_mock_micro, features = receptors, group.by = 'seurat_clusters', scale = FALSE)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_color_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white"),
                        values = c(1.0,0.7,0.4,0))+
  xlab('')+
  ylab('')+
  scale_size(range = c(0,8))

#look at receptores in knockout
mavs <- subset(ParseSeuratObj_int, Genotype == 'IPS1')
DotPlot(mavs, features = receptors, group.by = 'manualAnnotation', scale = FALSE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#CP subcluster receptors
ips_mock_cp <- subset(ParseSeuratObj_int, manualAnnotation == 'Choroid Plexus' & Genotype == 'IPS1')

ips_mock_cp <- prepSeuratObj(ips_mock_cp)
ElbowPlot(ips_mock_cp, ndims = 40)
ips_mock_cp <- prepUmapSeuratObj(ips_mock_cp, nDims = 15, reductionName = 'cp.umap', resolution_value = 0.6)

DimPlot(ips_mock_cp, reduction = 'cp.umap', label = FALSE, label.size = 6)

DotPlot(ips_mock_cp, features = receptors, group.by = 'seurat_clusters', scale = FALSE)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_color_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white"),
                        values = c(1.0,0.7,0.4,0))+
  xlab('')+
  ylab('')+
  scale_size(range = c(0,8))




