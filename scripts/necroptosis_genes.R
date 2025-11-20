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
ParseSeuratObj_int$treatment_celltype <- paste(ParseSeuratObj_int$Treatment, ParseSeuratObj_int$manualAnnotation, sep = '_')
#Gene lists
necroptosis <- c('Tnf', 'Tnfrsf1a', 'Ripk2', 'Mlkl', 'Ripk1', 'Ripk3')
pyroptosis <- c('Gsdmc', 'Nlrp3', 'Aim2', 'Gsdmd', 'Il18', 'Il1b', 'Casp9', 'Casp8', 'Casp6', 'Casp3', 'Casp4', 'Casp1')


DotPlot(ParseSeuratObj_int, features = c(necroptosis, pyroptosis), group.by = 'Timepoint', scale = FALSE)+
  coord_flip()

#Subset to same cells as in gal3 project
wt_cerebrum <- subset(ParseSeuratObj_int, Treatment %in% c('PBS', 'rLGTV') & Organ == 'Cerebrum' & Genotype == 'WT')

wt_cerebrum <- prepSeuratObj(wt_cerebrum)
ElbowPlot(wt_cerebrum, ndims = 40)
wt_cerebrum <- prepUmapSeuratObj(wt_cerebrum, nDims = 20, reductionName = 'wt.cerebrum.umap')

DimPlot(wt_cerebrum, reduction = 'wt.cerebrum.umap', label = TRUE)

#Want to do resident cells separate from infiltrating, so make lists here
resident_celltypes <- c('Astrocytes', 'Oligodendrocytes', 'Microglia', 'Endothelial', 'Choroid Plexus',
                        'Immature Neurons', 'Ependymal', 'Pericytes', 'Muscle cells', 'Neurons')

infiltrating_celltypes <- c("T cells", 'Macrophage/Monocytes', 'Nk cells', 'Granulocytes', 'B Cells')

###################Resident cell analysis###################
############################################################
#Split by treatment
wt_cerebrum_resident <- subset(wt_cerebrum, manualAnnotation %in% resident_celltypes)

#Not sure that this actually works, cannot visualize both celltype and gene well
celltype_treatment_necroptosis_scores <- list()
for(i in 1:length(necroptosis)){
  avg_data <- AverageExpression(wt_cerebrum_resident, necroptosis[i], group.by = 'treatment_celltype', assay = 'RNA', slot = 'data')$RNA
  celltype_treatment_necroptosis_scores[[i]] <- (t(avg_data))
}

wt_cerebrum_resident_nScores <- do.call(cbind, celltype_treatment_necroptosis_scores) %>% as.data.frame() %>% 
  tibble::rownames_to_column(var = "id") %>% 
  pivot_longer(!id, names_to = 'gene', values_to = 'avg_exp')
split_names <- str_split_fixed(wt_cerebrum_resident_nScores$id, "-", n = 2)
colnames(split_names) <- c('treatment', 'celltype')
wt_cerebrum_resident_nScores <- cbind(wt_cerebrum_resident_nScores, split_names)

#Differences between treatment and pbs avg expression for each gene
wt_cerebrum_resident_nScores %>%
  dplyr::group_by(gene, celltype) %>%
  dplyr::mutate(exp_change = avg_exp - first(avg_exp)) %>% 
  dplyr::filter(treatment == 'rLGTV') %>% 
  ggplot(aes(x = gene, y = celltype, fill = exp_change))+
  geom_tile()+
  scale_fill_gradientn(colors = c("purple", "white", "orange", "red"),
                       values = scales::rescale(c(-0.2, 0, 1, 2)))+
  ggtitle("Necroptosis gene LGTV - PBS difference")+
  geom_text(aes(label=round(exp_change, digits = 2)))

#check cell counts for each group
t(table(wt_cerebrum_resident$Treatment, wt_cerebrum_resident$manualAnnotation))

#Same thing for pyroptosis genes
celltype_treatment_pyroptosis_scores <- list()
for(i in 1:length(pyroptosis)){
  avg_data <- AverageExpression(wt_cerebrum_resident, pyroptosis[i], group.by = 'treatment_celltype', assay = 'RNA', slot = 'data')$RNA
  celltype_treatment_pyroptosis_scores[[i]] <- (t(avg_data))
}

wt_cerebrum_resident_pScores <- do.call(cbind, celltype_treatment_pyroptosis_scores) %>% as.data.frame() %>% 
  tibble::rownames_to_column(var = "id") %>% 
  pivot_longer(!id, names_to = 'gene', values_to = 'avg_exp')
split_names <- str_split_fixed(wt_cerebrum_resident_pScores$id, "-", n = 2)
colnames(split_names) <- c('treatment', 'celltype')
wt_cerebrum_resident_pScores <- cbind(wt_cerebrum_resident_pScores, split_names)

#Differences between treatment and pbs avg expression for each gene
wt_cerebrum_resident_pScores %>%
  dplyr::group_by(gene, celltype) %>%
  dplyr::mutate(exp_change = avg_exp - first(avg_exp)) %>% 
  dplyr::filter(treatment == 'rLGTV') %>% 
  ggplot(aes(x = gene, y = celltype, fill = exp_change))+
  geom_tile()+
  scale_fill_gradientn(colors = c("purple", "white", "orange", "red"),
                       values = scales::rescale(c(-0.2, 0, 1, 2)))+
  ggtitle("Pyroptosis gene LGTV - PBS difference")+
  geom_text(aes(label=round(exp_change, digits = 2)))

#Are any of these significant DEGs? Create pseudobulk object and check significance
wt_cerebrum_resident_bulk <- createPseudoBulk(wt_cerebrum_resident, c('Treatment', 'Timepoint'))
wt_cerebrum_resident_bulk <- DESeq(wt_cerebrum_resident_bulk)
wt_cerebrum_resident_bulk_res <- results(wt_cerebrum_resident_bulk, name = 'Treatment_rLGTV_vs_PBS')
wt_cerebrum_resident_bulk_res[necroptosis,]
wt_cerebrum_resident_bulk_res[pyroptosis,] %>% as.data.frame() %>% dplyr::filter(padj < 0.01)

###################Infiltrating cell analysis###################
################################################################
wt_cerebrum_infiltrating <- subset(wt_cerebrum, manualAnnotation %in% infiltrating_celltypes)

#Check cell counts for each group
t(table(wt_cerebrum_infiltrating$Treatment, wt_cerebrum_infiltrating$manualAnnotation))

#necroptosis genes
celltype_treatment_necroptosis_scores_infil <- list()
for(i in 1:length(necroptosis)){
  avg_data <- AverageExpression(wt_cerebrum_infiltrating, necroptosis[i], group.by = 'treatment_celltype', assay = 'RNA', slot = 'data')$RNA
  celltype_treatment_necroptosis_scores_infil[[i]] <- (t(avg_data))
}

wt_cerebrum_infil_nScores <- do.call(cbind, celltype_treatment_necroptosis_scores_infil) %>% as.data.frame() %>% 
  tibble::rownames_to_column(var = "id") %>% 
  pivot_longer(!id, names_to = 'gene', values_to = 'avg_exp')
split_names <- str_split_fixed(wt_cerebrum_infil_nScores$id, "-", n = 2)
colnames(split_names) <- c('treatment', 'celltype')
wt_cerebrum_infil_nScores <- cbind(wt_cerebrum_infil_nScores, split_names)

#Differences between treatment and pbs avg expression for each gene
wt_cerebrum_infil_nScores %>%
  dplyr::group_by(gene, celltype) %>%
  dplyr::mutate(exp_change = avg_exp - first(avg_exp)) %>% 
  dplyr::filter(treatment == 'rLGTV') %>% 
  ggplot(aes(x = gene, y = celltype, fill = exp_change))+
  geom_tile()+
  scale_fill_gradientn(colors = c("purple", "white", "orange", "red"),
                       values = scales::rescale(c(-0.2, 0, 1, 2)))+
  ggtitle("Necroptosis gene LGTV - PBS difference")+
  geom_text(aes(label=round(exp_change, digits = 2)))

#pyroptosis
#necroptosis genes
celltype_treatment_pyroptosis_scores_infil <- list()
for(i in 1:length(pyroptosis)){
  avg_data <- AverageExpression(wt_cerebrum_infiltrating, pyroptosis[i], group.by = 'treatment_celltype', assay = 'RNA', slot = 'data')$RNA
  celltype_treatment_pyroptosis_scores_infil[[i]] <- (t(avg_data))
}

wt_cerebrum_infil_pScores <- do.call(cbind, celltype_treatment_pyroptosis_scores_infil) %>% as.data.frame() %>% 
  tibble::rownames_to_column(var = "id") %>% 
  pivot_longer(!id, names_to = 'gene', values_to = 'avg_exp')
split_names <- str_split_fixed(wt_cerebrum_infil_pScores$id, "-", n = 2)
colnames(split_names) <- c('treatment', 'celltype')
wt_cerebrum_infil_pScores <- cbind(wt_cerebrum_infil_pScores, split_names)

#Differences between treatment and pbs avg expression for each gene
wt_cerebrum_infil_pScores %>%
  dplyr::group_by(gene, celltype) %>%
  dplyr::mutate(exp_change = avg_exp - first(avg_exp)) %>% 
  dplyr::filter(treatment == 'rLGTV') %>% 
  ggplot(aes(x = gene, y = celltype, fill = exp_change))+
  geom_tile()+
  scale_fill_gradientn(colors = c("purple", "white", "orange", "red"),
                       values = scales::rescale(c(-0.2, 0, 1, 2)))+
  ggtitle("Pyroptosis gene LGTV - PBS difference")+
  geom_text(aes(label=round(exp_change, digits = 2)))


#Are any of these significant DEGs? Create pseudobulk object and check significance
wt_cerebrum_infil_bulk <- createPseudoBulk(wt_cerebrum_infiltrating, c('Treatment', 'Timepoint'))
wt_cerebrum_infil_bulk <- DESeq(wt_cerebrum_infil_bulk)
wt_cerebrum_infil_bulk_res <- results(wt_cerebrum_infil_bulk, name = 'Treatment_rLGTV_vs_PBS')
wt_cerebrum_infil_bulk_res[necroptosis,]
wt_cerebrum_infil_bulk_res[pyroptosis,] %>% as.data.frame() %>% dplyr::filter(padj < 0.01)

#Don't plot celltypes with less than 50 cells
celltypes_for_plot <- names(table(wt_cerebrum_pbs$manualAnnotation)[table(wt_cerebrum_pbs$manualAnnotation) > 50])
wt_cerebrum_pbs_subset <- subset(wt_cerebrum_pbs, manualAnnotation %in% celltypes_for_plot)
pdf("~/Documents/ÖverbyLab/scPlots/necroptosis_genes_sc/wt_cerebrum_pbs_necroptosis_dot.pdf", width = 9, height = 6)
DotPlot(wt_cerebrum_pbs_subset, features = c(necroptosis, pyroptosis), group.by = 'manualAnnotation', scale = FALSE)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
  ggtitle("WT cerebrum PBS necroptosis and pyroptosis genes")
dev.off()

wt_cerebrum_lgtv <- subset(wt_cerebrum, Treatment == 'rLGTV')
celltypes_for_plot_lgtv <- names(table(wt_cerebrum_lgtv$manualAnnotation)[table(wt_cerebrum_lgtv$manualAnnotation) > 50])
wt_cerebrum_lgtv_subset <- subset(wt_cerebrum_lgtv, manualAnnotation %in% celltypes_for_plot_lgtv)
pdf("~/Documents/ÖverbyLab/scPlots/necroptosis_genes_sc/wt_cerebrum_LGTV_necroptosis_dot.pdf", width = 9, height = 6)
DotPlot(wt_cerebrum_lgtv, features = c(necroptosis, pyroptosis), group.by = 'manualAnnotation', scale = FALSE)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
  ggtitle("WT cerebrum LGTV necroptosis and pyroptosis genes")
dev.off()

mac_mono <- subset(wt_cerebrum, manualAnnotation == 'Macrophage/Monocytes')
DotPlot(mac_mono, features = c(necroptosis, pyroptosis), group.by = 'time_treatment', scale = FALSE)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
  ggtitle("Monocyte macrophage LGTV necroptosis and pyroptosis genes")
table(mac_mono$time_treatment)

#Split up data to look at gene expression across variables
