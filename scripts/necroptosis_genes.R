#Look at necroptosis/pyroptosis genes from email

#Load packages
library(Seurat)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(tidyr)
library(stringr)
library(gprofiler2)
source('~/Documents/ÖverbyLab/scripts/langatFunctions.R')

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

#Remove cells that appear to be falsely labelled macrophages - can be seen running through gal3Project.R script
#or load in from csv
false_macs_to_remove <- read.csv("~/Documents/ÖverbyLab/false_macs_to_remove.csv")[[2]]
#should be 409 cells
length(false_macs_to_remove)

ParseSeuratObj_int <- subset(ParseSeuratObj_int, cells = false_macs_to_remove, invert = TRUE)

#Subset to same cells as in gal3 project
wt_cerebrum_day5 <-  subset(ParseSeuratObj_int, Treatment %in% c('PBS', 'rLGTV') & Organ == 'Cerebrum' & Genotype == 'WT' & Timepoint == 'Day 5')

wt_cerebrum_day5 <- prepSeuratObj(wt_cerebrum_day5)
ElbowPlot(wt_cerebrum_day5, ndims = 40)
wt_cerebrum_day5 <- prepUmapSeuratObj(wt_cerebrum_day5, nDims = 20, reductionName = 'wt.cerebrum.umap')

DimPlot(wt_cerebrum_day5, reduction = 'wt.cerebrum.umap', label = TRUE)
DimPlot(wt_cerebrum_day5, reduction = 'wt.cerebrum.umap', group.by = 'manualAnnotation', cols = newCols)+
  ggtitle("WT Cerebrum Day 5")+
  xlab('Umap1')+
  ylab('Umap1')+
  theme(axis.ticks = element_blank(),
        axis.text=element_blank(),
        legend.text=element_text(size=17))

#Want to do resident cells separate from infiltrating, so make lists here
resident_celltypes <- c('Astrocytes', 'Oligodendrocytes', 'Microglia', 'Endothelial', 'Choroid Plexus',
                        'Immature Neurons', 'Ependymal', 'Pericytes', 'Muscle cells', 'Neurons')

infiltrating_celltypes <- c("T cells", 'Macrophage/Monocytes', 'Nk cells', 'Granulocytes', 'B Cells')

###################Resident cell analysis###################
############################################################
#Split by treatment
wt_cerebrum_day5_resident <- subset(wt_cerebrum_day5, manualAnnotation %in% resident_celltypes)

#Not sure that this actually works, cannot visualize both celltype and gene well
celltype_treatment_necroptosis_scores <- list()
for(i in 1:length(necroptosis)){
  avg_data <- AverageExpression(wt_cerebrum_day5_resident, necroptosis[i], group.by = 'treatment_celltype', assay = 'RNA', slot = 'data')$RNA
  celltype_treatment_necroptosis_scores[[i]] <- (t(avg_data))
}

wt_cerebrum_day5_resident_nScores <- do.call(cbind, celltype_treatment_necroptosis_scores) %>% as.data.frame() %>% 
  tibble::rownames_to_column(var = "id") %>% 
  tidyr::pivot_longer(!id, names_to = 'gene', values_to = 'avg_exp')
split_names <- stringr::str_split_fixed(wt_cerebrum_day5_resident_nScores$id, "-", n = 2)
colnames(split_names) <- c('treatment', 'celltype')
wt_cerebrum_day5_resident_nScores <- cbind(wt_cerebrum_day5_resident_nScores, split_names)

#Differences between treatment and pbs avg expression for each gene
wt_cerebrum_day5_resident_nScores %>%
  dplyr::group_by(gene, celltype) %>%
  dplyr::mutate(exp_change = avg_exp - dplyr::first(avg_exp)) %>% 
  dplyr::filter(treatment == 'rLGTV') %>% 
  ggplot(aes(x = gene, y = celltype, fill = exp_change))+
  geom_tile()+
  scale_fill_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white","blue"),
                       values = c(1, 0.7,0.2,0.05,-0.2),
                       limits = c(-0.45, 10))+
  ggtitle("Necroptosis gene LGTV - PBS difference")+
  geom_text(aes(label=round(exp_change, digits = 2)))

#Line plots showing difference
for(i in 1:length(necroptosis)){
  gene_plot <- wt_cerebrum_day5_resident_nScores %>% dplyr::filter(gene == necroptosis[i]) %>% 
    ggplot(aes(x = treatment, y = avg_exp, color = celltype, group = celltype))+
    geom_line(size = 1)+
    scale_color_manual(values=newCols[c(1,3,4,5,7,9,10,11,13,14)])+
    ggtitle(paste("Resident cell" , necroptosis[i],"change"))+
    theme(legend.text = element_text(size = 14))+
    ylim(0,11)
  print(gene_plot)
}


#check cell counts for each group
table(wt_cerebrum_day5_resident$Treatment, wt_cerebrum_day5_resident$manualAnnotation)%>% 
  as.data.frame() %>% dplyr::group_by(Var1) %>% dplyr::mutate(freq_props = Freq/sum(Freq))%>% 
  ggplot(aes(x = Var1, y = freq_props, fill = Var2))+
  geom_bar(stat = 'identity', position = 'stack', width = 0.6)+
  scale_fill_manual(values = newCols[c(1,3,4,5,7,9,10,11,13,14)])+
  theme_classic()

#Same thing for pyroptosis genes
celltype_treatment_pyroptosis_scores <- list()
for(i in 1:length(pyroptosis)){
  avg_data <- AverageExpression(wt_cerebrum_day5_resident, pyroptosis[i], group.by = 'treatment_celltype', assay = 'RNA', slot = 'data')$RNA
  celltype_treatment_pyroptosis_scores[[i]] <- (t(avg_data))
}

wt_cerebrum_day5_resident_pScores <- do.call(cbind, celltype_treatment_pyroptosis_scores) %>% as.data.frame() %>% 
  tibble::rownames_to_column(var = "id") %>% 
  pivot_longer(!id, names_to = 'gene', values_to = 'avg_exp')
split_names <- str_split_fixed(wt_cerebrum_day5_resident_pScores$id, "-", n = 2)
colnames(split_names) <- c('treatment', 'celltype')
wt_cerebrum_day5_resident_pScores <- cbind(wt_cerebrum_day5_resident_pScores, split_names)

#Differences between treatment and pbs avg expression for each gene
wt_cerebrum_day5_resident_pScores %>%
  dplyr::group_by(gene, celltype) %>%
  dplyr::mutate(exp_change = avg_exp - dplyr::first(avg_exp)) %>% 
  dplyr::filter(treatment == 'rLGTV') %>% 
  ggplot(aes(x = gene, y = celltype, fill = exp_change))+
  geom_tile()+
  scale_fill_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white","blue"),
                       values = c(1, 0.7,0.2,0.09,-0.1),
                       limits = c(-1, 10))+
  ggtitle("Pyroptosis gene LGTV - PBS difference")+
  geom_text(aes(label=round(exp_change, digits = 2)))

#Line plots showing difference
for(i in 1:length(pyroptosis)){
  gene_plot <- wt_cerebrum_day5_resident_pScores %>% dplyr::filter(gene == pyroptosis[i]) %>% 
    ggplot(aes(x = treatment, y = avg_exp, color = celltype, group = celltype))+
    geom_line(size = 1)+
    scale_color_manual(values=newCols[c(1,3,4,5,7,9,10,11,13,14)])+
    ggtitle(paste("Resident cell" , pyroptosis[i],"change"))+
    theme(legend.text = element_text(size = 14))+
    ylim(0,6)
  print(gene_plot)
}


#Use MAST to test for differences
treatment_markers <- FindAllMarkers(wt_cerebrum_day5_resident, group.by = 'Treatment', test.use = 'MAST', 
                                    only.pos = TRUE)
treatment_markers[necroptosis,]
treatment_markers[pyroptosis,]

############Resident pathway enrichment analysis############
############################################################

#Look at differentially expressed pathways from MAST results
#Upregulated in infection
upregulated_infection <- subset(treatment_markers, avg_log2FC > 1 & p_val_adj < 0.01 & cluster == 'rLGTV')

upregulated_infection_paths <- gprofiler2::gost(query = rownames(upregulated_infection), organism = 'mmusculus', evcodes = TRUE)
upregulated_infection_paths$result[upregulated_infection_paths$result$source == 'KEGG',]
upregulated_infection_paths$result[upregulated_infection_paths$result$source == 'GO:MF',]
upregulated_infection_paths$result[upregulated_infection_paths$result$source == 'GO:BP',]


#Are any of these significant DEGs? Create pseuTRUE#Are any of these significant DEGs? Create pseudobulk object and check significance
wt_cerebrum_day5_resident_bulk <- createPseudoBulk(wt_cerebrum_day5_resident, c('Treatment', 'manualAnnotation'))
wt_cerebrum_day5_resident_bulk <- DESeq(wt_cerebrum_day5_resident_bulk)
resultsNames(wt_cerebrum_day5_resident_bulk)
wt_cerebrum_day5_resident_bulk_res <- results(wt_cerebrum_day5_resident_bulk, name = 'Treatment_rLGTV_vs_PBS')
wt_cerebrum_day5_resident_bulk_res[necroptosis,]
wt_cerebrum_day5_resident_bulk_res[pyroptosis,] %>% as.data.frame() %>% dplyr::filter(padj < 0.01)

plotCounts(wt_cerebrum_day5_resident_bulk, gene = 'Tnf', intgroup = 'Treatment')

#Look at differentially expressed pathways
#Upregulated in infection
upregulated_infection <- subset(wt_cerebrum_day5_resident_bulk_res, log2FoldChange > 1 & padj < 0.01)
upregulated_infection_paths <- gprofiler2::gost(query = rownames(upregulated_infection), organism = 'mmusculus', evcodes = TRUE)
upregulated_infection_paths$result[upregulated_infection_paths$result$source == 'KEGG',]
upregulated_infection_paths$result[upregulated_infection_paths$result$source == 'GO:MF',]
upregulated_infection_paths$result[upregulated_infection_paths$result$source == 'GO:BP',]

downregulated_infection <- subset(wt_cerebrum_day5_resident_bulk_res, log2FoldChange < -1 & padj < 0.01)
downregulated_infection_paths <- gprofiler2::gost(query = rownames(downregulated_infection), organism = 'mmusculus', evcodes = TRUE)
downregulated_infection_paths$result[downregulated_infection_paths$result$source == 'KEGG',]
downregulated_infection_paths$result[downregulated_infection_paths$result$source == 'GO:MF',]
downregulated_infection_paths$result[downregulated_infection_paths$result$source == 'GO:BP',]

###################Infiltrating cell analysis###################
################################################################
wt_cerebrum_day5_infiltrating <- subset(wt_cerebrum_day5, manualAnnotation %in% infiltrating_celltypes)

#Compare infiltrating to resident pseudobulk
wt_cerebrum_day5[[]] <- wt_cerebrum_day5[[]] %>% dplyr::mutate(cell_class = case_when(manualAnnotation %in% resident_celltypes ~ 'residential',
                                                                                      manualAnnotation %in% infiltrating_celltypes ~ 'infiltrating',
                                                                                      .default = 'unknown'))
wt_cerebrum_day5_bulk <- createPseudoBulk(wt_cerebrum_day5, c('Treatment','cell_class'))
wt_cerebrum_day5_bulk <- DESeq(wt_cerebrum_day5_bulk)
resultsNames(wt_cerebrum_day5_bulk)
wt_cerebrum_day5_bulk_res <- results(wt_cerebrum_day5_bulk, name = 'cell_class_residential_vs_infiltrating')

#upregulated in residential
upregulated_residential <- subset(wt_cerebrum_day5_bulk_res, log2FoldChange > 1 & padj < 0.01)
upregulated_residential_paths <- gprofiler2::gost(query = rownames(upregulated_residential), organism = 'mmusculus', evcodes = TRUE)
upregulated_residential_paths$result

#upregulated in infiltrating
upregulated_infiltrating <- subset(wt_cerebrum_day5_bulk_res, log2FoldChange < -1 & padj < 0.01)
upregulated_infiltrating_paths <- gprofiler2::gost(query = rownames(upregulated_infiltrating), organism = 'mmusculus', evcodes = TRUE)
upregulated_infiltrating_paths$result

#Check cell counts for each group
t(table(wt_cerebrum_day5_infiltrating$Treatment, wt_cerebrum_day5_infiltrating$manualAnnotation))

#necroptosis genes
celltype_treatment_necroptosis_scores_infil <- list()
for(i in 1:length(necroptosis)){
  avg_data <- AverageExpression(wt_cerebrum_day5_infiltrating, necroptosis[i], group.by = 'treatment_celltype', assay = 'RNA', slot = 'data')$RNA
  celltype_treatment_necroptosis_scores_infil[[i]] <- (t(avg_data))
}

wt_cerebrum_day5_infil_nScores <- do.call(cbind, celltype_treatment_necroptosis_scores_infil) %>% as.data.frame() %>% 
  tibble::rownames_to_column(var = "id") %>% 
  pivot_longer(!id, names_to = 'gene', values_to = 'avg_exp')
split_names <- str_split_fixed(wt_cerebrum_day5_infil_nScores$id, "-", n = 2)
colnames(split_names) <- c('treatment', 'celltype')
wt_cerebrum_day5_infil_nScores <- cbind(wt_cerebrum_day5_infil_nScores, split_names)

#Differences between treatment and pbs avg expression for each gene
#So few pbs cells that this doesn't make sense, just report avg expression
wt_cerebrum_day5_infil_nScores %>%
  dplyr::group_by(gene, celltype) %>%
  dplyr::mutate(exp_change = avg_exp - first(avg_exp)) %>% 
  dplyr::filter(treatment == 'rLGTV') %>% 
  ggplot(aes(x = gene, y = celltype, fill = exp_change))+
  geom_tile()+
  scale_fill_gradientn(colors = c("purple", "white", "orange", "red"),
                       values = scales::rescale(c(-0.2, 0, 1, 2)))+
  ggtitle("Necroptosis gene LGTV - PBS difference")+
  geom_text(aes(label=round(exp_change, digits = 2)))

wt_cerebrum_day5_infil_nScores %>% 
  dplyr::filter(treatment == 'rLGTV') %>% 
  ggplot(aes(x = gene, y = celltype, fill = avg_exp))+
  geom_tile()+
  scale_fill_gradient2(low = "white", mid = "orange", high = "red", midpoint = 2.5)+
  ggtitle("Infiltrating necroptosis LGTV")+
  geom_text(aes(label=round(avg_exp, digits = 2)))

#pyroptosis
#necroptosis genes
celltype_treatment_pyroptosis_scores_infil <- list()
for(i in 1:length(pyroptosis)){
  avg_data <- AverageExpression(wt_cerebrum_day5_infiltrating, pyroptosis[i], group.by = 'treatment_celltype', assay = 'RNA', slot = 'data')$RNA
  celltype_treatment_pyroptosis_scores_infil[[i]] <- (t(avg_data))
}

wt_cerebrum_day5_infil_pScores <- do.call(cbind, celltype_treatment_pyroptosis_scores_infil) %>% as.data.frame() %>% 
  tibble::rownames_to_column(var = "id") %>% 
  pivot_longer(!id, names_to = 'gene', values_to = 'avg_exp')
split_names <- str_split_fixed(wt_cerebrum_day5_infil_pScores$id, "-", n = 2)
colnames(split_names) <- c('treatment', 'celltype')
wt_cerebrum_day5_infil_pScores <- cbind(wt_cerebrum_day5_infil_pScores, split_names)

#Differences between treatment and pbs avg expression for each gene
wt_cerebrum_day5_infil_pScores %>%
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
wt_cerebrum_day5_infil_bulk <- createPseudoBulk(wt_cerebrum_day5_infiltrating, c('Treatment', 'Timepoint'))
wt_cerebrum_day5_infil_bulk <- DESeq(wt_cerebrum_day5_infil_bulk)
wt_cerebrum_day5_infil_bulk_res <- results(wt_cerebrum_day5_infil_bulk, name = 'Treatment_rLGTV_vs_PBS')
wt_cerebrum_day5_infil_bulk_res[necroptosis,]
wt_cerebrum_day5_infil_bulk_res[pyroptosis,] %>% as.data.frame() %>% dplyr::filter(padj < 0.01)

#Don't plot celltypes with less than 50 cells
celltypes_for_plot <- names(table(wt_cerebrum_day5_pbs$manualAnnotation)[table(wt_cerebrum_day5_pbs$manualAnnotation) > 50])
wt_cerebrum_day5_pbs_subset <- subset(wt_cerebrum_day5_pbs, manualAnnotation %in% celltypes_for_plot)
pdf("~/Documents/ÖverbyLab/scPlots/necroptosis_genes_sc/wt_cerebrum_day5_pbs_necroptosis_dot.pdf", width = 9, height = 6)
DotPlot(wt_cerebrum_day5_pbs_subset, features = c(necroptosis, pyroptosis), group.by = 'manualAnnotation', scale = FALSE)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
  ggtitle("WT cerebrum PBS necroptosis and pyroptosis genes")
dev.off()

wt_cerebrum_day5_lgtv <- subset(wt_cerebrum_day5, Treatment == 'rLGTV')
celltypes_for_plot_lgtv <- names(table(wt_cerebrum_day5_lgtv$manualAnnotation)[table(wt_cerebrum_day5_lgtv$manualAnnotation) > 50])
wt_cerebrum_day5_lgtv_subset <- subset(wt_cerebrum_day5_lgtv, manualAnnotation %in% celltypes_for_plot_lgtv)
pdf("~/Documents/ÖverbyLab/scPlots/necroptosis_genes_sc/wt_cerebrum_day5_LGTV_necroptosis_dot.pdf", width = 9, height = 6)
DotPlot(wt_cerebrum_day5_lgtv, features = c(necroptosis, pyroptosis), group.by = 'manualAnnotation', scale = FALSE)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
  ggtitle("WT cerebrum LGTV necroptosis and pyroptosis genes")
dev.off()

mac_mono <- subset(wt_cerebrum_day5, manualAnnotation == 'Macrophage/Monocytes')
DotPlot(mac_mono, features = c(necroptosis, pyroptosis), group.by = 'time_treatment', scale = FALSE)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
  ggtitle("Monocyte macrophage LGTV necroptosis and pyroptosis genes")
table(mac_mono$time_treatment)

#Split up data to look at gene expression across variables

#Add module score to make vln plots 
wt_cerebrum_day5 <- AddModuleScore(wt_cerebrum_day5, features = list(necroptosis), name = 'necroptosis_score')
wt_cerebrum_day5 <- AddModuleScore(wt_cerebrum_day5, features = list(pyroptosis), name = 'pyroptosis_score')
wt_cerebrum_day5$cell_origin <- dplyr::case_when(wt_cerebrum_day5$manualAnnotation %in% infiltrating_celltypes ~ 'infiltrating',
                                                 wt_cerebrum_day5$manualAnnotation %in% resident_celltypes ~ 'resident',
                                                 .default = 'unknown')
wt_cerebrum_day5$cell_origin_treatment = paste(wt_cerebrum_day5$Treatment, wt_cerebrum_day5$cell_origin, sep = '_')
DotPlot(wt_cerebrum_day5, features = 'necroptosis_score1', group.by = 'cell_origin_treatment')
VlnPlot(subset(wt_cerebrum_day5, manualAnnotation != 'unknown'), features = 'necroptosis_score1', group.by = 'cell_origin_treatment', pt.size = 0)+
  theme(legend.position = 'none')

VlnPlot(subset(wt_cerebrum_day5, manualAnnotation != 'unknown'), features = 'pyroptosis_score1', group.by = 'cell_origin_treatment', pt.size = 0)+
  theme(legend.position = 'none')


