#Look at necroptosis/pyroptosis genes from email

#Load packages
library(Seurat)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(tidyr)
library(stringr)
library(gprofiler2)
library(msigdbr)
library(fgsea)
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
wt_cerebrum_day5 <-  subset(ParseSeuratObj_int, Treatment %in% c('PBS', 'rLGTV') & Organ == 'Cerebrum' & 
                              Genotype == 'WT' & (Timepoint == 'Day 5' | Treatment == 'PBS'))

#View top isgs to regress out when plotting
isg_markers <- FindMarkers(wt_cerebrum_day5, group.by = 'Treatment', only.pos = TRUE, test.use = 'MAST',
                           ident.1 = 'rLGTV')
head(isg_markers, n = 30)

isg_markers_to_regress <- c('Ifi209', 'Ifi204', 'Ifit2', 'Rsad2',
'Slfn4', 'Parp14', 'Ifi213', 'Slfn8',
'Ifih1', 'Nlrc5', 'Ifit3', 'Cxcl10', 'Ifi207',
'Ifi211', 'Irf7', 'Isg15')

wt_cerebrum_day5 <- AddModuleScore(wt_cerebrum_day5, features = list(isg_markers_to_regress), name = 'top_isgs_exp')

#No regressing out variables yet
wt_cerebrum_day5 <- prepSeuratObj(wt_cerebrum_day5, regress = FALSE)
ElbowPlot(wt_cerebrum_day5, ndims = 40)
wt_cerebrum_day5 <- prepUmapSeuratObj(wt_cerebrum_day5, nDims = 20, reductionName = 'wt.cerebrum.umap', num_neighbors = 30L)

#regress out ISGs for plotting. Using top degs from previous (not genes that are not found to be isgs elsewhere)
wt_cerebrum_day5_isgs_regressed <- prepSeuratObj(wt_cerebrum_day5, regress = TRUE, regressVars = 'top_isgs_exp1')
ElbowPlot(wt_cerebrum_day5_isgs_regressed, ndims = 40)
wt_cerebrum_day5_isgs_regressed <- prepUmapSeuratObj(wt_cerebrum_day5_isgs_regressed, nDims = 20, reductionName = 'wt.cerebrum.umap', num_neighbors = 20L)

#UMAP without isgs regressed
DimPlot(wt_cerebrum_day5, reduction = 'wt.cerebrum.umap', label = TRUE)

pdf(file = '~/Documents/ÖverbyLab/scPlots/galectin3_proj/gal3_project_umap_byTreatment.pdf',
    width = 7, height = 5)
DimPlot(wt_cerebrum_day5, reduction = 'wt.cerebrum.umap', group.by = 'Treatment')+
  xlab('Umap1')+
  ylab('Umap1')+
  theme(axis.ticks = element_blank(),
        axis.text=element_blank(),
        legend.text=element_text(size=17))+
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank())
dev.off()

pdf(file = '~/Documents/ÖverbyLab/scPlots/galectin3_proj/gal3_project_umap.pdf',
    width = 7, height = 5)
DimPlot(wt_cerebrum_day5, reduction = 'wt.cerebrum.umap', group.by = 'manualAnnotation', cols = newCols)+
  ggtitle("WT Cerebrum Day 5")+
  xlab('Umap1')+
  ylab('Umap1')+
  theme(axis.ticks = element_blank(),
        axis.text=element_blank(),
        legend.text=element_text(size=17))+
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank())
dev.off()


#UMAP with top isgs regressed
DimPlot(wt_cerebrum_day5_isgs_regressed, reduction = 'wt.cerebrum.umap', label = TRUE)

pdf(file = '~/Documents/ÖverbyLab/scPlots/galectin3_proj/gal3_project_umap_byTreatment_isgs_regressed.pdf',
    width = 7, height = 5)
DimPlot(wt_cerebrum_day5_isgs_regressed, reduction = 'wt.cerebrum.umap', group.by = 'Treatment')+
  ggtitle("Treatment - some ISGs regressed")+
  xlab('Umap1')+
  ylab('Umap1')+
  theme(axis.ticks = element_blank(),
        axis.text=element_blank(),
        legend.text=element_text(size=17))+
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank())
dev.off()

pdf(file = '~/Documents/ÖverbyLab/scPlots/galectin3_proj/gal3_project_isgs_regressed_umap.pdf',
    width = 7, height = 5)
DimPlot(wt_cerebrum_day5_isgs_regressed, reduction = 'wt.cerebrum.umap', group.by = 'manualAnnotation', cols = newCols)+
  ggtitle("WT Cerebrum Day 5 - some ISGs regressed out")+
  xlab('Umap1')+
  ylab('Umap1')+
  theme(axis.ticks = element_blank(),
        axis.text=element_blank(),
        legend.text=element_text(size=17))+
  xlab('Umap 1')+
  ylab('Umap 2')+  
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank())
dev.off()

#Want to do resident cells separate from infiltrating, so make lists here
resident_celltypes <- c('Astrocytes', 'Oligodendrocytes', 'Microglia', 'Endothelial', 'Choroid Plexus',
                        'Immature Neurons', 'Ependymal', 'Pericytes', 'Muscle cells', 'Neurons')

infiltrating_celltypes <- c("T cells", 'Macrophage/Monocytes', 'Nk cells', 'Granulocytes', 'B Cells')

###################Resident cell analysis###################
############################################################
#Split by treatment
wt_cerebrum_day5_resident <- subset(wt_cerebrum_day5, manualAnnotation %in% resident_celltypes)
wt_cerebrum_day5_infil <- subset(wt_cerebrum_day5, manualAnnotation %in% infiltrating_celltypes)

#Necroptosis gene expressions across celltype and gene
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

#Only plot celltypes with enough cells in infected
relevant_celltypes <- table(wt_cerebrum_day5_resident$Treatment, wt_cerebrum_day5_resident$manualAnnotation) %>% 
  as.data.frame() %>% 
  dplyr::filter(Var1 == 'rLGTV' & Freq > 30) 

#Get list of celltypes with enough cells
celltypes_over_30 <- relevant_celltypes$Var2

#Differences between treatment and pbs avg expression for each gene
pdf('~/Documents/ÖverbyLab/scPlots/galectin3_proj/cell_death_pathways/resident_necroptosis_scores.pdf', height = 5, width = 8)
wt_cerebrum_day5_resident_nScores %>%
  dplyr::group_by(gene, celltype) %>%
  dplyr::mutate(exp_change = avg_exp - dplyr::first(avg_exp)) %>% 
  dplyr::filter(treatment == 'rLGTV' & celltype %in% celltypes_over_30) %>% 
  ggplot(aes(x = gene, y = celltype, fill = exp_change))+
  geom_tile()+
  scale_fill_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white","blue"),
                       values = c(1, 0.7,0.2,0.05,-0.2),
                       limits = c(-0.45, 10))+
  ggtitle("Necroptosis gene LGTV - PBS difference")+
  geom_text(aes(label=round(exp_change, digits = 2)), size = 6) +
  theme(axis.text = element_text(size = 16),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 16),
        panel.grid= element_blank(),
        panel.background = element_rect(fill = 'white', colour = 'white'),
        axis.ticks = element_blank())+
  ylab('')+
  xlab('')
dev.off()

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
pdf('~/Documents/ÖverbyLab/scPlots/galectin3_proj/cell_death_pathways/wt_resident_celltype_props.pdf', height = 5, width = 8)
table(wt_cerebrum_day5_resident$Treatment, wt_cerebrum_day5_resident$manualAnnotation)%>% 
  as.data.frame() %>% dplyr::group_by(Var1) %>% dplyr::mutate(freq_props = Freq/sum(Freq))%>% 
  ggplot(aes(x = Var1, y = freq_props, fill = Var2))+
  geom_bar(stat = 'identity', position = 'stack', width = 0.6)+
  scale_fill_manual(values = newCols[c(1,3,4,5,7,9,10,11,13,14)])+
  theme_classic()+
  ylab('Cell type proportions')+
  xlab('') +
  theme(axis.text = element_text(size= 16),
        axis.title = element_text(size= 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16))+
  guides(fill=guide_legend(title="cell type"))
dev.off()

#total counts, rather than proportions
pdf('~/Documents/ÖverbyLab/scPlots/galectin3_proj/cell_death_pathways/wt_resident_celltype_counts.pdf', height = 5, width = 8)
table(wt_cerebrum_day5_resident$Treatment, wt_cerebrum_day5_resident$manualAnnotation)%>% 
  as.data.frame() %>% 
  ggplot(aes(x = Var1, y = Freq, fill = Var2))+
  geom_bar(stat = 'identity', position = 'stack', width = 0.6)+
  scale_fill_manual(values = newCols[c(1,3,4,5,7,9,10,11,13,14)])+
  theme_classic()+
  ylab('Cell type counts')+
  xlab('') +
  theme(axis.text = element_text(size= 16),
        axis.title = element_text(size= 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16))+
  guides(fill=guide_legend(title="cell type"))
dev.off()

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
pdf('~/Documents/ÖverbyLab/scPlots/galectin3_proj/cell_death_pathways/resident_pyroptosis_scores.pdf', height = 5, width = 8)
wt_cerebrum_day5_resident_pScores %>%
  dplyr::group_by(gene, celltype) %>%
  dplyr::mutate(exp_change = avg_exp - dplyr::first(avg_exp)) %>% 
  dplyr::filter(treatment == 'rLGTV' & celltype %in% celltypes_over_30) %>% 
  ggplot(aes(x = gene, y = celltype, fill = exp_change))+
  geom_tile()+
  scale_fill_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white","blue"),
                       values = c(1, 0.7,0.2,0.09,-0.1),
                       limits = c(-1, 10))+
  ggtitle("Pyroptosis gene LGTV - PBS difference")+
  geom_text(aes(label=round(exp_change, digits = 2)), size = 3) +
  theme(axis.text = element_text(size = 16),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 16),
        axis.text.x = element_text(angle = 45, vjust = 0.5),
        panel.grid= element_blank(),
        panel.background = element_rect(fill = 'white', colour = 'white'),
        axis.ticks = element_blank())+
  ylab('')+
  xlab('')
dev.off()

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

#Infiltrating cell proportions
pdf('~/Documents/ÖverbyLab/scPlots/galectin3_proj/cell_death_pathways/wt_infiltrating_celltype_props.pdf', height = 5, width = 8)
table(wt_cerebrum_day5_infil$Treatment, wt_cerebrum_day5_infil$manualAnnotation)%>% 
  as.data.frame() %>% dplyr::group_by(Var1) %>% dplyr::mutate(freq_props = Freq/sum(Freq))%>% 
  ggplot(aes(x = Var1, y = freq_props, fill = Var2))+
  geom_bar(stat = 'identity', position = 'stack', width = 0.6)+
  scale_fill_manual(values = newCols[c(2, 6, 8, 12, 15)])+
  theme_classic()
dev.off()

#Total counts
pdf('~/Documents/ÖverbyLab/scPlots/galectin3_proj/cell_death_pathways/wt_infiltrating_celltype_counts.pdf', height = 5, width = 8)
table(wt_cerebrum_day5_infil$Treatment, wt_cerebrum_day5_infil$manualAnnotation)%>% 
  as.data.frame() %>% 
  ggplot(aes(x = Var1, y = Freq, fill = Var2))+
  geom_bar(stat = 'identity', position = 'stack', width = 0.6)+
  scale_fill_manual(values = newCols[c(2, 6, 8, 12, 15)])+
  theme_classic()+
  ylab('Cell type counts')+
  xlab('') +
  theme(axis.text = element_text(size= 16),
        axis.title = element_text(size= 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16))+
  guides(fill=guide_legend(title="cell type"))
dev.off()

#Necroptosis and pyroptosis scores by celltype and treatment (addModuleScore then just dotplot)
wt_cerebrum_day5_resident <- AddModuleScore(wt_cerebrum_day5_resident, features = list(necroptosis), name = 'necroptosis')
wt_cerebrum_day5_resident <- AddModuleScore(wt_cerebrum_day5_resident, features = list(pyroptosis), name = 'pyroptosis')

necroptosis_heat_dat <- wt_cerebrum_day5_resident[[]] %>% dplyr::group_by(treatment_celltype) %>% 
  dplyr::summarise(avg_exp = mean(necroptosis1))

pdf('~/Documents/ÖverbyLab/scPlots/galectin3_proj/cell_death_pathways/celltype_necroptosis_scores.pdf', height = 5, width = 8)
necroptosis_heat_dat %>% tidyr::extract(col = treatment_celltype, regex = '(.+)_(.+)', into = c('treatment', 'celltype')) %>% 
  dplyr::filter(celltype %in% celltypes_over_30) %>% 
  ggplot(aes(x = treatment, y = celltype, fill = avg_exp))+
  geom_tile()+
  ggtitle('Necroptosis score')+
  scale_fill_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white","blue"),
                       values = c(1, 0.7,0.2,0.05,-0.2),
                       limits = c(-0.15, 0.5))+
  geom_text(aes(label=round(avg_exp, digits = 2)), size = 8)+
  theme(axis.text = element_text(size = 16),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 16),
        plot.title = element_text(size = 18),
        panel.grid= element_blank(),
        panel.background = element_rect(fill = 'white', colour = 'white'),
        axis.ticks = element_blank())+
  ylab('')+
  xlab('')
dev.off()

#Only use heat_dat, dot_dat actually uses mean of expm1 values as default rather than just mean. Leaving here as reminder
pyroptosis_dot_dat <- DotPlot(wt_cerebrum_day5_resident, features = 'pyroptosis1', group.by = 'treatment_celltype')$data
pyroptosis_heat_dat <- wt_cerebrum_day5_resident[[]] %>% dplyr::group_by(treatment_celltype) %>% 
  dplyr::summarise(avg_exp = mean(pyroptosis1))


pdf('~/Documents/ÖverbyLab/scPlots/galectin3_proj/cell_death_pathways/celltype_pyroptosis_scores.pdf', height = 5, width = 8)
pyroptosis_heat_dat %>% tidyr::extract(col = treatment_celltype, regex = '(.+)_(.+)', into = c('treatment', 'celltype')) %>% 
  dplyr::filter(celltype %in% celltypes_over_30) %>% 
  ggplot(aes(x = treatment, y = celltype, fill = avg_exp))+
  geom_tile()+
  ggtitle('Pyroptosis score')+
  scale_fill_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white","blue"),
                       values = c(1, 0.7,0.2,0.05,-0.2),
                       limits = c(-0.15, 0.5))+
  geom_text(aes(label=round(avg_exp, digits = 2)), size = 8)+
  theme(axis.text = element_text(size = 16),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 16),
        plot.title = element_text(size = 18),
        panel.grid= element_blank(),
        panel.background = element_rect(fill = 'white', colour = 'white'),
        axis.ticks = element_blank())+
  ylab('')+
  xlab('')
dev.off()

#Use MAST to test for differences
treatment_markers <- FindMarkers(wt_cerebrum_day5_resident, group.by = 'Treatment', test.use = 'MAST', ident.1 = 'rLGTV',
                                    only.pos = TRUE)
treatment_markers[necroptosis,]
treatment_markers[pyroptosis,]

############Resident pathway enrichment analysis############
############################################################

#Look at differentially expressed pathways from MAST results
#Upregulated in infection

upregulated_infection <- subset(treatment_markers, avg_log2FC > 1 & p_val_adj < 0.01)

#Set high p value threshold to see all pathways
upregulated_infection_paths <- gprofiler2::gost(query = rownames(upregulated_infection), organism = 'mmusculus', evcodes = TRUE,
                                                user_threshold = 1)

upregulated_infection_paths$result[upregulated_infection_paths$result$source == 'KEGG',]
upregulated_infection_paths$result[upregulated_infection_paths$result$source == 'GO:MF',]
upregulated_infection_paths$result[upregulated_infection_paths$result$source == 'GO:BP',]

#beginning of possible cell death pathway names, doing partial names to match pathways
paths_to_match <- c("apopt", "pyrop", "necrop", "ferrop", 'autopha', 'cuprop')
upregulated_infection_paths$result[grep(paste(paths_to_match,collapse="|"), 
     upregulated_infection_paths$result$term_name, ignore.case = TRUE),]$term_name %>% sort()

sig_pathways <- upregulated_infection_paths$result[upregulated_infection_paths$result$p_value < 0.05,]
sig_pathways[grep(paste(paths_to_match,collapse="|"), 
                  sig_pathways$term_name, ignore.case = TRUE),]$term_name %>% sort()

#Barplot of cell death pathway p values
pathways_to_plot <- upregulated_infection_paths$result[upregulated_infection_paths$result$term_name %in%
                                                         c('Apoptosis', 'autophagy', 'Necroptosis', 
                                                           'pyroptotic inflammatory response', 'Pyroptosis',
                                                           'ferroptosis', 'Ferroptosis', 'necroptotic process',
                                                           'positive regulation of apoptotic process', 'cell death'),]
pathways_to_plot <- pathways_to_plot[pathways_to_plot$source %in% c('GO:BP', 'KEGG'),]

pdf('~/Documents/ÖverbyLab/scPlots/galectin3_proj/cell_death_pathways/gprofiler_path_barplot.pdf', height = 5, width = 8)
ggplot(pathways_to_plot, aes(x = term_name, y = -log10(p_value), fill = source, color = source))+
  geom_bar(stat = 'identity', position = 'dodge')+
  geom_hline(yintercept = -log10(0.01), linetype= 2)+
  coord_flip()
dev.off()

#Look at relevant pathways by celltype
pathway_dotplot_dat <- pathways_to_plot[,c('term_name', 'intersection')]
pathway_dotplot_dat_usable <- list()
for(i in 1:nrow(pathway_dotplot_dat)){
  current_path <- pathway_dotplot_dat[i,]
  gene_list = stringr::str_split(current_path$intersection, pattern = ',') 
  pathway_dotplot_dat_usable[i] = gene_list
  names(pathway_dotplot_dat_usable)[i] = current_path$term_name
}

#Add module scores for plotting to seurat object
for(i in 1:length(pathway_dotplot_dat_usable)){
  wt_cerebrum_day5_resident <- AddModuleScore(wt_cerebrum_day5_resident, 
                                                              features = list(pathway_dotplot_dat_usable[[i]]),
                                              name = gsub(" ", "_", names(pathway_dotplot_dat_usable)[i]))
}

celltypes_with_150 <- table(wt_cerebrum_day5_resident$manualAnnotation) %>% as.data.frame() %>% 
  dplyr::filter(Freq > 150)

pdf('~/Documents/ÖverbyLab/scPlots/galectin3_proj/cell_death_pathways/gprofiler_paths_byCelltype.pdf', height = 5, width = 8)
wt_cerebrum_day5_resident[[]] %>% dplyr::filter(manualAnnotation %in% celltypes_with_150$Var1) %>% 
  dplyr::group_by(manualAnnotation, Treatment) %>% 
  dplyr:: summarise(across(cell_death1:Necroptosis1, ~ mean(.x, na.rm = TRUE))) %>% 
  tidyr::pivot_longer(cols = c(!c(manualAnnotation, Treatment)), names_to = 'pathway', values_to = 'meanExp') %>% 
  dplyr::group_by(manualAnnotation, pathway) %>% 
  dplyr::mutate(path_diff = meanExp - meanExp[Treatment == 'PBS']) %>% 
  dplyr::filter(Treatment == 'rLGTV' & !pathway %in% c('ferroptosis1', 'Ferroptosis1', 'necroptotic_process1')) %>% 
  dplyr::mutate(pathway = case_when(pathway == 'cell_death1'~ 'cell death',
                                    pathway == 'Apoptosis1'~ 'apoptosis',
                                    pathway == 'positive_regulation_of_apoptotic_process1'~ 'positive regulation apoptosis',
                                    pathway == 'pyroptotic_inflammatory_response1'~ 'pyroptotic inflammation',
                                    pathway == 'autophagy1'~ 'autophagy',
                                    pathway == 'Necroptosis1'~ 'necroptosis',
                                    pathway == 'necroptotic_process1'~ 'necroptotic process',
                                    .default = pathway)) %>% 
  ggplot(aes(x = pathway, y = manualAnnotation, fill = path_diff))+
  geom_tile()+
  scale_fill_gradient2(low = '#3DB9FF', mid = 'white', high = 'red', midpoint = 0)+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.6, size = 12),
        axis.text.y = element_text(size = 16))+
  xlab("")+
  ylab("")+
  ggtitle('LGTV - PBS average score')
dev.off()


#Are any of these significant DEGs? Create pseudobulk object and check significance
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

#Also look at pathways with fgsea, maybe best and can include custom lists.
#But think I won't go with gsea as we have plenty of highly significant degs
#for enrichment analyses

#Molecular signatures database
all_gene_sets <- msigdbr(species = "Mus musculus", db_species = "MM")

mouse_gene_sets <- all_gene_sets %>%
  split(x = .$gene_symbol, f = .$gs_name)

#Gene lists from email
anti_apoptosis <- c("Bcl2", "Bcl2l1", 'Bcl2l2', 'Bag1', 'Bag2', 'Bag3', 'Bag4')
pro_apoptosis <- c('Apaf1', 'Casp9', 'Casp8', 'Bax', 'Bak',
                   'Bid', 'Bad', 'Bim', 'Bcl10', 'Bik',
                   'Blk', 'Fas', 'Fasl', 'Tnfrsf1a', 'Tnf', 'Tyro3',
                   'Axl', 'Mertk', 'Tnfsf10', 'Tnfrsf10b', 'Casp3', 'Casp6', 'Casp7')

#Also lists Ferritin, Fth1? Ftl1? Think the map genes are right, should double check though
ferroptosis <- c('Gpx4', 'Acsl4', 'Ptgs2', 'Slc39a14', 'Prnp', 'Steap3', 'Vdac2', 'Vdac3', 'Alox15', 'Atf3')
autophagy <- c('Atg3', 'Atg5', 'Atg7', 'Atg10', 'Atg12', 'Atg13', 'Atg14', 'Ulk1', 'Becn1',
               'Ambra1', 'Map1lc3a', 'Map1lc3a')
cuproptosis <- c('Fdx1', 'Lias', 'Lipt1', 'Dld', 'Dlat', 'Pdha1', 'Pdhb', 'Mtf1',
                 'Gls', 'Cdkn2a', 'Atp7b', 'Slc31a1', 'Atp7a', 'Dlst', 'Dbt', 'Gcsh')

mouse_gene_sets$custom_necroptosis <- necroptosis
mouse_gene_sets$custom_pyrooptosis <- pyroptosis
mouse_gene_sets$custom_pro_apoptosis <- pro_apoptosis
mouse_gene_sets$custom_pro_ferroptosis <- ferroptosis
mouse_gene_sets$custom_autophagy <- autophagy
mouse_gene_sets$custom_cuproptosis <- cuproptosis

treatment_markers_sorted <- treatment_markers %>% dplyr::filter(cluster == 'rLGTV') %>% dplyr::arrange(desc(avg_log2FC)) 
treatment_markers_sorted_logfcs <- treatment_markers_sorted$avg_log2FC
names(treatment_markers_sorted_logfcs) <- rownames(treatment_markers_sorted)

GSEAres <- fgsea(pathways = mouse_gene_sets, # List of gene sets to check
                 stats = treatment_markers_sorted_logfcs,
                 scoreType = 'pos', # in this case we have both pos and neg rankings. if only pos or neg, set to 'pos', 'neg'
                 minSize = 5,
                 maxSize = 500) 

GSEAres[grep('custom', GSEAres$pathway, ignore.case = TRUE),] %>% 
  dplyr::arrange(padj)




###################Infiltrating cell analysis###################
################################################################

#Will subtract average pbs resident score from infiltrating scores
wt_cerebrum_day5_infiltrating <- subset(wt_cerebrum_day5, manualAnnotation %in% infiltrating_celltypes | Treatment == 'PBS')

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

#subtract avg pbs resident from infiltrating
pbs_resident_avg <- dplyr::filter(wt_cerebrum_day5_infil_nScores, treatment == 'PBS' & celltype %in% resident_celltypes)
pbs_resident_avg <- pbs_resident_avg %>% dplyr::group_by(gene) %>% dplyr::summarise(mean_pbs_exp = mean(avg_exp))

#Check which cells have high enough expression to keep in
cells_with_enough <- table(wt_cerebrum_day5_infiltrating$Treatment, wt_cerebrum_day5_infiltrating$manualAnnotation) %>% 
  as.data.frame() %>% 
  dplyr::filter(Var1 == 'rLGTV' & Freq > 30) 
celltypes_with_enough <- cells_with_enough$Var2

#Differences between treatment and pbs avg Var2#Differences between treatment and pbs avg expression for each gene
#So few pbs cells that this doesn't make sense, just report avg expression
pdf('~/Documents/ÖverbyLab/scPlots/galectin3_proj/cell_death_pathways/infiltrating_necroptosis_scores.pdf', height = 5, width = 8)
wt_cerebrum_day5_infil_nScores %>%
  dplyr::group_by(gene, celltype) %>%
  dplyr::left_join(pbs_resident_avg, by = 'gene') %>% 
  dplyr::mutate(exp_change = avg_exp - mean_pbs_exp) %>% 
  dplyr::filter(treatment == 'rLGTV' & celltype %in% celltypes_with_enough) %>% 
  ggplot(aes(x = gene, y = celltype, fill = exp_change))+
  geom_tile()+
  scale_fill_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white"),
                                   values = c(1, 0.7,0.2,0),
                                   limits = c(-0.01, 5))+
  ggtitle("Necroptosis gene LGTV - resident PBS average")+
  geom_text(aes(label=round(exp_change, digits = 2)), size = 6) +
  theme(axis.text = element_text(size = 16),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 16),
        panel.grid= element_blank(),
        panel.background = element_rect(fill = 'white', colour = 'white'),
        axis.ticks = element_blank())+
  ylab('')+
  xlab('')
dev.off()

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

#subtract avg pbs resident from infiltrating
pbs_resident_avg_pyrop <- dplyr::filter(wt_cerebrum_day5_infil_pScores, treatment == 'PBS' & celltype %in% resident_celltypes)
pbs_resident_avg_pyrop <- pbs_resident_avg_pyrop %>% dplyr::group_by(gene) %>% dplyr::summarise(mean_pbs_exp = mean(avg_exp))


#Differences between treatment and pbs avg expression for each gene
pdf('~/Documents/ÖverbyLab/scPlots/galectin3_proj/cell_death_pathways/infiltrating_pyroptosis_scores.pdf', height = 5, width = 10)
wt_cerebrum_day5_infil_pScores %>%
  dplyr::group_by(gene, celltype) %>%
  dplyr::left_join(pbs_resident_avg_pyrop, by = 'gene') %>% 
  dplyr::mutate(exp_change = avg_exp - mean_pbs_exp) %>% 
  dplyr::filter(treatment == 'rLGTV' & celltype %in% celltypes_with_enough) %>% 
  ggplot(aes(x = gene, y = celltype, fill = log(exp_change + 1)))+
  geom_tile()+
  scale_fill_gradientn(colours = rev(c("#F03C0C","#F57456","#FFB975", 'white', "blue")),
                      # values = c(1, 0.7,0.1, 0.001,  -0.01),
                      rescaler = ~ scales::rescale_mid(.x, mid = 1.4),
                       limits = c(-0.9, 5))+
  ggtitle("Pyroptosis gene LGTV - resident PBS average")+
  geom_text(aes(label=round(log(exp_change+1), digits = 2)), size = 4) +
  theme(axis.text = element_text(size = 16),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 16),
        panel.grid= element_blank(),
        panel.background = element_rect(fill = 'white', colour = 'white'),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 0.8))+
  ylab('')+
  xlab('')
dev.off()

wt_cerebrum_day5_infil_pScores %>% 
  dplyr::filter(treatment == 'rLGTV') %>% 
  ggplot(aes(x = gene, y = celltype, fill = log(avg_exp+1)))+
  geom_tile()+
  scale_fill_gradient2(low = "white", mid = "orange", high = "red", midpoint = 2)+
  ggtitle("Infiltrating pyroptosis LGTV (natural log scale)")+
  geom_text(aes(label=round(log(avg_exp+1), digits = 2)))

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


