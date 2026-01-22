#Load libraries 
library(dplyr)
library(Seurat)
library(ggplot2)
library(SingleR)
library(celldex)
library(scran)
library(RColorBrewer)

### SUBSET TO JUST WILDTYPE ###
#Source function
source('~/Documents/ÖverbyLab/scripts/langatFunctions.R')

#This is where the 10x data is
setwd('~/Documents/ÖverbyLab/single_nuclei_proj/')

#Read in processed data
sn_integrated_dat <- LoadSeuratRds('~/Documents/ÖverbyLab/single_nuclei_proj/LGTVscCombined.rds')

#Custom pathways
#Gene lists
necroptosis <- c('Tnf', 'Tnfrsf1a', 'Ripk2', 'Mlkl', 'Ripk1', 'Ripk3')
pyroptosis <- c('Gsdmc', 'Nlrp3', 'Aim2', 'Gsdmd', 'Il18', 'Il1b', 'Casp9', 'Casp8', 'Casp6', 'Casp3', 'Casp4', 'Casp1')

#Add column for easier plotting later
#Infected column matches up with new_inf column so good to use 
sn_integrated_dat$treatment_celltype <- paste(sn_integrated_dat$infected, sn_integrated_dat$manualAnnotation)
sn_integrated_dat$genotype_celltype <- paste(sn_integrated_dat$new_genotype, sn_integrated_dat$manualAnnotation)

#Check data
newCols <-  c(brewer.pal(12, 'Paired'), '#99FFE6', '#CE99FF', '#18662E','#737272',  '#FF8AEF')
DimPlot(sn_integrated_dat, group.by =   'manualAnnotation', cols = newCols, reduction = 'umap.integrated')

#Look at infection and see if labels make sense. Why not more infection in ifnar knockout?
FeaturePlot(sn_integrated_dat, features = 'LGTV', reduction = 'umap.integrated', split.by = 'new_inf')
FeaturePlot(sn_integrated_dat, features = 'rna_LGTV', reduction = 'umap.integrated',  split.by = 'new_genotype')
table(sn_integrated_dat$new_genotype, sn_integrated_dat$new_inf)
DotPlot(sn_integrated_dat, features = 'rna_LGTV', group.by = 'treatment_celltype')
DotPlot(sn_integrated_dat, features = 'rna_LGTV', group.by = 'genotype_celltype')
DotPlot(sn_integrated_dat, features = 'Rsad2', group.by = 'genotype_celltype')
DotPlot(sn_integrated_dat, features = 'Ifng', group.by = 'genotype_celltype')

#subset to wt
sn_integrated_dat_wt <- subset(sn_integrated_dat, new_genotype %in% c('wt', 'wt (same)'))

sn_integrated_dat_wt <- prepSeuratObj(sn_integrated_dat_wt)
ElbowPlot(sn_integrated_dat_wt) 
sn_integrated_dat_wt <- prepUmapSeuratObj(sn_integrated_dat_wt, nDims = 20,num_neighbors = 30L, 
                                          reductionName = 'wt.umap.integrated',
                                          resolution_value = 0.8)

pdf("~/Documents/ÖverbyLab/single_nuclei_proj/sn_plots/wt_celltype_umap.pdf", width = 8, height = 6)
DimPlot(sn_integrated_dat_wt, group.by =  'manualAnnotation', cols = newCols, reduction = 'wt.umap.integrated')+
  xlab('Umap1')+
  ylab('Umap1')+
  theme(axis.ticks = element_blank(),
        axis.text=element_blank(),
        legend.text=element_text(size=17))+
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank())+
  ggtitle('single nuclei WT')
dev.off()

pdf("~/Documents/ÖverbyLab/single_nuclei_proj/sn_plots/wt_celltype_umapbyInfection.pdf", width = 8, height = 6)
DimPlot(sn_integrated_dat_wt, group.by =  'infected', cols = newCols, reduction = 'wt.umap.integrated')+
  xlab('Umap1')+
  ylab('Umap1')+
  theme(axis.ticks = element_blank(),
        axis.text=element_blank(),
        legend.text=element_text(size=17))+
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank())+
  ggtitle('single nuclei WT by infection')
dev.off()

FeaturePlot(sn_integrated_dat_wt, features = 'rna_LGTV', reduction = 'wt.umap.integrated')

#regress out isgs for umap
lgtv_markers <- FindMarkers(sn_integrated_dat_wt, group.by = 'infected', ident.1 = 'TRUE', test.use = 'MAST',
                            only.pos = TRUE)
top_isgs <- head(rownames(lgtv_markers), n = 30)
sn_integrated_dat_wt <- AddModuleScore(sn_integrated_dat_wt, features = list(top_isgs), name = 'top_isg_score')

#reprocess data with isg regression
#I think a good option oculd be classify isgs by type, add multiple mole scores to seurat object, then regress on those scores
sn_integrated_dat_wt_isgs_regressed <- prepSeuratObj(sn_integrated_dat_wt, regress = TRUE, regressVars = top_isgs[!grepl('H2', top_isgs)])
ElbowPlot(sn_integrated_dat_wt_isgs_regressed, ndims = 40) 
sn_integrated_dat_wt_isgs_regressed <- prepUmapSeuratObj(sn_integrated_dat_wt_isgs_regressed, nDims = 25,num_neighbors = 30L, 
                                          reductionName = 'wt.umap.integrated.noisgs',
                                          resolution_value = 0.8)

DimPlot(sn_integrated_dat_wt_isgs_regressed, group.by =  'manualAnnotation', cols = newCols, reduction = 'wt.umap.integrated.noisgs')
DimPlot(sn_integrated_dat_wt_isgs_regressed, group.by =  'infected', cols = newCols, reduction = 'wt.umap.integrated.noisgs')

############ Begin Pathway Analyses ############ 
################################################ 

#necroptosis and pyroptosis plots

#Gene lists
necroptosis <- c('Tnf', 'Tnfrsf1a', 'Ripk2', 'Mlkl', 'Ripk1', 'Ripk3')
pyroptosis <- c('Gsdmc', 'Nlrp3', 'Aim2', 'Gsdmd', 'Il18', 'Il1b', 'Casp9', 'Casp8', 'Casp6', 'Casp3', 'Casp4', 'Casp1')

#Necroptosis gene expressions across celltype and gene
celltype_treatment_necroptosis_scores <- list()
for(i in 1:length(necroptosis)){
  avg_data <- AverageExpression(sn_integrated_dat_wt, necroptosis[i], group.by = 'treatment_celltype', assay = 'RNA', slot = 'data')$RNA
  celltype_treatment_necroptosis_scores[[i]] <- (t(avg_data))
}

sn_integrated_dat_wt_nScores <- do.call(cbind, celltype_treatment_necroptosis_scores) %>% as.data.frame() %>% 
  tibble::rownames_to_column(var = "id") %>% 
  tidyr::pivot_longer(!id, names_to = 'gene', values_to = 'avg_exp')
split_names <- stringr::str_split_fixed(sn_integrated_dat_wt_nScores$id, " ", n = 2)
colnames(split_names) <- c('treatment', 'celltype')
sn_integrated_dat_wt_nScores <- cbind(sn_integrated_dat_wt_nScores, split_names)

#Differences between treatment and pbs avg expression for each gene
pdf('~/Documents/ÖverbyLab/scPlots/galectin3_proj/cell_death_pathways/single_nuclei_necroptosis_scores.pdf', height = 5, width = 8)
sn_integrated_dat_wt_nScores %>%
  dplyr::group_by(gene, celltype) %>%
  dplyr::mutate(exp_change = avg_exp - dplyr::first(avg_exp)) %>% 
  dplyr::filter(treatment == TRUE) %>% 
  ggplot(aes(x = gene, y = celltype, fill = exp_change))+
  geom_tile()+
  scale_fill_gradientn(colours = rev(c("#F03C0C","#F57456","#FFB975", 'white', "lightblue")),
                       # values = c(1, 0.7,0.1, 0.001,  -0.01),
                       rescaler = ~ scales::rescale_mid(.x, mid = 1),
                       limits = c(-0.55, 3))+
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

#Pyroptosis gene expressions across celltype and gene
celltype_treatment_pyroptosis_scores <- list()
for(i in 1:length(pyroptosis)){
  avg_data <- AverageExpression(sn_integrated_dat_wt, pyroptosis[i], group.by = 'treatment_celltype', assay = 'RNA', slot = 'data')$RNA
  if(!is.null(avg_data)){
    celltype_treatment_pyroptosis_scores[[i]] <- (t(avg_data))
  }
}

sn_integrated_dat_wt_pScores <- do.call(cbind, celltype_treatment_pyroptosis_scores) %>% as.data.frame() %>% 
  tibble::rownames_to_column(var = "id") %>% 
  tidyr::pivot_longer(!id, names_to = 'gene', values_to = 'avg_exp')
split_names <- stringr::str_split_fixed(sn_integrated_dat_wt_pScores$id, " ", n = 2)
colnames(split_names) <- c('treatment', 'celltype')
sn_integrated_dat_wt_pScores <- cbind(sn_integrated_dat_wt_pScores, split_names)

#Differences between treatment and pbs avg expression for each gene
pdf('~/Documents/ÖverbyLab/scPlots/galectin3_proj/cell_death_pathways/single_nuclei_pyroptosis_scores.pdf', height = 5, width = 8)
sn_integrated_dat_wt_pScores %>%
  dplyr::group_by(gene, celltype) %>%
  dplyr::mutate(exp_change = avg_exp - dplyr::first(avg_exp)) %>% 
  dplyr::filter(treatment == TRUE) %>% 
  ggplot(aes(x = gene, y = celltype, fill = exp_change))+
  geom_tile()+
  scale_fill_gradientn(colours = rev(c("#F03C0C","#F57456","#FFB975", 'white', "lightblue")),
                       # values = c(1, 0.7,0.1, 0.001,  -0.01),
                       rescaler = ~ scales::rescale_mid(.x, mid = 1),
                       limits = c(-1, 3))+
  ggtitle("Pyroptosis gene LGTV - PBS difference")+
  geom_text(aes(label=round(exp_change, digits = 2)), size = 4) +
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


############Infected pathway enrichment analysis############
############################################################

#Use MAST to test for differences
treatment_markers <- FindMarkers(sn_integrated_dat_wt, group.by = 'infected', test.use = 'MAST', 
                                 only.pos = TRUE, ident.1 = TRUE)
treatment_markers[necroptosis,]
treatment_markers[pyroptosis,]

#Look at differentially expressed pathways from MAST results
#Upregulated in infection

upregulated_infection <- subset(treatment_markers, avg_log2FC > 1 & p_val_adj < 0.01)

#Set high p value threshold to see all pathways
upregulated_infection_paths <- gprofiler2::gost(query = rownames(upregulated_infection), organism = 'mmusculus', evcodes = TRUE,
                                                user_threshold = 1)

upregulated_infection_paths$result[upregulated_infection_paths$result$source == 'KEGG',]
upregulated_infection_paths$result[upregulated_infection_paths$result$source == 'GO:BP',]


#Barplot of cell death pathway p values
pathways_to_plot <- upregulated_infection_paths$result[upregulated_infection_paths$result$term_name %in%
                                                         c('Apoptosis', 'autophagy', 'Necroptosis', 
                                                           'pyroptotic inflammatory response', 'Pyroptosis',
                                                           'ferroptosis', 'Ferroptosis', 'necroptotic process',
                                                           'positive regulation of apoptotic process', 'cell death'),]
pathways_to_plot <- pathways_to_plot[pathways_to_plot$source %in% c('GO:BP', 'KEGG'),]

pdf('~/Documents/ÖverbyLab/scPlots/galectin3_proj/cell_death_pathways/gprofiler_path_barplot_single_nuclei.pdf', height = 5, width = 8)
ggplot(pathways_to_plot, aes(x = term_name, y = -log10(p_value), fill = source, color = source))+
  geom_bar(stat = 'identity', position = 'dodge')+
  geom_hline(yintercept = -log10(0.01), linetype= 2)+
  coord_flip()
dev.off()

#Chemokine expression
#Copying previously written function faceted_geno_heatmap in ISG_comps_heatmap.R script

ccl_chemokines <- rownames(sn_integrated_dat@assays$RNA$data)[grep('Ccl', rownames(sn_integrated_dat@assays$RNA$data))]
cxcl_chemokines <- rownames(sn_integrated_dat@assays$RNA$data)[grep('Cxcl', rownames(sn_integrated_dat@assays$RNA$data))]

cyto_data <- FetchData(object = sn_integrated_dat, vars = c(ccl_chemokines, cxcl_chemokines), layer = "data") %>% 
  tibble::rownames_to_column(var = 'cell_id') %>% 
  tidyr::pivot_longer(!cell_id, names_to = 'gene', values_to = 'expression')

#Create metadata dataframe
cell_metadata <- sn_integrated_dat[[]] %>% rownames_to_column(var = 'cell_id')

gene_plot_data <- dplyr::left_join(dplyr::select(cell_metadata, c(cell_id, manualAnnotation, new_inf)), 
                                   cyto_data, by = 'cell_id')

#Doing mean of exponential (expm1) because values are logged. Also how Seurat's DotPlot function does it
final_plot_data <- gene_plot_data %>% dplyr::group_by(manualAnnotation, new_inf, gene) %>% 
  dplyr::summarise(avg_expression = mean(expm1(x = expression))) %>% 
  dplyr::group_by(gene) %>% 
  dplyr::mutate(scaled_expression = scale(log1p(avg_expression))[,1]) %>% 
  dplyr::arrange(gene) %>% 
  replace(is.na(.), 0) 

final_plot_data %>% 
  ggplot(aes(x = new_inf, y = factor(gene), fill = scaled_expression))+
  geom_tile()+
  facet_wrap(~manualAnnotation, nrow = 1)+
  theme(strip.background = element_blank(), strip.placement = "outside",
        strip.text.x = element_text(angle = 80),
        strip.text = element_text(size=20),
        axis.text.y = element_text(size=16))+
  scale_x_discrete(position = "top") +
  scale_fill_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white","lightblue"),
                       values = c(1.0,0.7,0.4,0.35,-0.1),
                       limits = c(-2.5, 5))+
  ylab('')+
  xlab('')+
  ggtitle('Single Nuclei WT Chemokine Expression')

#Look at just astrocytes
sn_astrocytes <- subset(sn_integrated_dat, manualAnnotation == 'Astrocytes' & new_genotype %in% c('wt', 'wt (same)'))

cyto_data_astro <- FetchData(object = sn_astrocytes, vars = c(ccl_chemokines, cxcl_chemokines), layer = "data") %>% 
  tibble::rownames_to_column(var = 'cell_id') %>% 
  tidyr::pivot_longer(!cell_id, names_to = 'gene', values_to = 'expression')

sn_astrocytes <- prepSeuratObj(sn_astrocytes)
ElbowPlot(sn_astrocytes, ndims = 40)
sn_astrocytes <- prepUmapSeuratObj(sn_astrocytes, nDims = 15, reductionName = 'sn.astrocyte.umap')

DimPlot(sn_astrocytes, reduction = 'sn.astrocyte.umap', group.by = 'new_inf')

FeaturePlot(sn_astrocytes, reduction = 'sn.astrocyte.umap', features = 'Cxcl14')

#Chemotaxis genes dotplot

FeaturePlot(sn_integrated_dat_wt, features = ccl_chemokines, reduction = 'wt.umap.integrated')
FeaturePlot(sn_integrated_dat_wt, features = cxcl_chemokines, reduction = 'wt.umap.integrated')
FeaturePlot(sn_integrated_dat_wt, features = 'Cxcl10', reduction = 'wt.umap.integrated')


DotPlot(sn_integrated_dat_wt, features = ccl_chemokines, group.by = 'treatment_celltype', scale = TRUE)+
  theme(axis.text.x = element_text(angle = 90))
DotPlot(sn_integrated_dat_wt, features = cxcl_chemokines, group.by = 'treatment_celltype', scale = TRUE)+
  theme(axis.text.x = element_text(angle = 90))

#Combine cxcl10 and select cl for plotting
pdf('~/Documents/ÖverbyLab/single_nuclei_proj/sn_plots/select_chemokines_dotplot.pdf', width = 8, height = 6)
DotPlot(sn_integrated_dat_wt, features = c('Ccl2',  'Ccl3', 'Ccl5', 'Ccl7', 'Ccl11', 'Ccl12', 'Cxcl10'), group.by = 'treatment_celltype', scale = TRUE)+
  theme(axis.text.x = element_text(angle = 90))+
  ylab('')+
  xlab('')
dev.off()




