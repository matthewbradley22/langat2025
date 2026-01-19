library(Seurat)
library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(tidyr)
library(tibble)

#Load data
ParseSeuratObj_int <- LoadSeuratRds("~/Documents/ÖverbyLab/data/FilteredRpcaIntegratedDatNoDoublets.rds") 
ParseSeuratObj_int$hasVirus = ifelse(ParseSeuratObj_int$virusCountPAdj >= 10, 1, 0)
ParseSeuratObj_int$manualAnnotation[ParseSeuratObj_int$manualAnnotation == 'Macrophage/Monocytes'] = 'Macro/Mono'

#Load in single nuclei data
sn_integrated_dat <- LoadSeuratRds('~/Documents/ÖverbyLab/single_nuclei_proj/LGTVscCombined.rds')
sn_integrated_dat_wt <- subset(sn_integrated_dat, new_genotype %in% c('wt', 'wt (same)')) #only wt for gal 3 project
sn_integrated_dat_wt$genotype_same = 'WT'
table(sn_integrated_dat_wt$new_genotype, sn_integrated_dat_wt$genotype_same)

#Split sn data by infection for plotting
sn_integrated_dat_wt_infected <- subset(sn_integrated_dat_wt, infected == TRUE)
sn_integrated_dat_wt_mock_none <- subset(sn_integrated_dat_wt, infected == FALSE)

#Check data
newCols <-  c(brewer.pal(12, 'Paired'), '#99FFE6', '#CE99FF', '#18662E','#737272',  '#FF8AEF')
newCols[11] =  '#FF8AEF'

DimPlot(ParseSeuratObj_int, label = FALSE, group.by = 'manualAnnotation', reduction = 'umap.integrated',
        cols = newCols)+
  theme(axis.ticks = element_blank(),
        axis.text=element_blank(),
        legend.text=element_text(size=17))+
  xlab('Umap1')+
  ylab('Umap2')+
  guides(color=guide_legend(override.aes=list(size=8)))+
  ggtitle('')

#Look at total rna captured from celltypes
VlnPlot(ParseSeuratObj_int, features = 'nCount_RNA', group.by = 'manualAnnotation', pt.size = 0)

#subset data (no LGTV)
chimeric_mock <- subset(ParseSeuratObj_int, Treatment != 'rLGTV')
lgtv <- subset(ParseSeuratObj_int, Treatment == 'rLGTV')

#Split by treatment for plotting
chimeric_mock_infected <- subset(ParseSeuratObj_int, Treatment == 'rChLGTV')
chimeric_mock_mock <- subset(ParseSeuratObj_int, Treatment == 'PBS')

#Split by genotype
chimeric_mock_infected_wt <- subset(chimeric_mock_infected, Genotype == 'WT')
chimeric_mock_infected_ips <- subset(chimeric_mock_infected, Genotype == 'IPS1')

chimeric_mock_mock_wt <- subset(chimeric_mock_mock, Genotype == 'WT')
chimeric_mock_mock_ips <- subset(chimeric_mock_mock, Genotype == 'IPS1')

lgtv_wt <- subset(lgtv, Genotype == 'WT')
lgtv_ips <- subset(lgtv, Genotype == 'IPS1')

#Subset to cell types with at least 100 cells
infected_cells_meeting_min <- table(chimeric_mock_infected$manualAnnotation) %>% as.data.frame() %>% dplyr::filter(Freq > 100) %>% 
  dplyr::pull(Var1)

#Don't need to remove any cell types from infection samples
length(infected_cells_meeting_min)

#Subset to cell types with at least 100 cells
infected_wt_cells_meeting_min <- table(chimeric_mock_infected_wt$manualAnnotation) %>% as.data.frame() %>% dplyr::filter(Freq > 100) %>% 
  dplyr::pull(Var1)

infected_ips_cells_meeting_min <- table(chimeric_mock_infected_ips$manualAnnotation) %>% as.data.frame() %>% dplyr::filter(Freq > 100) %>% 
  dplyr::pull(Var1)

mock_wt_cells_meeting_min <- table(chimeric_mock_mock_wt$manualAnnotation) %>% as.data.frame() %>% dplyr::filter(Freq > 100) %>% 
  dplyr::pull(Var1)

mock_ips_cells_meeting_min <- table(chimeric_mock_mock_ips$manualAnnotation) %>% as.data.frame() %>% dplyr::filter(Freq > 100) %>% 
  dplyr::pull(Var1)

lgtv_wt_cells_meeting_min <- table(lgtv_wt$manualAnnotation) %>% as.data.frame() %>% dplyr::filter(Freq > 100) %>% 
  dplyr::pull(Var1)

lgtv_ips_cells_meeting_min <- table(lgtv_ips$manualAnnotation) %>% as.data.frame() %>% dplyr::filter(Freq > 100) %>% 
  dplyr::pull(Var1)

chimeric_mock_infected_wt <- subset(chimeric_mock_infected_wt, manualAnnotation %in% infected_wt_cells_meeting_min)
chimeric_mock_infected_ips <- subset(chimeric_mock_infected_ips, manualAnnotation %in% infected_ips_cells_meeting_min)
chimeric_mock_mock_wt <- subset(chimeric_mock_mock_wt, manualAnnotation %in% mock_wt_cells_meeting_min)
chimeric_mock_mock_ips <- subset(chimeric_mock_mock_ips, manualAnnotation %in% mock_ips_cells_meeting_min)
lgtv_wt <- subset(lgtv_wt, manualAnnotation %in% lgtv_wt_cells_meeting_min)
lgtv_ips <- subset(lgtv_ips, manualAnnotation %in% lgtv_wt_cells_meeting_min)

#Genes of interest for heatmap
chemokines <- factor(c('Ccl2', 'Ccl3', 'Ccl4',  'Ccl5', 'Ccl7', 'Ccl12',  'Ccr1',
                       'Ccr2', 'Ccr3', 'Ccr5', 'Cxcl9', 
                'Cxcl10',  'Cxcl11', 'Cxcr3',  'Cx3cr1'),
                levels = c('Ccl2', 'Ccl3', 'Ccl4',  'Ccl5', 'Ccl7', 'Ccl12',  'Ccr1',
                           'Ccr2', 'Ccr3', 'Ccr5', 'Cxcl9', 
                           'Cxcl10',  'Cxcl11', 'Cxcr3',  'Cx3cr1'))

interesting_genes_astrocytes <- factor(c('Cx3cl1', 'Il1b', 'Tnf', 'Il6', 'Il12', 'Vegf',
                                  'Bdnf', 'Gdnf', 'Fgf2', 'Mmp2', 'Mmp9',
                                  'Pge2', 'S1p', 'Il10', 'Tgfb1'),
                                  levels = c('Cx3cl1', 'Il1b', 'Tnf', 'Il6', 'Il12', 'Vegf',
                                             'Bdnf', 'Gdnf', 'Fgf2', 'Mmp2', 'Mmp9',
                                             'Pge2', 'S1p', 'Il10', 'Tgfb1'))

#Astro sig markers from script 'scRNA_astrocyte_analysis.R'
highest_upregulated_astrocytes = factor(rownames(astro_sig_markers)[1:15])

#Genes of interest for astrocytes

#Create function for genotype by gene heatmaps faceted by cell type
#Copying seurat scaling strategy - find mean first then scale log1p. Found in seurat visualization github code
#Split by infiltrating and resident?
faceted_geno_heatmap <- function(dat, genes, main = '', geno_column = NULL, returnData = FALSE){
  #Get expression data
  gene_data <- FetchData(object = dat, vars = genes, layer = "data") %>% 
    tibble::rownames_to_column(var = 'cell_id') %>% 
    tidyr::pivot_longer(!cell_id, names_to = 'gene', values_to = 'expression')
  
  #Create metadata dataframe
  cell_metadata <- dat[[]] %>% rownames_to_column(var = 'cell_id')
  
  gene_plot_data <- dplyr::left_join(dplyr::select(cell_metadata, c(cell_id, !!sym(geno_column), manualAnnotation)), 
                                     gene_data, by = 'cell_id')
  
  #Doing mean of exponential (expm1) because values are logged. Also how Seurat's DotPlot function does it
  final_plot_data <- gene_plot_data %>% dplyr::group_by(manualAnnotation, !!sym(geno_column), gene) %>% 
    dplyr::summarise(avg_expression = mean(expm1(x = expression))) %>% 
    dplyr::group_by(gene) %>% 
    dplyr::mutate(scaled_expression = scale(log1p(avg_expression))[,1]) %>% 
    dplyr::arrange(gene) %>% 
    replace(is.na(.), 0) 
  
  if(returnData){
    return(final_plot_data)
  }
  #Create faceted heatmap
  gene_plot <- final_plot_data %>% 
    ggplot(aes(x = !!sym(geno_column), y = factor(gene,levels = rev(levels(genes))), fill = scaled_expression))+
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
    ggtitle(main)
  print(gene_plot)
}

#Make chemokine gene plots
#rChLGTV
pdf('~/Documents/ÖverbyLab/scPlots/gene_heatmaps/chemokine_heatmap_infected_wt.pdf', width = 14, height = 7)
faceted_geno_heatmap(dat = chimeric_mock_infected_wt, genes = chemokines, geno_column = "Genotype",
                     main = 'Single cell: infected WT')
dev.off()

pdf('~/Documents/ÖverbyLab/scPlots/gene_heatmaps/chemokine_heatmap_infected_ips.pdf', width = 14, height = 7)
faceted_geno_heatmap(dat = chimeric_mock_infected_ips, genes = chemokines, geno_column = "Genotype",
                     main = 'Single cell: infected IPS')
dev.off()

DotPlot(chimeric_mock_mock_ips, features = 'Ccl2', group.by = 'manualAnnotation')$data

#Mock 
pdf('~/Documents/ÖverbyLab/scPlots/gene_heatmaps/chemokine_heatmap_mock_wt.pdf', width = 14, height = 7)
faceted_geno_heatmap(dat = chimeric_mock_mock_wt, genes = chemokines,  geno_column = 'Genotype',
                     main = 'Single cell: mock WT')
dev.off()

pdf('~/Documents/ÖverbyLab/scPlots/gene_heatmaps/chemokine_heatmap_mock_ips.pdf', width = 14, height = 7)
faceted_geno_heatmap(dat = chimeric_mock_mock_ips, genes = chemokines,  geno_column = 'Genotype',
                     main = 'Single cell: mock IPS')
dev.off()

#LGTV
pdf('~/Documents/ÖverbyLab/scPlots/gene_heatmaps/chemokine_heatmap_lgtv_wt.pdf', width = 14, height = 7)
faceted_geno_heatmap(dat = lgtv_wt, genes = chemokines, geno_column = "Genotype",
                     main = 'Single cell: LGTV WT')
dev.off()

pdf('~/Documents/ÖverbyLab/scPlots/gene_heatmaps/chemokine_heatmap_lgtv_ips.pdf', width = 14, height = 7)
faceted_geno_heatmap(dat = lgtv_ips, genes = chemokines, geno_column = "Genotype",
                     main = 'Single cell: LGTV WT')
dev.off()

#Single nuclei
pdf('~/Documents/ÖverbyLab/scPlots/gene_heatmaps/chemokine_heatmap_single_nuclei_infected.pdf', width = 14, height = 7)
faceted_geno_heatmap(dat = sn_integrated_dat_wt_infected, genes = chemokines,  geno_column = 'genotype_same',
                     main = 'Single nuclei: infected')
dev.off()

pdf('~/Documents/ÖverbyLab/scPlots/gene_heatmaps/chemokine_heatmap_single_nuclei_mock.pdf', width = 14, height = 7)
faceted_geno_heatmap(dat = sn_integrated_dat_wt_mock_none, genes = chemokines,  geno_column = 'genotype_same',
                     main = 'Single nuclei: mock')
dev.off()



#Make astrocyte gene plots
pdf('~/Documents/ÖverbyLab/scPlots/gene_heatmaps/astrocyte_heatmap_infected_wt.pdf', width = 14, height = 7)
faceted_geno_heatmap(dat = chimeric_mock_infected_wt, genes = interesting_genes_astrocytes, geno_column = "Genotype",
                     main = 'Single cell: infected WT')
dev.off()

pdf('~/Documents/ÖverbyLab/scPlots/gene_heatmaps/astrocyte_heatmap_infected_ips.pdf', width = 14, height = 7)
faceted_geno_heatmap(dat = chimeric_mock_infected_ips, genes = interesting_genes_astrocytes, geno_column = "Genotype",
                     main = 'Single cell: infected IPS')
dev.off()

pdf('~/Documents/ÖverbyLab/scPlots/gene_heatmaps/astrocyte_heatmap_mock_wt.pdf', width = 14, height = 7)
faceted_geno_heatmap(dat = chimeric_mock_mock_wt, genes = interesting_genes_astrocytes,  geno_column = 'Genotype',
                     main = 'Single cell: mock WT')
dev.off()

pdf('~/Documents/ÖverbyLab/scPlots/gene_heatmaps/astrocyte_heatmap_mock_ips.pdf', width = 14, height = 7)
faceted_geno_heatmap(dat = chimeric_mock_mock_ips, genes = interesting_genes_astrocytes,  geno_column = 'Genotype',
                     main = 'Single cell: mock IPS')
dev.off()

#Plots of most significantly upregulated genes in astrocytes
pdf('~/Documents/ÖverbyLab/scPlots/gene_heatmaps/astrocyteGenes_heatmap_infected_wt.pdf', width = 14, height = 7)
faceted_geno_heatmap(dat = chimeric_mock_infected_wt, genes = highest_upregulated_astrocytes, geno_column = "Genotype",
                     main = 'Single cell: infected WT')
dev.off()

pdf('~/Documents/ÖverbyLab/scPlots/gene_heatmaps/astrocyteGenes_heatmap_infected_ips.pdf', width = 14, height = 7)
faceted_geno_heatmap(dat = chimeric_mock_infected_ips, genes = highest_upregulated_astrocytes, geno_column = "Genotype",
                     main = 'Single cell: infected IPS')
dev.off()

pdf('~/Documents/ÖverbyLab/scPlots/gene_heatmaps/astrocyteGenes_heatmap_lgtv_wt.pdf', width = 14, height = 7)
faceted_geno_heatmap(dat = lgtv_wt, genes = highest_upregulated_astrocytes, geno_column = "Genotype",
                     main = 'Single cell: LGTV WT')
dev.off()

#Return data from function to compare directly celltypes
chLGTV_ips_data <- faceted_geno_heatmap(dat = chimeric_mock_infected_ips, genes = highest_upregulated_astrocytes, geno_column = "Genotype",
                     main = 'Single cell: infected WT', returnData = TRUE)
mock_ips_data <- faceted_geno_heatmap(dat = chimeric_mock_mock_ips, genes = highest_upregulated_astrocytes,  geno_column = 'Genotype',
                     main = 'Single cell: mock WT', returnData = TRUE)

left_join(chLGTV_ips_data, mock_ips_data, by = c('manualAnnotation', 'Genotype', 'gene'), suffix = c('_infected', '_mock')) %>% 
  dplyr::mutate(infected_vs_mock = avg_expression_infected - avg_expression_mock) %>% 
  dplyr::filter(!is.na(infected_vs_mock)) %>% 
  ggplot(aes(x = Genotype, y = gene, fill = infected_vs_mock))+
  geom_tile()+
  facet_wrap(~manualAnnotation, nrow = 1)+
  theme(strip.background = element_blank(), strip.placement = "outside",
        strip.text.x = element_text(angle = 80),
        strip.text = element_text(size=20),
        axis.text.y = element_text(size=16))+
  scale_x_discrete(position = "top") +
  scale_fill_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white","lightblue"),
                       values = c(1.0,0.7,0.6,0.48,-0.3),
                       limits = c(-13, 15))+
  ylab('')+
  xlab('')+
  ggtitle('Infected vs mock wt')

#Look at day 5 heatmaps
day5_cells = subset(ParseSeuratObj_int, Timepoint == 'Day 5' & Treatment == 'rChLGTV')
faceted_geno_heatmap(dat = day5_cells, genes = chemokines,  geno_column = 'Genotype',
                     main = '')

#Testing scale funciton to copy it
test_scaling <- DotPlot(chimeric_mock, features  = chemokines, group.by = 'manualAnnotation')$data
test_scaling %>% dplyr::arrange(features.plot) %>% 
  dplyr::group_by(features.plot) %>% dplyr::mutate(scaled_exp_testing = (scale(log1p(avg.exp))[,1])) %>% 
  dplyr::group_by(features.plot) %>% 
  dplyr::summarise(gene_mean = mean(avg.exp))


View(test_scaling)



