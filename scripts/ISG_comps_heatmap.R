library(Seurat)
library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(tidyr)

#Load data
ParseSeuratObj_int <- LoadSeuratRds("~/Documents/ÖverbyLab/data/FilteredRpcaIntegratedDatNoDoublets.rds") 
ParseSeuratObj_int$hasVirus = ifelse(ParseSeuratObj_int$virusCountPAdj >= 10, 1, 0)

#Load in single nuclei data
sn_integrated_dat <- LoadSeuratRds('~/Documents/ÖverbyLab/single_nuclei_proj/LGTVscCombined.rds')
sn_integrated_dat_wt <- subset(sn_integrated_dat, new_genotype %in% c('wt', 'wt (same)')) #only wt for gal 3 project

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

#subset data (no LGTV)
chimeric_mock <- subset(ParseSeuratObj_int, Treatment != 'rLGTV')

#Split by treatment for plotting
chimeric_mock_infected <- subset(ParseSeuratObj_int, Treatment == 'rChLGTV')
chimeric_mock_mock <- subset(ParseSeuratObj_int, Treatment == 'PBS')

#Split sn data by infection for plotting

#Subset to cell types with at least 100 cells
infected_cells_meeting_min <- table(chimeric_mock_infected$manualAnnotation) %>% as.data.frame() %>% dplyr::filter(Freq > 100) %>% 
  dplyr::pull(Var1)
#Don't need to remove any cell types from infection samples
length(infected_cells_meeting_min)


mock_cells_meeting_min <- table(chimeric_mock_mock$manualAnnotation) %>% as.data.frame() %>% dplyr::filter(Freq > 100) %>% 
  dplyr::pull(Var1)
chimeric_mock_mock <- subset(chimeric_mock_mock, manualAnnotation %in% mock_cells_meeting_min)

#Genes of interest for heatmap
chemokines <- c('Ccl7', 'Ccl2', 'Ccl12', 'Ccl4', 'Ccl3', 'Ccl5',
                'Cxcl10', 'Cxcl9', 'Cxcr3', 'Ccr1', 'Ccr2', 'Ccr5', 'Ccr3', 'Cxcl11', 'Cx3cr1')

#Create function for genotype by gene heatmaps faceted by cell type
#Copying seurat scaling strategy - find mean first then scale log1p. Found in seurat visualization github code
#Split by infiltrating and resident?
faceted_geno_heatmap <- function(dat, genes, main = '', geno_column = NULL){
  #Get expression data
  gene_data <- FetchData(object = dat, vars = genes, layer = "data") %>% 
    rownames_to_column(var = 'cell_id') %>% 
    tidyr::pivot_longer(!cell_id, names_to = 'gene', values_to = 'expression')
  
  #Create metadata dataframe
  cell_metadata <- dat[[]] %>% rownames_to_column(var = 'cell_id')
  
  gene_plot_data <- dplyr::left_join(dplyr::select(cell_metadata, c(cell_id, !!sym(geno_column), manualAnnotation, Treatment)), 
                                     gene_data, by = 'cell_id')
  #Create faceted heatmap
  gene_plot <- gene_plot_data %>% dplyr::group_by(manualAnnotation, !!sym(geno_column), gene) %>% 
    dplyr::summarise(avg_expression = mean(expression)) %>% 
    dplyr::group_by(gene) %>% 
    dplyr::mutate(scaled_expression = scale(log1p(avg_expression))[,1]) %>% 
    dplyr::arrange(gene) %>% 
    ggplot(aes(x = !!sym(geno_column), y = gene, fill = scaled_expression))+
    geom_tile()+
    facet_wrap(~manualAnnotation, nrow = 1)+
    theme(strip.background = element_blank(), strip.placement = "outside",
          strip.text.x = element_text(angle = 80),
          strip.text = element_text(size=16),
          axis.text.y = element_text(size=16))+
    scale_x_discrete(position = "top") +
    scale_fill_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white","lightblue"),
                         values = c(1.0,0.7,0.4,0.3,-0.1),
                         limits = c(-2, 5))+
    ylab('')+
    ggtitle(main)
  print(gene_plot)
}

pdf('~/Documents/ÖverbyLab/scPlots/chemokine_heatmap_infected.pdf', width = 14, height = 7)
faceted_geno_heatmap(dat = chimeric_mock_infected, genes = chemokines, geno_column = "Genotype")
dev.off()

pdf('~/Documents/ÖverbyLab/scPlots/chemokine_heatmap_mock.pdf', width = 14, height = 7)
faceted_geno_heatmap(chimeric_mock_mock, genes = chemokines,  geno_column = 'Genotype')
dev.off()


#Testing scale funciton to copy it
test_scaling <- DotPlot(chimeric_mock, features  = chemokines, group.by = 'manualAnnotation')$data
test_scaling %>% dplyr::arrange(features.plot) %>% 
  dplyr::group_by(features.plot) %>% dplyr::mutate(scaled_exp_testing = (scale(log1p(avg.exp))[,1])) %>% 
  dplyr::group_by(features.plot) %>% 
  dplyr::summarise(gene_mean = mean(avg.exp))


View(test_scaling)



