library(Seurat)
library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(tidyr)

#Load data
ParseSeuratObj_int <- LoadSeuratRds("~/Documents/ÖverbyLab/data/FilteredRpcaIntegratedDatNoDoublets.rds") 
ParseSeuratObj_int$hasVirus = ifelse(ParseSeuratObj_int$virusCountPAdj >= 10, 1, 0)


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

#Genes of interest for heatmap
chemokines <- c('Ccl7', 'Ccl2', 'Ccl12', 'Ccl4', 'Ccl3', 'Ccl5',
                'Cxcl10', 'Cxcl9', 'Cxcr3', 'Ccr1', 'Ccr2', 'Ccr5', 'Ccr3', 'Cxcl11', 'Cx3cr1')

chemokine_expression <- FetchData(object = chimeric_mock, vars = chemokines, layer = "data") %>% 
  rownames_to_column(var = 'cell_id') %>% 
  tidyr::pivot_longer(!cell_id, names_to = 'gene', values_to = 'expression')

heatmap_meta <- chimeric_mock[[]] %>% rownames_to_column(var = 'cell_id')

chemokine_expression_plot_dat <- dplyr::left_join(dplyr::select(heatmap_meta, c(cell_id, Genotype, manualAnnotation, Treatment)), 
                                                  chemokine_expression, by = 'cell_id')

#Copying seurat scaling strategy - find mean first then scale log1p. Found in seurat visualization github code
chemokine_expression_plot_dat %>% dplyr::group_by(manualAnnotation, Genotype, gene) %>% 
  dplyr::summarise(avg_expression = mean(expression)) %>% 
  dplyr::mutate(scaled_expression = scale(log1p(avg_expression))[,1]) %>% 
  dplyr::arrange(gene) %>% 
  ggplot(aes(x = Genotype, y = gene, fill = scaled_expression))+
  geom_tile()+
  facet_wrap(~manualAnnotation, nrow = 1)+
  theme(strip.background = element_blank(), strip.placement = "outside",
        strip.text.x = element_text(angle = 80))+
  scale_x_discrete(position = "top") +
  scale_fill_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white","lightblue"),
                       values = c(1.0,0.7,0.4,0.3,-0.1),
                       limits = c(-1.5, 3.5))

#Arrange separate plots of use faceting to combine heatmaps?


#Testing scale funciton to copy it
test_scaling <- DotPlot(chimeric_mock, features  = chemokines, group.by = 'manualAnnotation')$data
test_scaling %>% dplyr::arrange(features.plot) %>% 
  dplyr::group_by(features.plot) %>% dplyr::mutate(scaled_exp_testing = (scale(log1p(avg.exp))))






