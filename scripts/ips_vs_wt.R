library(Seurat)
library(RColorBrewer)
library(dplyr)
library(UpSetR)
library(ggplot2)
library(tidyr)
source('~/Documents/ÖverbyLab//scripts/langatFunctions.R')

#Load data
ParseSeuratObj_int <- LoadSeuratRds("~/Documents/ÖverbyLab/data/FilteredRpcaIntegratedDatNoDoublets.rds") 

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

#Compare infected ips vs infected wt over all celltypes
chimeric_infected <- subset(ParseSeuratObj_int, Treatment == 'rChLGTV')

chimeric_infected <- prepSeuratObj(chimeric_infected)
ElbowPlot(chimeric_infected, ndims = 30)
chimeric_infected <- prepUmapSeuratObj(chimeric_infected, nDims = 20, reductionName = 'ch_umap')
DimPlot(chimeric_infected, group.by = 'Genotype', reduction = 'ch_umap')

#Markers upregulated in wt
wt_up <- FindMarkers(chimeric_infected, group.by = 'Genotype', ident.1 = 'WT', only.pos = TRUE, test.use = 'MAST')
wt_up_sig <- wt_up %>% as.data.frame() %>% dplyr::filter(p_val_adj < 0.01 & avg_log2FC > 1)

#Compare mavs by orig.ident to confirm it lines up with barcoding plate csv we have, 09 - 24 is ips
mavs_dat <- DotPlot(ParseSeuratObj_int, features = 'Mavs', group.by = 'orig.ident', scale = FALSE)$data
mavs_dat$id <- mavs_dat$id %>% as.character() %>% as.numeric()
mavs_dat <- mavs_dat %>% dplyr::mutate(genotype = ifelse(id <= 24, yes = 'IPS1', no = 'WT'))

ggplot(mavs_dat, aes(x = features.plot, y = id, size = pct.exp, color = avg.exp.scaled))+
  facet_wrap(~genotype)+
  geom_point()+ 
  theme(panel.background = element_rect(fill = 'white', colour = 'white'))+
  scale_color_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white"), 
                       values = c(1.0,0.7,0.4,0))+
  ylab('plate id')+
  xlab('')+
  ggtitle('MAVS expression')


