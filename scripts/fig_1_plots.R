library(Seurat)
library(ggplot2)
library(dplyr)

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

chimeric_mock <- subset(ParseSeuratObj_int, Treatment != 'rLGTV' & manualAnnotation != 'unknown')

pdf("~/Documents/ÖverbyLab/single_cell_ISG_figures/fig_1_plots/chimeric_mock_umap.pdf", width = 9, height = 7)
DimPlot(chimeric_mock, label = FALSE, group.by = 'manualAnnotation', reduction = 'umap.integrated',
        cols = newCols)+
  theme(axis.ticks = element_blank(),
        axis.text=element_blank(),
        legend.text=element_text(size=17))+
  xlab('Umap1')+
  ylab('Umap2')+
  guides(color=guide_legend(override.aes=list(size=8)))+
  ggtitle('')
dev.off()

#Split data for plotting
chLgtv_ips <- subset(chimeric_mock, Treatment == 'rChLGTV' & Genotype == 'IPS1')
chLgtv_wt <- subset(chimeric_mock, Treatment == 'rChLGTV' & Genotype == 'WT')
mock_ips <- subset(chimeric_mock, Treatment == 'PBS'  & Genotype == 'IPS1')
mock_wt <- subset(chimeric_mock, Treatment == 'PBS'  & Genotype == 'WT')


#IPS chLGTV
pdf("~/Documents/ÖverbyLab/single_cell_ISG_figures/fig_1_plots/ChLgtv_ips_cellPopBars.pdf", width = 5, height = 8)
table(chLgtv_ips$Timepoint, chLgtv_ips$manualAnnotation, chLgtv_ips$Organ) %>% 
  as.data.frame() %>% dplyr::group_by(Var1) %>% dplyr::mutate(freq_props = Freq/sum(Freq))%>% 
  ggplot(aes(x = Var1, y = freq_props, fill = Var2))+
  geom_bar(stat = 'identity', position = 'stack', width = 0.6)+
  scale_fill_manual(values = newCols)+
  theme_classic()+
  theme(legend.position = 'none',
        axis.text.x = element_text(size = 22),
        axis.text.y = element_text(size = 22),
        axis.title=element_text(size=22),
        plot.title = element_text(size = 22))+
  xlab('')+
  ggtitle('ChLGTV IPS1')+
  ylab('Proportion of cells')
dev.off()

#WT chLGTV
pdf("~/Documents/ÖverbyLab/single_cell_ISG_figures/fig_1_plots/ChLgtv_wt_cellPopBars.pdf", width = 5, height = 8)
table(chLgtv_wt$Timepoint, chLgtv_wt$manualAnnotation, chLgtv_wt$Organ) %>% 
  as.data.frame() %>% dplyr::group_by(Var1) %>% dplyr::mutate(freq_props = Freq/sum(Freq))%>% 
  ggplot(aes(x = Var1, y = freq_props, fill = Var2))+
  geom_bar(stat = 'identity', position = 'stack', width = 0.6)+
  scale_fill_manual(values = newCols)+
  theme_classic()+
  theme(legend.position = 'none',
        axis.text.x = element_text(size = 22),
        axis.text.y = element_text(size = 22),
        axis.title=element_text(size=22),
        plot.title = element_text(size = 22))+
  xlab('')+
  ggtitle('ChLGTV WT')+
  ylab('Proportion of cells')
dev.off()

#Mock IPS
pdf("~/Documents/ÖverbyLab/single_cell_ISG_figures/fig_1_plots/mock_ips_cellPopBars.pdf", width = 5, height = 8)
table(mock_ips$Timepoint, mock_ips$manualAnnotation) %>% 
  as.data.frame() %>% dplyr::group_by(Var1) %>% dplyr::mutate(freq_props = Freq/sum(Freq))%>% 
  ggplot(aes(x = Var1, y = freq_props, fill = Var2))+
  geom_bar(stat = 'identity', position = 'stack', width = 0.6)+
  scale_fill_manual(values = newCols)+
  theme_classic()+
  theme(legend.position = 'none',
        axis.text.x = element_text(size = 22),
        axis.text.y = element_text(size = 22),
        axis.title=element_text(size=22),
        plot.title = element_text(size = 22))+
  xlab('')+
  ggtitle('Mock IPS1')+
  ylab('Proportion of cells')
dev.off()

#Mock WT
pdf("~/Documents/ÖverbyLab/single_cell_ISG_figures/fig_1_plots/mock_wt_cellPopBars.pdf", width = 5, height = 8)
table(mock_wt$Timepoint, mock_wt$manualAnnotation) %>% 
  as.data.frame() %>% dplyr::group_by(Var1) %>% dplyr::mutate(freq_props = Freq/sum(Freq))%>% 
  ggplot(aes(x = Var1, y = freq_props, fill = Var2))+
  geom_bar(stat = 'identity', position = 'stack', width = 0.6)+
  scale_fill_manual(values = newCols)+
  theme_classic()+
  theme(legend.position = 'none',
        axis.text.x = element_text(size = 22),
        axis.text.y = element_text(size = 22),
        axis.title=element_text(size=22),
        plot.title = element_text(size = 22))+
  xlab('')+
  ggtitle('Mock WT')+
  ylab('Proportion of cells')
dev.off()

chimeric_mock_organ_counts <- chimeric_mock[[]] %>% 
  group_by(Genotype, Treatment, Organ) %>% 
  dplyr::summarise(cell_count = n()) 

pdf("~/Documents/ÖverbyLab/single_cell_ISG_figures/fig_1_plots/cerebrum_cellPopBars.pdf", width = 5, height = 8)
chimeric_mock_organ_counts %>% 
  dplyr::filter(Organ == 'Cerebrum') %>% 
  ggplot(aes(x = Treatment, y = cell_count, fill = Genotype))+
  geom_bar(stat = 'identity', position = 'dodge')+
  xlab('')+
  ylab('cell count')+
  scale_fill_manual(values = c(newCols[4], newCols[7]))+
  theme_classic()+
  theme(legend.position = 'none',
        axis.text.x = element_text(size = 22),
        axis.text.y = element_text(size = 22),
        axis.title=element_text(size=22),
        plot.title = element_text(size = 22))+
  ggtitle('Cerebrum')+
  ylim(c(0, 26000))
dev.off()

pdf("~/Documents/ÖverbyLab/single_cell_ISG_figures/fig_1_plots/cerebellum_cellPopBars.pdf", width = 5, height = 8)
chimeric_mock_organ_counts %>% 
  dplyr::filter(Organ == 'Cerebellum') %>% 
  ggplot(aes(x = Treatment, y = cell_count, fill = Genotype))+
  geom_bar(stat = 'identity', position = 'dodge')+
  xlab('')+
  ylab('cell count')+
  scale_fill_manual(values = c(newCols[4], newCols[7]))+
  theme_classic()+
  theme(axis.text.x = element_text(size = 22),
        axis.text.y = element_text(size = 22),
        axis.title=element_text(size=22),
        plot.title = element_text(size = 22))+
  ggtitle('Cerebellum')+
  ylim(c(0, 26000))
dev.off()

chimeric_mock$manualAnnotation <- factor(chimeric_mock$manualAnnotation, 
                                              levels = c( 'unknown',  'T cells',   'Nk cells', 'Macrophage/Monocytes', 
                                                          'Granulocytes', 'B Cells',  'Pericytes', 'Oligodendrocytes','Neurons',
                                                          'Muscle cells', 'Microglia', 'Immature Neurons', 'Ependymal','Endothelial', 
                                                          'Choroid Plexus', 'Astrocytes'))


pdf('~/Documents/ÖverbyLab/single_cell_ISG_figures/fig_1_plots/mock_chimeric_cell_type_dotplot.pdf', height = 6, width = 12)
DotPlot(chimeric_mock, features = c('Aqp4', 'Fgfr3','Gfap', 'Ttr','Kl', 
                                       'Flt1', 'Pecam1','Vwf',
                                       'Cfap54' ,'Nnat', 'Mia',
                                       'Sox11', 'Celf4','Csf1r', 'Cx3cr1', 
                                       'Tmem119', 'Acta2', 'Tagln',
                                       'Snap25', 'Syt1',  
                                       'Mag', 'Mog', 'Abcc9', 'Vtn', 
                                       'Cd19', 'Ms4a1', 'S100a9', 'Il1r2', 
                                       'Klra2', 'Ccr2', 'Lyz2', 
                                       'Clnk','Nkg7','Cd3e', 'Cd3d'), 
        group.by = 'manualAnnotation', assay = 'RNA')+
  theme(axis.text.x = element_text(angle = 75, vjust = 0.5),
        axis.text = element_text(size = 16))+
  scale_color_gradient2(low = 'blue', mid = 'white', high = 'orange', labels = c(-1, 0, 1,2),
                        breaks = c(-1, 0, 1,2), , limits = c(-2,2.5))+
  # scale_color_gradientn(colours = c('lightblue','white', '#FFD991', '#FF4024'), 
  #              values = c(0, 0.35, 0.5, 1),
  #             name = 'Average Expression')+
  scale_size_continuous(range = c(0,8), limits = c(0,100))+
  theme(legend.position = "bottom",
        legend.justification = "center",
        legend.direction = "horizontal",
        legend.title = element_text(hjust = 0.3),
        legend.spacing.x = unit(2, "cm"))+
  guides(size = guide_legend(title.position = "top", title = 'Percent Expressed', theme = theme(legend.key.width  = unit(1, "cm"))),
         color = guide_colorbar(title.position = "top", title = 'Average Scaled Expression', theme = theme(legend.key.width  = unit(4, "cm"))))
dev.off()
geno_treatment
#Mavs expression
chimeric_mock$geno_treatment <- paste(chimeric_mock$Genotype, chimeric_mock$Treatment, sep = '_')
mavs_dat <- DotPlot(object = chimeric_mock, features = c("Mavs"), group.by = 'geno_treatment', scale = FALSE)$data

mavs_dat_meta <- str_split_fixed(mavs_dat$id, "_", 2)
colnames(mavs_dat_meta) <- c('genotype', 'treatment')
mavs_dat <- cbind(mavs_dat, mavs_dat_meta)

#Could split by time
pdf('~/Documents/ÖverbyLab/single_cell_ISG_figures/fig_1_plots/mavs_expression.pdf', height = 9, width = 7)
ggplot(mavs_dat, aes(x = genotype, y = treatment, color = avg.exp.scaled, size = pct.exp))+
  geom_point()+
  scale_color_gradient2(low = 'white', mid = 'orange', high = 'red', midpoint = 0.06)+
  theme_classic()+
  theme(axis.text = element_text(size = 24),
        legend.text = element_text(size = 24),
        legend.title = element_text(size = 24),
        plot.title = element_text(size =30))+
  xlab('')+
  ylab('')+
  ggtitle('Mavs Expression')+
  scale_size_continuous(range = c(3,9))
dev.off()



