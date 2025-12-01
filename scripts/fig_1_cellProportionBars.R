
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

#Split data for plotting
lgtv_ips <- subset(ParseSeuratObj_int, Treatment == 'rLGTV' & Genotype == 'IPS1')
lgtv_wt <- subset(ParseSeuratObj_int, Treatment == 'rLGTV' & Genotype == 'WT')
chLgtv_ips <- subset(ParseSeuratObj_int, Treatment == 'rChLGTV' & Genotype == 'IPS1')
chLgtv_wt <- subset(ParseSeuratObj_int, Treatment == 'rChLGTV' & Genotype == 'WT')
mock_ips <- subset(ParseSeuratObj_int, Treatment == 'PBS'  & Genotype == 'IPS1')
mock_wt <- subset(ParseSeuratObj_int, Treatment == 'PBS'  & Genotype == 'WT')

#IPS LGTV
pdf("~/Documents/ÖverbyLab/single_cell_ISG_figures/fig_1_plots/lgtv_ips_cellPopBars.pdf", width = 5, height = 8)
table(lgtv_ips$Timepoint, lgtv_ips$manualAnnotation) %>% 
  as.data.frame() %>% dplyr::group_by(Var1) %>% dplyr::mutate(freq_props = Freq/sum(Freq))%>% 
  ggplot(aes(x = Var1, y = freq_props, fill = Var2))+
  geom_bar(stat = 'identity', position = 'stack', width = 0.6)+
  scale_fill_manual(values = newCols)+
  theme_classic()+
  theme(legend.position = 'none',
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title=element_text(size=16),
        plot.title = element_text(size = 20))+
  xlab('')+
  ggtitle('LGTV IPS1')+
  ylab('Proportion of cells')
dev.off()

#WT LGTV
pdf("~/Documents/ÖverbyLab/single_cell_ISG_figures/fig_1_plots/lgtv_wt_cellPopBars.pdf", width = 5, height = 8)
table(lgtv_wt$Timepoint, lgtv_wt$manualAnnotation) %>% 
  as.data.frame() %>% dplyr::group_by(Var1) %>% dplyr::mutate(freq_props = Freq/sum(Freq))%>% 
  ggplot(aes(x = Var1, y = freq_props, fill = Var2))+
  geom_bar(stat = 'identity', position = 'stack', width = 0.6)+
  scale_fill_manual(values = newCols)+
  theme_classic()+
  theme(legend.position = 'none',
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title=element_text(size=16),
        plot.title = element_text(size = 20))+
  xlab('')+
  ggtitle('LGTV WT')+
  ylab('Proportion of cells')
dev.off()

#IPS chLGTV
pdf("~/Documents/ÖverbyLab/single_cell_ISG_figures/fig_1_plots/ChLgtv_ips_cellPopBars.pdf", width = 5, height = 8)
table(chLgtv_ips$Timepoint, chLgtv_ips$manualAnnotation) %>% 
  as.data.frame() %>% dplyr::group_by(Var1) %>% dplyr::mutate(freq_props = Freq/sum(Freq))%>% 
  ggplot(aes(x = Var1, y = freq_props, fill = Var2))+
  geom_bar(stat = 'identity', position = 'stack', width = 0.6)+
  scale_fill_manual(values = newCols)+
  theme_classic()+
  theme(legend.position = 'none',
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title=element_text(size=16),
        plot.title = element_text(size = 20))+
  xlab('')+
  ggtitle('ChLGTV IPS1')+
  ylab('Proportion of cells')
dev.off()

#WT chLGTV
pdf("~/Documents/ÖverbyLab/single_cell_ISG_figures/fig_1_plots/ChLgtv_wt_cellPopBars.pdf", width = 5, height = 8)
table(chLgtv_wt$Timepoint, chLgtv_wt$manualAnnotation) %>% 
  as.data.frame() %>% dplyr::group_by(Var1) %>% dplyr::mutate(freq_props = Freq/sum(Freq))%>% 
  ggplot(aes(x = Var1, y = freq_props, fill = Var2))+
  geom_bar(stat = 'identity', position = 'stack', width = 0.6)+
  scale_fill_manual(values = newCols)+
  theme_classic()+
  theme(legend.position = 'none',
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title=element_text(size=16),
        plot.title = element_text(size = 20))+
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
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title=element_text(size=16),
        plot.title = element_text(size = 20))+
  xlab('')+
  ggtitle('ChLGTV IPS1')+
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
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title=element_text(size=16),
        plot.title = element_text(size = 20))+
  xlab('')+
  ggtitle('ChLGTV IPS1')+
  ylab('Proportion of cells')
dev.off()






## Looking at data with lgtv removed
#IPS infected
pdf("~/Documents/ÖverbyLab/single_cell_ISG_figures/fig_1_plots/chimeric_ips_cellPopBars.pdf", width = 5, height = 7)
table(subset(chimeric_mock_ips, Treatment == 'rChLGTV')$Timepoint, subset(chimeric_mock_ips, Treatment == 'rChLGTV')$manualAnnotation) %>% 
  as.data.frame() %>% dplyr::group_by(Var1) %>% dplyr::mutate(freq_props = Freq/sum(Freq))%>% 
  ggplot(aes(x = Var1, y = freq_props, fill = Var2))+
  geom_bar(stat = 'identity', position = 'stack', width = 0.6)+
  scale_fill_manual(values = newCols)+
  theme_classic()+
  theme(legend.position = 'none',
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title=element_text(size=16),
        plot.title = element_text(size = 20))+
  xlab('')+
  ggtitle('IPS infected')+
  ylab('Proportion of cells')
dev.off()

#IPS mock
pdf("~/Documents/ÖverbyLab/single_cell_ISG_figures/fig_1_plots/chimeric_ips_mock_cellPopBars.pdf", width = 5, height = 7)
table(subset(chimeric_mock_ips, Treatment == 'PBS')$Timepoint, subset(chimeric_mock_ips, Treatment == 'PBS')$manualAnnotation) %>% 
  as.data.frame() %>% dplyr::group_by(Var1) %>% dplyr::mutate(freq_props = Freq/sum(Freq))%>% 
  ggplot(aes(x = Var1, y = freq_props, fill = Var2))+
  geom_bar(stat = 'identity', position = 'stack', width = 0.6)+
  scale_fill_manual(values = newCols)+
  theme_classic()+
  theme(legend.position = 'none',
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title=element_text(size=16),
        plot.title = element_text(size = 20))+
  xlab('')+
  ggtitle('IPS mock')+
  ylab('Proportion of cells')
dev.off()

#WT infected
pdf("~/Documents/ÖverbyLab/single_cell_ISG_figures/fig_1_plots/wt_ips_cellPopBars.pdf", width = 5, height = 7)
table(subset(chimeric_mock_wt, Treatment == 'rChLGTV')$Timepoint, subset(chimeric_mock_wt, Treatment == 'rChLGTV')$manualAnnotation) %>% 
  as.data.frame() %>% dplyr::group_by(Var1) %>% dplyr::mutate(freq_props = Freq/sum(Freq))%>% 
  ggplot(aes(x = Var1, y = freq_props, fill = Var2))+
  geom_bar(stat = 'identity', position = 'stack', width = 0.6)+
  scale_fill_manual(values = newCols)+
  theme_classic()+
  theme(legend.position = 'none',
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title=element_text(size=16),
        plot.title = element_text(size = 20))+
  xlab('')+
  ggtitle('WT Infected')+
  ylab('Proportion of cells')
dev.off()

#WT mock
pdf("~/Documents/ÖverbyLab/single_cell_ISG_figures/fig_1_plots/wt_ips_mock_cellPopBars.pdf", width = 5, height = 7)
table(subset(chimeric_mock_wt, Treatment == 'PBS')$Timepoint, subset(chimeric_mock_wt, Treatment == 'PBS')$manualAnnotation) %>% 
  as.data.frame() %>% dplyr::group_by(Var1) %>% dplyr::mutate(freq_props = Freq/sum(Freq))%>% 
  ggplot(aes(x = Var1, y = freq_props, fill = Var2))+
  geom_bar(stat = 'identity', position = 'stack', width = 0.6)+
  scale_fill_manual(values = newCols)+
  theme_classic()+
  theme(legend.position = 'none',
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title=element_text(size=16),
        plot.title = element_text(size = 20))+
  xlab('')+
  ggtitle('WT mock')+
  ylab('Proportion of cells')
dev.off()
