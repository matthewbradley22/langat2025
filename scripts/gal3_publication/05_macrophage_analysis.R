library(Seurat)
library(tibble)
library(UCell)
source('./langatFunctions.R')

#Load data
ParseSeuratObj_int <- LoadSeuratRds("~/Documents/ÖverbyLab/data/FilteredRpcaIntegratedDatNoDoublets.rds") 

#Subset to monocytes/macrophages in infected samples
macrophages_wt_infected <- subset(ParseSeuratObj_int, Treatment == 'rLGTV'& Organ == 'Cerebrum' & 
                                    Genotype == 'WT' & manualAnnotation == 'Macrophage/Monocytes')  

#Prepare umap for infected samples
macrophages_wt_infected <- prepSeuratObj(macrophages_wt_infected)
ElbowPlot(macrophages_wt_infected, ndims = 40)
macrophages_wt_infected <- prepUmapSeuratObj(macrophages_wt_infected, nDims = 20, reductionName = 'wt.infected.mac.umap',
                                             resolution_value = 0.8)

DimPlot(macrophages_wt_infected, reduction = 'wt.infected.mac.umap', label = FALSE, group.by = 'Timepoint',
        label.size = 8)+
  ggtitle('WT Cerebrum Macrophages')+
  xlab('UMAP 1')+
  ylab('UMAP 2')+
  theme(legend.text=element_text(size=16),
        plot.title = element_text(size = 22))+
  guides(colour = guide_legend(override.aes = list(size=8)))+
  scale_colour_manual(values = c("#292270","#166DF0", "#6DC3F8"))

#Proportions of cells at 
macrophages_wt_infected[[]] %>% dplyr::group_by(Timepoint) %>% 
  dplyr::summarise(count = n()) %>% 
  dplyr::mutate(prop = count / sum(count)) %>% 
  ggplot(aes(y = prop, x = '', fill = Timepoint))+
  geom_bar(position = 'stack', stat = 'identity')+
  theme_classic()+
  scale_fill_manual(values = c("#292270","#166DF0", "#6DC3F8"))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#### Differentially expressed genes between early and late timepoints ####
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


#Examine what separates day 5 from other timepoints
macrophages_wt_infected$time_collapsed <- ifelse(macrophages_wt_infected$Timepoint == 'Day 5', yes = 'Day 5',no = 'Day 3/4')

macro_markers <- FindMarkers(macrophages_wt_infected, group.by = 'time_collapsed', test.use = 'MAST',
                                  ident.1 = 'Day 5')

#Volcano plot of markers
macro_markers %>% dplyr::mutate(direction = case_when(p_val_adj < 0.01 & (avg_log2FC) > 1 ~ 'Day 5',
                                                           p_val_adj < 0.01 & (avg_log2FC) < -1 ~ 'Day 3/4',
                                                           .default = 'not_sig')) %>%
  ggplot(aes(x = avg_log2FC, y = -log10(p_val_adj), color = direction))+
  geom_point(alpha = 1)+
  scale_color_manual(values = c("#166DF0", "#6DC3F8", 'lightgrey'))+
  theme_classic()+
  geom_hline(yintercept = -log10(0.01))+
  geom_vline(xintercept = -1)+
  geom_vline(xintercept = 1)

sum(macro_markers$avg_log2FC > 1 & macro_markers$p_val_adj < 0.01)
sum(macro_markers$avg_log2FC < -1 & macro_markers$p_val_adj < 0.01)

#Plot top up and top down genes
top_up <- dplyr::filter(macro_markers, avg_log2FC > 1 & p_val_adj < 0.01) %>% head(n = 10)
top_down <- dplyr::filter(macro_markers, avg_log2FC < -1 & p_val_adj < 0.01) %>% head(n = 10)

deg_dotplot_dat <- DotPlot(macrophages_wt_infected, features = c(rownames(top_up), 
                           rownames(top_down)), 
                           group.by = 'Timepoint', scale = FALSE)$data

ggplot(deg_dotplot_dat, aes(x = features.plot, y = id))+
  geom_point(aes(fill = avg.exp.scaled, size = pct.exp), pch=21)+
  theme(axis.text.x = element_text(angle = 90, size = 16),
        axis.text.y = element_text(size = 16),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))+
  scale_size(range = c(0.01, 12), name = 'Percent expressed') +
  scale_fill_gradientn(colours = c('white', '#FFD991', '#FF7530', '#FF4024'), 
                       values = c(0, 0.3, 0.6, 1),
                       name = 'Average Expression')+
  ylab('')+
  xlab('')


# - - - - - - - - - - - - - - - - - - - - - - 
#### Monocyte/Macrophage identity scores ####
# - - - - - - - - - - - - - - - - - - - - - - 

#Macrophage gene sets from Fig. 2 https://link.springer.com/article/10.1186/s12974-016-0581-z

macrophage_subset_markers <- list(M1_signature = c('Tnf', 'Il1b', 'Il6', 'Il12a', 'Il23a', 'Il27', 'Nos2', 'Ido1', 'Irf5',
                                                   'Socs1', 'Marco', 'Cd80', 'Cd86', 'Cd274', 'C1d', 'C1qa', 'C1qb', 'C1qbp', 'C1qc',
                                                   'C1qtnf1', 'C3', 'C3ar1', 'C5ar1', 'Itgb2', 'Ccl2', 'Ccl3', 'Ccl4', 'Ccl5', 'Cxcl9',
                                                   'Cxcl10', 'Cxcl16', 'Ccr7'),
                                  Il4_alt_signature = c('Arg1', 'Mrc1', 'Chil3', 'Tgm2', 'Il1rn', 'Msr1', 'Pdcd1lg2', 'Cd209a', 'Clec10a', 'Retnla', 'Alox15',
                                                        'Socs2', 'Irf4', 'Ppard', 'Pparg', 'Ccl17', 'Ccl24'),
                                  Il10_M2c_signature = c('Il10', 'Il4ra', 'Tgfb1', 'Nfil3', 'Sbno2', 'Socs3', 'Fcgr1', 'Fcgr2b', 'Fcgr3', 'Marco',
                                                         'Cxcl13'),
                                  mhc2_sig = c('H2-Aa', 'H2-Ab1', 'H2-DMa', 'H2-DMb1', 'H2-DMb2', 'H2-Eb1'))

#Use UCell to calculate cell scores for gene signatures

#Calculate median gene count for maxRank parameter per the UCell docs
median_genes <- round(median(macrophages_wt_infected$nFeature_RNA))

#Macrophage scores
macrophages_wt_infected <- AddModuleScore_UCell(macrophages_wt_infected, features = macrophage_subset_markers, maxRank=median_genes)

#M1 sgnature
VlnPlot(macrophages_wt_infected, features = 'M1_signature_UCell', group.by = 'Timepoint', pt.size = 0, cols = c("#292270","#166DF0", "#6DC3F8"))+
  theme(legend.position = 'none')+
  ggtitle('Classical Activation Score')+
  geom_boxplot(width = 0.1, position = position_dodge(0.9), alpha = 0.5,
               fill = 'white')+
  ylim(c(-0.01, 0.4))+
  xlab('')

#Statistical comparison
kruskal.test(M1_signature_UCell ~ Timepoint, data = macrophages_wt_infected[[]])
dunn_test(macrophages_wt_infected[[]], M1_signature_UCell ~ Timepoint, p.adjust.method = "holm", detailed = TRUE)

#Alternative activation scores
VlnPlot(macrophages_wt_infected, features = 'Il4_alt_signature_UCell', group.by = 'Timepoint', pt.size = 0, cols = c("#292270","#166DF0", "#6DC3F8"))+
  theme(legend.position = 'none')+
  geom_boxplot(width = 0.1, position = position_dodge(0.9), alpha = 0.5,
               fill = 'white')+
  ylim(c(-0.01, 0.4))+
  ggtitle('Alternate Activation Score')

#Statistical comparison
kruskal.test(Il4_alt_signature_UCell ~ Timepoint, data = macrophages_wt_infected[[]])
dunn_test(macrophages_wt_infected[[]], Il4_alt_signature_UCell ~ Timepoint, p.adjust.method = "holm", detailed = FALSE)

#M2c score
VlnPlot(macrophages_wt_infected, features = 'Il10_M2c_signature_UCell', group.by = 'Timepoint', pt.size = 0,  cols = c("#292270","#166DF0", "#6DC3F8"))+
  theme(legend.position = 'none')+
  geom_boxplot(width = 0.1, position = position_dodge(0.9), alpha = 0.5,
               fill = 'white')+
  ylim(c(-0.01, 0.4))+
  ggtitle('M2c Score')

#Statistical comparison
kruskal.test(Il10_M2c_signature_UCell ~ Timepoint, data = macrophages_wt_infected[[]])
dunn_test(macrophages_wt_infected[[]], Il10_M2c_signature_UCell ~ Timepoint, p.adjust.method = "holm", detailed = FALSE)

#Compare slopes of M1 and M2 signatures
ggplot(macrophages_wt_infected[[]], aes (x = M1_signature_UCell, y = Il4_alt_signature_UCell))+
  geom_point(alpha = 0.8, aes(color = Timepoint))+
  geom_smooth(method = 'lm', se = FALSE, color = 'black')+
  xlab('M1 Signature')+
  ylab('M2 Signature')+
  xlim(c(0,0.38))+
  ylim(c(0,0.38))+
  scale_colour_manual(values =c("#292270","#166DF0", "#6DC3F8"))


