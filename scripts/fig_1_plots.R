library(Seurat)
library(ggplot2)
library(dplyr)
library(networkD3)

#Load data
ParseSeuratObj_int <- LoadSeuratRds("~/Documents/ÖverbyLab/data/FilteredRpcaIntegratedDatNoDoublets.rds") 
ParseSeuratObj_int$hasVirus = ifelse(ParseSeuratObj_int$virusCountPAdj >= 10, 1, 0)

#Check data
umap_color_list <- c( "#7047A1", "#B370AE","#292270",  "#166DF0","#6D92F8",  "#6DC3F8", "#8a0000","#F76363", "#FF96A2", "#D6644B", 
                      "#F08C3A", "#fdc087","#074F00", "#208d1f","#7bcd79", 
                        "gray")

DimPlot(ParseSeuratObj_int, label = FALSE, group.by = 'manualAnnotation', reduction = 'umap.integrated',
        cols = umap_color_list)+
  theme(axis.ticks = element_blank(),
        axis.text=element_blank(),
        legend.text=element_text(size=17))+
  xlab('Umap1')+
  ylab('Umap2')+
  guides(color=guide_legend(override.aes=list(size=8)))+
  ggtitle('')

chimeric_mock <- subset(ParseSeuratObj_int, Treatment != 'rLGTV' & manualAnnotation != 'unknown')

chimeric_mock$manualAnnotation <- factor(chimeric_mock$manualAnnotation, 
                                         levels = rev(c( 'unknown',  'T cells',   'Nk cells', 'Macrophage/Monocytes', 
                                                     'Granulocytes', 'B Cells',  'Pericytes', 'Oligodendrocytes','Neurons',
                                                     'Muscle cells', 'Microglia', 'Immature Neurons', 'Ependymal','Endothelial', 
                                                     'Choroid Plexus', 'Astrocytes')))

png("~/Documents/ÖverbyLab/single_cell_ISG_figures/fig_1_plots/chimeric_mock_umap.png", width = 900, height = 700)
DimPlot(chimeric_mock, label = FALSE, group.by = 'manualAnnotation', reduction = 'umap.integrated',
        cols = umap_color_list)+
  theme(axis.ticks = element_blank(),
        axis.text=element_blank(),
        legend.text=element_text(size=24))+
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
table(chLgtv_ips$Timepoint, chLgtv_ips$manualAnnotation) %>% 
  as.data.frame() %>% dplyr::group_by() %>% dplyr::mutate(freq_props = Freq/sum(Freq))%>% 
  ggplot(aes(x = Var1, y = freq_props, fill = Var2))+
  geom_bar(stat = 'identity', position = 'stack', width = 0.6)+
  scale_fill_manual(values = umap_color_list)+
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
table(chLgtv_wt$Timepoint, chLgtv_wt$manualAnnotation) %>% 
  as.data.frame() %>% dplyr::group_by(Var1) %>% dplyr::mutate(freq_props = Freq/sum(Freq))%>% 
  ggplot(aes(x = Var1, y = freq_props, fill = Var2))+
  geom_bar(stat = 'identity', position = 'stack', width = 0.6)+
  scale_fill_manual(values = umap_color_list)+
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
  scale_fill_manual(values = umap_color_list)+
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
  scale_fill_manual(values = umap_color_list)+
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
  scale_fill_manual(values = c(umap_color_list[4], umap_color_list[7]))+
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
  scale_fill_manual(values = c(umap_color_list[4], umap_color_list[7]))+
  theme_classic()+
  theme(axis.text.x = element_text(size = 22),
        axis.text.y = element_text(size = 22),
        axis.title=element_text(size=22),
        plot.title = element_text(size = 22))+
  ggtitle('Cerebellum')+
  ylim(c(0, 26000))
dev.off()


chimeric_mock$manualAnnotation <- factor(chimeric_mock$manualAnnotation, 
                                         levels = rev(levels(chimeric_mock$manualAnnotation)))

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
  scale_color_gradient2(low = '#9CCCFF', mid = 'white', high = '#9f0000', labels = c(-1, 0, 1,2),
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

#Mavs expression
chimeric_mock$geno_treatment_time_celltype <- paste(chimeric_mock$Genotype, chimeric_mock$Treatment, 
                                                    chimeric_mock$Timepoint, chimeric_mock$manualAnnotation, sep = '_')

mavs_dat <- DotPlot(object = chimeric_mock, features = c("Mavs"), group.by = 'geno_treatment_time_celltype', scale = FALSE)$data

mavs_dat_meta <- str_split_fixed(mavs_dat$id, "_", 4)
colnames(mavs_dat_meta) <- c('genotype', 'treatment', 'time', 'celltype')
mavs_dat <- cbind(mavs_dat, mavs_dat_meta)

#Could split by time
pdf('~/Documents/ÖverbyLab/single_cell_ISG_figures/fig_1_plots/mavs_expression.pdf', height = 9, width = 7)
ggplot(mavs_dat, aes(x = time, y = genotype, color = avg.exp.scaled, size = pct.exp))+
  facet_wrap(~treatment)+
  geom_point()+
  scale_color_gradient2(low = 'white', mid = 'orange', high = 'red', midpoint = 0.06)+
  theme_classic()+
  theme(axis.text = element_text(size = 24),
        legend.text = element_text(size = 24),
        legend.title = element_text(size = 24),
        plot.title = element_text(size =30),
        axis.text.x = element_text(angle = 90),
        strip.text = element_text(size = 24))+
  xlab('')+
  ylab('')+
  ggtitle('Mavs Expression')+
  scale_size_continuous(range = c(3,9))
dev.off()

#Split by celltype as well, not splitting by time for the sake of plotting
chimeric_mock$geno_treatment_celltype <- paste(chimeric_mock$Genotype, chimeric_mock$Treatment, 
                                                    chimeric_mock$manualAnnotation, sep = '_')
mavs_dat <- DotPlot(object = chimeric_mock, features = c("Mavs"), group.by = 'geno_treatment_celltype', scale = FALSE)$data
mavs_dat_meta <- str_split_fixed(mavs_dat$id, "_", 3)
colnames(mavs_dat_meta) <- c('genotype', 'treatment', 'celltype')
mavs_dat <- cbind(mavs_dat, mavs_dat_meta)

#Remove groups with low cell counts
table(chimeric_mock$Genotype, chimeric_mock$manualAnnotation, chimeric_mock$Treatment)
mavs_dat <- mavs_dat[!mavs_dat$celltype %in% c('Granulocytes', 'Nk cells') | !mavs_dat$treatment == 'PBS',]

pdf('~/Documents/ÖverbyLab/single_cell_ISG_figures/fig_1_plots/mavs_expression_celltype.pdf', height = 13, width = 13)
ggplot(mavs_dat, aes(x = genotype, y = avg.exp.scaled, color = celltype))+
  facet_wrap(~treatment)+
  geom_point(aes(size = pct.exp))+
  geom_line(aes(group = celltype))+
  theme_classic()+
  theme(axis.text = element_text(size = 24),
        legend.text = element_text(size = 24),
        legend.title = element_text(size = 24),
        plot.title = element_text(size =30),
        axis.text.x = element_text(angle = 90),
        strip.text = element_text(size = 24))+
  xlab('')+
  ylab('')+
  ggtitle('Mavs Expression')+
  scale_size_continuous(range = c(3,9))
dev.off()

#Plot of neo resistance gene reads 
#Data created on hpc2n at /home/m/mahogny/mystore/dataset/250729_scflavi/neo_read_counts from bam files
#Also reading in original barcoding data
final_neo_count <- read_table("~/Documents/ÖverbyLab/data/final_neo_count.txt",  col_names = FALSE)
colnames(final_neo_count) = c('total_reads', 'orig_id')
final_neo_count$orig_id = as.numeric(final_neo_count$orig_id)

ParseBarcodePlate <- read_csv("~/Documents/ÖverbyLab/ParseBarcodePlate.csv")
ParseBarcodePlate$orig_id = 1:48

neo_data_to_plot <- left_join(ParseBarcodePlate, final_neo_count, by = 'orig_id')
neo_data_to_plot$total_reads[is.na(neo_data_to_plot$total_reads)] <- 0
neo_data_to_plot_chimeric <- dplyr::filter(neo_data_to_plot, Treatment != 'rLGTV')

pdf('~/Documents/ÖverbyLab/single_cell_ISG_figures/fig_1_plots/neomycin_levels.pdf', height = 6, width = 8)
ggplot(neo_data_to_plot_chimeric, aes(x = Genotype, y = total_reads, fill = Treatment))+
  geom_boxplot()+
  ggtitle('Neomycin-resistance gene')+
  ylab('Count')+
  theme_classic()+
  theme(text = element_text(size = 24))+
  scale_fill_manual(values = c(umap_color_list[4], umap_color_list[10]))
dev.off()

#Infiltrating cell counts over time, normalized to day 3
wt_chimeric_infiltrating <- subset(chimeric_mock, manualAnnotation %in% c("T cells", "Macrophage/Monocytes", "Nk cells",
                                                                        "Granulocytes", "B Cells"))
levels(wt_chimeric_infiltrating$manualAnnotation) <- c(levels(wt_chimeric_infiltrating$manualAnnotation), "Macro")

wt_chimeric_infiltrating$manualAnnotation[wt_chimeric_infiltrating$manualAnnotation == 'Macrophage/Monocytes'] = 'Macro'

baseline_values = wt_chimeric_infiltrating[[]] %>% dplyr::group_by(Timepoint, manualAnnotation, Genotype) %>% 
  dplyr::summarise(baseline = n()) %>% dplyr::filter(Timepoint == 'Day 3') %>% 
  dplyr::ungroup() %>% 
  dplyr::select(manualAnnotation, Genotype, baseline)
  
pdf('~/Documents/ÖverbyLab/single_cell_ISG_figures/fig_1_plots/infiltrating_counts.pdf', height = 6, width = 11)
wt_chimeric_infiltrating[[]] %>% dplyr::group_by(Timepoint, manualAnnotation, Genotype) %>% 
  dplyr::summarise(cell_count = n()) %>% 
  dplyr::left_join(baseline_values, by = c('manualAnnotation', 'Genotype')) %>% 
  dplyr::mutate(change = cell_count - baseline,
                pct.change = change / baseline,
                Genotype = factor(Genotype, levels = c('WT', 'IPS1'))) %>% 
  dplyr::filter(Timepoint != 'Day 3') %>% 
  dplyr::mutate(manualAnnotation = factor(manualAnnotation, levels = c('B Cells', 'Granulocytes', 
                                                                       'Macro', 'Nk cells',
                                                                       'T cells'))) %>% 
  ggplot(aes(x = Genotype, y = pct.change, fill = Timepoint))+
  geom_bar(stat = 'identity', position = 'dodge')+
  facet_wrap(~manualAnnotation, nrow = 1) +
  theme_classic()+
  theme(text = element_text(size = 24))+
  ylab('Percent change from day 3')+
  xlab('')+
  scale_fill_manual(values = c(umap_color_list[4], umap_color_list[10]))

dev.off()


