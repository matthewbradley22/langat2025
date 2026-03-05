library(Seurat)
library(dplyr)
library(UCell)

#Load in data
yang_data <- LoadSeuratRds("~/Documents/ÖverbyLab/yang_public_data/yang_rpca_integrated_obj.RDS")
newCols <-  c(brewer.pal(12, 'Paired'), '#99FFE6', '#CE99FF', '#18662E','#737272',  '#FF8AEF')
newCols[11] =  '#FF8AEF'

#Plot data
DimPlot(yang_data,reduction = "umap.rpca", label = FALSE, group.by = 'manualAnnotation', cols = newCols)
DimPlot(yang_data,reduction = "umap.rpca", label = FALSE, group.by = 'treatment', split.by = 'treatment' ,cols = newCols)

#Label by treatment group
yang_sample <- substring(yang_data$id, 1, 1) 
yang_sample_labels <- case_when(yang_sample == 1 ~ 'healthy',
          yang_sample %in% c(2,3) ~'mild',
          yang_sample %in% c(4) ~ 'moderate',
          yang_sample %in% c(5,6) ~'severe')
yang_data$treatment <- yang_sample_labels

DimPlot(yang_data,reduction = "umap.rpca", split.by = 'treatment', group.by = 'treatment')
table(yang_data$manualAnnotation, yang_data$treatment)

#Look at genes of interest
FeaturePlot(yang_data, features = 'Nos2', reduction = 'umap.rpca')
DotPlot(yang_data, features = c('Nos2'), group.by = 'manualAnnotation', scale = FALSE)

#Look at macrophage specific expression
macs <- subset(yang_data, manualAnnotation == 'Macro/Mono')

DotPlot(macs, features = c('Nos2'), group.by = 'treatment', scale = FALSE)+
  ggtitle('Macrophage Nos2 Expression')

#MHC-2 markers https://link.springer.com/article/10.1186/s12974-016-0581-z
DotPlot(macs, features = c('H2-Aa', 'H2-Ab1', 'H2-DMa', 'H2-DMb1', 'H2-DMb2', 'H2-Eb1'), group.by = 'treatment', scale = FALSE)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
  ggtitle('mhc-2 mac markers (Yang)')+
  coord_flip()+
  xlab('')+
  ylab('')+
  scale_color_gradientn(colours = c('white', '#FFD991', '#FF7530', '#FF4024'), 
                        values = c(0, 0.3, 0.6, 1),
                        name = 'Average Expression')

DotPlot(yang_data, features = c('Ccl2', 'Ccl5', 'Ccl7'), group.by = 'manualAnnotation', scale = FALSE)
FeaturePlot(yang_data, features = c('Ccl2'))

#Trafficking genes
yang_endo <- subset(yang_data, manualAnnotation == 'Endothelial')
endo_trafficking_dot_dat <- DotPlot(yang_endo, features = c('Vcam1', 'Icam1', 'Icam2', 'Sele', 'Selp', 'Mcam', 'F11r', 'Jam2', 'Pecam1', 'Pvr', 'Cd99l2',
                                                                 'Cdh5'), group.by = 'treatment', scale = FALSE)$data
ggplot(endo_trafficking_dot_dat, aes(x = id, y = features.plot))+
  geom_point(aes(size = pct.exp, fill = avg.exp.scaled), pch = 21)+
  #coord_flip()+
  scale_size_continuous(range = c(0.5,6), limits = c(0,100))+
  scale_fill_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white"), 
                       values = c(1.0,0.7,0.4,0),
                       limits = c(0,2.7))+
  ggtitle('Endothelial trafficking genes Yang data')+
  theme_bw()+
  theme(panel.grid = element_blank(),
        #panel.border = element_rect(colour = "white", fill = NA),
        panel.spacing = unit(0, "line"),
        strip.background = element_blank())+
  scale_size(range = c(0,9), limits = c(0, 100))+
  ylab('')+
  xlab('')

#Monocyte trafficking genes
#Still deciding which cells will be labelled monocytes, as some cluster closely with microglia
yang_mono <- subset(yang_data, manualAnnotation == 'Macro/Mono')
mono_trafficking_dot_dat <- DotPlot(yang_mono, features = c('Pecam1', 'Pvr', 'Cd99l2', 'Epha1', 'Ephb1', 'Itga4', 'Itgb1', 'Ccr1', 'Ccr2', 'Ccr5'), 
                                    group.by = 'treatment', scale = FALSE)$data

ggplot(mono_trafficking_dot_dat, aes(x = id, y = features.plot))+
  geom_point(aes(size = pct.exp, fill = avg.exp.scaled), pch = 21)+
  #coord_flip()+
  scale_size_continuous(range = c(0.5,6), limits = c(0,100))+
  scale_fill_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white"), 
                       values = c(1.0,0.7,0.4,0),
                       limits = c(0,1.7))+
  ggtitle('Monocyte trafficking genes Yang data')+
  theme_bw()+
  theme(panel.grid = element_blank(),
        #panel.border = element_rect(colour = "white", fill = NA),
        panel.spacing = unit(0, "line"),
        strip.background = element_blank())+
  scale_size(range = c(0,9), limits = c(0, 100))+
  ylab('')+
  xlab('')

#Monocyte/Macropahge polarization
macrophage_subset_markers <- list(M1_signature = c('Tnf', 'Il1b', 'Il6', 'Il12a', 'Il23a', 'Il27', 'Nos2', 'Ido1', 'Irf5', 'Nfkblz',
                                                   'Socs1', 'Marco', 'Cd80', 'Cd86', 'Cd274', 'C1d', 'C1qa', 'C1qb', 'C1qbp', 'C1qc',
                                                   'C1qtnf1', 'C3', 'C3ar1', 'C5ar1', 'Itgb2', 'Ccl2', 'Ccl3', 'Ccl4', 'Ccl5', 'Cxcl9',
                                                   'Cxcl10', 'Cxcl11', 'Cxcl16', 'Ccr7'),
                                  Il4_alt_signature = c('Arg1', 'Mrc1', 'Chil3', 'Tgm2', 'Il1rn', 'Msr1', 'Pdcd1lg2', 'Cd209a', 'Clec10a', 'Retnla', 'Alox15',
                                                        'Socs2', 'Irf4', 'Ppard', 'Pparg', 'Ccl17', 'Ccl24'),
                                  Il10_M2c_signature = c('Il10', 'Il4ra', 'Tgfb1', 'Nfil3', 'Sbno2', 'Socs3', 'Fcgr1', 'Fcgr2b', 'Fcgr3', 'Marco',
                                                         'Cxcl13'),
                                  mhc2_sig = c('H2-Aa', 'H2-Ab1', 'H2-DMa', 'H2-DMb1', 'H2-DMb2', 'H2-Eb1'))

median_genes <- round(median(yang_mono$nFeature_RNA))
yang_mono <- AddModuleScore_UCell(yang_mono, features = macrophage_subset_markers, maxRank=median_genes)

#M1 Signature
VlnPlot(yang_mono, features = 'M1_signature_UCell', group.by = 'treatment', pt.size = 0)+
  geom_boxplot(width = 0.1, position = position_dodge(0.9), alpha = 0.5,fill = 'white')+
  ggtitle('M1 signature Yang')+
  ylim(c(-0.01,0.65))

#M2 signature
VlnPlot(yang_mono, features = 'Il4_alt_signature_UCell', group.by = 'treatment', pt.size = 0)+
  geom_boxplot(width = 0.1, position = position_dodge(0.9), alpha = 0.5,fill = 'white')+
  ggtitle('M2a signature Yang')+
  ylim(c(-0.01,0.65))

#M2c signature
VlnPlot(yang_mono, features = 'Il10_M2c_signature_UCell', group.by = 'treatment', pt.size = 0)+
  geom_boxplot(width = 0.1, position = position_dodge(0.9), alpha = 0.5,fill = 'white')+
  ggtitle('M2c signature Yang')+
  ylim(c(-0.01,0.65))

#Mhc-II signature
VlnPlot(yang_mono, features = 'mhc2_sig_UCell', group.by = 'treatment', pt.size = 0)+
  geom_boxplot(width = 0.1, position = position_dodge(0.9), alpha = 0.5,fill = 'white')

#Chemokines
yang_data$treatment_celltype <- paste(yang_data$treatment, yang_data$manualAnnotation, sep = '_')
ccl_genes <- c('Ccl2', 'Ccl3', 'Ccl5', 'Ccl7')
ccl_exp <- DotPlot(yang_data, features = ccl_genes,
        group.by = 'treatment_celltype', scale = FALSE)$data
ccl_meta <- stringr::str_split_fixed(ccl_exp$id, "_", 2)
colnames(ccl_meta) = c('treatment', 'celltype')
ccl_dat <- cbind(ccl_exp, ccl_meta)

ggplot(ccl_dat, aes(x = celltype, y = features.plot))+
  geom_point(aes(size = pct.exp, fill = avg.exp.scaled), pch = 21)+
  facet_grid(cols = vars(treatment), scales = "free")+
  coord_flip()+
  scale_size_continuous(range = c(0.5,6), limits = c(0,100))+
  scale_fill_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white"), 
                       values = c(1.0,0.7,0.4,0),
                       limits = c(0,6))+
  ggtitle('Ccl genes Yang data')+
  theme_bw()+
  theme(panel.grid = element_blank(),
        #panel.border = element_rect(colour = "white", fill = NA),
        panel.spacing = unit(0, "line"),
        strip.background = element_blank())+
  scale_size(range = c(0,9), limits = c(0, 100))+
  ylab('')+
  xlab('')

#Plot ifna genes
yang_data <- JoinLayers(yang_data)
DotPlot(yang_data, features = c(paste0('Ifna', 1:16), 'Ifnb1', 'Ifng'), group.by = 'treatment', scale = FALSE)+
  coord_flip()+
  theme(axis.text.x = element_text(angle = 90))+
  xlab('')+
  ylab('')+
  scale_color_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white"), 
                       values = c(1.0,0.7,0.4,0),
                       limits = c(0,1.3))+
  ggtitle('Interferon expression Yang data')

table(yang_data$manualAnnotation, yang_data$treatment) %>% 
  as.data.frame() %>% 
  dplyr::group_by(Var2) %>% 
  dplyr::mutate(prop = Freq / sum(Freq)) %>% 
  ggplot(aes(x = Var2, y = prop, fill = Var1))+
  geom_bar(stat = 'identity')+
  scale_fill_manual(values=newCols)+
  xlab('')+
  ylab('')

