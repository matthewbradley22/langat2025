library(Seurat)
library(dplyr)

#Load in data
yang_data <- LoadSeuratRds("~/Documents/ÖverbyLab/yang_public_data/yang_rpca_integrated_obj.RDS")
newCols <-  c(brewer.pal(12, 'Paired'), '#99FFE6', '#CE99FF', '#18662E','#737272',  '#FF8AEF')
newCols[11] =  '#FF8AEF'
DimPlot(yang_data,reduction = "umap.rpca", label = FALSE, group.by = 'manualAnnotation', cols = newCols)

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

