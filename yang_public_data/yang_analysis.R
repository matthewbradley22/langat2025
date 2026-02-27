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
endo_trafficking_dot_dat <- DotPlot(yang_data, features = c('Vcam1', 'Icam1', 'Icam2', 'Sele', 'Selp', 'Mcam', 'F11r', 'Jam2', 'Pecam1', 'Pvr', 'Cd99l2',
                                                                 'Cdh5'), group.by = 'treatment', scale = FALSE)$data

