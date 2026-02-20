library(Seurat)
library(dplyr)

#Load in data
yang_data <- LoadSeuratRds("~/Documents/ÖverbyLab/yang_public_data/yang_rpca_integrated_obj.RDS")
DimPlot(yang_data,reduction = "umap.rpca", label = TRUE)+
  theme(legend.position = 'none')

#Label by treatment group
yang_sample <- substring(yang_data$id, 1, 1) 
yang_sample_labels <- case_when(yang_sample == 1 ~ 'healthy',
          yang_sample %in% c(2,3) ~'mild',
          yang_sample %in% c(4) ~ 'moderate',
          yang_sample %in% c(5,6) ~'severe')
yang_data$treatment <- yang_sample_labels

DimPlot(yang_data,reduction = "umap.rpca", group.by = 'treatment')
table(yang_data$manualAnnotation, yang_data$treatment)

#Look at genes of interest
FeaturePlot(yang_data, features = 'Nos2', reduction = 'umap.rpca')

FeaturePlot(yang_data, features = 'H2-Aa', reduction = 'umap.rpca')

DotPlot(yang_data, features = c('H2-Aa', 'H2-Ab1', 'H2-DMa'), group.by = 'treatment', scale = FALSE)
