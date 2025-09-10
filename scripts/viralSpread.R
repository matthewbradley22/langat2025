#Script analyzing viral spread

#Load packages
library(Seurat)

#Load in data
ParseSeuratObj_int <- LoadSeuratRds("./data/FilteredRpcaIntegratedDat.rds")

#For now, manual labels assigned in LangatCellAnnotations.R and doublets removed in manualDoubletCheck.R
DimPlot(ParseSeuratObj_int, label = FALSE, group.by = 'manualAnnotation', 
        reduction = 'umap.integrated',
        cols = newCols)

ParseSeuratObj_int$hasVirus = ifelse(ParseSeuratObj_int$virusCountPAdj > 10, 1, 0)

#Where is virus localized
FeaturePlot(ParseSeuratObj_int, 'hasVirus',
            reduction = 'umap.integrated')+
  theme(legend.position = 'NONE')

#How many infected cells per cell type
ParseSeuratObj_int[[]] %>% subset(hasVirus == 1) %>% group_by(manualAnnotation) %>% 
   dplyr::summarise(infectedCells = n()) %>% arrange(desc(infectedCells))

ParseSeuratObj_int[[]] %>% subset(hasVirus == 1) %>% group_by(manualAnnotation, Timepoint) %>% 
  dplyr::summarise(infectedCells = n()) %>% arrange(desc(infectedCells))

#Look at proportion of infected cells per cell type
ParseSeuratObj_int[[]] %>%  group_by(manualAnnotation, hasVirus) %>% 
  dplyr::summarise(cellCount = n()) %>% 
  ggplot(aes(x = manualAnnotation, y = cellCount, fill = hasVirus))+
  geom_bar(stat = 'identity')+
  theme(axis.text.x = element_text(angle = 90))



#Create gene counts table with virus
#This is not perfect as virus counts have not been corrected for fragmentation 
geneCounts <- ParseSeuratObj_int[['RNA']]$counts
sparseVirus = sparseMatrix(i = which(ParseSeuratObj_int$virusCountPAdj > 0), j = rep(1, 40368), 
                           x = ParseSeuratObj_int$virusCountPAdj[which(ParseSeuratObj_int$virusCountPAdj > 0)], 
                           dims = c(119145, 1))
sparseVirus <- t(sparseVirus)
dimnames(sparseVirus)[[1]] = 'LGTV'
datWithVirus <- rbind(geneCounts, sparseVirus)

seuObjWithVirus <- CreateSeuratObject(datWithVirus, project = 'datWithVirus', assay = 'RNA')
seuObjWithVirus <- NormalizeData(seuObjWithVirus)
#These variable features are slightly different than the variable features in the integrated seurat
#Object, I think because I removed doublets after integration
seuObjWithVirus <- FindVariableFeatures(seuObjWithVirus)
seuObjWithVirus <- ScaleData(seuObjWithVirus)

twoVarDotPlot <- function(){
  exp_mat<- datWithVirusScaled['LGTV',]
  meta<- ParseSeuratObj_int[[]]
  as.matrix(exp_mat) %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column(var="cell") %>% 
    mutate('cell1' = substr(cell, start = 9, stop = 16)) %>% 
    mutate('sublib' = substr(cell, start = 1, stop = 7)) 
}

