#Script analyzing viral spread

#Load packages
library(Seurat)
library(ComplexHeatmap)

#Load in data
ParseSeuratObj_int <- LoadSeuratRds("./data/FilteredRpcaIntegratedDat.rds")

#For now, manual labels assigned in LangatCellAnnotations.R and doublets removed in manualDoubletCheck.R
DimPlot(ParseSeuratObj_int, label = FALSE, group.by = 'manualAnnotation', 
        reduction = 'umap.integrated',
        cols = newCols)

ParseSeuratObj_int$hasVirus = ifelse(ParseSeuratObj_int$virusCountPAdj >= 10, 1, 0)

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

createCountsWithVirus <- function(countDataObj, viralCounts){
  geneCounts <- countDataObj[['RNA']]$counts
  totalInfected <- sum(countDataObj$virusCountPAdj > 0)
  numSamples = dim(geneCounts)[2]
  sparseVirus = sparseMatrix(i = which(countDataObj$virusCountPAdj > 0), j = rep(1, totalInfected), 
                             x = countDataObj$virusCountPAdj[which(countDataObj$virusCountPAdj > 0)], 
                             dims = c(numSamples, 1))
  sparseVirus <- t(sparseVirus)
  dimnames(sparseVirus)[[1]] = 'LGTV'
  datWithVirus <- rbind(geneCounts, sparseVirus)
  
  seuObjWithVirus <- CreateSeuratObject(datWithVirus, project = 'datWithVirus', assay = 'RNA')
  seuObjWithVirus <- NormalizeData(seuObjWithVirus)
  #These variable features are slightly different than the variable features in the integrated seurat
  #Object, I think because I removed doublets after integration
  seuObjWithVirus <- FindVariableFeatures(seuObjWithVirus)
  seuObjWithVirus <- ScaleData(seuObjWithVirus)
  seuObjWithVirus
}

lgtvSamples <- subset(ParseSeuratObj_int, Treatment == 'rLGTV')


#### Function to prep data for two variable dot plot ####
#Taken from article here https://divingintogeneticsandgenomics.com/post/how-to-make-a-multi-group-dotplot-for-single-cell-rnaseq-data/

ParseSeuratObj_int$timeGenotype <- factor(paste(ParseSeuratObj_int$Timepoint, ParseSeuratObj_int$Genotype))

#Not even sure this should be a function, keep running through it line by line, we'll see
twoVarDotPlot <- function(datObj, metaDatObj, feature, virus = 'rLGTV'){
  exp_mat<- datObj[['RNA']]$data[feature,, drop = FALSE]
  count_mat<- datObj[['RNA']]$counts[feature,,drop=FALSE ]
  meta<- metaDatObj[[]] 
  
    
  #Creates a matrix of average expression across whichever groups are chosen in 
  #group_by
  exp_df <-  as.matrix(exp_mat) %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column(var="gene") %>%
    tidyr::pivot_longer(!gene, names_to = "cell", values_to = "expression") %>%
    mutate('cell1' = substr(cell, start = 9, stop = 16)) %>% 
    mutate('sublib' = substr(cell, start = 1, stop = 7)) %>% 
    left_join(meta, by = c('cell1' = 'cell', 'sublib' = 'subLib')) %>%
    group_by(gene, manualAnnotation, timeGenotype) %>% 
    dplyr::summarise(average_expression = mean(expression)) %>%
    tidyr::pivot_wider(names_from = manualAnnotation, 
                       values_from= average_expression) 
  
  #Remove character solumns (gene name and Genotype) from df to make a matrix, add back
  #Genotype as row name
  exp_mat<- exp_df[, -c(1, 2)] %>% as.matrix()
  rownames(exp_mat)<- exp_df %>% pull(timeGenotype)
  
  # get the percentage positive cell matrix
  count_df<- as.matrix(count_mat) %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column(var="gene") %>% 
    tidyr::pivot_longer(!gene, names_to = "cell", values_to = "count") %>% 
    mutate('cell1' = substr(cell, start = 9, stop = 16)) %>% 
    mutate('sublib' = substr(cell, start = 1, stop = 7)) %>% 
    left_join(meta, by = c('cell1' = 'cell', 'sublib' = 'subLib')) %>% 
    group_by(gene, manualAnnotation, timeGenotype) %>% 
    
    #Need to change anything about expression matrix to match that this is > 10
    dplyr::summarise(percentage = 100*mean(count >= 10)) %>% 
    tidyr::pivot_wider(names_from = manualAnnotation, 
                       values_from= percentage) 
  
  percent_mat<- count_df[, -c(1,2)] %>% as.matrix()
  rownames(percent_mat)<- count_df %>% pull(timeGenotype)
  
  identical(dim(exp_mat), dim(percent_mat))
  all.equal(colnames(exp_mat), colnames(percent_mat))
  all.equal(rownames(exp_mat), rownames(percent_mat))
  
  dat <- data.frame(timeGenotype = rep(rownames(percent_mat), each = ncol(percent_mat)), 
             cellType = rep(colnames(percent_mat),6),
             exp = as.vector(t(exp_mat)),
             percent = as.vector(t(percent_mat)))
  return(dat)
}

LGTV_dat <- twoVarDotPlot(seuObjWithVirus, ParseSeuratObj_int, feature = 'LGTV')

#Percent is % of cells with over 10 reads

#Need to split this by treatment, right now day 4 dominating because no PBS samples
ggplot(LGTV_dat, aes(x = timeGenotype, y = cellType, col = exp, size = percent))+
  geom_point()+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  scale_colour_gradient2 (low = "#2389A8" ,mid = "white", high = "#A82323", midpoint = 0.6,
                          name = "Average scaled expression")+
  theme(legend.title = element_text(size = 12, angle = 90))+
  guides( colour = guide_colourbar(theme = theme(
           legend.key.width  = unit(1, "lines"),
           legend.key.height = unit(10, "lines")),
           title.position = "left"))

ParseSeuratObj_int[[]] %>% mutate(virusHigh = ifelse(virusCountPAdj >= 10, 1, 0)) %>% 
  group_by(manualAnnotation, Genotype) %>% dplyr::summarise(virusHighProp = mean(virusHigh))

table(ParseSeuratObj_int$timeGenotype, ParseSeuratObj_int$Treatment)

#Which virus spreads more
#caveat here is that both viruses have higher expression in cerebrum
table(ParseSeuratObj_int$Treatment, ParseSeuratObj_int$Genotype)
ParseSeuratObj_int[[]] %>% filter(virusCountPAdj >= 10) %>% group_by(Treatment, Genotype) %>% 
  dplyr::summarise(highVirus = n())

#Make a model to see which variables have largest effect on viral load
#start linear
linearModel <- lm(virusCountPAdj~ Genotype + Treatment + Timepoint + Organ, #+ manualAnnotation,
                  data = ParseSeuratObj_int[[]])

summary(linearModel)



