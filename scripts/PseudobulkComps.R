library(Seurat)
library(ggrepel)
library(pheatmap)
library(VennDiagram)

#Source functions
source('~/Documents/ÖverbyLab/scripts/langatFunctions.R')

#Load data
ParseSeuratObj_int <- LoadSeuratRds("~/Documents/ÖverbyLab/data/FilteredRpcaIntegratedDatNoDoublets.rds") 
ParseSeuratObj_int$hasVirus = ifelse(ParseSeuratObj_int$virusCountPAdj >= 10, 1, 0)

#Check data
newCols <-  c(brewer.pal(12, 'Paired'), '#99FFE6', '#CE99FF', '#18662E','#737272',  '#FF8AEF')
newCols[11] =  '#FF8AEF'
DimPlot(ParseSeuratObj_int, label = FALSE, group.by = 'manualAnnotation', reduction = 'umap.integrated',
        cols = newCols)
ParseSeuratObj_int$time_treatment <- paste(ParseSeuratObj_int$Timepoint, ParseSeuratObj_int$Treatment, sep = '_')
ParseSeuratObj_int$treatment_celltype <- paste(ParseSeuratObj_int$Treatment, ParseSeuratObj_int$manualAnnotation, sep = '_')

#Subset to just chimeric and resident cells 
resident_celltypes <- c('Astrocytes', 'Oligodendrocytes', 'Microglia', 'Endothelial', 'Choroid Plexus',
                        'Immature Neurons', 'Ependymal', 'Pericytes', 'Muscle cells', 'Neurons')


resident_no_lgtv <- subset(ParseSeuratObj_int, Treatment != 'rLGTV' & manualAnnotation %in% resident_celltypes)
table(resident_no_lgtv$Treatment)
table(resident_no_lgtv$manualAnnotation)

#Pseudobulk comparison between treatments and of timepoints
resident_no_lgtv_bulk <- createPseudoBulk(resident_no_lgtv, variables = c('Genotype', 'Treatment', 'Timepoint', 'Organ'))
resident_no_lgtv_bulk <- DESeq(resident_no_lgtv_bulk)
resultsNames(resident_no_lgtv_bulk)

#Day 4 vs day 3

#Not many upregulated at day 4
#grpofiler path function defined at bottom of script
day4_vs_day3_up_paths <- get_gprofiler_paths(resident_no_lgtv_bulk, compName = 'Timepoint_Day.4_vs_Day.3', direction = 'up')

#Lots of synapse and ion pathways showing up
day4_vs_day3_down_paths <- get_gprofiler_paths(resident_no_lgtv_bulk, compName = 'Timepoint_Day.4_vs_Day.3', direction = 'down')
day4_vs_day3_down_paths$result[day4_vs_day3_down_paths$result$source == 'KEGG',]
day4_vs_day3_down_paths$result[day4_vs_day3_down_paths$result$source == 'GO:MF',]
day4_vs_day3_down_paths$result[day4_vs_day3_down_paths$result$source == 'GO:BP',]

#Day 5 vs day 3
#Cytokine and tnf signalling showing upregulation now
day5_vs_day3_up_paths <- get_gprofiler_paths(resident_no_lgtv_bulk, compName = 'Timepoint_Day.5_vs_Day.3', direction = 'up')
day5_vs_day3_up_paths$result[day5_vs_day3_up_paths$result$source == 'KEGG',]
day5_vs_day3_up_paths$result[day5_vs_day3_up_paths$result$source == 'GO:MF',]
day5_vs_day3_up_paths$result[day5_vs_day3_up_paths$result$source == 'GO:BP',]

day5_vs_day3_down_paths <- get_gprofiler_paths(resident_no_lgtv_bulk, compName = 'Timepoint_Day.5_vs_Day.3', direction = 'down')
day5_vs_day3_down_paths$result[day5_vs_day3_down_paths$result$source == 'KEGG',]
day5_vs_day3_down_paths$result[day5_vs_day3_down_paths$result$source == 'GO:MF',]
day5_vs_day3_down_paths$result[day5_vs_day3_down_paths$result$source == 'GO:BP',]

##########GAL 3 Project Samples########
#######################################
wt_cerebrum_day5 <- subset(ParseSeuratObj_int, Treatment %in% c('PBS', 'rLGTV') & Organ == 'Cerebrum' & Genotype == 'WT' & Timepoint == 'Day 5')
wt_cerebrum_day5_bulk <- createPseudoBulk(wt_cerebrum_day5, variables = c('Treatment', 'manualAnnotation'))
wt_cerebrum_day5_bulk <- DESeq(wt_cerebrum_day5_bulk)
resultsNames(wt_cerebrum_day5_bulk)
results(wt_cerebrum_day5_bulk, name = "Treatment_rLGTV_vs_PBS")
#Upregulated LGTV
wt_cerebrum_day5_paths_up <- get_gprofiler_paths(wt_cerebrum_day5_bulk, compName = "Treatment_rLGTV_vs_PBS", direction = 'up')
wt_cerebrum_day5_paths_up$result[wt_cerebrum_day5_paths_up$result$source == 'KEGG',]$term_name
wt_cerebrum_day5_paths_up$result[wt_cerebrum_day5_paths_up$result$source == 'GO:MF',]$term_name
wt_cerebrum_day5_paths_up$result[wt_cerebrum_day5_paths_up$result$source == 'GO:BP',]$term_name

#Select pathways to plot
wt_cerebrum_day5_paths_up$result[wt_cerebrum_day5_paths_up$result$term_name %in% 
                                   c('apoptotic process', 'TNF signaling pathway', 'Apoptosis', 'Chemokine signaling pathway',
                                     'Cytokine-cytokine receptor interaction', 'Necroptosis',
                                     'Efferocytosis', 'CCR chemokine receptor binding'),] %>% 
  ggplot(aes(x = -log10(p_value), y = term_name, fill = source))+
  geom_bar(stat = 'identity', position = 'dodge')+
  ggtitle('WT Cerebrum Day 5 Upregulated in LGTV')

#Compare genes across similar paths to see if they are similar or different
#returns one long string of genes so need to split by comma
apoptotic_process <- unlist(strsplit(wt_cerebrum_day5_paths_up$result[wt_cerebrum_day5_paths_up$result$term_name %in% 
                                                   c('apoptotic process'),]$intersection, ","))
apoptosis <- unlist(strsplit(wt_cerebrum_day5_paths_up$result[wt_cerebrum_day5_paths_up$result$term_name %in% 
                                                                       c('Apoptosis'),]$intersection, ","))
#pro apotosis list from anna email
pro_apoptosis <- c('Apaf1', 'Casp9', 'Casp8', 'Bax', 'Bak',
                   'Bid', 'Bad', 'Bim', 'Bcl10', 'Bik',
                   'Blk', 'Fas', 'Fasl', 'Tnfrsf1a', 'Tnf', 'Tyro3',
                   'Axl', 'Mertk', 'Tnfsf10', 'Tnfrsf10b', 'Casp3', 'Casp6', 'Casp7')
venn.diagram(
  x = list(apoptotic_process, apoptosis, pro_apoptosis),
  category.names = c("apoptotic_process" , "apoptosis" , "pro_apoptosis"),
  filename = '~/Documents/ÖverbyLab/scPlots/apoptosis_venn_diagram.png',
  output=TRUE,
  fill = c('red', 'blue', 'yellow'),
  cat.dist = c(0.055, 0.055, 0.055),
  cat.cex = 1.2
)

#Same thing as above apoptosis genes but with necroptosis genes
necroptosis_kegg <- unlist(strsplit(wt_cerebrum_day5_paths_up$result[wt_cerebrum_day5_paths_up$result$term_name %in% 
                                                                        c('Necroptosis'),]$intersection, ","))
necroptosis <- c('Tnf', 'Tnfrsf1a', 'Ripk2', 'Mlkl', 'Ripk1', 'Ripk3')


draw.pairwise.venn(length(necroptosis_kegg),  length(necroptosis), length(intersect(necroptosis_kegg, necroptosis)),
                   fill = c('blue', 'yellow'))


#Downregulated LGTV
wt_cerebrum_day5_paths_down <- get_gprofiler_paths(wt_cerebrum_day5_bulk, compName = "Treatment_rLGTV_vs_PBS", direction = 'down')
wt_cerebrum_day5_paths_down$result[wt_cerebrum_day5_paths_down$result$source == 'KEGG',]$term_name
wt_cerebrum_day5_paths_down$result[wt_cerebrum_day5_paths_down$result$source == 'GO:MF',]$term_name
wt_cerebrum_day5_paths_down$result[wt_cerebrum_day5_paths_down$result$source == 'GO:BP',]$term_name


#Write functions to
#Run fgsea
#Load msgigdbr data
all_gene_sets <- msigdbr(species = "Mus musculus", db_species = "MM")
mouse_gene_sets <- all_gene_sets %>%
  split(x = .$gene_symbol, f = .$gs_name)

run_fgsea <- function(dat, compName){
  dat_comp <- results(dat, name=compName)
  ranks <- dat_comp %>% as.data.frame() %>% dplyr::arrange(desc(log2FoldChange)) %>% rownames()
  fgseaRes <- fgsea(pathways = mouse_gene_sets, 
                    stats    = ranks,
                    minSize  = 15,
                    maxSize  = 500)
}

#Run gprofiler2
get_gprofiler_paths <- function(dat, compName, direction){
  dat_comp <- results(dat, name=compName)
  if(direction == 'up'){
    dat_comp <- dat_comp %>% as.data.frame() %>% dplyr::filter(padj < 0.05 & log2FoldChange > 1) %>% 
      dplyr::arrange(padj)
  }
  if(direction == 'down'){
    dat_comp <- dat_comp %>% as.data.frame() %>% dplyr::filter(padj < 0.05, log2FoldChange < -1) %>% 
      dplyr::arrange(padj)
  }
  if(direction == 'all'){
    dat_comp <- dat_comp %>% as.data.frame() %>% dplyr::filter(padj < 0.05 & abs(log2FoldChange) > 1) %>% 
      dplyr::arrange(padj)
  }
  dat_comp_paths <- gprofiler2::gost(query = rownames(dat_comp), organism = 'mmusculus', evcodes = TRUE)
  dat_comp_paths
}




