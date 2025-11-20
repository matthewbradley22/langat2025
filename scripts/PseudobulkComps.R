library(Seurat)
library(ggrepel)
library(pheatmap)

#Source functions
source('./scripts/langatFunctions.R')

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
day4_vs_day3_up_paths <- get_gprofiler_paths(resident_no_lgtv_bulk, compName = 'Timepoint_Day.4_vs_Day.3', direction = 'up')

#Lots of synapse and ion pathways showing up
day4_vs_day3_down_paths <- get_gprofiler_paths(resident_no_lgtv_bulk, compName = 'Timepoint_Day.4_vs_Day.3', direction = 'down')
day4_vs_day3_down_paths$result[day4_vs_day3_paths$result$source == 'KEGG',]
day4_vs_day3_down_paths$result[day4_vs_day3_paths$result$source == 'GO:MF',]
day4_vs_day3_down_paths$result[day4_vs_day3_paths$result$source == 'GO:BP',]

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




