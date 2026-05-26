library(Seurat)
library(ggplot2)
library(dplyr)
library(networkD3)
library(ggpubr)
library(SeuratExtend)
library(irGSEA)
library(clusterProfiler)
library(org.Mm.eg.db)
source('~/Documents/ÖverbyLab//scripts/langatFunctions.R')

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

#Get rid of lgtv and unknown samples
chimeric_mock <- subset(ParseSeuratObj_int, Treatment != 'rLGTV' & manualAnnotation != 'unknown')

astros <- subset(chimeric_mock, manualAnnotation == 'Astrocytes')

#Test gsea with day 3 wt vs mock
astros_day3 <- subset(astros, Timepoint == 'Day 3' | Treatment == 'PBS')
astros_day3_wt <- subset(astros, (Genotype == 'WT' & Timepoint == 'Day 3') | Treatment == 'PBS')

#GSEA
#Calculate enrichment scores. Tried with 4 cores which froze computer. 1 seems to work fine
irgsea_score_fun <- function(dat){
  return(irGSEA.score(object = dat, assay = "RNA", 
               slot = "data", seeds = 123, ncores = 1,
               min.cells = 3, min.feature = 0,
               custom = F, geneset = NULL, msigdb = T, 
               species = "Mus musculus", category = "H",  
               subcategory = NULL, geneid = "symbol",
               method = c("AUCell", "UCell", "singscore", 
                          "ssgsea", "JASMINE", "viper"),
               aucell.MaxRank = NULL, ucell.MaxRank = NULL, 
               kcdf = 'Gaussian'))
}

astros_day3 <- irgsea_score_fun(astros_day3)
astros_day3_wt <- irgsea_score_fun(astros_day3_wt)

#Integrate differential gene sets
irgsea_result_fun <- function(dat){
  return(irGSEA.integrate(object = dat, 
                                 group.by = "Treatment",
                                 metadata = NULL, col.name = NULL,
                                 method = c("AUCell","UCell","singscore",
                                            "ssgsea", "JASMINE", "viper")))
}

result_dge_3 <- irgsea_result_fun(astros_day3)
result_dge_wt_3 <- irgsea_result_fun(astros_day3_wt)

#Heatmap of sig pathways
irGSEA.heatmap(object = result_dge_wt_3, 
               method = "RRA",
               top = 50, 
               show.geneset = NULL)


irGSEA.heatmap(object = result_dge_3,
               method = "RRA",
               top = 50,
               show.geneset = NULL)

#UMAP astrocytes and do density scatter plots
prep_astro_umap <- function(astro_dat, num_dims = 20, returnElbow = FALSE, reduc_name = NULL){
  astro_dat <- prepSeuratObj(astro_dat, use_all_genes = FALSE)
  if(returnElbow){
    print(ElbowPlot(astro_dat, ndims = 40))
    return()
  }
  astro_dat <- prepUmapSeuratObj(astro_dat, nDims = num_dims, reductionName = reduc_name, resolution_value = 0.8)
  return(astro_dat)
}

astros_day3 <- prep_astro_umap(astros_day3, num_dims = 20, reduc_name = 'day_3_astro')
astros_day3_wt <- prep_astro_umap(astros_day3_wt, num_dims = 20, reduc_name = 'day_wt_3_astro')

#density scatter
irGSEA.density.scatterplot(object = astros_day3,
                           method = "UCell",
                           show.geneset = "HALLMARK-INTERFERON-ALPHA-RESPONSE",
                           reduction = "day_3_astro")

#Switch to clusterprofiler for multi enrichment plots
#Make ranked gene list for cluterprofiler
day3_wt_markers <- FindMarkers(astros_day3_wt, group.by = 'Treatment', ident.1 = 'rChLGTV', only.pos = TRUE, logfc.threshold = 0, 
                               min.pct = 0)
day3_wt_markers <- day3_wt_markers %>% as.data.frame() %>% dplyr::arrange(desc(avg_log2FC)) 
day3_wt_markers_list <- as.list(day3_wt_markers$avg_log2FC)
names(day3_wt_markers_list) <- rownames(day3_wt_markers)
#Probably a better way to do this but the above creates a list of lists, so get rid of that now
day3_wt_markers_list <- unlist(day3_wt_markers_list)

day3_wt_gsea<- gseGO(geneList     = day3_wt_markers_list,
              OrgDb        = org.Mm.eg.db,
              ont          = "BP",
              minGSSize    = 100,
              maxGSSize    = 600,
              pvalueCutoff = 0.05,
              verbose      = FALSE)
