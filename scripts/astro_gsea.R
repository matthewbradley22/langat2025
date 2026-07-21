library(Seurat)
library(ggplot2)
library(dplyr)
library(networkD3)
library(ggpubr)
library(SeuratExtend)
library(irGSEA)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(DOSE)
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
astros_day4 <- subset(astros, Timepoint == 'Day 4' | Treatment == 'PBS')
astros_day5 <- subset(astros, Timepoint == 'Day 5' | Treatment == 'PBS')

astros_day3_wt <- subset(astros, (Timepoint == 'Day 3' & Treatment == 'rChLGTV' & Genotype == 'WT') | Treatment == 'PBS')
astros_day5_wt <- subset(astros, (Timepoint == 'Day 5' & Treatment == 'rChLGTV' & Genotype == 'WT')  | Treatment == 'PBS')
astros_day3_ips <- subset(astros, (Timepoint == 'Day 3' & Treatment == 'rChLGTV' & Genotype == 'IPS1')  | Treatment == 'PBS')
astros_day5_ips <- subset(astros, (Timepoint == 'Day 5' & Treatment == 'rChLGTV' & Genotype == 'IPS1')  | Treatment == 'PBS')

#GSEA
#Only run on specific go lists
defense_response_to_virus <- read.table("~/Documents/ÖverbyLab/gene_lists_for_gsea/defense_response_to_virus.txt", quote="\"", comment.char="")
pos_reg_type_1_inf_production <- read.table("~/Documents/ÖverbyLab/gene_lists_for_gsea/pos_reg_type_1_int_production.txt", quote="\"", comment.char="")
type_1_inf_production <- read.table("~/Documents/ÖverbyLab/gene_lists_for_gsea/type_1_inf_production.txt", quote="\"", comment.char="")

defense_response_to_virus <- defense_response_to_virus$V1
pos_reg_type_1_inf_production <- pos_reg_type_1_inf_production$V1
type_1_inf_production <- type_1_inf_production$V1

custom_lists = list('defense_response_to_virus' = defense_response_to_virus, 
                    'type_1_inf_production' = type_1_inf_production)

#Calculate enrichment scores. Tried with 4 cores which froze computer. 1 seems to work fine
irgsea_score_fun <- function(dat){
  return(irGSEA.score(object = dat, assay = "RNA", 
               slot = "data", seeds = 123, ncores = 1,
               min.cells = 3, min.feature = 0,
               custom = TRUE, geneset = custom_lists, msigdb = T, 
               species = "Mus musculus", category = "C5",  
               subcategory = 'GO:BP', geneid = "symbol",
               method = c("UCell", "singscore", 
                          "ssgsea", "JASMINE", "viper"), #,"AUCell", ),
               aucell.MaxRank = NULL, ucell.MaxRank = NULL, 
               kcdf = 'Gaussian'))
}

astros_day3 <- irgsea_score_fun(astros_day3)
astros_day4 <- irgsea_score_fun(astros_day4)
astros_day5 <- irgsea_score_fun(astros_day5)

#Integrate differential gene sets
irgsea_result_fun <- function(dat){
  return(irGSEA.integrate(object = dat, 
                                 group.by = "Treatment",
                                 metadata = NULL, col.name = NULL,
                                 method = c("UCell","singscore",
                                            "ssgsea", "JASMINE", "viper")))
}

result_dge_3 <- irgsea_result_fun(astros_day3)
result_dge_4 <- irgsea_result_fun(astros_day4)
result_dge_5 <- irgsea_result_fun(astros_day5)

#Heatmap of sig pathways
irGSEA.heatmap(object = result_dge_3,
               method = "RRA",
               top = 50,
               show.geneset = NULL)

irGSEA.heatmap(object = result_dge_4, 
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
astros_day4 <- prep_astro_umap(astros_day4, num_dims = 20, reduc_name = 'day_4_astro')
astros_day5 <- prep_astro_umap(astros_day5, num_dims = 20, reduc_name = 'day_5_astro')

#density scatter
pdf('~/Documents/ÖverbyLab/single_cell_ISG_figures/astrocytes_fig/density_feature_plots/ucell_day3_defenseVirus.pdf', height = 6, width = 6)
irGSEA.density.scatterplot(object = astros_day3,
                           method = "UCell",
                           show.geneset = "defense-response-to-virus",
                           reduction = "day_3_astro")+ 
  scale_color_viridis_c(limits = c(0,0.0011))
  
dev.off()

pdf('~/Documents/ÖverbyLab/single_cell_ISG_figures/astrocytes_fig/density_feature_plots/ucell_day4_defenseVirus.pdf', height = 6, width = 6)
irGSEA.density.scatterplot(object = astros_day4,
                           method = "UCell",
                           show.geneset = "defense-response-to-virus",
                           reduction = "day_4_astro")+ 
  scale_color_viridis_c(limits = c(0,0.0011))

dev.off()

pdf('~/Documents/ÖverbyLab/single_cell_ISG_figures/astrocytes_fig/density_feature_plots/ucell_day5_defenseVirus.pdf', height = 6, width = 6)
irGSEA.density.scatterplot(object = astros_day5,
                           method = "UCell",
                           show.geneset = "defense-response-to-virus",
                           reduction = "day_5_astro")+ 
  scale_color_viridis_c(limits = c(0,0.0011))
dev.off()

pdf('~/Documents/ÖverbyLab/single_cell_ISG_figures/astrocytes_fig/density_feature_plots/ucell_day3_type1Inf.pdf', height = 6, width = 6)
irGSEA.density.scatterplot(object = astros_day3,
                           method = "UCell",
                           show.geneset = "type-1-inf-production",
                           reduction = "day_3_astro")+
  scale_color_viridis_c(limits = c(0,0.0006))
dev.off()

pdf('~/Documents/ÖverbyLab/single_cell_ISG_figures/astrocytes_fig/density_feature_plots/ucell_day4_type1Inf.pdf', height = 6, width = 6)
irGSEA.density.scatterplot(object = astros_day4,
                           method = "UCell",
                           show.geneset = "type-1-inf-production",
                           reduction = "day_4_astro")+
  scale_color_viridis_c(limits = c(0,0.0006))
dev.off()

pdf('~/Documents/ÖverbyLab/single_cell_ISG_figures/astrocytes_fig/density_feature_plots/ucell_day5_type1Inf.pdf', height = 6, width = 6)
irGSEA.density.scatterplot(object = astros_day5,
                           method = "UCell",
                           show.geneset = "type-1-inf-production",
                           reduction = "day_5_astro")+
  scale_color_viridis_c(limits = c(0,0.0006))
dev.off()

# - - - - - - - - - - - - - - - - - - 
#### GSEA erichment line plots ####
# - - - - - - - - - - - - - - - - - - 

#Switch to clusterprofiler for multi enrichment plots
#Make ranked gene list for cluterprofiler
create_gsea_object <- function(dat, min.p.val = 0.05){
  dat_markers <- FindMarkers(dat, group.by = 'Treatment', ident.1 = 'rChLGTV', logfc.threshold = 0, 
                                 min.pct = 0)
  dat_markers <- dat_markers %>% as.data.frame() %>% dplyr::arrange(desc(avg_log2FC)) 
  dat_markers_list <- as.list(dat_markers$avg_log2FC)
  names(dat_markers_list) <- rownames(dat_markers)
  #Probably a better way to do this but the above creates a list of lists, so get rid of that now
  dat_markers_list <- unlist(dat_markers_list)
  
  dat_gsea<- gseGO(geneList = dat_markers_list,
                       OrgDb = org.Mm.eg.db,
                       ont   = "BP",
                       minGSSize = 100,
                       maxGSSize = 600,
                       pvalueCutoff = min.p.val,
                       verbose = FALSE, 
                       keyType = 'SYMBOL')
}

#Visualization options for gsea results
day3_wt_gsea <- create_gsea_object(astros_day3_wt)
day3_ips_gsea <- create_gsea_object(astros_day3_ips, min.p.val = 1)
day5_wt_gsea <- create_gsea_object(astros_day5_wt)
day5_ips_gsea <- create_gsea_object(astros_day5_ips)


if(FALSE){
  create_emap <- function(gsea_dat){
    print(emapplot(pairwise_termsim(gsea_dat), showCategory = 7, node_label="none" )+
            theme_classic() +
            scale_color_gradientn(colours = c("#F03C0C","#F57456","#FFB975","#FFD1A6"), 
                                  values = c(0,0.4,0.7,1.0),
                                  limits = c(1, 2*10^-8))+
            theme(axis.text = element_blank(),
                  axis.ticks = element_blank())+
            ylab('')+
            xlab(''))+
      ggrepel::geom_text_repel(aes(x = x, y = y, label = name), size = 5)
  }
  
  pdf('~/Documents/ÖverbyLab/single_cell_ISG_figures/astrocytes_fig/wt_day3_emap.pdf', width = 7, height = 6)
  create_emap(day3_wt_gsea) 
  dev.off()
  
  pdf('~/Documents/ÖverbyLab/single_cell_ISG_figures/astrocytes_fig/ips_day3_emap.pdf', width = 7, height = 6)
  create_emap(day3_ips_gsea)
  dev.off()
  
  pdf('~/Documents/ÖverbyLab/single_cell_ISG_figures/astrocytes_fig/wt_day5_emap.pdf', width = 7, height = 6)
  create_emap(day5_wt_gsea)
  dev.off()
  
  pdf('~/Documents/ÖverbyLab/single_cell_ISG_figures/astrocytes_fig/ips_day5_emap.pdf', width = 7, height = 6)
  create_emap(day5_ips_gsea)
  dev.off()
  
}

#Dot plots
pdf('~/Documents/ÖverbyLab/single_cell_ISG_figures/astrocytes_fig/wt_day3_dotplot.pdf', width = 5.5, height = 6)
enrichplot::dotplot(day3_wt_gsea, showCategory=7, x = 'NES')+
  theme_classic() +
  scale_fill_gradientn(colours = c("#F03C0C","#FFD1A6","#FFEEE0"), 
                       values = c(0,0.5,1.0),
                       limits = c(1, 2*10^-8))+
  scale_size(range = c(0, 9), limits = c(0, 200))
dev.off()

pdf('~/Documents/ÖverbyLab/single_cell_ISG_figures/astrocytes_fig/ips_day3_dotplot.pdf', width = 5.5, height = 6)
enrichplot::dotplot(day3_ips_gsea, showCategory=7, x = 'NES')+
  theme_classic() +
  scale_fill_gradientn(colours = c("#F03C0C","#FFD1A6","#FFEEE0"), 
                       values = c(0,0.5,1.0),
                       limits = c(1, 2*10^-8))+
  scale_size(range = c(0, 9), limits = c(0, 200))
dev.off()

pdf('~/Documents/ÖverbyLab/single_cell_ISG_figures/astrocytes_fig/wt_day5_dotplot.pdf', width = 5.5, height = 6)
enrichplot::dotplot(day5_wt_gsea, showCategory=7, x = 'NES')+
  theme_classic() +
  scale_fill_gradientn(colours = c("#F03C0C","#FFD1A6","#FFEEE0"), 
                       values = c(0,0.5,1.0),
                       limits = c(1, 2*10^-8))+
  scale_size(range = c(0, 9), limits = c(0, 200))
dev.off()

pdf('~/Documents/ÖverbyLab/single_cell_ISG_figures/astrocytes_fig/ips_day5_dotplot.pdf', width = 5.5, height = 6)
enrichplot::dotplot(day5_ips_gsea, showCategory=7, x = 'NES')+
  theme_classic() +
  scale_fill_gradientn(colours = c("#F03C0C","#FFD1A6","#FFEEE0"), 
                       values = c(0,0.5,1.0),
                       limits = c(1, 2*10^-8))+
  scale_size(range = c(0, 9), limits = c(0, 200))
dev.off()

gseaplot2(day3_wt_gsea, geneSetID = 1:5, pvalue_table = TRUE)

pdf('~/Documents/ÖverbyLab/single_cell_ISG_figures/astrocytes_fig/wt_day3_ridgeplot.pdf', width = 7, height = 6)
ridgeplot(day3_wt_gsea, showCategory = 7) + labs(x = "enrichment distribution") + 
  theme_classic()+
  scale_fill_gradientn(colours = c("#D9530B","#F57456","#FFD1A6","#FFEEE0"), 
                       values = c(0,0.4,0.6,1.0),
                       limits = c(1, 2*10^-8))+
  theme(text = element_text(size = 18))+
  xlim(0, 8.5)
dev.off()

pdf('~/Documents/ÖverbyLab/single_cell_ISG_figures/astrocytes_fig/ips_day3_ridgeplot.pdf', width = 7, height = 6)
ridgeplot(day3_ips_gsea, showCategory = 7) + labs(x = "enrichment distribution") + 
  theme_classic()+
  scale_fill_gradientn(colours = c("#D9530B","#F57456","#FFD1A6","#FFEEE0"), 
                       values = c(0,0.4,0.6,1.0),
                       limits = c(1, 2*10^-8))+
  xlim(0, 8.5)+
  theme(text = element_text(size = 18))
dev.off()

pdf('~/Documents/ÖverbyLab/single_cell_ISG_figures/astrocytes_fig/wt_day5_ridgeplot.pdf', width = 7, height = 6)
ridgeplot(day5_wt_gsea, showCategory = 7) + labs(x = "enrichment distribution") + 
  theme_classic()+
  scale_fill_gradientn(colours = c("#D9530B","#F57456","#FFD1A6","#FFEEE0"), 
                       values = c(0,0.4,0.6,1.0),
                       limits = c(1, 2*10^-8))+
  theme(text = element_text(size = 18))+
  xlim(0, 8.5)
dev.off()

pdf('~/Documents/ÖverbyLab/single_cell_ISG_figures/astrocytes_fig/ips_day5_ridgeplot.pdf', width = 7, height = 6)
ridgeplot(day5_ips_gsea, showCategory = 7) + labs(x = "enrichment distribution") + 
  theme_classic()+
  scale_fill_gradientn(colours = c("#D9530B","#F57456","#FFD1A6","#FFEEE0"), 
                       values = c(0,0.4,0.6,1.0),
                       limits = c(1, 2*10^-8))+
  theme(text = element_text(size = 18))+
  xlim(0, 8.5)
dev.off()

#Top pathway gsea running score plots
pdf('~/Documents/ÖverbyLab/single_cell_ISG_figures/astrocytes_fig/wt_day3_gsea2_plot.pdf', width = 9, height = 5)
gseaplot2(day3_wt_gsea, geneSetID = 1:6, subplots = 1) +
  ylim(c(-0.15,0.7))+
  theme(legend.position = 'right',
        text = element_text(size = 18))
dev.off()

#Paths to plot in ips day 3 
top_paths <- day3_wt_gsea@result[1:6, 'ID']

dotplot(day3_ips_gsea, showCategory=10)
ips_paths_to_plot <- which(day3_ips_gsea@result$ID %in% top_paths)
pdf('~/Documents/ÖverbyLab/single_cell_ISG_figures/astrocytes_fig/ips_day3_gsea2_plot.pdf', width = 9, height = 5)
gseaplot2(day3_ips_gsea, geneSetID = ips_paths_to_plot, subplots = 1) + 
  ylim(c(-0.15,0.7))+
  theme(legend.position = 'right',
        text = element_text(size = 18))
dev.off()

#Day 5 wt
dotplot(day5_wt_gsea, showCategory=10)
pdf('~/Documents/ÖverbyLab/single_cell_ISG_figures/astrocytes_fig/wt_day5_gsea2_plot.pdf', width = 10, height = 5)
gseaplot2(day5_wt_gsea, geneSetID = 1:6, subplots = 1) +
  ylim(c(-0.15,0.7))+
  theme(legend.position = 'right',
        text = element_text(size = 18))
dev.off()

#Paths to plot in ips day 5 
top_paths_5 <- day5_wt_gsea@result[1:6, 'ID']
ips_5_paths_to_plot <- which(day5_ips_gsea@result$ID %in% top_paths_5)

pdf('~/Documents/ÖverbyLab/single_cell_ISG_figures/astrocytes_fig/ips_day5_gsea2_plot.pdf', width = 9, height = 5)
gseaplot2(day5_ips_gsea, geneSetID = 1:6, subplots = 1) +
  ylim(c(-0.15,0.7))+
  theme(legend.position = 'right',
        text = element_text(size = 18))
dev.off()

#Look at genes in negative regulation of immune response pathway
neg_reg_genes <- day5_wt_gsea@result[4,]$core_enrichment
neg_reg_genes_split <- strsplit(neg_reg_genes, "/")[[1]]

#Most clear negative regulators from pathway based on research
neg_reg_genes_curated <- c('Usp18', 'Serping1', 'H2-T23', 'Parp14', 'Lgals9',
                           'Dhx58', 'Clec2d')
#'Slamf8', 'Serpinb9b', 'Clec12b'

astros$treatment_geno_time <- paste(astros$Treatment, astros$Genotype, astros$Timepoint, sep = '_')
astro_neg_dot_dat <- DotPlot(astros, features = neg_reg_genes_curated, scale = FALSE, group.by = 'treatment_geno_time')$data

pdf('~/Documents/ÖverbyLab/single_cell_ISG_figures/astrocytes_fig/inflammation_neg_reg.pdf', width = 6, height = 6)
astro_neg_dot_dat %>% tidyr::separate(id, into = c('treatment', 'geno', 'time'), sep = '_') %>% 
  dplyr::mutate(treatment_time = paste(treatment, time, sep = '_')) %>% 
  dplyr::mutate(geno = factor(geno, levels = c('WT', 'IPS1'))) %>% 
  ggplot(aes(x = treatment_time, y = features.plot, fill = avg.exp.scaled, size = pct.exp))+
  facet_wrap(~geno) + 
  geom_point(pch = 21)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 75, hjust = 1))+
  scale_fill_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white"), 
                       values = c(1.0,0.7,0.4,0))
dev.off()

#Custom gsea plots with multiple days on one plot

#First, defense response to virus pathway

timepoint_gsea_plot <- function(pathway_name){
  dat_list = list(day3_wt_gsea, day5_wt_gsea, day3_ips_gsea, day5_ips_gsea)
  plot_labels = c('wt day 3', 'wt day 5', 'ips day 3', 'ips day 5')
  
  dat_to_plot_list <- lapply(dat_list, FUN = function(x){
    pathway_row <- which(x@result$Description == pathway_name)
    gsInfo(x, pathway_row)
  })
  
  for(i in 1:4){
    dat_to_plot_list[[i]]$sample = plot_labels[i]
  }
  
  dat_to_plot <- do.call(rbind, dat_to_plot_list)
  
  length(unique(c(dat_to_plot$Description))) == 1  #Double checking same pathway is plotted
  
  dat_to_plot$sample = factor(dat_to_plot$sample, levels = plot_labels)
  
  ggplot(dat_to_plot, aes_(x = ~x)) + xlab(NULL) +
    theme_classic(11) +
    theme(panel.grid.major = element_line(colour = "grey92"),
          panel.grid.minor = element_line(colour = "grey92"),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          text = element_text(size = 24)) +
    scale_x_continuous(expand=c(0,0))+
    geom_line(aes_(y = ~runningScore, color= ~sample),
              size=1) +
    theme(legend.position = c(.8, .8), legend.title = element_blank(),
          legend.background = element_rect(fill = "transparent"))+
    ggtitle(pathway_name)
}

pdf('~/Documents/ÖverbyLab/single_cell_ISG_figures/astrocytes_fig/response_to_virus_gsea.pdf', width = 7, height = 6)
timepoint_gsea_plot('defense response to virus')
dev.off()

pdf('~/Documents/ÖverbyLab/single_cell_ISG_figures/astrocytes_fig/type_1_inf_production.pdf', width = 7, height = 6)
timepoint_gsea_plot('type I interferon production')
dev.off()
#timepoint_gsea_plot('negative regulation of innate immune response')

## extract gsea result of selected geneSet
#Copied from enrichplot docs https://rdrr.io/bioc/enrichplot/src/R/gseaplot.R

gseaScores <- getFromNamespace("gseaScores", "DOSE")

gsInfo <- function(object, geneSetID) {
  geneList <- object@geneList
  
  if (is.numeric(geneSetID))
    geneSetID <- object@result[geneSetID, "ID"]
  
  geneSet <- object@geneSets[[geneSetID]]
  exponent <- object@params[["exponent"]]
  df <- gseaScores(geneList, geneSet, exponent, fortify=TRUE)
  df$ymin <- 0
  df$ymax <- 0
  pos <- df$position == 1
  h <- diff(range(df$runningScore))/20
  df$ymin[pos] <- -h
  df$ymax[pos] <- h
  df$geneList <- geneList
  
  df$Description <- object@result[geneSetID, "Description"]
  return(df)
}
