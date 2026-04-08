library(Seurat)
library(RColorBrewer)
library(dplyr)
library(UpSetR)
library(ggplot2)
library(stringr)
library(tidyr)
source('~/Documents/ÖverbyLab//scripts/langatFunctions.R')

#Load data
ParseSeuratObj_int <- LoadSeuratRds("~/Documents/ÖverbyLab/data/FilteredRpcaIntegratedDatNoDoublets.rds") 
ParseSeuratObj_int$hasVirus = ifelse(ParseSeuratObj_int$virusCountPAdj >= 10, 1, 0)
ParseSeuratObj_int$manualAnnotation[ParseSeuratObj_int$manualAnnotation == 'Macrophage/Monocytes'] = 'Macro/Mono'

#Create new column for use in pseudobulk later
ParseSeuratObj_int$Treatment_celltype <- paste(ParseSeuratObj_int$Treatment, ParseSeuratObj_int$manualAnnotation, sep = '_')

#Add extra column for umap plot
ParseSeuratObj_int$geno_timepoint_treatment = paste(ParseSeuratObj_int$Genotype, ParseSeuratObj_int$Timepoint, ParseSeuratObj_int$Treatment)

#Check data
newCols <-  c(brewer.pal(12, 'Paired'), '#99FFE6', '#CE99FF', '#18662E','#737272',  '#FF8AEF')
newCols[11] =  '#FF8AEF'

DimPlot(ParseSeuratObj_int, label = FALSE, group.by = 'manualAnnotation', reduction = 'umap.integrated',
        cols = newCols)+
  theme(axis.ticks = element_blank(),
        axis.text=element_blank(),
        legend.text=element_text(size=17))+
  xlab('Umap1')+
  ylab('Umap2')+
  guides(color=guide_legend(override.aes=list(size=8)))+
  ggtitle('')

DimPlot(ParseSeuratObj_int, label = FALSE, group.by = 'Treatment', reduction = 'umap.integrated',
        cols = newCols)

#Start looking at overall astrocyte degs, then split by group
astrocytes <- subset(ParseSeuratObj_int, manualAnnotation == 'Astrocytes')

#Recluster astrocytes
#This was true before so can switch to get old umap
astrocytes <- prepSeuratObj(astrocytes, use_all_genes = FALSE)
ElbowPlot(astrocytes, ndims = 40)
astrocytes <- prepUmapSeuratObj(astrocytes, nDims = 20, reductionName = 'astrocytes_umap', resolution_value = 0.8)

#UMAPs
DimPlot(astrocytes, reduction = 'astrocytes_umap')
astrocytes$infection_group <- ifelse(astrocytes$Treatment %in% c('rChLGTV', 'rLGTV'), 'infected', 'uninfected')
DimPlot(astrocytes, reduction = 'astrocytes_umap', group.by = 'infection_group')
DimPlot(astrocytes, reduction = 'astrocytes_umap', group.by = 'Genotype')
DimPlot(astrocytes, reduction = 'astrocytes_umap', group.by = 'Organ')+
  xlab('')+
  ylab('')+
  ggtitle('Astrocyte UMAP')
DimPlot(astrocytes, reduction = 'astrocytes_umap', group.by = 'Timepoint')
DimPlot(astrocytes, reduction = 'astrocytes_umap', group.by = 'geno_timepoint_treatment', cols = newCols)
table(astrocytes$infection_group, astrocytes$Genotype)

#Astrocyte markers
astrocytes$time_treatment <- paste(astrocytes$Timepoint, astrocytes$Treatment, sep = '_')

#Change time_treatment so pbs is all combined for plotting and comparisons
astrocytes$time_treatment <- ifelse(astrocytes$time_treatment == 'Day 5_PBS' | astrocytes$time_treatment == 'Day 3_PBS', 'PBS', astrocytes$time_treatment)

ips_astrocytes <- subset(astrocytes, Genotype == 'IPS1' & Treatment != 'rLGTV')
ips_astrocytes_cerebrum <- subset(ips_astrocytes, Organ == 'Cerebrum')
wt_astrocytes <- subset(astrocytes, Genotype == 'WT' & Treatment != 'rLGTV')
wt_astrocytes_cerebrum <- subset(wt_astrocytes,  Organ == 'Cerebrum')

ips_astrocytes_cerebrum_markers <- FindAllMarkers(ips_astrocytes_cerebrum, only.pos = TRUE, test.use = 'MAST',
                                                group.by = 'time_treatment')
wt_astrocytes_cerebrum_markers <- FindAllMarkers(wt_astrocytes_cerebrum, only.pos = TRUE, test.use = 'MAST',
                                                  group.by = 'time_treatment')

#For now only plotting genes from scale.data, otherwise have to use data slot which ruins heatmaps
#can try later to manually scale from data slow
genes_to_choose_from_ips <- rownames(ips_astrocytes_cerebrum[['RNA']]$scale.data)
genes_to_plot_ips <- c()
for(i in 1:length(unique(ips_astrocytes_cerebrum_markers$cluster))){
  cur_dat <- ips_astrocytes_cerebrum_markers[ips_astrocytes_cerebrum_markers$cluster == unique(ips_astrocytes_cerebrum_markers$cluster)[i],]
  cur_dat <- dplyr::arrange(cur_dat, p_val_adj)
  cur_dat <- cur_dat[!cur_dat$gene %in% genes_to_plot_ips,]
  cur_dat <- cur_dat[cur_dat$gene %in% genes_to_choose_from_ips,]
  genes_to_plot_ips <- c(genes_to_plot_ips, head(cur_dat$gene, n = 15))
}

#Need to make it so genes do not repeat.
genes_to_choose_from_wt <- rownames(wt_astrocytes_cerebrum[['RNA']]$scale.data)
genes_to_plot_wt <- c()
for(i in 1:length(unique(wt_astrocytes_cerebrum_markers$cluster))){
  cur_dat <- wt_astrocytes_cerebrum_markers[wt_astrocytes_cerebrum_markers$cluster == unique(wt_astrocytes_cerebrum_markers$cluster)[i],]
  cur_dat <- dplyr::arrange(cur_dat, p_val_adj)
  cur_dat <- cur_dat[!cur_dat$gene %in% genes_to_plot_wt,]
  cur_dat <- cur_dat[cur_dat$gene %in% genes_to_choose_from_wt,]
  genes_to_plot_wt <- c(genes_to_plot_wt, head(cur_dat$gene, n = 15))
}

ips_astrocytes_cerebrum$time_treatment <- factor(ips_astrocytes_cerebrum$time_treatment, 
                                                 levels = c("PBS", "Day 3_rChLGTV", "Day 4_rChLGTV", "Day 5_rChLGTV"))
wt_astrocytes_cerebrum$time_treatment <- factor(wt_astrocytes_cerebrum$time_treatment, 
                                                 levels = c("PBS", "Day 3_rChLGTV", "Day 4_rChLGTV", "Day 5_rChLGTV"))

tiff('~/Documents/ÖverbyLab/scPlots/ips_astro_cerebrum_time_heatmap.tiff', width = 900, height = 1100, res = 100)
DoHeatmap(ips_astrocytes_cerebrum, features = genes_to_plot_ips, group.by = 'time_treatment')+
  scale_fill_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white","lightblue"),
                       values = c(1.0,0.9,0.7,0.5,0))+
  ggtitle('IPS Cerebrum Astrocytes')+
  theme(plot.title = element_text(size = 25, vjust = 3))
dev.off()

tiff('~/Documents/ÖverbyLab/scPlots/wt_astro_cerebrum_time_heatmap.tiff', width = 900, height = 1100, res = 100)
DoHeatmap(wt_astrocytes_cerebrum, features = genes_to_plot_wt, group.by = 'time_treatment')+
  scale_fill_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white","lightblue"),
                       values = c(1.0,0.7,0.7,0.5,0))+
  ggtitle('WT Cerebrum Astrocytes')+
  theme(plot.title = element_text(size = 25, vjust = 3))
dev.off()

#Loot at markers from above groups

#PBS markers
pbs_markers_ips <- ips_astrocytes_cerebrum_markers[ips_astrocytes_cerebrum_markers$cluster == 'PBS' & 
                                  ips_astrocytes_cerebrum_markers$avg_log2FC > 1 & 
                                    ips_astrocytes_cerebrum_markers$p_val_adj < 0.01,]
pbs_paths <- gprofiler2::gost(query = rownames(pbs_markers_ips), organism = 'mmusculus', evcodes = TRUE)

#Day 3 markers
day_3_markers_ips <- ips_astrocytes_cerebrum_markers[ips_astrocytes_cerebrum_markers$cluster == 'Day 3_rChLGTV' & 
                                                     ips_astrocytes_cerebrum_markers$avg_log2FC > 1 & 
                                                     ips_astrocytes_cerebrum_markers$p_val_adj < 0.01,]
day_3_paths <- gprofiler2::gost(query = rownames(day_3_markers_ips), organism = 'mmusculus', evcodes = TRUE)

#Day 4 markers
day_4_markers_ips <- ips_astrocytes_cerebrum_markers[ips_astrocytes_cerebrum_markers$cluster == 'Day 4_rChLGTV' & 
                                                       ips_astrocytes_cerebrum_markers$avg_log2FC > 1 & 
                                                       ips_astrocytes_cerebrum_markers$p_val_adj < 0.01,]
day_4_paths <- gprofiler2::gost(query = rownames(day_4_markers_ips), organism = 'mmusculus', evcodes = TRUE)


#Also look at pseudobulk markers
#Maybe use all astrocytes to create usable pseudobulk object
ips_astrocyte_cerebrum_bulk <- createPseudoBulk(ips_astrocytes_cerebrum, c('time_treatment'), intercept = FALSE)
ips_astrocyte_cerebrum_bulk <- DESeq(ips_astrocyte_cerebrum_bulk)
resultsNames(ips_astrocyte_cerebrum_bulk)
ips_astrocyte_bulk_treatment_markers <- results(ips_astrocyte_bulk, name = 'Treatment_rChLGTV_vs_PBS') %>% 
  as.data.frame() %>% dplyr::filter(padj < 0.01 & abs(log2FoldChange) > 1) %>% 
  rownames_to_column(var = 'gene') %>% 
  dplyr::arrange(padj) %>% dplyr::mutate(direction = ifelse(log2FoldChange > 0, 'up', 'down')) %>% 
  dplyr::group_by(direction) %>% 
  dplyr::slice_head(n = 10)



#Look at markers
astro_infected_markers <- FindMarkers(astrocytes, ident.1 = 'infected', group.by = 'infection_group',
                             test.use = 'MAST', only.pos = TRUE)
astro_sig_markers <- astro_infected_markers %>% dplyr::filter(avg_log2FC > 1 & p_val_adj < 0.01)
astro_marker_paths <- gprofiler2::gost(query = rownames(astro_sig_markers), organism = 'mmusculus', evcodes = TRUE)
astro_marker_paths$result

#Look into path 558 - seems to be evidence of signalling
astro_marker_paths$result[astro_marker_paths$result$source == 'KEGG',]

#Plot markers - overall, a lot of clear ISGs
#Irf1/8 paper https://pmc.ncbi.nlm.nih.gov/articles/PMC4821649/
FeaturePlot(astrocytes, features = 'Irf1', reduction = 'astrocytes_umap')
DotPlot(astrocytes, features = 'Irf1', scale = FALSE, group.by = 'Treatment')

#Rig 1 type 1 interferon response https://www.sciencedirect.com/science/article/pii/S2211124715004647
FeaturePlot(astrocytes, features = 'Ddx60', reduction = 'astrocytes_umap')
DotPlot(astrocytes, features = 'Ddx60', scale = FALSE, group.by = 'Treatment')

#MHC class 1 https://www.nature.com/articles/nri3339
FeaturePlot(astrocytes, features = 'Nlrc5', reduction = 'astrocytes_umap')
DotPlot(astrocytes, features = 'Nlrc5', scale = FALSE, group.by = 'Treatment')

DotPlot(astrocytes, features = c('Stat1', 'Stat2'), scale = FALSE, group.by = 'Treatment')

FeaturePlot(astrocytes, features = 'Ccl2', reduction = 'astrocytes_umap')
FeaturePlot(astrocytes, features = 'Cxcl10', reduction = 'astrocytes_umap')

#See which groups express Ccl2
gene_presence_group_bars <- function(dat, gene, main = ''){
  dat$gene_dat <- FetchData(dat, vars = gene, slot = 'data')
  gene_plot <- table(dat$geno_timepoint_treatment, dat$gene_dat > 0) %>% 
    as.data.frame() %>% 
    dplyr::group_by(Var1) %>% 
    dplyr::mutate(prop = Freq/sum(Freq)) %>% 
    ggplot(aes(x = Var1, y = prop, fill = Var2))+
    geom_bar(stat = 'identity', position = 'stack')+
    ylab(gene) +
    theme(axis.text.x = element_text(angle = 75, vjust = 0.5, size = 15),
          axis.title.y = element_text(size = 15))+
    ggtitle(main)
  print(gene_plot)
}

gene_presence_group_bars(astrocytes, 'Ccl2')
gene_presence_group_bars(ependymal, 'Ccl2')
gene_presence_group_bars(microglia, 'Ccl2')

gene_presence_group_bars(astrocytes, 'Cxcl10', main = 'Astrocytes')
gene_presence_group_bars(ependymal, 'Cxcl10',main = 'Ependymal')
gene_presence_group_bars(microglia, 'Cxcl10')

#'Cxcl10' should also be plotted but drowns out other features
DotPlot(astrocytes, features = c('Ccl2', 'Ccl7', 'Ccl12', 'Ccl4',
                                 'Ccl3', 'Ccl5', 'Cxcl9',
                                 'Cxcr3', 'Ccr1'), scale = FALSE,
        group.by = 'geno_timepoint_treatment')+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.8))

#Look at other celltypes to compare
DotPlot(ependymal, features = c('Ccl2', 'Ccl7', 'Ccl12', 'Ccl4',
                                 'Ccl3', 'Ccl5', 'Cxcl9',
                                 'Cxcr3', 'Ccr1'), scale = FALSE,
        group.by = 'geno_timepoint_treatment')+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.8))

table(astrocytes$geno_timepoint_treatment)
FeaturePlot(astrocytes, features = 'Ccl3', reduction = 'astrocytes_umap')

#Cxcl10
DotPlot(astrocytes, features = c('Cxcl10'), scale = FALSE,
        group.by = 'geno_timepoint_treatment')+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.8))

DotPlot(ependymal, features = c('Cxcl10'), scale = FALSE,
        group.by = 'geno_timepoint_treatment')+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.8))

######### Look at number of degs by celltype ############
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##


#Pseudobulk comparison of celltype number of degs between infected and uninfected - split between genotypes to reduce time
#Also check without pseudobulk to see if major changes
ParseSeuratObj_int_bulk <- createPseudoBulk(ParseSeuratObj_int, variables = c('Genotype', 'Treatment_celltype', 'Timepoint', 'Organ'))

ParseSeuratObj_int_bulk <- DESeq(ParseSeuratObj_int_bulk)

resultsNames(ParseSeuratObj_int_bulk)

#Non pseudobulk version
celltypes <- unique(ParseSeuratObj_int$manualAnnotation)
genotypes <- unique(ParseSeuratObj_int$Genotype)
num_markers_df <- data.frame(celltype = rep(celltypes,2), genotype = rep(genotypes, each = length(celltypes)),
                             degs_lgtv_pbs = 0, degs_chlgtv_pbs = 0)

#Need to split by genotype
for(i in 1:length(genotypes)){
  cur_genotype = genotypes[i]
  cur_genotype_dat <- subset(ParseSeuratObj_int, Genotype == cur_genotype)
  for(j in 1:length(celltypes)){
    cur_celltype = celltypes[j]
    marker_ident_1 <- paste('PBS', cur_celltype, sep = '_')
    marker_ident_2 <- paste('rLGTV', cur_celltype, sep = '_')
    marker_ident_3 <- paste('rChLGTV', cur_celltype, sep = '_')
    sig_markers_lgtv <- FindMarkers(cur_genotype_dat, group.by = 'Treatment_celltype', ident.1 = marker_ident_1, ident.2 = marker_ident_2, 
                               test.use = 'MAST') %>% 
      dplyr::filter(abs(avg_log2FC) > 1 & p_val_adj < 0.01)
    sig_markers_ch <- FindMarkers(cur_genotype_dat, group.by = 'Treatment_celltype', ident.1 = marker_ident_1, ident.2 = marker_ident_3,
                                  test.use = 'MAST') %>% 
      dplyr::filter(abs(avg_log2FC) > 1 & p_val_adj < 0.01)
    num_markers_df[num_markers_df$celltype == cur_celltype & num_markers_df$genotype == cur_genotype,]['degs_lgtv_pbs'] = nrow(sig_markers_lgtv)
    num_markers_df[num_markers_df$celltype == cur_celltype & num_markers_df$genotype == cur_genotype,]['degs_chlgtv_pbs'] = nrow(sig_markers_ch)
    print(paste("Done with", cur_celltype))
  }
}

write.csv(num_markers_df, '~/Documents/ÖverbyLab/single_cell_degs_per_celltype.csv', row.names = FALSE)


#Validate some of results to make sure function worked properly
#Seems to be working
wt_dat <- subset(ParseSeuratObj_int, Genotype == 'WT')
ips_dat <- subset(ParseSeuratObj_int, Genotype == 'IPS1')
#WT astrocytes
lgtv_astro_sig_markers <- FindMarkers(wt_dat, group.by = 'Treatment_celltype', ident.1 = 'PBS_Astrocytes', ident.2 = 'rLGTV_Astrocytes', 
            test.use = 'MAST')
lgtv_astro_sig_markers %>% dplyr::filter(abs(avg_log2FC) > 1 & p_val_adj < 0.01) %>% nrow() == num_markers_df[num_markers_df$celltype == 'Astrocytes' &
                                                                                                                num_markers_df$genotype == 'WT',]['degs_lgtv_pbs']

#IPS astrocytes
ips_chlgtv_astro_sig_markers <- FindMarkers(ips_dat, group.by = 'Treatment_celltype', ident.1 = 'PBS_Astrocytes', ident.2 = 'rChLGTV_Astrocytes', 
                                      test.use = 'MAST')
ips_chlgtv_astro_sig_markers %>% dplyr::filter(abs(avg_log2FC) > 1 & p_val_adj < 0.01) %>% nrow() == num_markers_df[num_markers_df$celltype == 'Astrocytes' &
                                                                                                                num_markers_df$genotype == 'IPS1',]['degs_chlgtv_pbs']

#IPS t cells
ips_chlgtv_tcell_sig_markers <- FindMarkers(ips_dat, group.by = 'Treatment_celltype', ident.1 = 'PBS_T cells', ident.2 = 'rChLGTV_T cells', 
                                            test.use = 'MAST')
ips_chlgtv_tcell_sig_markers %>% dplyr::filter(abs(avg_log2FC) > 1 & p_val_adj < 0.01) %>% nrow() == num_markers_df[num_markers_df$celltype == 'T cells' &
                                                                                                                      num_markers_df$genotype == 'IPS1',]['degs_chlgtv_pbs']


#Plot results
num_markers_df <- read.csv('~/Documents/ÖverbyLab/single_cell_degs_per_celltype.csv')
num_markers_df
num_markers_df %>% tidyr::pivot_longer(cols = tidyr::starts_with('degs'), names_to = 'comp', values_to = 'num_degs') %>% 
  ggplot(aes(x = celltype, y = num_degs, fill = comp))+
  geom_bar(stat = 'identity', position = 'dodge')+
  facet_wrap(~genotype)+
  coord_flip()+
  ggtitle('Total DEGs per celltype')


###### Infected astrocytes vs other cell types #######
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

#Can look at degs in astrocytes in mock mice, the infected,
#then compare which genes are upregulated in infection but not mock
pbs <- subset(ParseSeuratObj_int, Treatment == 'PBS')
chlgtv <- subset(ParseSeuratObj_int, Treatment == 'rChLGTV')

#pbs astro markers
pbs_astro_markers <- FindMarkers(pbs, group.by = 'manualAnnotation', 
                                 ident.1 = 'Astrocytes', test.use = 'MAST',
                                 only.pos = TRUE)

chlgtv_astro_markers <- FindMarkers(chlgtv, group.by = 'manualAnnotation', 
                                 ident.1 = 'Astrocytes', test.use = 'MAST',
                                 only.pos = TRUE)

pbs_sig_markers <- pbs_astro_markers %>% dplyr::filter(avg_log2FC > 1 & p_val_adj < 0.01)
chlgtv_sig_markers <- chlgtv_astro_markers %>% dplyr::filter(avg_log2FC > 1 & p_val_adj < 0.01)

#Majority (94%) of degs from pbs remain in chimeric 
length(intersect(rownames(pbs_sig_markers), rownames(chlgtv_sig_markers)))/nrow(pbs_sig_markers)
length(intersect(rownames(pbs_sig_markers), rownames(chlgtv_sig_markers)))/nrow(chlgtv_sig_markers) #Only 62% other way

markers_only_chlgtv <- chlgtv_sig_markers[!rownames(chlgtv_sig_markers) %in% rownames(pbs_sig_markers),]
markers_only_chlgtv_paths <- gprofiler2::gost(query = rownames(markers_only_chlgtv), organism = 'mmusculus', evcodes = TRUE)
markers_only_chlgtv_paths$result[21:30,]

#Make UMAPs by organ
astrocytes_cerebrum <- subset(astrocytes, Organ == 'Cerebrum')
astrocytes_cerebellum <- subset(astrocytes, Organ == 'Cerebellum')

#Recluster astrocytes
astrocytes_cerebrum <- prepSeuratObj(astrocytes_cerebrum)
ElbowPlot(astrocytes_cerebrum, ndims = 40)
astrocytes_cerebrum <- prepUmapSeuratObj(astrocytes_cerebrum, nDims = 20, reductionName = 'astrocytes_cerebrum_umap', resolution_value = 0.8)

DimPlot(astrocytes_cerebrum, reduction = 'astrocytes_cerebrum_umap', group.by = 'Genotype')+
  ylab('')+
  xlab('')+
  ggtitle('Cerebrum Astrocytes')

#Recluster astrocytes
astrocytes_cerebellum <- prepSeuratObj(astrocytes_cerebellum)
ElbowPlot(astrocytes_cerebellum, ndims = 40)
astrocytes_cerebellum <- prepUmapSeuratObj(astrocytes_cerebellum, nDims = 20, reductionName = 'astrocytes_cerebellum_umap', resolution_value = 0.8)

DimPlot(astrocytes_cerebellum, reduction = 'astrocytes_cerebellum_umap', group.by = 'Genotype')+
  xlab('')+
  ylab('')+
  ggtitle('Cerebellum Astrocytes')

######### Genes upregulated by celltype for upsetplot ##########
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## #

wt_chimeric <- subset(ParseSeuratObj_int, Genotype == 'WT' & Treatment %in% c('PBS', 'rChLGTV'))
ips_chimeric <- subset(ParseSeuratObj_int, Genotype == 'IPS1' & Treatment %in% c('PBS', 'rChLGTV'))

#Can skip for loop and read in data here
degs_by_celltype_wt <- readRDS('~/Documents/ÖverbyLab/cell_upsetPlots/wt_celltype_degs.rds')
degs_by_celltype_ips <- readRDS('~/Documents/ÖverbyLab/cell_upsetPlots/ips_celltype_degs.rds')

#Or run for loop
celltypes <- unique(wt_chimeric$manualAnnotation)
degs_by_celltype_wt <- list()
degs_by_celltype_ips <- list()

for(i in 1:length(celltypes)){
  print(paste('Beginning celltype', celltypes[i]))
  cur_dat_wt <- subset(wt_chimeric, manualAnnotation == celltypes[i])
  cur_dat_ips <- subset(ips_chimeric, manualAnnotation == celltypes[i])
  celltype_markers_wt <- FindMarkers(cur_dat_wt, group.by = 'Treatment', ident.1 = 'rChLGTV', ident.2 = 'PBS', 
              test.use = 'MAST')
  celltype_markers_ips <- FindMarkers(cur_dat_ips, group.by = 'Treatment', ident.1 = 'rChLGTV', ident.2 = 'PBS', 
                                     test.use = 'MAST')
  degs_by_celltype_wt[[i]] = celltype_markers_wt
  degs_by_celltype_ips[[i]] = celltype_markers_ips
  names(degs_by_celltype_wt)[[i]] = celltypes[i]
  names(degs_by_celltype_ips)[[i]] = celltypes[i]
  print(paste('Done with celltype', celltypes[i]))
}

#saveRDS(degs_by_celltype_wt, '~/Documents/ÖverbyLab/cell_upsetPlots/wt_celltype_degs.rds')
#saveRDS(degs_by_celltype_ips, '~/Documents/ÖverbyLab/cell_upsetPlots/ips_celltype_degs.rds')

celltype_markers_sig_wt <- lapply(degs_by_celltype_wt, FUN = function(x){
  celltype_dat <- as.data.frame(x)
  celltype_dat[celltype_dat$p_val_ad < 0.01 & celltype_dat$avg_log2FC > 1,]
})

celltype_markers_sig_ips <- lapply(degs_by_celltype_ips, FUN = function(x){
  celltype_dat <- as.data.frame(x)
  celltype_dat[celltype_dat$p_val_ad < 0.01 & celltype_dat$avg_log2FC > 1,]
})

num_degs_wt <- lapply(celltype_markers_sig_wt, nrow)
num_degs_wt_df <- data.frame(celltypes = names(celltype_markers_sig_wt), num_degs = unlist(num_degs_wt))
num_degs_wt_df

deg_names_wt <- lapply(celltype_markers_sig_wt, FUN = function(x){
  rownames(x)
})

deg_names_ips <- lapply(celltype_markers_sig_ips, FUN = function(x){
  rownames(x)
})

upset_dat_wt <- UpSetR::fromList(deg_names_wt)
upset_dat_ips <- UpSetR::fromList(deg_names_ips)

#Upset plots of upregulated degs
#Copying these with ~700 x 1400 dimensions
UpSetR::upset(upset_dat_wt, nsets = 10, order.by = 'freq', text.scale = 2)
UpSetR::upset(upset_dat_ips, nsets = 10, order.by = 'freq', text.scale = 2)

#Get commonly upregulated genes between celltypes in upset data
shared_wt_genes <- Reduce(intersect, deg_names_wt[c('Microglia', 'Choroid Plexus', 'Ependymal', 'Endothelial', 'Macro/Mono', 'T cells', 'Astrocytes')])
shared_ips_genes <- Reduce(intersect, deg_names_ips[c('Microglia', 'Choroid Plexus', 'Ependymal', 'Endothelial', 'Macro/Mono', 'T cells', 'Astrocytes')])


unique_wt <- length(shared_wt_genes[!shared_wt_genes %in% shared_ips_genes])
unique_ips <- length(shared_ips_genes[!shared_ips_genes %in% shared_wt_genes])
shared_genes <- length(intersect(x = shared_ips_genes, y = shared_wt_genes))

ips_vs_wt_up_genes <- data.frame('type' = factor(c('wt_only', 'ips_only', 'shared'), levels = c('wt_only', 'ips_only', 'shared')), 
                                 counts = c(unique_wt, unique_ips, shared_genes))

ggplot(ips_vs_wt_up_genes, aes(x = type, y = counts, fill = type))+
  geom_bar(stat = 'identity')+
  theme(legend.position = 'none')+
  scale_fill_manual(values = c("#A0AEF2", "#88BBBD",  "#AD88BD"))

#Function that, given a list of findmarkers output (dataframes with lfc and pvals etc...), named by celltype,
#returns a genes only upregulated in given celltype
unqiue_celltype_markers <- function(markers_list, celltype){
  all_marker_names <- lapply(markers_list, FUN = function(x){
    rownames(x)
  })
  non_celltype_markers <- unname(unlist(all_marker_names[names(all_marker_names) != celltype]))
  celltype_markers <- all_marker_names[[celltype]]
  unique_celltype_markers <- celltype_markers[!celltype_markers %in% non_celltype_markers]
  return(unique_celltype_markers)
}

#Get all genes upregulated in non astrocytes in infection
astro_genes_wt <- unqiue_celltype_markers(celltype_markers_sig_wt, 'Astrocytes')
micro_genes_wt <- unqiue_celltype_markers(celltype_markers_sig_wt, 'Microglia')

#No pathways
gprofiler2::gost(query = astro_genes_wt, organism = 'mmusculus', evcodes = TRUE, sources = c('GO:BP', 'KEGG'))

#Lots of cell cycle genes
micro_paths <- gprofiler2::gost(query = micro_genes_wt, organism = 'mmusculus', evcodes = TRUE, sources = c('GO:BP', 'KEGG'))
micro_paths$result

#Dotplot specific genes
wt_chimeric_known <- subset(wt_chimeric, manualAnnotation != 'unknown')

for(i in 1:length(astro_genes)){
  gene_dot <- DotPlot(wt_chimeric, features = c(astro_genes[[i]]), group.by = 'Treatment_celltype', scale = FALSE)$data
  gene_dot_meta <- str_split_fixed(gene_dot$id, pattern = '_', n = 2)
  colnames(gene_dot_meta) <- c('treatment', 'celltype')
  gene_dot <- cbind(gene_dot, gene_dot_meta)
  
  print(ggplot(gene_dot, aes(x = treatment, y = celltype, fill = avg.exp.scaled))+
    geom_tile()+
    scale_fill_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white"), 
                         values = c(1.0,0.7,0.4,0))+
    ggtitle(astro_genes[[i]]))
}

gene_dot <- DotPlot(wt_chimeric_known, features = c('Grm5'), group.by = 'Treatment_celltype', scale = FALSE)$data
gene_dot_meta <- str_split_fixed(gene_dot$id, pattern = '_', n = 2)
colnames(gene_dot_meta) <- c('treatment', 'celltype')
gene_dot <- cbind(gene_dot, gene_dot_meta)
ggplot(gene_dot, aes(x = treatment, y = celltype, fill = avg.exp.scaled))+
        geom_tile()+
        scale_fill_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white"), 
                             values = c(1.0,0.7,0.4,0))+
        ggtitle('Grm5')


########### UpSet plots by timepoint ########### 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

#Only day 3 and 5 have pbs and infected
wt_chimeric_three <- subset(wt_chimeric, Timepoint == 'Day 3' | Treatment == 'PBS')
wt_chimeric_five <- subset(wt_chimeric, Timepoint == 'Day 5' | Treatment == 'PBS')

ips_chimeric_three <- subset(ips_chimeric, Timepoint == 'Day 3' | Treatment == 'PBS')
ips_chimeric_five <- subset(ips_chimeric, Timepoint == 'Day 5' | Treatment == 'PBS')

#Can skip running loop and load here 
wt_three_degs <- readRDS(file = '~/Documents/ÖverbyLab/cell_upsetPlots/wt_three_celltype_upset.rds')
wt_five_degs <- readRDS(file = '~/Documents/ÖverbyLab/cell_upsetPlots/wt_five_celltype_upset.rds')
ips_three_degs <- readRDS(file = '~/Documents/ÖverbyLab/cell_upsetPlots/ips_three_celltype_upset.rds')
ips_five_degs <- readRDS(file = '~/Documents/ÖverbyLab/cell_upsetPlots/ips_five_celltype_upset.rds')

#Otherwise can run from beginning
wt_three_degs <- list()
wt_five_degs <- list()
ips_three_degs <- list()
ips_five_degs <- list()

for(i in 1:length(celltypes)){
  print(paste('Beginning celltype', celltypes[i]))
  
  cur_dat_three <- subset(wt_chimeric_three, manualAnnotation == celltypes[i])
  cur_dat_five <- subset(wt_chimeric_five, manualAnnotation == celltypes[i])
  cur_dat_three_ips <- subset(ips_chimeric_three, manualAnnotation == celltypes[i])
  cur_dat_five_ips <- subset(ips_chimeric_five, manualAnnotation == celltypes[i])
  
  celltype_markers_three <- FindMarkers(cur_dat_three, group.by = 'Treatment', ident.1 = 'rChLGTV', ident.2 = 'PBS', 
                                     test.use = 'MAST')
 
  celltype_markers_five <- FindMarkers(cur_dat_five, group.by = 'Treatment', ident.1 = 'rChLGTV', ident.2 = 'PBS', 
                                       test.use = 'MAST')
  
  celltype_markers_three_ips <- FindMarkers(cur_dat_three_ips, group.by = 'Treatment', ident.1 = 'rChLGTV', ident.2 = 'PBS', 
                                        test.use = 'MAST')
  
  celltype_markers_five_ips <- FindMarkers(cur_dat_five_ips, group.by = 'Treatment', ident.1 = 'rChLGTV', ident.2 = 'PBS',
                                            test.use = 'MAST')
  
  wt_three_degs[[i]] = celltype_markers_three
  wt_five_degs[[i]] = celltype_markers_five
  ips_three_degs[[i]] = celltype_markers_three_ips
  ips_five_degs[[i]] = celltype_markers_five_ips
  
  names(wt_three_degs)[[i]] = celltypes[i]
  names(wt_five_degs)[[i]] = celltypes[i]
  names(ips_three_degs)[[i]] = celltypes[i]
  names(ips_five_degs)[[i]] = celltypes[i]
  
  print(paste('Done with celltype', celltypes[i]))
}

#saveRDS(wt_three_degs, file = '~/Documents/ÖverbyLab/cell_upsetPlots/wt_three_celltype_upset.rds')
#saveRDS(wt_five_degs, file = '~/Documents/ÖverbyLab/cell_upsetPlots/wt_five_celltype_upset.rds')
#saveRDS(ips_three_degs, file = '~/Documents/ÖverbyLab/cell_upsetPlots/ips_three_celltype_upset.rds')
#saveRDS(ips_five_degs, file = '~/Documents/ÖverbyLab/cell_upsetPlots/ips_five_celltype_upset.rds')

#function to generate upset data from list of degs
generate_upset_dat <- function(deg_list_dat){
  sig_dat <- lapply(deg_list_dat, FUN = function(x){
    celltype_dat <- as.data.frame(x)
    celltype_dat[celltype_dat$p_val_ad < 0.01 & celltype_dat$avg_log2FC > 1,]
  })
  deg_names <- lapply(sig_dat, FUN = function(x){
    rownames(x)
  })
  upset_dat<- UpSetR::fromList(deg_names)
}

upset_dat_wt_three <- generate_upset_dat(wt_three_degs)
upset_dat_wt_five <- generate_upset_dat(wt_five_degs)
upset_dat_ips_three <- generate_upset_dat(ips_three_degs)
upset_dat_ips_five <- generate_upset_dat(ips_five_degs)

pdf('~/Documents/ÖverbyLab/single_cell_ISG_figures/sc_celltype_fig_plots/wt_three_upset.pdf', height = 8, width = 13, onefile = FALSE)
UpSetR::upset(upset_dat_wt_three, nsets = 5, order.by = 'freq', text.scale = 3, point.size = 4, show.numbers = FALSE)
dev.off()

pdf('~/Documents/ÖverbyLab/single_cell_ISG_figures/sc_celltype_fig_plots/wt_five_upset.pdf', height = 8, width = 13, onefile = FALSE)
UpSetR::upset(upset_dat_wt_five, nsets = 5, order.by = 'freq', text.scale = 3, point.size = 4, show.numbers = FALSE)
dev.off()

pdf('~/Documents/ÖverbyLab/single_cell_ISG_figures/sc_celltype_fig_plots/ips_three_upset.pdf', height = 8, width = 13, onefile = FALSE)
UpSetR::upset(upset_dat_ips_three, nsets = 5, order.by = 'freq', text.scale = 3, point.size = 4, show.numbers = FALSE)
dev.off()

pdf('~/Documents/ÖverbyLab/single_cell_ISG_figures/sc_celltype_fig_plots/ips_five_upset.pdf', height = 8, width = 13, onefile = FALSE)
UpSetR::upset(upset_dat_ips_five, nsets = 5, order.by = 'freq', text.scale = 3, point.size = 4, show.numbers = FALSE)
dev.off()

#Look at microglia genes specifically
micro_genes_wt_three <- unqiue_celltype_markers(wt_three_sig, 'Microglia')
micro_genes_wt_five <- unqiue_celltype_markers(wt_five_sig, 'Microglia')

micro_wt_three_paths <- gprofiler2::gost(query = micro_genes_wt_three, organism = 'mmusculus', evcodes = TRUE, sources = c('GO:BP', 'KEGG'))
micro_wt_five_paths <- gprofiler2::gost(query = micro_genes_wt_five, organism = 'mmusculus', evcodes = TRUE, sources = c('GO:BP', 'KEGG'))

ggplot(head(micro_wt_five_paths$result), aes(x = -log10(p_value), y = term_name))+
  geom_bar(stat = 'identity')

#Number of DEGs
get_sig_degs <- function(deg_dat, direction = 'up'){
  sig_dat <- lapply(deg_dat, FUN = function(x){
    celltype_dat <- as.data.frame(x)
    if(direction == 'up'){
      sig_cell_dat <- celltype_dat[celltype_dat$p_val_ad < 0.01 & celltype_dat$avg_log2FC > 1,]
    }
    if(direction == 'down'){
      sig_cell_dat <- celltype_dat[celltype_dat$p_val_ad < 0.01 & celltype_dat$avg_log2FC < -1,]
    }
    sig_cell_dat
  })
}

wt_three_sig <- get_sig_degs(wt_three_degs, direction = 'up')
wt_five_sig <- get_sig_degs(wt_five_degs, direction = 'up')
ips_three_sig <- get_sig_degs(ips_three_degs, direction = 'up')
ips_five_sig <- get_sig_degs(ips_five_degs, direction = 'up')

#Get downregulated genes too
wt_three_sig_down <- get_sig_degs(wt_three_degs, direction = 'down')
wt_five_sig_down <- get_sig_degs(wt_five_degs, direction = 'down')
ips_three_sig_down <- get_sig_degs(ips_three_degs, direction = 'down')
ips_five_sig_down <- get_sig_degs(ips_five_degs, direction = 'down')


get_deg_count_df <- function(deg_sig_dat, label){
  deg_counts <- lapply(deg_sig_dat, FUN = function(x){
    nrow(x)
  })
  deg_counts_df <- as.data.frame(t(as.data.frame(deg_counts)))
  colnames(deg_counts_df) <- label
  deg_counts_df
}

#Counts of upregulated degs
wt_3_deg_counts <- get_deg_count_df(wt_three_sig, label = 'wt_3')
wt_5_deg_counts <- get_deg_count_df(wt_five_sig, label = 'wt_5')
ips_3_deg_counts <- get_deg_count_df(ips_three_sig, label = 'ips_3')
ips_5_deg_counts <- get_deg_count_df(ips_five_sig, label = 'ips_5')

#Counts of downregulated degs
wt_3_deg_count_down <- get_deg_count_df(wt_three_sig_down, label = 'wt_3')
wt_5_deg_count_down <- get_deg_count_df(wt_five_sig_down, label = 'wt_5')
ips_3_deg_count_down <- get_deg_count_df(ips_three_sig_down, label = 'ips_3')
ips_5_deg_count_down <- get_deg_count_df(ips_five_sig_down, label = 'ips_5')

deg_count_df <- cbind(wt_3_deg_counts, wt_5_deg_counts, ips_3_deg_counts, ips_5_deg_counts)
deg_count_down_df <- cbind(wt_3_deg_count_down, wt_5_deg_count_down, ips_3_deg_count_down, ips_5_deg_count_down)

deg_count_df <- deg_count_df %>% rownames_to_column(var = 'celltype') %>%   
  pivot_longer(cols = c(wt_3, wt_5, ips_3, ips_5), 
               names_to = 'comp', values_to = 'count')
deg_count_down_df <- deg_count_down_df %>% rownames_to_column(var = 'celltype') %>%   
  pivot_longer(cols = c(wt_3, wt_5, ips_3, ips_5), 
               names_to = 'comp', values_to = 'count')

#Prep data for plotting
deg_count_df$genotype <- substr(deg_count_df$comp, start = 1, stop = 2)
deg_count_df[deg_count_df$genotype == 'ip',]$genotype = 'ips'
deg_count_df$day <- str_sub(deg_count_df$comp,-1,-1)
deg_count_df$day <- paste('Day', deg_count_df$day)

deg_count_down_df$genotype <- substr(deg_count_down_df$comp, start = 1, stop = 2)
deg_count_down_df[deg_count_down_df$genotype == 'ip',]$genotype = 'ips'
deg_count_down_df$day <- str_sub(deg_count_down_df$comp,-1,-1)
deg_count_down_df$day <- paste('Day', deg_count_down_df$day)

ggplot(deg_count_df, aes(x = count, y = celltype, fill = genotype))+
  facet_wrap(~day)+
  geom_bar(stat = 'identity', position = 'dodge')+
  ggtitle('upregulated genes')

ggplot(deg_count_down_df, aes(x = count, y = celltype, fill = genotype))+
  facet_wrap(~day)+
  geom_bar(stat = 'identity', position = 'dodge')+
  ggtitle('downregulated genes')

#Pathways for significant defs
ips_chimeric_five_resident <- subset(ips_chimeric_five, manualAnnotation %in% c('Astrocytes', 'Choroid Plexus', 'Endothelial', 'Ependymal', 
                                                                                'Immature Neurons', 'Microglia', 'Muscle cells', 'Neurons',
                                                                                'Oligodendrocytes', 'Pericytes'))
wt_chimeric_five_resident <- subset(wt_chimeric_five, manualAnnotation %in% c('Astrocytes', 'Choroid Plexus', 'Endothelial', 'Ependymal', 
                                                                                'Immature Neurons', 'Microglia', 'Muscle cells', 'Neurons',
                                                                                'Oligodendrocytes', 'Pericytes'))

ips_5_resident_degs_up <- FindMarkers(ips_chimeric_five_resident, group.by = 'Treatment', ident.1 = 'rChLGTV', ident.2 = 'PBS', only.pos = TRUE,
                             test.use = 'MAST')
wt_5_resident_degs_up <- FindMarkers(wt_chimeric_five_resident, group.by = 'Treatment', ident.1 = 'rChLGTV', ident.2 = 'PBS', only.pos = TRUE,
                                      test.use = 'MAST')

ips_5_degs_up_sig <- ips_5_resident_degs_up[ips_5_resident_degs_up$p_val_ad < 0.01 & ips_5_resident_degs_up$avg_log2FC > 1,]
wt_5_degs_up_sig <- wt_5_resident_degs_up[wt_5_resident_degs_up$p_val_ad < 0.01 & wt_5_resident_degs_up$avg_log2FC > 1,]

ips_5_paths <- gprofiler2::gost(query = rownames(ips_5_degs_up_sig), organism = 'mmusculus', evcodes = TRUE, sources = c('GO:BP', 'KEGG'))
ips_5_paths$result

DotPlot(wt_chimeric, group.by = 'Treatment', features = 'Ifit3', scale = FALSE)

#Can make venn diagram of genes up in wt and ips at day 3 and 5
sum(rownames(wt_5_degs_up_sig) %in% rownames(ips_5_degs_up_sig)) / nrow(ips_5_degs_up_sig)

