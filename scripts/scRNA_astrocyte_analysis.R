library(Seurat)
library(RColorBrewer)
library(dplyr)
library(ggplot2)
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

#Get some other celltpyes for some deg comparisons
endothelial <- subset(ParseSeuratObj_int, manualAnnotation == 'Endothelial')
ependymal <- subset(ParseSeuratObj_int, manualAnnotation == 'Ependymal')
microglia <- subset(ParseSeuratObj_int, manualAnnotation == 'Microglia')

#Recluster astrocytes
astrocytes <- prepSeuratObj(astrocytes, use_all_genes = TRUE)
ElbowPlot(astrocytes, ndims = 40)
astrocytes <- prepUmapSeuratObj(astrocytes, nDims = 20, reductionName = 'astrocytes_umap', resolution_value = 0.8)

#UMAPs
DimPlot(astrocytes, reduction = 'astrocytes_umap')
astrocytes$infection_group <- ifelse(astrocytes$Treatment %in% c('rChLGTV', 'rLGTV'), 'infected', 'uninfected')
DimPlot(astrocytes, reduction = 'astrocytes_umap', group.by = 'infection_group')
DimPlot(astrocytes, reduction = 'astrocytes_umap', group.by = 'Genotype')
DimPlot(astrocytes, reduction = 'astrocytes_umap', group.by = 'Organ')
DimPlot(astrocytes, reduction = 'astrocytes_umap', group.by = 'Timepoint')
DimPlot(astrocytes, reduction = 'astrocytes_umap', group.by = 'geno_timepoint_treatment', cols = newCols)
table(astrocytes$infection_group, astrocytes$Genotype)

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

#########Look at number of degs by celltype############
#######################################################


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
######################################################

#Can look at degs in astrocytes in mock mice, the infected,
#then compare which genes are upregulated in infection but not mock
pbs <- subset(ParseSeuratObj_int, Treatment == 'PBS')
chlgtv <- subset(ParseSeuratObj_int, Treatment == 'rChLGTV')

#pbs astro markers
pbs_astro_markers <- FindMarkers(pbs, group.by = 'manualAnnotation', 
                                 ident.1 = 'Astrocytes')





