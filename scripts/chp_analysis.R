library(Seurat)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(tibble)
library(msigdbr)

source('~/Documents/ÖverbyLab//scripts/langatFunctions.R')

#Load data
ParseSeuratObj_int <- LoadSeuratRds("~/Documents/ÖverbyLab/data/FilteredRpcaIntegratedDatNoDoublets.rds") 

#Check data
newCols <-  c(brewer.pal(12, 'Paired'), '#99FFE6', '#CE99FF', '#18662E','#737272',  '#FF8AEF')
newCols[11] =  '#FF8AEF'
DimPlot(ParseSeuratObj_int, label = FALSE, group.by = 'manualAnnotation', reduction = 'umap.integrated',
        cols = newCols)

#Subset to only ChP cells
chp <- subset(ParseSeuratObj_int, manualAnnotation == 'Choroid Plexus')

chp <- prepSeuratObj(chp)
ElbowPlot(chp, ndims = 40)
chp <- prepUmapSeuratObj(chp, nDims = 15, reductionName = 'chp.umap', resolution_value = 0.5)

DimPlot(chp, label = FALSE, reduction = 'chp.umap',
        cols = newCols)

chp_clust_markers <- FindAllMarkers(chp, group.by = 'seurat_clusters', only.pos = TRUE, 
                                    test.use = 'MAST')


chp_clust_markers %>% dplyr::filter(p_val_adj < 0.01) %>% 
  dplyr::group_by(cluster) %>% 
  dplyr::slice_head(n = 10)


# - - - - - - - - - - - - 
#### Marker genes ####
# - - - - - - - - - - - - 


#This paper has adult marker genes
#https://www.cell.com/cell/fulltext/S0092-8674(21)00438-4

adult_markers <- c('Slc1a3', 'Ly86', 'P2ry6', 'Hmha1', 'Ptprc',
                   'Pecam1', 'Plvap', 'Flt1', 'Pear1', 'Cxcl12',
                   'Dcn', 'Adam12', 'Emilin1', 'Lum', 'Col3a1', 'Col1a2',
                   'Col1a1', 'Aqp1', 'Kcnj13', 'Enpp2', 'Htr2c')

DotPlot(chp, features = adult_markers, scale = FALSE)$data %>% 
  ggplot(aes(x = id, y = features.plot, fill = avg.exp.scaled, size = pct.exp))+
  geom_point(pch = 21)+
  theme_classic() +
  scale_fill_gradientn(colours = c('white', '#FFD991', '#FF7530', '#FF4024'), 
                       values = c(0, 0.3, 0.6, 1))+
  xlab('ChP Cluster')+
  ylab('')
  
  
#Markers from paper
#https://advanced.onlinelibrary.wiley.com/doi/full/10.1002/advs.202415984
gao_markers <- c('Ttr', 'Htr2c', 'Kenj13', 'Aqp1',
                 'Flt1', 'Plvap', 'Pecam1',
                 'Col1a1', 'Col1a2', 'Emilin1', 'Lum', 'Col18a1',
                 'Dcn', 'Cxcl12', 'Pdgfrb', 'Kenj8', 'Rgs5', 'Acta2', 'Tagln', 'Myh11', 'Ptprc', 'Ly86',
                 'Cd74')

DotPlot(chp, features = gao_markers, scale = FALSE)$data %>% 
  ggplot(aes(x = id, y = features.plot, fill = avg.exp.scaled, size = pct.exp))+
  geom_point(pch = 21)+
  theme_classic() +
  scale_fill_gradientn(colours = c('white', '#FFD991', '#FF7530', '#FF4024'), 
                       values = c(0, 0.3, 0.6, 1))+
  xlab('ChP Cluster')+
  ylab('')

#General markers from panglao db
epi_markers <- c('Krt14', 'Epcam', 'Cdh1')
endo_markers <- c('Cd93', 'Vwf', 'Emcn', 'Egfl7')

#Very low expression in general
DotPlot(chp, features = c(epi_markers, endo_markers), scale = FALSE)$data %>% 
  ggplot(aes(x = id, y = features.plot, fill = avg.exp.scaled, size = pct.exp))+
  geom_point(pch = 21)+
  theme_classic() +
  scale_fill_gradientn(colours = c('white', '#FFD991', '#FF7530', '#FF4024'), 
                       values = c(0, 0.3, 0.6, 1))

#Tight junction and adherens junction proteins
tj <- c('Cldn1', 'Cldn2', 'Cldn3', 'Cldn11', 'Tjp1', 'Ocln')

chp$genotype_treatment <- paste(chp$Genotype, chp$Treatment, sep = '_')
DotPlot(chp, features = tj, scale = FALSE, group.by = 'genotype_treatment')$data %>% 
  ggplot(aes(x = id, y = features.plot, fill = avg.exp.scaled, size = pct.exp))+
  geom_point(pch = 21)+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 75, hjust = 1))+
  scale_fill_gradientn(colours = c('white', '#FFD991', '#FF7530', '#FF4024'), 
                       values = c(0, 0.3, 0.6, 1))

aj <- c('Cdh1', 'Cdh2', 'Cdh3')

DotPlot(chp, features = aj, scale = FALSE, group.by = 'genotype_treatment')$data %>% 
  ggplot(aes(x = id, y = features.plot, fill = avg.exp.scaled, size = pct.exp))+
  geom_point(pch = 21)+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 75, hjust = 1))+
  scale_fill_gradientn(colours = c('white', '#FFD991', '#FF7530', '#FF4024'), 
                       values = c(0, 0.3, 0.6, 1))

DotPlot(chp, features = 'Jup', scale = FALSE, group.by = 'genotype_treatment')$data %>% 
  ggplot(aes(x = id, y = features.plot, fill = avg.exp.scaled, size = pct.exp))+
  geom_point(pch = 21)+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 75, hjust = 1))+
  scale_fill_gradientn(colours = c('white', '#FFD991', '#FF7530', '#FF4024'), 
                       values = c(0, 0.3, 0.6, 1))


# - - - - - - - - - - - - - - - -
#### Single nuclei ChP genes ####
# - - - - - - - - - - - - - - - 

#Read in processed data
sn_integrated_dat <- LoadSeuratRds('~/Documents/ÖverbyLab/single_nuclei_proj/LGTVscCombined.rds')

chp_sn <- subset(sn_integrated_dat, manualAnnotation == 'ChP')

chp_sn <- prepSeuratObj(chp_sn)
ElbowPlot(chp_sn, ndims = 40)
chp_sn <- prepUmapSeuratObj(chp_sn, nDims = 10, reductionName = 'chp.sn.umap', resolution_value = 0.5)
DimPlot(chp_sn, label = FALSE, reduction = 'chp.sn.umap',
        cols = newCols)

DotPlot(chp_sn, features = 'rna_LGTV', scale = FALSE)$data %>% 
  ggplot(aes(x = id, y = features.plot, fill = avg.exp.scaled, size = pct.exp))+
  geom_point(pch = 21)+
  theme_classic() +
  scale_fill_gradientn(colours = c('white', '#FFD991', '#FF7530', '#FF4024'), 
                       values = c(0, 0.3, 0.6, 1))+
  xlab('ChP Cluster')+
  ylab('')


# - - - - - - - - - - - - - - - - - - -
#### Try integrating the datasets ####
# - - - - - - - - - - - - - - - - - -

#only use wt for sc data
chp_wt <- subset(chp, Genotype == 'WT')

merged_chp <- merge(chp_wt, y = chp_sn, add.cell.ids = c("sc", "sn"))

merged_chp <- prepSeuratObj(merged_chp)
ElbowPlot(merged_chp, ndims = 40)
merged_chp <- prepUmapSeuratObj(merged_chp, nDims = 15, reductionName = 'chp.umap', resolution_value = 0.5)

sc_sn_int <- IntegrateLayers(object = merged_chp, method = CCAIntegration, orig.reduction = "pca",
                       new.reduction = "integrated.cca", verbose = FALSE, k.weight = 50)

sc_sn_int <- FindNeighbors(sc_sn_int, reduction = "integrated.cca", dims = 1:30)
sc_sn_int <- FindClusters(sc_sn_int)
sc_sn_int <- RunUMAP(sc_sn_int, reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.cca")
sc_sn_int$project <- stringr::str_sub(colnames(sc_sn_int), start = 1, end = 2) 

DimPlot(sc_sn_int, reduction = "umap.cca", group.by = 'project', label.size = 2)+
  xlab('')+
  ylab('')+
  scale_color_manual(values = c('#F09D46', '#3992F7'))+
  ggtitle('Single cell vs single nuclei ChP')


# - - - - - - - - - - - - - - - -
#### sc ChP cells vs others ####
# - - - - - - - - - - - - - - - -

#Object created in ips_vs_wt.R
deg_counts_by_time <- readRDS("~/Documents/ÖverbyLab/single_cell_ISG_figures/isg_fig_plots/deg_counts_time_celltype.rds")

deg_counts_select <- deg_counts_by_time[grep('(Astro)|(Endo)|(Choroid)', names(deg_counts_by_time))]
deg_counts_select <- deg_counts_select[grep('wt_up', names(deg_counts_select))]

lapply(deg_counts_select, nrow) %>% 
  as.data.frame() %>% t() %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  tidyr::separate(rowname, into = c('celltype', 'time', 'geno', 'direction'), sep = '_') %>% 
  ggplot(aes(x = celltype, fill = time, y = V1))+
  geom_bar(stat = 'identity', position = 'dodge', color = 'black')+
  theme_classic()+
  ylab('Number upregulated DEGs')+
  xlab('')+
  scale_fill_manual(values= c('#BBCEFC', '#7CA2FC', '#346AEB'))

#Plot counts of type 1 interferons
#Load in ISGs
#Molecular signatures database
#Should compare using db_species vs not using it
all_gene_sets <- msigdbr(species = "Mus musculus", db_species = "MM")

mouse_gene_sets <- all_gene_sets %>%
  split(x = .$gene_symbol, f = .$gs_name)

#Get total number of positive DEGs as above, but subset to ISGs
ifnA_response <- mouse_gene_sets$HALLMARK_INTERFERON_ALPHA_RESPONSE
ifnA_GOBP_response <- mouse_gene_sets$GOBP_RESPONSE_TO_INTERFERON_ALPHA
type1_response <- mouse_gene_sets$GOBP_RESPONSE_TO_TYPE_I_INTERFERON

all_ISGs_type1 = unique(c(ifnA_response, ifnA_GOBP_response, type1_response))

lapply(deg_counts_select, FUN = function(x){
  cur_dat <- x %>% as.data.frame() %>% 
    rownames_to_column(var = 'gene') %>% 
    dplyr::filter(gene %in% all_ISGs_type1)
  nrow(cur_dat)
}) %>% 
  as.data.frame() %>% t() %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  tidyr::separate(rowname, into = c('celltype', 'time', 'geno', 'direction'), sep = '_') %>% 
  ggplot(aes(x = celltype, fill = time, y = V1))+
  geom_bar(stat = 'identity', position = 'dodge', color = 'black')+
  theme_classic()+
  ylab('Number upregulated ISGs')+
  xlab('')+
  scale_fill_manual(values= c('#BBCEFC', '#7CA2FC', '#346AEB'))+
  ggtitle('Type 1 ISGs')

#Which DEGs are upregulated especially much in
#ChP cells at each timepoint compared to other cell types

deg_counts_select_gene <- lapply(deg_counts_select, FUN = function(x){
  x %>% as.data.frame() %>% tibble::rownames_to_column(var = 'gene')
})

deg_counts_select_df <- bind_rows(deg_counts_select_gene, .id = "cell_time") %>% 
  tidyr::separate(cell_time, into = c('celltype', 'time', 'genotype', 'direction'), sep = '_')

#Average log2FC of non chp cells, then subtract from
#ChP log2FC and see which genes have largest difference
deg_counts_select_df <- deg_counts_select_df %>% 
  dplyr::mutate(cell_broad = ifelse(celltype == 'Choroid Plexus',
                                    yes = 'ChP', no = 'other'))

times = c('Day 3', 'Day 4', 'Day 5')

gene_fc_diffs <- lapply(times, FUN = function(x){
  cur_dat = dplyr::filter(deg_counts_select_df, time == x)
  
  cur_dat %>% dplyr::group_by(cell_broad, gene) %>% 
    dplyr::summarise(mean_log2FC = mean(avg_log2FC)) %>% 
    tidyr::pivot_wider(names_from = cell_broad, values_from = mean_log2FC) %>% 
    dplyr::mutate(cell_diff = ChP - other) %>% 
    dplyr::arrange(desc(cell_diff))
})


#Which genes are only upregulated in ChP at a given time
day4_chp = deg_counts_select_df[deg_counts_select_df$cell_broad == 'ChP' &
                                  deg_counts_select_df$time == 'Day 4',]$gene 

day4_other= deg_counts_select_df[deg_counts_select_df$cell_broad == 'other' &
                                  deg_counts_select_df$time == 'Day 4',]$gene 

day4_chp_unique <- day4_chp[!day4_chp %in% day4_other]
day4_chp_unique_paths <- gprofiler2::gost(query = day4_chp_unique, 
                                          organism = 'mmusculus', 
                                          evcodes = TRUE, sources = c('GO:BP', 'KEGG'))
day4_chp_unique_paths$result[c('term_name', 'intersection')]

#Plvap look
FeaturePlot(ParseSeuratObj_int, features = 'Plvap', reduction = 'umap.integrated') 
VlnPlot(ParseSeuratObj_int, features = 'Plvap', group.by = 'manualAnnotation', pt.size = 0) 


# - - - - - - - - - - - - - - - - 
#### Check endothelial cells ####
# - - - - - - - - - - - - - - - - 

#Do any labelled endothlial cells have chp markers single-cell data
endothelial <- subset(ParseSeuratObj_int, manualAnnotation == 'Endothelial')

endothelial <- prepSeuratObj(endothelial)
ElbowPlot(endothelial, ndims = 40)
endothelial <- prepUmapSeuratObj(endothelial, nDims = 15, reductionName = 'endo.umap', resolution_value = 0.5)

DimPlot(endothelial, label = FALSE, reduction = 'endo.umap',
        cols = newCols)

DotPlot(endothelial, features = adult_markers, scale = FALSE)$data %>% 
  ggplot(aes(x = id, y = features.plot, fill = avg.exp.scaled, size = pct.exp))+
  geom_point(pch = 21)+
  theme_classic() +
  scale_fill_gradientn(colours = c('white', '#FFD991', '#FF7530', '#FF4024'), 
                       values = c(0, 0.3, 0.6, 1))+
  xlab('ChP Cluster')+
  ylab('')

#Chp associated markers
DotPlot(endothelial, features = c('Plvap','Ramp3', 'Plpp1', 'Esm1'), scale = FALSE)$data %>% 
  ggplot(aes(x = id, y = features.plot, fill = avg.exp.scaled, size = pct.exp))+
  geom_point(pch = 21)+
  theme_classic() +
  scale_fill_gradientn(colours = c('white', '#FFD991', '#FF7530', '#FF4024'), 
                       values = c(0, 0.3, 0.6, 1))+
  xlab('ChP Cluster')+
  ylab('')

#Non Chp associated markers
DotPlot(endothelial, features = c('Slc2a1', 'Mfsd2a', 'Ly6a', 'Cldn5'), scale = FALSE)$data %>% 
  ggplot(aes(x = id, y = features.plot, fill = avg.exp.scaled, size = pct.exp))+
  geom_point(pch = 21)+
  theme_classic() +
  scale_fill_gradientn(colours = c('white', '#FFD991', '#FF7530', '#FF4024'), 
                       values = c(0, 0.3, 0.6, 1))+
  xlab('ChP Cluster')+
  ylab('')

#Plvap featureplot
FeaturePlot(endothelial, features = 'Plvap',reduction = 'endo.umap')
