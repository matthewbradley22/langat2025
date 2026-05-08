#Look at genes upregulated in astrocytes upon infection across datasets
#Start with single cell data

#Packages and functions
library(Seurat)
library(UCell)
library(RColorBrewer)
library(VennDiagram)
source('~/Documents/ÖverbyLab//scripts/langatFunctions.R')

############ Loading in all data (single-cell, bulk, single-nuclei) ############ 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

## Load single-cell data
ParseSeuratObj_int <- LoadSeuratRds("~/Documents/ÖverbyLab/data/FilteredRpcaIntegratedDatNoDoublets.rds") 
ParseSeuratObj_int$hasVirus = ifelse(ParseSeuratObj_int$virusCountPAdj >= 10, 1, 0)

#Check data
newCols <-  c(brewer.pal(12, 'Paired'), '#99FFE6', '#CE99FF', '#18662E','#737272',  '#FF8AEF')
newCols[11] =  '#FF8AEF'
DimPlot(ParseSeuratObj_int, label = FALSE, group.by = 'manualAnnotation', reduction = 'umap.integrated',
        cols = newCols)

#Setting treatment to NOT equal lgtv, so we keep mock and chimeric
sc_astros <- subset(ParseSeuratObj_int, manualAnnotation == 'Astrocytes' & Genotype == 'WT' & Treatment != 'rLGTV')

#Going to split by timepoint for comparison as well, comparing day 3 wt infected to all pbs (same for ips)
sc_astros_day3 <- subset(ParseSeuratObj_int, manualAnnotation == 'Astrocytes' & Treatment != 'rLGTV' &
                           (Timepoint == 'Day 3' & Genotype == 'WT'| Treatment == 'PBS'))

sc_astros_day4 <- subset(ParseSeuratObj_int, manualAnnotation == 'Astrocytes' & Treatment != 'rLGTV' &
                            ((Timepoint == 'Day 4') & Genotype == 'WT'| Treatment == 'PBS'))

sc_astros_ips <- subset(ParseSeuratObj_int, manualAnnotation == 'Astrocytes' & Genotype == 'IPS1' & Treatment != 'rLGTV')
sc_astros_ips_day3 <- subset(ParseSeuratObj_int, manualAnnotation == 'Astrocytes' & Treatment != 'rLGTV' &
                           (Timepoint == 'Day 3' & Genotype == 'IPS1'| Treatment == 'PBS'))
sc_astros_ips_day4 <- subset(ParseSeuratObj_int, manualAnnotation == 'Astrocytes' & Treatment != 'rLGTV' &
                            ((Timepoint == 'Day 4') & Genotype == 'IPS1'| Treatment == 'PBS'))

VlnPlot(sc_astros, features = 'Socs1', group.by = 'Treatment')
VlnPlot(sc_astros_ips, features = 'Stat1', group.by = 'Treatment')

## Loading bulk data takes quite a few lines, so will just run bulk_astrocyte_heatmaps.R to load both datasets
#Just check that they're loaded
dds_mavs
dds_wt
plotPCA(vsd_mavs, intgroup=c("treatment_time"))
plotPCA(vsd_wt, intgroup=c("treatment_time"))

########## Get lists of DEGs for respective datasets ########## 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

#SC markers
sc_astro_markers <- FindAllMarkers(sc_astros, group.by = 'Treatment', test.use = 'MAST', only.pos = TRUE)
sc_astro_markers_chlgtv <- subset(sc_astro_markers, cluster == 'rChLGTV' & p_val_adj < 0.01 & avg_log2FC > 1)

sc_astro_markers_mavs <- FindAllMarkers(sc_astros_ips, group.by = 'Treatment',test.use = 'MAST', only.pos = TRUE)
sc_astro_markers_mavs_chlgtv <- subset(sc_astro_markers_mavs, cluster == 'rChLGTV' & p_val_adj < 0.01 & avg_log2FC > 1)

#Day 3 markers
sc_astro_wt_3_markers <- FindAllMarkers(sc_astros_day3, group.by = 'Treatment', test.use = 'MAST', only.pos = TRUE)
sc_astro_wt_3_markers_sig <- subset(sc_astro_wt_3_markers, cluster == 'rChLGTV' & p_val_adj < 0.01 & avg_log2FC > 1)

sc_astro_ips_3_markers <- FindAllMarkers(sc_astros_ips_day3, group.by = 'Treatment', test.use = 'MAST', only.pos = TRUE)
sc_astro_ips_3_markers_sig <- subset(sc_astro_ips_3_markers, cluster == 'rChLGTV' & p_val_adj < 0.01 & avg_log2FC > 1)

#Day 4 markers
sc_astro_wt_4_markers <- FindAllMarkers(sc_astros_day4, group.by = 'Treatment', test.use = 'MAST', only.pos = TRUE)
sc_astro_wt_4_markers_sig <- subset(sc_astro_wt_4_markers, cluster == 'rChLGTV' & p_val_adj < 0.01 & avg_log2FC > 1)

sc_astro_ips_4_markers <- FindAllMarkers(sc_astros_ips_day4, group.by = 'Treatment', test.use = 'MAST', only.pos = TRUE)
sc_astro_ips_4_markers_sig <- subset(sc_astro_ips_4_markers, cluster == 'rChLGTV' & p_val_adj < 0.01 & avg_log2FC > 1)

#Bulk markers

#Compare mock to average of others
#Changing dds formula for this comparison, just compare between treatment groups 
dds_wt <- DESeqDataSetFromTximport(txi_wt, metadata_wt, ~Treatment)
dds_wt <- DESeq(dds_wt)
resultsNames(dds_wt)

dds_wt_res <- results(dds_wt, name = "Treatment_chLGTV_vs_Mock")
dds_wt_res_sig <- dds_wt_res %>% as.data.frame() %>% dplyr::filter(padj < 0.01 & log2FoldChange > 1) %>% 
  dplyr::arrange(padj)
plotCounts(dds_wt, gene = 'ENSMUSG00000040033', intgroup = 'Treatment')
dds_wt_res_sig$GENEID = rownames(dds_wt_res_sig)

#Mavs
dds_mavs <- DESeqDataSetFromTximport(txi_mavs, metadata_mavs, ~Treatment)
dds_mavs <- DESeq(dds_mavs)
resultsNames(dds_mavs)

dds_mavs_res <- results(dds_mavs, name = "Treatment_chLGTV_vs_Mock")
dds_mavs_res_sig <- dds_mavs_res %>% as.data.frame() %>% dplyr::filter(padj < 0.01 & log2FoldChange > 1) %>% 
  dplyr::arrange(padj)
plotCounts(dds_mavs, gene = 'ENSMUSG00000054160', intgroup = 'Treatment')
dds_mavs_res_sig$GENEID = rownames(dds_mavs_res_sig)

########## Venn diagrams comparing single cell and bulk degs ########## 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

#Venn diagram, turn off log file saving
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

VennDiagram::venn.diagram(x = list(wt = rownames(dds_wt_res_sig), ips = rownames(dds_mavs_res_sig)),
                          category.names = c('WT', 'IPS'),
                          filename = '~/Documents/ÖverbyLab/single_cell_ISG_figures/sc_vs_bulk/bulk_ips_vs_wt.png',
                          resolution = 300,
                          fill = c("#B3E2CD", "#FDCDAC"),
                          cat.cex = 3,
                          cex = 3)

bulk_wt_vs_ips <- VennDiagram::venn.diagram(list(wt = rownames(dds_wt_res_sig), ips = rownames(dds_mavs_res_sig)), filename = NULL,
                                        fill = brewer.pal(3, "Pastel2")[1:2], cex = 1.5, cat.cex = 1.5)
grid::grid.draw(bulk_wt_vs_ips)

########## Comapre up and downregulated genes between datasets ########## 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

#Convert bulk deg ensembl names to symbols
bulk_sig_genes <- rownames(dds_wt_res_sig)
geneConversion <- ensembldb::select(EnsDb.Mmusculus.v79, keys= bulk_sig_genes, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
bulk_sig_gene_symbols <- geneConversion$SYMBOL
dds_wt_res_sig <- left_join(dds_wt_res_sig, geneConversion, by = c("GENEID"))

bulk_sig_genes_mavs <- rownames(dds_mavs_res_sig)
geneConversion_mavs <- ensembldb::select(EnsDb.Mmusculus.v79, keys= bulk_sig_genes_mavs, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
bulk_sig_genes_mavs_symbols <- geneConversion_mavs$SYMBOL
dds_mavs_res_sig <- left_join(dds_mavs_res_sig, geneConversion_mavs, by = c("GENEID"))

#Proportion of overlap between sc and bulk
total_gene_overlap <- sum(sc_astro_wt_3_markers_sig$gene %in% bulk_sig_gene_symbols)
total_gene_overlap/nrow(sc_astro_wt_3_markers_sig) #45% sc genes in bulk
total_gene_overlap/length(bulk_sig_gene_symbols) #14% bulk genes in sc


sc_day_3_vs_bulk_wt <- VennDiagram::venn.diagram(list(single_cell = sc_astro_wt_3_markers_sig$gene, bulk = bulk_sig_gene_symbols), filename = NULL,
                                            fill = brewer.pal(3, "Pastel2")[1:2], cex = 1.5, cat.cex = 1.5)
grid::grid.draw(sc_day_3_vs_bulk_wt)

#Save plot
VennDiagram::venn.diagram(x = list(sc_astro_wt_3_markers_sig$gene, bulk_sig_gene_symbols),
                          category.names = c('SC', 'bulk'),
                          filename = '~/Documents/ÖverbyLab/single_cell_ISG_figures/sc_vs_bulk/bulk_sc_astro_wt_3_venn.png',
                          resolution = 300,
                          fill = c("#B3E2CD", "#FDCDAC"),
                          cat.cex = 3,
                          cex = 3)

#Day 4 and 5 wt
sc_day_45_vs_bulk_wt <- VennDiagram::venn.diagram(list(single_cell = sc_astro_wt_45_markers_sig$gene, bulk = bulk_sig_gene_symbols), filename = NULL,
                                                 fill = brewer.pal(3, "Pastel2")[1:2], cex = 1.5, cat.cex = 1.5)
grid::grid.draw(sc_day_45_vs_bulk_wt)

VennDiagram::venn.diagram(x = list(sc_astro_wt_45_markers_sig$gene, bulk_sig_gene_symbols),
                          category.names = c('SC', 'bulk'),
                          filename = '~/Documents/ÖverbyLab/single_cell_ISG_figures/sc_vs_bulk/bulk_sc_astro_wt_45_venn.png',
                          resolution = 300,
                          fill = c("#B3E2CD", "#FDCDAC"),
                          cat.cex = 3,
                          cex = 3)

#Mavs
total_gene_overlap_mavs <- sum(sc_astro_ips_3_markers_sig$gene %in% bulk_sig_genes_mavs_symbols)
total_gene_overlap_mavs/nrow(sc_astro_markers_mavs_chlgtv) #3.3% sc genes in bulk
total_gene_overlap_mavs/length(bulk_sig_genes_mavs_symbols) #3.4% bulk genes in sc

#Mavs day 3
VennDiagram::venn.diagram(x = list(sc_astro_ips_3_markers_sig$gene, bulk_sig_genes_mavs_symbols),
                          category.names = c('SC', 'bulk'),
                          filename = '~/Documents/ÖverbyLab/single_cell_ISG_figures/sc_vs_bulk/bulk_sc_astro_ips_3_venn.png',
                          resolution = 300,
                          fill = c("#B3E2CD", "#FDCDAC"),
                          cat.cex = 3,
                          cex = 3)

#Mavs day 4/5

VennDiagram::venn.diagram(x = list(sc_astro_ips_45_markers_sig$gene, bulk_sig_genes_mavs_symbols),
                          category.names = c('SC', 'bulk'),
                          filename = '~/Documents/ÖverbyLab/single_cell_ISG_figures/sc_vs_bulk/bulk_sc_astro_ips_45_venn.png',
                          resolution = 300,
                          fill = c("#B3E2CD", "#FDCDAC"),
                          cat.cex = 3,
                          cex = 3)

#Upset plot to get multiple groups
#All day 3
all_3_upset <- UpSetR::fromList(list('sc_wt_3' = sc_astro_wt_3_markers_sig$gene, 'bulk_wt' = bulk_sig_gene_symbols, 
                                      'sc_ips_3' = sc_astro_ips_3_markers_sig$gene, 'bulk_ips' = bulk_sig_genes_mavs_symbols))

UpSetR::upset(all_3_upset, nsets = 4, order.by = 'freq', text.scale = 2, nintersects = 10)

#All day 4/5
all_45_upset <- UpSetR::fromList(list('sc_wt_45' = sc_astro_wt_45_markers_sig$gene, 'bulk_wt' = bulk_sig_gene_symbols, 
                                      'sc_ips_45' = sc_astro_ips_45_markers_sig$gene, 'bulk_ips' = bulk_sig_genes_mavs_symbols))

UpSetR::upset(all_45_upset, nsets = 4, order.by = 'freq', text.scale = 2, nintersects = 8)


#Look at pathways between shared genes
shared_genes <- sc_astro_markers_chlgtv[sc_astro_markers_chlgtv$gene %in% bulk_sig_gene_symbols,]$gene
shared_genes_paths <- gprofiler2::gost(shared_genes, organism = 'mmusculus', evcodes = TRUE, sources = c('GO:BP', 'KEGG'))
shared_genes_paths$result

shared_genes_paths$result[37,]

paths_for_barplot <- shared_genes_paths$result[c(2, 7, 12, 25, 41, 44, 49, 50),]
paths_for_barplot$term_name <- factor(paths_for_barplot$term_name, levels = paths_for_barplot$term_name[order(paths_for_barplot$p_value)])
ggplot(paths_for_barplot, aes(x = -log10(p_value), y = factor(term_name), fill = -log10(p_value)))+
  geom_bar(stat = 'identity')+
  xlab('')+
  ylab('')+
  geom_vline(xintercept=-log10(0.01), linetype="dotted", col ='white')+
  theme(axis.text = element_text(size = 16))

#Rank genes by avg p value across datasets
shared_genes_whole_df <- dplyr::inner_join(sc_astro_markers_chlgtv, dds_wt_res_sig, by = join_by(gene == SYMBOL))
#For some reason mean() returning weird results, just adding and dividing for now
shared_genes_whole_df <-  shared_genes_whole_df %>% dplyr::mutate(avg_p_value = (p_val_adj+padj)/2) %>% 
  dplyr::arrange(avg_p_value)
head(shared_genes_whole_df, n = 30)

#Same for MAVS
shared_genes_mavs <- sc_astro_markers_mavs_chlgtv[sc_astro_markers_mavs_chlgtv$gene %in% bulk_sig_genes_mavs_symbols,]$gene
shared_genes_mavs_paths <- gprofiler2::gost(shared_genes_mavs, organism = 'mmusculus', evcodes = TRUE, sources = c('GO:BP', 'KEGG'))
shared_genes_mavs_paths$result

shared_genes_mavs_paths$result$term_name <- factor(shared_genes_mavs_paths$result$term_name, 
                                                   levels = shared_genes_mavs_paths$result$term_name[order(shared_genes_mavs_paths$result$p_value)])
#Only plotting significant paths
ggplot(shared_genes_mavs_paths$result[1:3,], aes(x = -log10(p_value), y = factor(term_name), fill = -log10(p_value)))+
  geom_bar(stat = 'identity')+
  xlab('')+
  ylab('')+
  geom_vline(xintercept=-log10(0.01), linetype="dotted", col ='white')+
  theme(axis.text = element_text(size = 16))

########## Comapre MAVS and WT in bulk and sc respectively ########## 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
sum(dds_wt_res_sig$SYMBOL %in% dds_mavs_res_sig$SYMBOL)

VennDiagram::venn.diagram(x = list(bulk_sig_gene_symbols, bulk_sig_genes_mavs_symbols),
                          category.names = c('', ''),
                          filename = '~/Documents/ÖverbyLab/scPlots/bulk_mavs_vs_wt_venn.png',
                          resolution = 50,
                          fill = c("#B3E2CD", "#FDCDAC"),
                          cat.cex = 10,
                          cex = 15,width = 3000)

wt_mavs_overlap <- dds_wt_res_sig$SYMBOL[dds_wt_res_sig$SYMBOL %in% dds_mavs_res_sig$SYMBOL]
wt_mavs_overlap_paths <- gprofiler2::gost(wt_mavs_overlap, organism = 'mmusculus', evcodes = TRUE, sources = c('GO:BP', 'KEGG'))
wt_mavs_overlap_paths$result

#Plot select paths
wt_mavs_overlap_paths$result$term_name <- factor(wt_mavs_overlap_paths$result$term_name, 
                                                   levels = wt_mavs_overlap_paths$result$term_name[order(wt_mavs_overlap_paths$result$p_value)])

ggplot(wt_mavs_overlap_paths$result[c(1,3,5,6,9,11,12),], aes(x = -log10(p_value), y = factor(term_name), fill = -log10(p_value)))+
  geom_bar(stat = 'identity')+
  xlab('')+
  ylab('')+
  geom_vline(xintercept=-log10(0.01), linetype="dotted", col ='white')+
  theme(axis.text = element_text(size = 16))

########## Compare ISG expression between sc and bulk ########## 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## #

#Look at comparisons of just isgs
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

#single cell ISG degs
sc_wt_3_isgs <- sc_astro_wt_3_markers_sig[rownames(sc_astro_wt_3_markers_sig) %in% all_ISGs_type1,]
sc_ips_3_isgs <- sc_astro_ips_3_markers_sig[rownames(sc_astro_ips_3_markers_sig) %in% all_ISGs_type1,]
sc_wt_45_isgs <- sc_astro_wt_45_markers_sig[rownames(sc_astro_wt_45_markers_sig) %in% all_ISGs_type1,]
sc_ips_45_isgs <- sc_astro_ips_45_markers_sig[rownames(sc_astro_ips_45_markers_sig) %in% all_ISGs_type1,]

#Bulk ISG degs
wt_bulk_isgs <- bulk_sig_gene_symbols[bulk_sig_gene_symbols %in% all_ISGs_type1]
ips_bulk_isgs <- bulk_sig_genes_mavs_symbols[bulk_sig_genes_mavs_symbols %in% all_ISGs_type1]

#sc vs bulk isg expression
#Wt day 3
VennDiagram::venn.diagram(x = list(sc_wt_3_isgs$gene, wt_bulk_isgs),
                          category.names = c('SC', 'bulk'),
                          filename = '~/Documents/ÖverbyLab/single_cell_ISG_figures/sc_vs_bulk/bulk_sc_wt_3_isgs.png',
                          resolution = 300,
                          fill = c("#B3E2CD", "#FDCDAC"),
                          cat.cex = 3,
                          cex = 3)

#WT day 4/5
VennDiagram::venn.diagram(x = list(sc_wt_45_isgs$gene, wt_bulk_isgs),
                          category.names = c('SC', 'bulk'),
                          filename = '~/Documents/ÖverbyLab/single_cell_ISG_figures/sc_vs_bulk/bulk_sc_wt_45_isgs.png',
                          resolution = 300,
                          fill = c("#B3E2CD", "#FDCDAC"),
                          cat.cex = 3,
                          cex = 3)

#IPS day 3
VennDiagram::venn.diagram(x = list(sc_ips_3_isgs$gene, ips_bulk_isgs),
                          category.names = c('SC', 'bulk'),
                          filename = '~/Documents/ÖverbyLab/single_cell_ISG_figures/sc_vs_bulk/bulk_sc_ips_3_isgs.png',
                          resolution = 300,
                          fill = c("#B3E2CD", "#FDCDAC"),
                          cat.cex = 3,
                          cex = 3)

#IPSday 4/5
VennDiagram::venn.diagram(x = list(sc_ips_45_isgs$gene, ips_bulk_isgs),
                          category.names = c('SC', 'bulk'),
                          filename = '~/Documents/ÖverbyLab/single_cell_ISG_figures/sc_vs_bulk/bulk_sc_ips_45_isgs.png',
                          resolution = 300,
                          fill = c("#B3E2CD", "#FDCDAC"),
                          cat.cex = 3,
                          cex = 3)

#Do isgs go up over time in bulk?
wt_gene_names <- ensembldb::select(EnsDb.Mmusculus.v79, keys= rownames(dds_wt), keytype = "GENEID", columns = c("SYMBOL","GENEID"))
wt_gene_names_isgs <- wt_gene_names[wt_gene_names$SYMBOL %in% all_ISGs_type1,]
dds_wt$Timepoint <- factor(dds_wt$Timepoint)
dds_mavs$Timepoint <- factor(dds_mavs$Timepoint)

plotCounts(dds_wt, gene = 'ENSMUSG00000037523', intgroup = 'Timepoint') 
plotCounts(dds_mavs, gene = 'ENSMUSG00000031601', intgroup = 'Timepoint') 

#Compare sc to bulk at each timepoint
#Set dds equation to account for treatment and time
dds_wt <- DESeqDataSetFromTximport(txi_wt, metadata_wt, ~treatment_time)
dds_mavs <- DESeqDataSetFromTximport(txi_mavs, metadata_mavs, ~treatment_time)

#Want to compare everything to mock
dds_wt$treatment_time <- relevel(dds_wt$treatment_time, ref = 'Mock_0')
dds_mavs$treatment_time <- relevel(dds_mavs$treatment_time, ref = 'Mock_0')

dds_wt <- DESeq(dds_wt)
dds_mavs <- DESeq(dds_mavs)

resultsNames(dds_wt)
resultsNames(dds_mavs)

res_time_list <- list()
res_time_list_mavs <- list()

for(i in 2:4){
  cur_res <- results(dds_wt, name = resultsNames(dds_wt)[i])
  cur_res_mavs <- results(dds_mavs, name = resultsNames(dds_mavs)[i])
  
  cur_res_sig <- cur_res %>% as.data.frame() %>% dplyr::filter(padj < 0.01 & log2FoldChange > 1) %>% 
    dplyr::arrange(padj) %>% rownames_to_column(var = 'GENEID')
  cur_res_mavs_sig <- cur_res_mavs %>% as.data.frame() %>% dplyr::filter(padj < 0.01 & log2FoldChange > 1) %>% 
    dplyr::arrange(padj) %>% rownames_to_column(var = 'GENEID')
  
  geneConversion_cur_res <- ensembldb::select(EnsDb.Mmusculus.v79, keys= cur_res_sig$GENEID, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
  geneConversion_cur_res_mavs <- ensembldb::select(EnsDb.Mmusculus.v79, keys= cur_res_mavs_sig$GENEID, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
  
  cur_res_sig <- left_join(cur_res_sig, geneConversion_cur_res, by = c("GENEID"))
  cur_res_mavs_sig <- left_join(cur_res_mavs_sig, geneConversion_cur_res_mavs, by = c("GENEID"))
  
  res_time_list[[length(res_time_list) + 1]] = cur_res_sig
  res_time_list_mavs[[length(res_time_list_mavs) + 1]] = cur_res_mavs_sig
  
  names(res_time_list)[length(res_time_list)] = resultsNames(dds_wt)[i]
  names(res_time_list_mavs)[length(res_time_list_mavs)] = resultsNames(dds_mavs)[i]
}

lapply(res_time_list, nrow)
lapply(res_time_list_mavs, nrow)


wt_isg_by_time <- lapply(res_time_list, FUN = function(x){
  sig_isgs <- dplyr::filter(x, SYMBOL %in% all_ISGs_type1)
})
ips_isg_by_time <- lapply(res_time_list_mavs, FUN = function(x){
  sig_isgs <- dplyr::filter(x, SYMBOL %in% all_ISGs_type1)
})

lapply(wt_isg_by_time, nrow)
lapply(ips_isg_by_time, nrow)
