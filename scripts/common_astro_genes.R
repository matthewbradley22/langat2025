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
sc_astros_ips <- subset(ParseSeuratObj_int, manualAnnotation == 'Astrocytes' & Genotype == 'IPS1' & Treatment != 'rLGTV')
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
sc_astro_markers <- FindAllMarkers(sc_astros, group.by = 'Treatment',test.use = 'MAST', only.pos = TRUE)
sc_astro_markers_chlgtv <- subset(sc_astro_markers, cluster == 'rChLGTV' & p_val_adj < 0.01 & avg_log2FC > 1)

sc_astro_markers_mavs <- FindAllMarkers(sc_astros_ips, group.by = 'Treatment',test.use = 'MAST', only.pos = TRUE)
sc_astro_markers_mavs_chlgtv <- subset(sc_astro_markers_mavs, cluster == 'rChLGTV' & p_val_adj < 0.01 & avg_log2FC > 1)

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
total_gene_overlap <- sum(sc_astro_markers_chlgtv$gene %in% bulk_sig_gene_symbols)
total_gene_overlap/nrow(sc_astro_markers_chlgtv) #45% sc genes in bulk
total_gene_overlap/length(bulk_sig_gene_symbols) #17% bulk genes in sc

VennDiagram::venn.diagram(x = list(sc_astro_markers_chlgtv$gene, bulk_sig_gene_symbols),
                          category.names = c('SC', 'bulk'),
                          filename = '~/Documents/ÖverbyLab/scPlots/bulk_sc_astro_venn.png',
                          resolution = 300,
                          fill = c("#B3E2CD", "#FDCDAC"),
                          cat.cex = 3,
                          cex = 3)

#Mavs
total_gene_overlap_mavs <- sum(sc_astro_markers_mavs_chlgtv$gene %in% bulk_sig_genes_mavs_symbols)
total_gene_overlap_mavs/nrow(sc_astro_markers_mavs_chlgtv) #45% sc genes in bulk
total_gene_overlap_mavs/length(bulk_sig_genes_mavs_symbols) #17% bulk genes in sc

VennDiagram::venn.diagram(x = list(sc_astro_markers_mavs_chlgtv$gene, bulk_sig_genes_mavs_symbols),
                          category.names = c('SC', 'bulk'),
                          filename = '~/Documents/ÖverbyLab/scPlots/bulk_sc_astro_mavs_venn.png',
                          resolution = 300,
                          fill = c("#B3E2CD", "#FDCDAC"),
                          cat.cex = 3,
                          cex = 3)

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


