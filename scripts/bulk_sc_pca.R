#What happens if we combine pseudobulk single cell and bulk data on one pca plot?

library(AnnotationDbi)
library(org.Hs.eg.db)

#Load bulk data in ifnTreatedAnalysis.R
resultsNames(dds)

plotPCA(vsd, intgroup=c("Treatment"))+
  theme_classic()+
  theme(text = element_text(size = 20))

#Load sc data
ParseSeuratObj_int <- LoadSeuratRds("~/Documents/ÖverbyLab/data/FilteredRpcaIntegratedDatNoDoublets.rds") 
astrocytes <- subset(ParseSeuratObj_int, manualAnnotation == 'Astrocytes' & Treatment != 'rLGTV')

#Create pseudobulk object
parse_pb <- AggregateExpression(astrocytes, assays = "RNA", return.seurat = T, group.by = c("orig.ident", 'Genotype', 'Treatment', 'Timepoint', 'Organ'))

pb_counts <- parse_pb[['RNA']]$counts
meta <- parse_pb[[]]

pb_dds <- DESeqDataSetFromMatrix(countData = pb_counts,
                              colData = meta,
                              design = ~ Genotype + Treatment + Timepoint)

pb_dds <- DESeq(pb_dds)
pb_vsd <- vst(pb_dds)

plotPCA(pb_vsd, intgroup=c("Timepoint"))+
  theme_classic()+
  theme(text = element_text(size = 20))

#Combine bulk data with pseudobulk data
bulk_counts = dds@assays@data@listData$counts

gene_symbols <- mapIds(
  org.Mm.eg.db,
  keys = rownames(bulk_counts),
  column = "SYMBOL",
  keytype = "ENSEMBL"
)

gene_symbols <- data.frame(ensembl = names(gene_symbols), symbol = gene_symbols)

bulk_counts_symbol <- bulk_counts %>% as.data.frame() %>% rownames_to_column(var = 'ensembl') %>% 
  dplyr::left_join(gene_symbols, by = 'ensembl')
genes_with_multiple <- names(table(bulk_counts_symbol$symbol)[table(bulk_counts_symbol$symbol) > 1])

bulk_counts_symbol <- bulk_counts_symbol %>% dplyr::filter(!symbol %in% genes_with_multiple & !is.na(symbol)) %>% dplyr::select(!ensembl)

rownames(bulk_counts_symbol) <- bulk_counts_symbol$symbol
bulk_counts_symbol$symbol = NULL

#Match pseudobulk rows w bulk data
#Genes in both
genes_in_both <- rownames(bulk_counts_symbol)[rownames(bulk_counts_symbol) %in% rownames(pb_counts)]
bulk_counts_symbol <- bulk_counts_symbol[genes_in_both,]
pb_counts <- pb_counts[genes_in_both,]

combined_counts <- cbind(bulk_counts_symbol, pb_counts)

dds$experiment <- 'bulk'
pb_dds$experiment <- 'single_cell'

#Fill in columns and prepare for rbind
colData(pb_dds) <- colData(pb_dds)[, c("Genotype", "Treatment", "Timepoint", "experiment")]
dds$Timepoint = factor('filler')
dds$Genotype = factor(dds$Genotype)
dds$Treatment = factor(dds$Treatment)
colData(dds) <- colData(dds)[, c("Genotype", "Treatment", "Timepoint", "experiment")]


combined_meta <- rbind(colData(dds), colData(pb_dds))

#Combined dds object
combined_dds <- DESeqDataSetFromMatrix(countData = combined_counts,
                                 colData = combined_meta,
                                 design = ~ Genotype + Treatment)


combined_dds <- DESeq(combined_dds)
combined_vsd <- vst(combined_dds)

plotPCA(combined_vsd, intgroup=c("experiment"))+
  theme_classic()+
  theme(text = element_text(size = 20))


#Take genes separating ifnb from other ifns and look at expression patters in single cell
#Get pc2 loadings
rv <- rowVars(assay(vsd))
select_500 <- order(rv, decreasing=TRUE)[seq_len(min(500, length(rv)))]
pca <- prcomp(t(assay(vsd)[select_500,]))
#Get top genes associated with pc1
top_genes <- pca$rotation %>%as.data.frame() %>%
  dplyr::arrange(desc(abs(PC2))) %>%
  rownames() %>% 
  head(n = 25)

top_genes_symbol <- ensembldb::select(EnsDb.Mmusculus.v79, keys= top_genes, keytype = "GENEID", columns = c("SYMBOL","GENEID")) 
vsd_min = min(assay(vsd))
treatment_levels = c('Mock', 'IFNa1', 'IFNa4', 'IFNa5', 'IFNa6', 'IFNa9', 'IFNa11', 'IFNB')

#Try with z score
z_scores_for_heatmap <- t(scale(t(assay(vsd)[top_genes,]), scale = FALSE))

z_scores_for_heatmap %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'GENEID') %>% 
  tidyr::pivot_longer(cols = !GENEID, names_to = 'sample', values_to = 'exp.level') %>% 
  dplyr::left_join(top_genes_symbol, by = c('GENEID')) %>% 
  tidyr::separate(col = sample, into = c('sample', 'sample_2', 'genotype', 'treatment', 'num'), sep = '-') %>% 
  # dplyr::mutate(exp.level = exp.level - vsd_min) %>% 
  dplyr::group_by(SYMBOL, treatment) %>% 
  dplyr::summarise(mean_exp = mean(exp.level)) %>% 
  dplyr::mutate(SYMBOL = factor(SYMBOL, levels = rev(top_genes_symbol$SYMBOL))) %>% 
  dplyr::mutate(treatment = factor(treatment, levels = treatment_levels)) %>% 
  ggplot(aes(x = treatment, y = SYMBOL, fill = mean_exp))+
  geom_tile()+
  scale_fill_gradientn(colours = c("#F57456","#FFB975","#FFEDB8","#9CCCFF","#4272C2"),
                       values = c(1.0,0.84,0.75, 0.3, 0))+
  theme_classic()+
  theme(text = element_text(size = 22),
        axis.text.x = element_text(angle = 75, hjust = 1))


#Using degList from ifnTreatedAnalysis.R which compared each group to mock
degList_no_rownames <- lapply(degList, FUN = function(x){
  x %>% rownames_to_column(var = 'GENEID')
})

deg_df <- do.call(rbind, degList_no_rownames) 
deg_df_top_genes <- dplyr::arrange(deg_df, compGroup, padj, desc(log2FoldChange)) %>% 
  dplyr::group_by(compGroup) %>% 
  dplyr::slice_head(n = 5)

ifn_b_genes <- deg_df %>% dplyr::filter(compGroup == 'Treatment_IFNB_vs_Mock') %>% 
  dplyr::pull(GENEID)

deg_counts <- table(deg_df$gene_id) 

#Get ifn beta up genes that are not up in any other groups

gene_symbols_ifnb <- mapIds(
  org.Mm.eg.db,
  keys = as.character(deg_df_top_genes$GENEID),
  column = "SYMBOL",
  keytype = "ENSEMBL"
)

gene_symbols_ifnb_df <- data.frame(GENEID = names(gene_symbols_ifnb), symbol = gene_symbols_ifnb) %>% distinct()

ifnb_z_score <- t(scale(t(assay(vsd)[gene_symbols_ifnb_df$GENEID,]), scale = FALSE))

ifnb_z_score %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'GENEID') %>% 
  dplyr::filter(GENEID %in% deg_df_top_genes$GENEID) %>% 
  tidyr::pivot_longer(cols = !GENEID, names_to = 'sample', values_to = 'exp.level') %>% 
  dplyr::left_join(gene_symbols_ifnb_df, by = c('GENEID')) %>% 
  tidyr::separate(col = sample, into = c('sample', 'sample_2', 'genotype', 'treatment', 'num'), sep = '-') %>% 
  # dplyr::mutate(exp.level = exp.level - vsd_min) %>% 
  dplyr::group_by(symbol, treatment) %>% 
  dplyr::summarise(mean_exp = mean(exp.level)) %>% 
  dplyr::mutate(treatment = factor(treatment, levels = treatment_levels)) %>% 
  ggplot(aes(x = treatment, y = symbol, fill = mean_exp))+
  geom_tile()+
  scale_fill_gradientn(colours = c("#F57456","#FFB975","#FFEDB8","#9CCCFF","#4272C2"),
                       values = c(1.0,0.84,0.75, 0.3, 0))+
  theme_classic()+
  theme(text = element_text(size = 22),
        axis.text.x = element_text(angle = 75, hjust = 1))




#Do sc cells show a pattern regarding Ifnb genes
DotPlot(astrocytes, features = c('Acod1', 'Cd69', 'Cxcl10', 'Gm12185', 'Ifit1bl1', 'Rsad2',
                                 'Tnfsf10'), group.by = 'Treatment', scale = FALSE)+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.7))



