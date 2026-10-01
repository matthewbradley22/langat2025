library(stringr)
library(readr)
library(gprofiler2)
library(gridExtra)
library(tidyr)
library(ggrepel)
library(dplyr)
library(DESeq2)
library(tidyverse)

setwd('~/Documents/ÖverbyLab/ebba_24_data/')


# - - - - - - - - - - - - - 
#### Data preprocessing ####
# - - - - - - - - - - - - - 


# Problem with quotations in csv file, so need to read line by line and fix quotation marks
lines <- readLines("nunya_count.csv", encoding = "UTF-8")

#Remove any starting and ending quotes, and turn double quotes to single, to correct excel export effects
fixed <- ifelse(
  startsWith(lines, '"') & endsWith(lines, '"'),
  gsub('""', '"', substr(lines, 2, nchar(lines) - 1)),
  lines
)

#Rewrite fixed version of csv
writeLines(fixed, "nunya_count_fixed.csv", useBytes = TRUE)

nunya_count_fixed <- as.data.frame(read_table("nunya_count_fixed.csv",   skip = 1))
rownames(nunya_count_fixed) = nunya_count_fixed$Geneid

#Only interested in gene counts
nunya_count_fixed_counts = dplyr::select(nunya_count_fixed, !(Geneid:Length))
nunya_count_fixed_counts[is.na(nunya_count_fixed_counts)]  = 0

#Read in metadata
nunya_cond <- read_csv("nunya_cond.csv")

#Read in metadata to sample mapping, and combine map with metadata
map_nunyaid_to_sampleid <- read.csv("~/Documents/ÖverbyLab/ebba_24_data/map_nunyaid_to_sampleid.csv", row.names=NULL)
mapped_metadata <- left_join(nunya_cond, map_nunyaid_to_sampleid, by = 'nunyaid')
mapped_metadata <- dplyr::filter(mapped_metadata, !is.na(sampleid))

#Only keep columns present in metadata
cols_to_keep <- which(substr(colnames(nunya_count_fixed_counts), start = 1, stop = 7) %in% map_nunyaid_to_sampleid$sampleid)
nunya_counts_subset <- nunya_count_fixed_counts[,cols_to_keep]

#Reorder metadata to match column order
sample_order <- substr(colnames(nunya_counts_subset), start = 1, stop = 7)

#Which order should metadata be in
meta_order <- unlist(lapply(sample_order, FUN = function(x){
  which(mapped_metadata$sampleid == x)
}))

meta_correct_order <- mapped_metadata[meta_order,]

#Subset data by celltype for further analysis (when unsplit pc1 explains 83% of variance)
cell_split_metadata <- lapply(unique(meta_correct_order$celltype), FUN = function(x){
  celltype_meta <- dplyr::filter(meta_correct_order, celltype == x)
  celltype_meta$treatment_genotype = paste(celltype_meta$treatment, celltype_meta$genotype, sep = '_')
  celltype_meta
})

cell_split_counts <- lapply(cell_split_metadata, FUN = function(x){
  celltype_counts <- nunya_counts_subset[,substr(colnames(nunya_counts_subset), start = 1, stop = 7) %in% x$sampleid]
})

#Create count and metadata for each celltype
neuron_counts = cell_split_counts[[1]]
neuron_meta = cell_split_metadata[[1]]

astro_counts <- cell_split_counts[[2]]
astro_meta = cell_split_metadata[[2]]

micro_counts <- cell_split_counts[[3]]
micro_meta <- cell_split_metadata[[3]]

count_meta_groups <- list(neurons = list(neuron_counts, neuron_meta),
                          astrocytes = list(astro_counts, astro_meta),
                          microglia = list(micro_counts, micro_meta))


# - - - - - - - - - - - - -
#### Deseq  analysis ####
# - - - - - - - - - - - - - 


#Create dds objects

celltpye_dds <- lapply(count_meta_groups, FUN = function(x){
  print(unique(x[[2]]$celltype))
  
  dds_celltype <- DESeqDataSetFromMatrix(countData = x[[1]],
                                        colData = x[[2]], 
                                        design = ~treatment + genotype)
  dds <- DESeq(dds_celltype)
})

# Perform variance stabilizing transformation
vsd_neurons <- vst(celltpye_dds$neurons, blind = FALSE)
vsd_astros <- vst(celltpye_dds$astrocytes, blind = FALSE)
vsd_micro <- vst(celltpye_dds$microglia, blind = FALSE)

plotPCA(vsd_neurons, intgroup=c('treatment', 'genotype'))+
  scale_color_manual(values = c('#FFDEA8', '#FFC175', '#BAD49F', '#76B535',
                                '#C9EAFF', '#99D7FF', '#D9A0D5', '#B244AA'))
plotPCA(vsd_astros, intgroup=c('treatment', 'genotype'))+
  scale_color_manual(values = c('#FFDEA8', '#FFC175', '#BAD49F', '#76B535',
                                '#C9EAFF', '#99D7FF', '#D9A0D5', '#B244AA'))
plotPCA(vsd_micro, intgroup=c('treatment', 'genotype'))+
  scale_color_manual(values = c('#FFDEA8', '#FFC175', '#BAD49F', '#76B535',
                                '#C9EAFF', '#99D7FF', '#D9A0D5', '#B244AA'))

#Receptors of interest
receptors <- c('Axl', 'Havcr1', 'Havcr2', 'Timd4', 'Tyro3', 'Lrp8', 'Lrp1',
               'Lrp4', 'Cd209a', 'Cd209b', 'Cd209c', 'Itgb4', 'Ncam1', 'Hspa1a', 'Vim', 'Itgav',
               'Itgb3', 'Cldn1', 'Clec5a', 'Mrc1', 'Mer', 'Scarb1', 'Dag1') #'Hspa5',

#Convert gene ids to names
gene_names <-  gprofiler2::gconvert(nunya_count_fixed$Geneid, organism = 'mmusculus', target = "ENSG")
gene_names <- gene_names[c('input', 'name')]
receptor_names <- dplyr::filter(gene_names, name %in% receptors)
colnames(receptor_names) <- c('gene', 'gene_name')

#Look at receptor levels in each celltype
vsd_neurons

plot_list <- list(list(vsd_neurons, 'Bulk neurons'), list(vsd_micro, 'Bulk microglia'),
                  list(vsd_astros, 'Bulk astrocytes'))

receptor_heatmaps <- lapply(plot_list, FUN = function(x){
  assay(x[[1]])[receptor_names$gene,] %>% 
    as.data.frame() %>% 
    rownames_to_column(var = 'gene') %>% 
    tidyr::pivot_longer(starts_with('samp')) %>% 
    dplyr::left_join(receptor_names, by = 'gene') %>% 
    dplyr::mutate(name = substr(name, start = 1, stop = 7)) %>% 
    ggplot(aes(x = name, y = reorder(gene_name, value), fill = value))+
    geom_tile() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 12))+
    scale_fill_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white"),
                         values = c(1.0,0.7,0.4,0))+
    ylab('')+
    xlab('')+
    ggtitle(x[[2]])
})

receptor_heatmaps[[1]]
receptor_heatmaps[[2]]
receptor_heatmaps[[3]]

#Mock vs mock receptor expression by celltype
dds_all <- DESeqDataSetFromMatrix(nunya_counts_subset,
                                  colData = meta_correct_order, 
                                  design = ~treatment + genotype + celltype)

dds_all <- DESeq(dds_all)

meta_mock <- dplyr::filter(meta_correct_order, treatment == 'MOCK')

#Only keep mock columns
mock_cols_to_keep <- which(substr(colnames(nunya_count_fixed_counts), start = 1, stop = 7) %in% meta_mock$sampleid)
mock_counts <- nunya_count_fixed_counts[,mock_cols_to_keep]

#Make sure metadata is in right order
sample_order_mock <- substr(colnames(mock_counts), start = 1, stop = 7)

mock_meta_order <- unlist(lapply(sample_order_mock, FUN = function(x){
  which(meta_mock$sampleid == x)
}))

mock_meta_correct_order <- meta_mock[mock_meta_order,]

dds_mock <- DESeqDataSetFromMatrix(mock_counts,
                                  colData = mock_meta_correct_order, 
                                  design = ~genotype + celltype)

dds_mock <- DESeq(dds_mock)

vsd_mock <- vst(dds_mock, blind = FALSE)

#Plot receptors between mock of all celltypes
assay(vsd_mock)[receptor_names$gene,] %>% as.data.frame() %>% 
  rownames_to_column(var = 'gene') %>% tidyr::pivot_longer(starts_with('samp')) %>% 
  dplyr::left_join(receptor_names, by = 'gene') %>% 
  dplyr::mutate(sampleid = substr(name, start = 1, stop = 7)) %>% 
  dplyr::left_join(mock_meta_correct_order, by = 'sampleid') %>% 
  dplyr::group_by(celltype, gene_name) %>% 
  dplyr::summarise(mean_exp = mean(value)) %>% 
  ggplot(aes(x = celltype, y = reorder(gene_name, mean_exp), fill = mean_exp))+
  geom_tile()+
  scale_fill_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white"),
                       values = c(1.0,0.7,0.4,0))+
  theme_classic()+
  ylab('')+
  xlab('')+
  theme(axis.text = element_text(size = 14))

# #Any receptors separating bystander and LGTV?
# #Compare receptors up in infected lgtv samples with those in tbev
# lgtv_vs_bystander <- lapply(celltpye_dds, FUN = function(x){
#   lgtv_bystander <-  results(x, name = 'treatment_LGTV_vs_BYSTANDER')
#   
#   lgtv_bystander[receptor_names$input,] %>% as.data.frame() %>% 
#     dplyr::arrange(padj) %>% 
#     rownames_to_column(var = 'input') %>% 
#     dplyr::left_join(gene_names, by = 'input') %>% 
#     dplyr::mutate(comp = 'lgtv')
# })
# 
# tbe_vs_bystander <- lapply(celltpye_dds, FUN = function(x){
#   tbe_bystander <- results(x, name = 'treatment_TBEV_vs_BYSTANDER')
#   
#   tbe_bystander[receptor_names$input,] %>% as.data.frame() %>% 
#     dplyr::arrange(padj) %>% 
#     rownames_to_column(var = 'input') %>% 
#     dplyr::left_join(gene_names, by = 'input')%>% 
#     dplyr::mutate(comp = 'tbev')
# })
# 
# lgtv_vs_tbe_receptor_plot <- function(lgtv_dat, tbe_dat, main){
#   rbind(lgtv_dat, tbe_dat) %>% 
#     dplyr::filter(!is.na(padj)) %>% 
#     dplyr::select(padj, name, comp) %>% 
#     tidyr::pivot_wider(names_from = comp, values_from = padj) %>% 
#     ggplot(aes(x = -log10(lgtv), y = -log10(tbev), label = name))+
#     geom_point()+
#     geom_text_repel(color = 'red')+
#     ylab('TBEV -log10 pvalue')+
#     xlab('LGTV -log10 pvalue')+
#     xlim(0, 30)+
#     theme_classic()+
#     geom_vline(xintercept = -log10(0.01), linetype = 'dotted')+
#     geom_hline(yintercept = -log10(0.01), linetype = 'dotted') +
#     ggtitle(main)
# }
# 
# pdf('~/Documents/ÖverbyLab/ebba_24_data/plots/neuron_virus_vs_bystander.pdf', height = 6, width = 7)
# lgtv_vs_tbe_receptor_plot(lgtv_vs_bystander$neurons, tbe_vs_bystander$neurons, main = 'Neurons')
# dev.off()
# 
# pdf('~/Documents/ÖverbyLab/ebba_24_data/plots/astro_virus_vs_bystander.pdf', height = 6, width = 7)
# lgtv_vs_tbe_receptor_plot(lgtv_vs_bystander$astrocytes, tbe_vs_bystander$astrocytes, main = 'astrocytes')
# dev.off()
# 
# pdf('~/Documents/ÖverbyLab/ebba_24_data/plots/micro_virus_vs_bystander.pdf', height = 6, width = 7)
# lgtv_vs_tbe_receptor_plot(lgtv_vs_bystander$microglia, tbe_vs_bystander$microglia, main = 'microglia')
# dev.off()
# 

# - - - - - - - - - - - - - - - - - - - - - 
#### Deseq comparison between celltypes ####
# - - - - - - - - - - - - - - - - - - - - - 

pdf('~/Documents/ÖverbyLab/ebba_24_data/plots/astrocyte_vs_neurons.pdf', height = 6, width = 7)
results(dds_all, name = 'celltype_Neuron_vs_Astrocyte') %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'input') %>% 
  dplyr::left_join(gene_names, by = 'input') %>% 
  dplyr::filter(name %in% receptors & !is.na(padj)) %>% 
  dplyr::arrange(padj) %>% 
  ggplot(aes(x = -log10(padj), y = log2FoldChange, label = name))+
  geom_point()+
  geom_text_repel(color = 'red')+
  theme_classic() +
  geom_hline(yintercept = 0, linetype = 'dashed')+
  ggtitle('Neurons vs Astrocytes')
dev.off()

plotCounts(dds_all, 'ENSMUSG00000039542', intgroup = 'celltype')

pdf('~/Documents/ÖverbyLab/ebba_24_data/plots/neurons_vs_microglia.pdf', height = 6, width = 7)
results(dds_all, contrast = c('celltype', 'Neuron', 'Microglia')) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'input') %>% 
  dplyr::left_join(gene_names, by = 'input') %>% 
  dplyr::filter(name %in% receptors & !is.na(padj)) %>% 
  dplyr::arrange(padj) %>% 
  ggplot(aes(x = -log10(padj), y = log2FoldChange, label = name))+
  geom_point()+
  geom_text_repel(color = 'red')+
  theme_classic() +
  geom_hline(yintercept = 0, linetype = 'dashed')+
  ggtitle('Neurons vs Microglia')
dev.off()


