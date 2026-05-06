library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)
library(msigdbr)
source('~/Documents/ÖverbyLab//scripts/langatFunctions.R')

#Load data
ParseSeuratObj_int <- LoadSeuratRds("~/Documents/ÖverbyLab/data/FilteredRpcaIntegratedDatNoDoublets.rds") 

umap_color_list <- c( "#7047A1", "#B370AE","#292270",  "#166DF0","#6D92F8",  "#6DC3F8", "#8a0000","#F76363", "#FF96A2", "#D6644B", 
                      "#F08C3A", "#fdc087","#074F00", "#208d1f","#7bcd79", 
                      "gray")

#Not interested in langat samples right now
chimeric_mock <- subset(ParseSeuratObj_int, Treatment != 'rLGTV')
DimPlot(chimeric_mock, reduction = 'umap.integrated', cols = umap_color_list, group.by = 'manualAnnotation')

#Load ISGs for comparison to deg results
all_gene_sets <- msigdbr(species = "Mus musculus", db_species = "MM")

mouse_gene_sets <- all_gene_sets %>%
  split(x = .$gene_symbol, f = .$gs_name)

#Get total number of positive DEGs as above, but subset to ISGs
ifnA_response <- mouse_gene_sets$HALLMARK_INTERFERON_ALPHA_RESPONSE
ifnA_GOBP_response <- mouse_gene_sets$GOBP_RESPONSE_TO_INTERFERON_ALPHA
type1_response <- mouse_gene_sets$GOBP_RESPONSE_TO_TYPE_I_INTERFERON
all_ISGs_type1 = unique(c(ifnA_response, ifnA_GOBP_response, type1_response))

#Pseudobulk comparison by timepoint, celltype, and genotype
#Aggregate counts by group
bulk <- AggregateExpression(
  chimeric_mock,
  return.seurat = TRUE,
  assays = "RNA",
  group.by =  c("Genotype", "Treatment", "Timepoint", "manualAnnotation", "Well")
)

#Number of cells per group
n_cells <- chimeric_mock@meta.data %>% 
  dplyr::count(Genotype, Treatment, Timepoint, manualAnnotation, Well) %>% 
  rename("n"="n_cells")
meta_bulk <- left_join(bulk@meta.data, n_cells)
rownames(meta_bulk) <- meta_bulk$orig.ident
bulk@meta.data <- meta_bulk

#All cell, genotype, time combinations
all_var_combos <- expand.grid(unique(bulk$Genotype), unique(bulk$manualAnnotation), unique(bulk$Timepoint))
colnames(all_var_combos) <- c('genotype', 'celltype', 'timepoint')

#Create list to fill with for loop of degs from bulk comps
bulk_comps = list()

#For loop to get bulk degs from each celltype, timepoint, genotype combo (infected vs uninfected)
for(i in 1:nrow(all_var_combos)){
  cur_geno = all_var_combos[i,1]
  cur_celltype = all_var_combos[i,2]
  cur_time = all_var_combos[i,3]
  
  #Subset bulk data to current vars of interest
  cur_bulk_sample <- subset(bulk, subset= (manualAnnotation == cur_celltype & Timepoint == cur_time & Genotype == cur_geno))
  
  # Get count matrix
  cluster_counts <- FetchData(cur_bulk_sample, layer="counts", vars=rownames(cur_bulk_sample))
  
  # Create DESeq2 object
  # transpose it to get genes as rows
  dds <- DESeqDataSetFromMatrix(t(cluster_counts),
                                colData = cur_bulk_sample@meta.data,
                                design = ~ Treatment)
  
  dds <- DESeq(dds)
  
  
  # Results of Wald test
  res <- results(dds, 
                 name = "Treatment_rChLGTV_vs_PBS",
                 alpha = 0.05)
  
  bulk_comps[[length(bulk_comps) + 1]] = res
  
  #Do not retain any spaces in names, to make data wrangling easier later
  names(bulk_comps)[length(bulk_comps)] = paste(gsub(' ', '', cur_geno), 
                                                gsub(' ', '', cur_celltype), 
                                                gsub(' ', '', cur_time), sep = '_')
  print(paste('Done with', cur_geno, cur_celltype, cur_time))
}

sig_bulk_degs <- lapply(bulk_comps, FUN = function(x){
  x %>% as.data.frame() %>% dplyr::filter(padj < 0.01 & log2FoldChange > 0)
})

deg_count_by_geno <- data.frame(genotype = character(), celltype = character(), 
                                timepoint = character(), count = integer())

for(i in seq(1, length(sig_bulk_degs), 2)){
  #which samples are being looked at
  current_ips_sample <- names(sig_bulk_degs)[i]
  current_wt_sample <- names(sig_bulk_degs)[i+1]
  
  #Make sure celltypes match
  cur_ips_celltype <- stringr::str_split(current_ips_sample, '_', n = 3)[[1]][2]
  cur_wt_celltype <- stringr::str_split(current_wt_sample, '_', n = 3)[[1]][2]
  
  if(cur_ips_celltype != cur_wt_celltype){
    print('Celltypes do not match')
    break
  }
  
  cur_ips_time <- stringr::str_split(current_ips_sample, '_', n = 3)[[1]][3]
  cur_wt_time <- stringr::str_split(current_wt_sample, '_', n = 3)[[1]][3]
  
  if(cur_ips_time != cur_wt_time){
    print('times do not match')
    break
  }
  
  #Get names of degs in each sample
  ips_celltype_sample_degs <- rownames(sig_bulk_degs[[i]])
  wt_celltype_sample_degs <- rownames(sig_bulk_degs[[i+1]])
  
  #Which degs are in both or just one
  ips_only <- ips_celltype_sample_degs[!ips_celltype_sample_degs %in% wt_celltype_sample_degs]
  wt_only <- wt_celltype_sample_degs[!wt_celltype_sample_degs %in% ips_celltype_sample_degs]
  both <- ips_celltype_sample_degs[ips_celltype_sample_degs %in% wt_celltype_sample_degs]
  
  count_df <- data.frame(genotype = c('ips', 'wt', 'both'), celltype = cur_ips_celltype,
                         timepoint = cur_ips_time,
                         count = c(length(ips_only), length(wt_only),
                                   length(both)))
  
  deg_count_by_geno <- rbind(deg_count_by_geno, count_df)
}

deg_count_by_geno %>% dplyr::mutate(genotype = factor(genotype, levels = c('wt', 'both', 'ips'))) %>% 
ggplot(aes(x = genotype, y = celltype, fill = count))+
  geom_tile()+
  facet_wrap(~timepoint)+
  geom_text(aes(label=count), size = 4)+
  scale_fill_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white"), 
                       values = c(1.0,0.7,0.4,0))

