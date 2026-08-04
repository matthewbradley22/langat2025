#Functions script to source from
#packages
library(tibble)
library(dplyr)
library(DESeq2)
#Pseudobulk info found here
#https://hbctraining.github.io/Pseudobulk-for-scRNAseq/schedule/self-learning.html 

createPseudoBulk <- function(data, variables, intercept = TRUE){
  
  meta_columns <- variables
  meta <- data[[]] %>%
    dplyr::select(all_of(meta_columns)) %>%
    unique() %>%
    remove_rownames()
  
  bulk <- AggregateExpression(
    data,
    return.seurat = TRUE,
    assays = "RNA",
    group.by = meta_columns
  )
  
  n_cells <- data[[]] %>% 
    dplyr::group_by(across(all_of(meta_columns))) %>% 
    dplyr::summarise('n_cells' = n())
  
  #Need to make sure any underscores become dashes to count properly (think this is automatically done in AggregateExpression)
  n_cells <- data.frame(apply(n_cells, 2, function(x) if(!is.numeric(x)) gsub("_", "-", x)))
  meta_bulk <- left_join(bulk[[]], n_cells)
  rownames(meta_bulk) <- meta_bulk$orig.ident
  bulk[[]] <- meta_bulk
  
  # Get count matrix
  cluster_counts <- FetchData(bulk, layer="counts", vars=rownames(bulk))
  
  if(intercept){
    designForm <- reformulate(variables)
  }
  else{
    designForm <- reformulate(c('0',variables))
  }
  
  dds <- DESeqDataSetFromMatrix(t(cluster_counts),
                                colData = bulk@meta.data,
                                design = designForm)
}

#Reprocess a subset seurat object (rerun normalizing, scaling etc...)
prepSeuratObj <- function(obj, regress = FALSE, regressVars = NULL, use_all_genes = FALSE){
  #Rerun through data processing and visualization
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj)
  if(!regress){
    all.genes <- rownames(obj)
    if(use_all_genes){
      obj <- ScaleData(obj, features = all.genes)
    }
    else{
      obj <- ScaleData(obj)
    }
  }
  if(regress){
    obj <- ScaleData(obj, vars.to.regress = c(regressVars), features = rownames(obj))
  }
  obj <- RunPCA(obj)
  obj
}

prepUmapSeuratObj <- function(obj, nDims, reductionName, num_neighbors = 30L, resolution_value = 2){
  obj <- FindNeighbors(obj, dims = 1:nDims, reduction = "pca")
  obj <- FindClusters(obj, resolution = resolution_value)  
  obj <- RunUMAP(obj, dims = 1:nDims, reduction = "pca", reduction.name = reductionName, n.neighbors = num_neighbors)
  obj
}

#Create a dotplot of specific gene(s) in style used for figures
create_dot_plot <- function(dat, gene, main_title, facet = NULL, celltypes_to_plot = NULL, x_var = 'time', combine_pbs = FALSE,
                            flip_coords = FALSE){
  
  #Create dotplot data for the 4 variables we commonly plot
  dat$geno_treat_time_celltype <- paste(dat$Genotype,
                                        dat$Treatment,
                                        dat$Timepoint,
                                        dat$manualAnnotation,
                                        sep = '_')
  
  gene_exp <- DotPlot(dat, features = gene, group.by = 'geno_treat_time_celltype', scale = FALSE)$data
  
  #Separate id column so each variable is separated
  gene_exp <- gene_exp %>% tidyr::separate(col = id, into = c('geno', 'treatment', 'time', 'celltype'), sep = '_')
  gene_exp <- dplyr::filter(gene_exp, celltype != 'unknown')
  gene_exp$geno_treatment = paste(gene_exp$geno, gene_exp$treatment, sep = '_')
  
  #If we are grouping by treatment it can be helpful to average pbs day 3 and 5 
  if(combine_pbs){
    #If there is only one timepoint for infected, we could just use this
    average_expressions = gene_exp %>%  dplyr::group_by(treatment, celltype) %>% 
      dplyr::summarise(mean_exp = mean(avg.exp.scaled))
    
    #Making sure it works with multiple timepoints (don't want to average day3 and day5 infected for example)
    #We can just filter out one of the pbs groups at the end since they now have equal expressions
    gene_exp <- gene_exp %>% dplyr::left_join(average_expressions, by = c('treatment', 'celltype')) %>% 
       dplyr::mutate(avg.exp.scaled = case_when(treatment != 'PBS' ~ avg.exp.scaled,
                                                   treatment == 'PBS' ~ mean_exp)) %>% 
      dplyr::arrange(celltype) %>% 
      dplyr::filter(!(treatment == 'PBS' & time == 'Day 3'))
  }
  
  #Select certain celltypes
  if(!is.null(celltypes_to_plot)){
    gene_exp <- dplyr::filter(gene_exp, celltype %in% celltypes_to_plot)
  }
  
  if(length(unique(gene_exp$celltype)) == 1){
    custom_dot <- ggplot(gene_exp, aes(x = !!sym(x_var), y = treatment, size = pct.exp, fill = avg.exp.scaled))+
      geom_point(pch = 21)+
      theme_classic()+
      scale_fill_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white"), 
                           values = c(1.0,0.7,0.4,0))+
      ylab('')+
      xlab('')+
      ggtitle(main_title)
    
  } else{
    custom_dot <- ggplot(gene_exp, aes(x = !!sym(x_var), y = celltype, size = pct.exp, fill = avg.exp.scaled))+
      geom_point(pch = 21)+
      theme_classic()+
      scale_fill_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white"), 
                           values = c(1.0,0.7,0.4,0))+
      ylab('')+
      xlab('')+
      ggtitle(main_title)
  }
 
  if(!is.null(facet)){
    custom_dot = custom_dot +
      facet_wrap(vars(!!sym(facet)), scales = 'free')
  }
  
  if(flip_coords){
    custom_dot = custom_dot +
      coord_flip()+
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
  
  print(custom_dot)
}




