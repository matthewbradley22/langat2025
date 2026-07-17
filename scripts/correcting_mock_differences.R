#Reading in mock degs generated in ips_vs_wt.R

mock_wt_vs_ips <- readRDS(file = "~/Documents/ÖverbyLab/mock_wt_vs_ips_dat.rds")

#Read in celltype degs by time and genotype
treatment_degs <-  readRDS("~/Documents/ÖverbyLab/single_cell_ISG_figures/isg_fig_plots/deg_counts_time_celltype.rds")

#Load in data for celltypes
ParseSeuratObj_int <- LoadSeuratRds("~/Documents/ÖverbyLab/data/FilteredRpcaIntegratedDatNoDoublets.rds") 

shared_genes = (lapply(unique(ParseSeuratObj_int$manualAnnotation), FUN = function(x){
  
  cur_mock_samples <- mock_wt_vs_ips[grep(x, names(mock_wt_vs_ips))]
  cur_treatment_samples <- treatment_degs[grep(x, names(treatment_degs))]
  
  subset_common_degs <- lapply(cur_treatment_samples, FUN = function(y){
    y[rownames(y) %in% unique(unname(unlist(lapply(cur_mock_samples, rownames)))),]
  })
  
  subset_common_degs
  
}))

names(shared_genes) = unique(ParseSeuratObj_int$manualAnnotation)
#Check if comparisons stay significant if only comparing to respsective mock 
#(so compare wt infected astros to wt mock if it is up in wt mock vs ips mock)
#We have, up to now, compared infected astros to combined mock samples

#Subset data for comparisons
chimeric_mock <- subset(ParseSeuratObj_int, Treatment != 'rLGTV')

#Change this from ips_vs_wt script to only include wt pbs for wt and opposite for ips
chimeric_mock_wt <- subset(chimeric_mock, Genotype == 'WT')
chimeric_mock_ips <- subset(chimeric_mock, Genotype == 'IPS1')

#Modifying funciton from ips_vs_wt.R
celltypes_of_interest = c('Astrocytes')
times = c('Day 3', 'Day 4', 'Day 5')

deg_counts_time_celltype_geno <- list()

for(i in 1:length(celltypes_of_interest)){
  print(paste("starting cell", celltypes_of_interest[i]))
  current_celltype <- celltypes_of_interest[[i]]
  current_cell_data_wt <- subset(chimeric_mock_wt, manualAnnotation == current_celltype)
  current_cell_data_ips <- subset(chimeric_mock_ips, manualAnnotation == current_celltype)
  for(j in 1:length(times)){
    print(paste("starting time", times[j]))
    cur_timepoint_wt <- subset(current_cell_data_wt, Timepoint == times[[j]] | Treatment == 'PBS')
    cur_timepoint_ips <- subset(current_cell_data_ips, Timepoint == times[[j]] | Treatment == 'PBS')
    
    #Need at least 3 of each condition to compare
    if(all(table(cur_timepoint_wt$Treatment) > 3)){
      wt_markers <- FindMarkers(cur_timepoint_wt, test.use = 'MAST', group.by = 'Treatment', ident.1 = 'rChLGTV')
      wt_markers_up_sig <- wt_markers[wt_markers$p_val_adj < 0.01 & (wt_markers$avg_log2FC) > 1,]
      wt_markers_down_sig <- wt_markers[wt_markers$p_val_adj < 0.01 & (wt_markers$avg_log2FC) < -1,]
      deg_counts_time_celltype_geno[[length(deg_counts_time_celltype_geno) + 1]] <- wt_markers_up_sig
      names(deg_counts_time_celltype_geno)[length(deg_counts_time_celltype_geno)] = paste(current_celltype, times[j], 'wt_up', sep = '_')
      deg_counts_time_celltype_geno[[length(deg_counts_time_celltype_geno) + 1]] <- wt_markers_down_sig
      names(deg_counts_time_celltype_geno)[length(deg_counts_time_celltype_geno)] = paste(current_celltype, times[j], 'wt_down', sep = '_')
    }
    
    
    if(all(table(cur_timepoint_ips$Treatment) > 3 )){
      ips_markers <- FindMarkers(cur_timepoint_ips, test.use = 'MAST', group.by = 'Treatment', ident.1 = 'rChLGTV')
      ips_markers_up_sig <- ips_markers[ips_markers$p_val_adj < 0.01 & (ips_markers$avg_log2FC) > 1,]
      ips_markers_down_sig <- ips_markers[ips_markers$p_val_adj < 0.01 & (ips_markers$avg_log2FC) < -1,]
      deg_counts_time_celltype_geno[[length(deg_counts_time_celltype_geno) + 1]] <- ips_markers_up_sig
      names(deg_counts_time_celltype_geno)[length(deg_counts_time_celltype_geno)] = paste(current_celltype, times[j], 'ips_up', sep = '_')
      deg_counts_time_celltype_geno[[length(deg_counts_time_celltype_geno) + 1]] <- ips_markers_down_sig
      names(deg_counts_time_celltype_geno)[length(deg_counts_time_celltype_geno)] = paste(current_celltype, times[j], 'ips_down', sep = '_')
    }
  }
}


#Check if degs remain after only comparing to same genotype
#So just compare shared genes from above and see if they still exist in deg_counts_time_celltype_geno
sig_gene_comp_astro <- lapply(names(shared_genes$Astrocytes), FUN = function(x){
  cur_shared_genes = rownames(shared_genes$Astrocytes[[x]])
  shared_genes_still_present = rownames(deg_counts_time_celltype_geno[[x]])[rownames(deg_counts_time_celltype_geno[[x]]) %in% 
                                                                              cur_shared_genes]
  
  no_longer_significant = sum(!cur_shared_genes %in% shared_genes_still_present)
  still_significant = length(shared_genes_still_present)
  
  data.frame('not_sig_now' = no_longer_significant, 'still_sig' = still_significant)
})

names(sig_gene_comp_astro) = names(shared_genes$Astrocytes)
sig_gene_comp_astro[grep('up', names(sig_gene_comp_astro))]

sig_gene_comp_cp <- lapply(names(shared_genes$`Choroid Plexus`), FUN = function(x){
  cur_shared_genes = rownames(shared_genes$`Choroid Plexus`[[x]])
  shared_genes_still_present = rownames(deg_counts_time_celltype_geno[[x]])[rownames(deg_counts_time_celltype_geno[[x]]) %in% 
                                                                              cur_shared_genes]
  
  no_longer_significant = sum(!cur_shared_genes %in% shared_genes_still_present)
  still_significant = length(shared_genes_still_present)
  
  data.frame('not_sig_now' = no_longer_significant, 'still_sig' = still_significant)
})

names(sig_gene_comp_cp) = names(shared_genes$`Choroid Plexus`)
sig_gene_comp_cp[grep('up', names(sig_gene_comp_cp))]

