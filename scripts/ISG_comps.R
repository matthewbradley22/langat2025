library(Seurat)
library(RColorBrewer)
library(msigdbr)
library(dplyr)
library(matrixStats)
library(ggplot2)
library(gt)
library(ggrepel)
library(purrr)
library(forcats)

#Load data
ParseSeuratObj_int <- LoadSeuratRds("~/Documents/ÖverbyLab/data/FilteredRpcaIntegratedDatNoDoublets.rds") 
ParseSeuratObj_int$hasVirus = ifelse(ParseSeuratObj_int$virusCountPAdj >= 10, 1, 0)

#Create new columns for plotting and functions later later where we want things grouped
ParseSeuratObj_int$time_treatment <- paste(ParseSeuratObj_int$Treatment, ParseSeuratObj_int$Timepoint, sep = '_')
ParseSeuratObj_int$time_celltype <-  paste(ParseSeuratObj_int$Timepoint, ParseSeuratObj_int$manualAnnotation, sep = '_')
ParseSeuratObj_int$celltype_time <-  paste(ParseSeuratObj_int$manualAnnotation, ParseSeuratObj_int$Timepoint, sep = '_')
ParseSeuratObj_int$time_treatment_celltype <-  paste(ParseSeuratObj_int$Timepoint, ParseSeuratObj_int$Treatment, ParseSeuratObj_int$manualAnnotation, 
                                                     sep = '_')
 
#Check data
newCols <-  c(brewer.pal(12, 'Paired'), '#99FFE6', '#CE99FF', '#18662E','#737272',  '#FF8AEF')
newCols[11] =  '#FF8AEF'

pdf("~/Documents/ÖverbyLab/single_cell_ISG_figures/fig_1_plots/main_umap.pdf", height = 6, width = 9)
DimPlot(ParseSeuratObj_int, label = FALSE, group.by = 'manualAnnotation', reduction = 'umap.integrated',
        cols = newCols)+
  theme(axis.ticks = element_blank(),
        axis.text=element_blank(),
        legend.text=element_text(size=17))+
  xlab('Umap1')+
  ylab('Umap2')+
  guides(color=guide_legend(override.aes=list(size=8)))+
  ggtitle('')
dev.off()

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
#reactome_ifn_antiviral <- mouse_gene_sets$REACTOME_ANTIVIRAL_MECHANISM_BY_IFN_STIMULATED_GENES
all_ISGs_type1 = unique(c(ifnA_response, ifnA_GOBP_response, type1_response))

#Look at select isg in cells
ParseSeuratObj_int = AddModuleScore(ParseSeuratObj_int, features = list(all_ISGs_type1), name = 'ISG_score',
                                    slot = 'data', assay = 'RNA') 
chimeric_mock <- subset(ParseSeuratObj_int, Treatment != 'rLGTV')
chimeric_mock$genotype_treatment <- paste(chimeric_mock$Genotype, chimeric_mock$Treatment, sep = '_')

#Subset by genotype for plotting later
ips <- subset(ParseSeuratObj_int, Genotype == 'IPS1')
wt <- subset(ParseSeuratObj_int, Genotype == 'WT')
wt_cerebrum <- subset(wt, Organ == 'Cerebrum')
wt_cerebellum <- subset(wt, Organ == 'Cerebellum')

#Violin plot separated by genotype and treatment
chimeric_mock$genotype_treatment <- factor(chimeric_mock$genotype_treatment, levels = c('WT_PBS', 'WT_rChLGTV', 'IPS1_PBS', 'IPS1_rChLGTV'))
pdf('~/Documents/ÖverbyLab/scPlots/ISG_chimeric_vln_ips.pdf', width = 7, height = 4)
VlnPlot(chimeric_mock, features = 'ISG_score1', group.by = 'genotype_treatment', split.by = 'Timepoint',
        pt.size = 0)+
  ggtitle('ISG Module Scores')
dev.off()


#Plot ISGs
pdf('~/Documents/ÖverbyLab/scPlots/ISG_vln_ips.pdf', width = 7, height = 4)
VlnPlot(ips, features = 'ISG_score1', group.by = 'Treatment', split.by = 'Timepoint',
        pt.size = 0)+
  ggtitle('ISG Module Scores IPS1')
dev.off()

isg_data <- ParseSeuratObj_int[['RNA']]$data[which(rownames(ParseSeuratObj_int[['RNA']]$data) %in% all_ISGs_type1),]
isg_variances <- rowVars(as.matrix(isg_data))
high_var_ISGs <- rev(names(tail(sort(isg_variances), n = 20)))

pdf('~/Documents/ÖverbyLab/scPlots/ISG_dotplot.pdf', width = 9, height = 4)
DotPlot(ParseSeuratObj_int, group.by = 'time_treatment', features = high_var_ISGs, assay = 'RNA',
        scale = FALSE)+
  theme(axis.text.x = element_text(angle = 90))+
  ggtitle('High variance ISGs')
dev.off()

pdf('~/Documents/ÖverbyLab/scPlots/ISG_dotplot_ips.pdf', width = 9, height = 4)
DotPlot(ips, group.by = 'time_treatment', features = high_var_ISGs, assay = 'RNA',
        scale = FALSE)+
  theme(axis.text.x = element_text(angle = 90))+
  ggtitle('High variance ISGs - IPS1')
dev.off()

pdf('~/Documents/ÖverbyLab/scPlots/ISG_dotplot_wt.pdf', width = 9, height = 4)
DotPlot(wt, group.by = 'time_treatment', features = high_var_ISGs, assay = 'RNA',
        scale = FALSE)+
  theme(axis.text.x = element_text(angle = 90))+
  ggtitle('High variance ISGs - WT')
dev.off()

pdf('~/Documents/ÖverbyLab/scPlots/ISG_Umap.pdf', width = 7, height = 5)
FeaturePlot(ParseSeuratObj_int, features = 'ISG_score1', reduction = 'umap.integrated')
dev.off()

#Compare wt and IPS
ParseSeuratObj_int$time_treatment_genotype <- paste(ParseSeuratObj_int$time_treatment, ParseSeuratObj_int$Genotype, sep = '_')
DotPlot(subset(ParseSeuratObj_int, Treatment != 'PBS'), group.by = 'time_treatment_genotype', features = high_var_ISGs, assay = 'RNA',
        scale = FALSE)+
  theme(axis.text.x = element_text(angle = 90))+
  ggtitle('High variance ISGs')

#################################Groups#################################
#Split everything by group, look at celltypes in each group as well ie UMAPs/dotplots

#Create a table with y axis cell types, x axis days, and values as avg.expression of isg scores. split by virus and organ and genotype?
#Turn into function at some point
lgtv_cerebrum_wt <- subset(ParseSeuratObj_int, Organ == 'Cerebrum' & Treatment == 'rLGTV' & Genotype == 'WT')
lgtv_cerebrum_wt[[]] %>% dplyr::group_by(manualAnnotation, Timepoint) %>% 
  dplyr::summarise(mean_isg_score = mean(ISG_score1)) %>%
  pivot_wider(names_from = Timepoint, values_from = mean_isg_score) %>% dplyr::ungroup(manualAnnotation) %>% 
  gt(rowname_col = 'manualAnnotation') %>% 
  gtsave(filename = '~/Documents/ÖverbyLab/isg_table_test.pdf')

#Create dot plots of isg avg.expression
#But since its score of 200 genes, percent expressed isn't so important. Can do heatmap too

pdf("~/Documents/ÖverbyLab/scPlots/ISG_plots_byOrgan/wt_cerebrum_ISG_violin.pdf", width = 7, height = 4)
VlnPlot(wt_cerebrum, features = 'ISG_score1', group.by = 'Treatment', split.by = 'Timepoint',
        pt.size = 0)+
  ggtitle('ISG Module Scores WT Cerebrum')+
  scale_y_continuous(breaks = seq(0, 1, by = 0.2))
dev.off()

#Looks less intense than cerebrum on violin plot
pdf("~/Documents/ÖverbyLab/scPlots/ISG_plots_byOrgan/wt_cerebrellum_ISG_violin.pdf", width = 7, height = 4)
VlnPlot(wt_cerebellum, features = 'ISG_score1', group.by = 'Treatment', split.by = 'Timepoint',
        pt.size = 0)+
  ggtitle('ISG Module Scores WT Cerebellum')+
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(-0.25, 1))
dev.off()

#rank celltypes
rank_isg_scores_by_celltype <- function(dat){
  treatments <- c(unique(dat[['Treatment']]))[['Treatment']]
  #We don't care about pbs here, maybe as control later
  treatments <- treatments[treatments != 'PBS']
  for(i in 1:length(treatments)){
    dat_subset <- subset(dat, Treatment == treatments[i])
    dat_subset_dotplot <- DotPlot(dat_subset, features = 'ISG_score1', group.by = 'manualAnnotation', scale = FALSE)
    dat_subset_dotplot_data <- arrange(dat_subset_dotplot$data, desc(avg.exp))
    print(paste0("For treatment ", treatments[i], ":"))
    print(dat_subset_dotplot_data['id'])
  }
}

rank_isg_scores_by_celltype(wt_cerebrum)
rank_isg_scores_by_celltype(wt_cerebellum)
rank_isg_scores_by_celltype(ips_cerebrum)
rank_isg_scores_by_celltype(ips_cerebellum)

#Look at difference between wt organs (dotplot of subtraction)
left_join(wt_cerebrum_isg_plot_data, wt_cerebellum_isg_plot_data, by = 'id', suffix = c('.cerebrum',
                                                                                        '.cerebellum')) %>% 
  dplyr::mutate(organ_difference = avg.exp.cerebrum - avg.exp.cerebellum) %>% 
  tidyr::separate(id, c('ignore_that' ,'time', 'celltype')) %>% 
  mutate(celltype = replace(celltype, celltype == "Macrophage", "Macrophage/Monocytes")) %>% 
  dplyr::filter(celltype != 'Immature') %>% 
  ggplot(aes(x = time, y = celltype, size = 3, color = organ_difference))+
  geom_point()+
  xlab('Day')+
  ylab('Cell type')+
  ggtitle('WT Cerebrum Cerebellum DIfference ISG Scores')+
  scale_colour_gradient2(low = "blue", high = "red")+
  guides(size = "none")

data.frame(cerebrum_exp = wt_cerebrum_isg_plot_data$avg.exp, cerebellum_exp = wt_cerebellum_isg_plot_data$avg.exp)

#IPS groups
ips_cerebrum <- subset(ips, Organ == 'Cerebrum')

pdf("~/Documents/ÖverbyLab/scPlots/ISG_plots_byOrgan/ips_cerebrum_ISG_violin.pdf", width = 7, height = 4)
VlnPlot(ips_cerebrum, features = 'ISG_score1', group.by = 'Treatment', split.by = 'Timepoint',
        pt.size = 0)+
  ggtitle('ISG Module Scores IPS1 Cerebrum')+
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(-0.25, 1))
dev.off()

ips_cerebrum_isg_plot <- DotPlot(ips_cerebrum, features = 'ISG_score1', group.by = 'time_celltype', scale = FALSE) #Can split this by day too at some point, or make heatmap but that's annoying.
ips_cerebrum_isg_plot_data <- ips_cerebrum_isg_plot$data

ips_cerebellum <- subset(ips, Organ == 'Cerebellum')

pdf("~/Documents/ÖverbyLab/scPlots/ISG_plots_byOrgan/ips_cerebellum_ISG_violin.pdf", width = 7, height = 4)
VlnPlot(ips_cerebellum, features = 'ISG_score1', group.by = 'Treatment', split.by = 'Timepoint',
        pt.size = 0)+
  ggtitle('ISG Module Scores IPS1 Cerebellum')+
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(-0.25, 1))
dev.off()

ips_cerebellum_isg_plot <- DotPlot(ips_cerebellum, features = 'ISG_score1', group.by = 'time_celltype', scale = FALSE) #Can split this by day too at some point, or make heatmap but that's annoying.
ips_cerebellum_isg_plot_data <- ips_cerebellum_isg_plot$data

##########################ISG Dotplots##########################

#Will order heatmap by resident and infiltrating
resident_celltypes <- rev(c('Astrocytes', 'Choroid Plexus', 'Endothelial', 'Ependymal', 'Immature Neurons',
                        'Microglia', 'Muscle cells','Neurons', 'Oligodendrocytes', 'Pericytes'))


infiltrating_celltypes <- rev(c( 'B cells',  'Granulocytes', 'Macrophage/Monocytes', 'Nk cells',  "T cells"))

#Need a function to plug in data , return dotplot splot by day and celltype for each virus type (3 dot plots)
#Color and size based on isg scores
isg_heatmap_create <- function(dat , main, file_name_start = NULL, returnData = FALSE, infected_vs_mock = FALSE,
                               resident_only = FALSE){
  
  #There must be a better way to extract name of each treatment group from a seurat object, but 
  #for now this works fine, just hard to read
  treatments <- c(unique(dat[['Treatment']]))[['Treatment']]
  
  if(returnData){
    data_list = list()
  }
  
  if(!infected_vs_mock){
  #Create dotplot for each treatment
  #2 to 3 plots total 
  for(i in 1:length(treatments)){
    treatment_subset_dat <-subset(dat, Treatment == treatments[i])
    isg_plot <- DotPlot(treatment_subset_dat, features = 'ISG_score1', group.by = 'time_celltype', scale = FALSE)
    isg_plot_data <- isg_plot$data
    if(returnData){
      data_list[[i]] = isg_plot_data
      next
    }
    pdf(paste0("~/Documents/ÖverbyLab/scPlots/ISG_heatmaps_by_genotype/", file_name_start, treatments[i], '_ISG_heatmap.pdf'), width = 6, height = 5)
    plot = isg_plot_data %>% tidyr::separate(id, c('ignore_that' ,'time', 'celltype')) %>% 
      mutate(celltype = case_when(celltype == "Macrophage" ~ "Macrophage/Monocytes",
                                  celltype == "B" ~ "B cells",
                                  celltype == "T" ~ "T cells",
                                  celltype == "Nk" ~ "Nk cells",
                                  celltype == "Immature" ~ "Immature Neurons",
                                  celltype == "Choroid" ~ "Choroid Plexus",
                                  celltype == "Muscle" ~ "Muscle cells",
                                  .default = celltype))  %>% 
      mutate(celltype=forcats::fct_relevel(celltype,c(resident_celltypes, infiltrating_celltypes))) %>% 
      ggplot(aes(x = time, y = celltype, fill = avg.exp))+
      geom_tile()+
      xlab('Day')+
      ylab('Cell type')+
      ggtitle(paste0(main, ' ', treatments[i]))+
      #Make limits argument for function
      scale_fill_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white","lightblue"),
                            values = c(1.0,0.7,0.4,0.2,0),
                           limits = c(-0.15, 0.95))+
      geom_text(hjust=0.5, vjust=0.5, aes(label = round(avg.exp, digits = 2)), color = 'black', size = 4)+
      theme_classic()+
      theme(axis.line=element_blank())
    print(plot)
    dev.off()
    
    #Can replace heatmap with dotplot
    # isg_plot_data %>% tidyr::separate(id, c('ignore_that' ,'time', 'celltype')) %>% 
    #   mutate(celltype = replace(celltype, celltype == "Macrophage", "Macrophage/Monocytes"))  %>% 
    #   ggplot(aes(x = time, y = celltype, size = pct.exp, color = avg.exp))+
    #   geom_point()+
    #   xlab('Day')+
    #   ylab('Cell type')+
    #   ggtitle(paste0(main, ' ', treatments[i]))+
    #   #Make limits argument for function
    #   scale_colour_gradient2(low = "blue", high = "red", limits = c(-0.1, 0.95))+
    #   geom_text(hjust=0.5, vjust=0.5, aes(label = round(avg.exp, digits = 2)), color = 'black', size = 4)
    
    }
  }
  
  if(returnData){
    names(data_list) = treatments
    data_list
  }
  
  if(infected_vs_mock){
    #Get average expression for isgs of both groups then subtract
    infected_subset_dat <-subset(dat, Treatment == 'rChLGTV')
    mock_subset_dat <- subset(dat, Treatment == 'PBS')
    
    infected_isg_plot <- DotPlot(infected_subset_dat, features = 'ISG_score1', group.by = 'time_celltype', scale = FALSE)
    mock_isg_plot <- DotPlot(mock_subset_dat, features = 'ISG_score1', group.by = 'manualAnnotation', scale = FALSE)
    
    #subtract average celltype pbs isg score from all infected, rather then going day by day (no longer day 4 - day 4, just day 4 infected
    # - average mock)
    infected_isg_plot_dat <- infected_isg_plot$data
    mock_isg_plot_dat <- mock_isg_plot$data
    colnames(mock_isg_plot_dat)[4] = 'celltype'
    
    #Extract celltype from infected data
    infected_isg_plot_dat <- tidyr::extract(infected_isg_plot_dat, col = id, regex = "(.+)_(.+)", into = c('time', 'celltype'))
    
    combined_isg_scores <- dplyr::left_join(mock_isg_plot_dat, infected_isg_plot_dat, by = 'celltype', suffix = c('_mock', '_infected'))
    combined_isg_scores$inf_diff <- combined_isg_scores$avg.exp_infected - combined_isg_scores$avg.exp_mock
    
    plot_infected_vs_mock_dat = combined_isg_scores %>% 
      mutate(celltype = case_when(celltype == "B Cells" ~ "B cells",
                                  .default = celltype))  %>% 
      mutate(celltype=forcats::fct_relevel(celltype,c(resident_celltypes)))#, infiltrating_celltypes)))
    
    if(resident_only){
      plot_infected_vs_mock_dat <- plot_infected_vs_mock_dat[plot_infected_vs_mock_dat$celltype %in% resident_celltypes,]
    }
    
    pdf(paste0("~/Documents/ÖverbyLab/scPlots/ISG_heatmaps_by_genotype/", file_name_start, 'infected_vs_mock_ISG_heatmap.pdf'), width = 8, height = 5)
    plot_infected_vs_mock <- 
      ggplot(plot_infected_vs_mock_dat, aes(x = time, y = celltype, fill = inf_diff))+
      geom_tile()+
      xlab('Day')+
      ylab('Cell type')+
      ggtitle(main)+
      #Make limits argument for function
      scale_fill_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white"),
                           values = c(1.0,0.7,0.4,0),
                           limits = c(-0.02, 0.95))+
      geom_text(hjust=0.5, vjust=0.5, aes(label = round(inf_diff, digits = 2)), color = 'black', size = 4)+
      theme_classic()+
      theme(axis.line=element_blank(),
            legend.title = element_text(size = 15),
            axis.text = element_text(size = 17),
            plot.title = element_text(size = 17))+
      guides(fill=guide_legend(title="isg exp difference"))
      print(plot_infected_vs_mock)
      dev.off()
    }
}

#Warnings expected here since we're splitting data columns roughly right now. WIll fix for final plotting
#wt cerebrum
isg_heatmap_create(wt_cerebrum, main = "WT Cerebrum ISG Scores", file_name_start = 'wt_cerebrum_', infected_vs_mock = TRUE,
                   resident_only = TRUE)

# organs combined
isg_heatmap_create(wt, main = "WT ISG Scores", file_name_start = 'wt_')
isg_heatmap_create(ips, main = "IPS ISG Scores", file_name_start = 'ips_')

#wrt cerebellum
isg_heatmap_create(wt_cerebellum, main = "WT Cerebellum ISG Scores", file_name_start = 'wt_cerebellum_')

#ips cerebrum
isg_heatmap_create(ips_cerebrum, main = "IPS Cerebrum ISG Scores", file_name_start = 'ips_cerebrum_')
table(subset(ips_cerebrum, Treatment == 'rLGTV')$Timepoint, subset(ips_cerebrum, Treatment == 'rLGTV')$manualAnnotation) %>% 
  as.data.frame() %>% pivot_wider(names_from = Var1, values_from = Freq) %>% map_df(rev)
table(subset(ips_cerebrum, Treatment == 'rChLGTV')$Timepoint, subset(ips_cerebrum, Treatment == 'rLGTV')$manualAnnotation) %>% 
  as.data.frame() %>% pivot_wider(names_from = Var1, values_from = Freq) %>% map_df(rev)

#ips cerebellum
isg_heatmap_create(ips_cerebellum, main = "IPS Cerebellum ISG Scores", file_name_start = 'ips_cerebellum_')
table(subset(ips_cerebellum, Treatment == 'rChLGTV')$Timepoint, subset(ips_cerebellum, Treatment == 'rChLGTV')$manualAnnotation) %>% 
  as.data.frame() %>% pivot_wider(names_from = Var1, values_from = Freq) %>% map_df(rev)

#wt both organs
#Annoying column names, but Var1 = day, Var2 = cell type
chimeric_mock_ips <- subset(chimeric_mock, Genotype == 'IPS1')
chimeric_mock_wt <- subset(chimeric_mock, Genotype == 'WT')
isg_heatmap_create(chimeric_mock_ips, main = "IPS ISG Score", file_name_start = 'ips_isg_')
isg_heatmap_create(chimeric_mock_wt, main = "WT ISG Score", file_name_start = 'wt_isg_')

#Both organs, isg score infected vs mock
isg_heatmap_create(chimeric_mock_ips, main = "IPS ISG Infected vs Mock", file_name_start = 'ips_', infected_vs_mock = TRUE, resident_only = TRUE)
isg_heatmap_create(chimeric_mock_wt, main = "WT ISG Infected vs Mock", file_name_start = 'wt_', infected_vs_mock = TRUE, resident_only = TRUE)

#Difference between ips organ groups
left_join(ips_cerebrum_isg_plot_data, ips_cerebellum_isg_plot_data, by = 'id', suffix = c('.cerebrum',
                                                                                        '.cerebellum')) %>% 
  dplyr::mutate(organ_difference = avg.exp.cerebrum - avg.exp.cerebellum) %>% 
  tidyr::separate(id, c('ignore_that' ,'time', 'celltype')) %>% 
  mutate(celltype = replace(celltype, celltype == "Macrophage", "Macrophage/Monocytes")) %>% 
  dplyr::filter(celltype != 'Immature') %>% 
  ggplot(aes(x = time, y = celltype, size = 3, color = organ_difference))+
  geom_point()+
  xlab('Day')+
  ylab('Cell type')+
  ggtitle('WT Cerebrum Cerebellum DIfference ISG Scores')+
  scale_colour_gradient2(low = "blue", high = "red")+
  guides(size = "none")

#Can try scatter plot with y axis isg score, x axis day, color by celltype. Only for top ~5 celltypes for viewability
#Do not have percent infected here
wt_cerebrum_day_celltype_ISGs <- isg_dotplot_create(wt_cerebrum, main = '', returnData = TRUE)

#2nd in list is lgtv, 3rd is chimeric
wt_cerebrum_day_celltype_ISGs[[3]] %>% tidyr::separate(id, c('ignore_that' ,'day', 'celltype')) %>% 
  mutate(celltype = replace(celltype, celltype == "Macrophage", "Macrophage/Monocytes")) %>% 
  ggplot(aes(x = day, y = avg.exp, fill = celltype, colour = celltype))+
  geom_point()+
  geom_line(aes(group = celltype))+ 
  scale_colour_manual(values = newCols)+
  ggtitle("rChLGTV ISG Scores")

#Can look at which ISGs most expressed by celltype
#Average expression of each isg per celltype and timepoint?
wt_cerebrum[['RNA']]$data
isg_evg_exp <- AverageExpression(wt_cerebrum, assays = 'RNA', features = all_ISGs_type1, group.by = 'manualAnnotation')
for(i in 1:length(colnames(isg_evg_exp$RNA))){
  celltype_data <- isg_evg_exp$RNA[,i]
  print(colnames(isg_evg_exp$RNA)[i])
  print(tail(sort(celltype_data)))
}
isg_evg_exp

#Look at a few isg markers without addModuleScore and see celtype differences in lgtv vs pbs
#Could make dotplot with gene on x axis, celltype y axis, and fill by difference of lgtv vs pbs
lgtv_pbs_gene_comparison <- function(dat, celltypes, gene){
  for(i in 1:length(celltypes)){
    current_celltype = subset(dat, manualAnnotation == celltypes[i])
    print(celltypes[i])
    print(AverageExpression(current_celltype, features = gene, group.by = 'Treatment'))
  }
}

isgs_to_plot <- c('Rsad2', 'Ifitm1', 'Ifitm3', 'Ifit1', 'Ifit2', 'Isg15')
for(i in 1:length(isgs_to_plot)){
  gene = isgs_to_plot[i]
  avg_exp <- AverageExpression(wt, features = gene, group.by = 'time_treatment_celltype')
  avg_exp_df <- as.data.frame(t(as.matrix(avg_exp$RNA)))
  avg_exp_df$id = rownames(avg_exp_df)
  avg_exp_df <- tidyr::extract(avg_exp_df, into = c('timepoint_treatment', 'celltype'), col = 'id', regex = '(.+)-(.+)')
  gene_heatmap <- ggplot(avg_exp_df, aes_string(x = 'timepoint_treatment', y = 'celltype', fill = gene))+
    geom_tile()+
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5))+
    scale_fill_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white"),
                          values = c(1.0,0.65,0.3,0))+
    ggtitle(gene)
  print(gene_heatmap)
  }

wt[['RNA']]$data[isgs_to_plot,]

lgtv_pbs_gene_comparison(wt, celltypes = c('Astrocytes', 'Microglia', 'Choroid Plexus'), gene = 'Ifit2')

wt_cerebrum_pbs <- subset(wt_cerebrum, Treatment == 'PBS')
wt_cerebrum_chlgtv <- subset(wt_cerebrum, Treatment == 'rChLGTV')


##############################Other Marker Genes#########################################
#########################################################################################

#IFN genes
ifn_genes <- rownames(ParseSeuratObj_int[['RNA']]$data)[grep('Ifn', rownames(ParseSeuratObj_int[['RNA']]$data))]
ifn_genes <- ifn_genes[grep('Ifna', ifn_genes)]

pdf('~/Documents/ÖverbyLab/scPlots/IfnDotPlot.pdf', width = 8, height = 4)
DotPlot(ParseSeuratObj_int, features = ifn_genes, scale = FALSE, group.by = 'Timepoint')
dev.off()

#Reln
dplot <- DotPlot(chimeric_and_mock, features = 'Reln', group.by = 'cellType_treatment') 
dplot_dat_Reln <- dplot$data
dplot_Reln_meta <- str_split_fixed(dplot_dat_Reln$id, "_", 2)
colnames(dplot_Reln_meta) = c('cellType', 'treatment')
dplot_dat_Reln <- cbind(dplot_dat_Reln, dplot_Reln_meta)
ggplot(dplot_dat_Reln, aes(x = treatment, y = cellType, color = avg.exp.scaled, size = pct.exp))+
  geom_point()+
  ggtitle('Reln expression single-cell')
VlnPlot(chimeric_and_mock, 'Reln', group.by = 'cellType_treatment') + theme(legend.position = 'None')

chi_mock_pseudo <- createPseudoBulk(chimeric_and_mock, c('Genotype', 'Timepoint', 'Treatment'))
chi_mock_pseudo <- DESeq(chi_mock_pseudo)
resultsNames(chi_mock_pseudo)
treatment_res <- results(chi_mock_pseudo, name = 'Treatment_rChLGTV_vs_PBS')
treatment_res['Reln',]

#LRP8
FeaturePlot(ParseSeuratObj_int, reduction = 'umap.integrated', features = 'Lrp8')
FeaturePlot(ParseSeuratObj_int, reduction = 'umap.integrated', features = 'Rsad2')

VlnPlot(ParseSeuratObj_int, features = 'Rsad2', group.by = 'hasVirus', pt.size = 0)
DotPlot(subset(ParseSeuratObj_int, Treatment == 'rChLGTV' | Treatment == 'PBS'), features = 'Lrp8', 
        group.by = 'manualAnnotation',split.by = 'Treatment', scale = FALSE)

chimeric_and_mock <- subset(ParseSeuratObj_int, Treatment == 'rChLGTV' | Treatment == 'PBS')
chimeric_and_mock$cellType_treatment <- paste(chimeric_and_mock$manualAnnotation, chimeric_and_mock$Treatment, sep = '_')
chimeric_and_mock_ips <- subset(chimeric_and_mock, Genotype == 'IPS1')
chimeric_and_mock_wt <- subset(chimeric_and_mock, Genotype == 'WT')
dplot <- DotPlot(chimeric_and_mock_wt, features = 'rna_Lrp8', 
                 group.by = 'cellType_treatment', scale = TRUE)

dplot_dat <- dplot$data
dplot_meta <- str_split_fixed(dplot_dat$id, "_", 2)
colnames(dplot_meta) = c('cellType', 'treatment')
dplot_dat <- cbind(dplot_dat, dplot_meta)
ggplot(dplot_dat, aes(x = treatment, y = cellType, color = avg.exp.scaled, size = pct.exp))+
  geom_point()+
  ggtitle('Lrp8 expression WT')+
  scale_color_gradient2(low = '#3DB9FF', mid = 'white', high = 'red', midpoint = 0)+
  scale_size(range = c(2, 10))+
  theme_bw()+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 


FeaturePlot(ParseSeuratObj_int, features = 'Lrp8', reduction = 'umap.integrated')


#Could ISG differences in cerebrum and cerebellum be due to cell type porportion diffs?
#Particularly granulocytes + macrophages

cerebellum <- subset(ParseSeuratObj_int, Organ == 'Cerebellum')
cerebrum <- subset(ParseSeuratObj_int, Organ == 'Cerebrum')

cerebellum_chLgtv <- subset(cerebellum, Treatment == 'rChLGTV')
cerebrum_chLgtv <- subset(cerebrum, Treatment == 'rChLGTV')

table(cerebrum_chLgtv$Timepoint, cerebrum_chLgtv$manualAnnotation)
table(cerebellum_chLgtv$Timepoint, cerebellum_chLgtv$manualAnnotation)

#Ifng could be driving things, look at levels across groups
#Look overall at ips
DotPlot(ips, features = 'Ifng', group.by = 'time_treatment', scale = F)$data
DotPlot(ips, features = 'Ifng', group.by = 'time_treatment', scale = F)+
  scale_size_continuous(range = c(2,6), breaks = seq(0.5, 3, 0.5))

#Look by celltype
wt_lgtv <-subset(wt, Treatment == 'rLGTV')
wt_chimeric_mock <-subset(wt, Treatment %in% c('PBS', 'rChLGTV'))

ips_lgtv <-subset(ips, Treatment == 'rLGTV')
ips_chimeric_mock <-subset(ips, Treatment %in% c('PBS', 'rChLGTV'))

time_celltype_gene_dot <- function(dat, title, gene, exp_max = 3.2, max_pct_exp = 100, maxDotSize = 6,
                                   plotTextSize = 4){
  #Get number of cells of each cell type at each timepoint
  cell_counts <- data.frame(id = names(table(dat$time_treatment_celltype)), counts = c(unname(table(dat$time_treatment_celltype))))
  #get plotting data
  dot_dat <- DotPlot(dat, features = gene, scale = FALSE, group.by = 'time_treatment_celltype')$data
  #combine plotting data with cell counts
  dot_dat <- left_join(dot_dat, cell_counts, by = 'id')
  #prep data for plotting
  #Use extract to split column by second underscore (keep time and treatment together for x axis)
  dot_dat <- tidyr::extract(dot_dat, col = id, regex = "(.+?_.+?)_(.+)", into = c('time_treatment', 'celltype'))
  dot_dat$celltype = factor(dot_dat$celltype, levels = c(resident_celltypes, infiltrating_celltypes, 'unknown'))
  
  dot <- ggplot(dot_dat, aes(x = time_treatment, y = celltype, color = avg.exp, size = pct.exp, label = counts))+
    geom_point()+
    geom_text(color = 'black', size = plotTextSize, hjust = -0.5)+
  
    ggtitle(title)+
    scale_color_gradientn(colours = c("#F03C0C","#F57456","#FFB975","lightgray"),
                          values = c(1.0,0.65,0.3,0),
                          limits = c(0,exp_max))+
    scale_size_continuous(range = c(1,maxDotSize), limits = c(0, max_pct_exp))+
    theme_bw()+
    theme(panel.border = element_blank(),
          axis.line = element_line(color = 'black'),
          axis.text = element_text(size = 21),
          axis.text.x = element_text(vjust = 0.5, angle = 30),
          plot.title = element_text(size = 20),
          legend.title = element_text(size = 15),
          legend.text=element_text(size=16))+
    xlab('')+
    ylab('')
  print(dot)
}

#WT Ifng and Ifngr plots
pdf("~/Documents/ÖverbyLab/single_cell_ISG_figures/isg_chemokine_fig_plots/wt_ifng_dot.pdf", width = 12, height = 6)
time_celltype_gene_dot(wt_chimeric_mock,  title = 'WT Ifng Expression', gene = 'Ifng', max_pct_exp = 55,
                       maxDotSize = 11, plotTextSize = 6)
dev.off()

pdf("~/Documents/ÖverbyLab/single_cell_ISG_figures/isg_chemokine_fig_plots/wt_ifngr1_dot.pdf", width = 12, height = 6)
time_celltype_gene_dot(wt_chimeric_mock,  title = 'WT Ifngr1 Expression', gene = 'Ifngr1', exp_max = 9.5,
                       maxDotSize = 11, plotTextSize = 6)
dev.off()

pdf("~/Documents/ÖverbyLab/single_cell_ISG_figures/isg_chemokine_fig_plots/wt_ifngr2_dot.pdf", width = 12, height = 6)
time_celltype_gene_dot(wt_chimeric_mock,  title = 'WT Ifngr2 Expression', gene = 'Ifngr2', exp_max = 8,
                       maxDotSize = 11, plotTextSize = 6)
dev.off()

#Ips Ifng and Ifngr plots
pdf("~/Documents/ÖverbyLab/single_cell_ISG_figures/isg_chemokine_fig_plots/ips_ifng_dot.pdf", width = 12, height = 6)
time_celltype_gene_dot(ips_chimeric_mock,  title = 'IPS Ifng Expression', gene = 'Ifng',  max_pct_exp = 55,
                       maxDotSize = 11, plotTextSize = 6)
dev.off()

pdf("~/Documents/ÖverbyLab/single_cell_ISG_figures/isg_chemokine_fig_plots/ips_ifngr1_dot.pdf", width = 12, height = 6)
time_celltype_gene_dot(ips_chimeric_mock,  title = 'IPS Ifngr1 Expression', gene = 'Ifngr1', exp_max = 9.5,
                       maxDotSize = 11, plotTextSize = 6)
dev.off()

pdf("~/Documents/ÖverbyLab/single_cell_ISG_figures/isg_chemokine_fig_plots/ips_ifngr2_dot.pdf", width = 12, height = 6)
time_celltype_gene_dot(ips_chimeric_mock,  title = 'IPS Ifngr2 Expression', gene = 'Ifngr2', exp_max = 8,
                       maxDotSize = 11, plotTextSize = 6)
dev.off()

#ips chimeric seems to have unique nk cell infiltration
time_celltype_gene_dot(ips_chLgtv, title = 'IPS ChLGTV Ifng Expression')

ips$time_treatment <- factor(ips$time_treatment)
ips_dds <- createPseudoBulk(ips, variables = c("time_treatment", "Organ"))
ips_dds <- DESeq(ips_dds)

#Choose ref level for comparison
colData(ips_dds)$time_treatment <- relevel(colData(ips_dds)$time_treatment, ref = "Day 3-rLGTV")
colData(ips_dds)$time_treatment <- relevel(colData(ips_dds)$time_treatment, ref = "rChLGTV-Day 3")
ips_dds <- DESeq(ips_dds)
resultsNames(ips_dds)

#Depending on ref choice
results(ips_dds, name ='time_treatment_Day.4.rLGTV_vs_Day.3.rLGTV')['Ifng',]
results(ips_dds, name ='time_treatment_rChLGTV.Day.4_vs_rChLGTV.Day.3')['Ifng',]
results(ips_dds, name ='time_treatment_rChLGTV.Day.5_vs_rChLGTV.Day.3')['Ifng',]

AverageExpression(ips, features = 'Ifng', group.by = 'time_treatment')

#Socs genes
chimeric_mock_wt_infected <- subset(chimeric_mock_wt, Treatment == 'rChLGTV')
chimeric_mock_wt_infected_resident <- subset(chimeric_mock_wt_infected, manualAnnotation %in% resident_celltypes)
chimeric_mock_wt_infected_infiltrating<- subset(chimeric_mock_wt_infected, manualAnnotation %in% infiltrating_celltypes)
DotPlot(chimeric_mock_wt_infected_resident, features = paste0('Socs', c('1', '2', '3', '4' ,'5', '6', '7')), 
                    group.by = 'celltype_time', scale = FALSE)+
  theme(panel.grid.major = element_line(colour = "lightgrey"),
        axis.text.x = element_text(angle = 45, vjust = 0.6))

DotPlot(chimeric_mock_wt_infected_infiltrating, features = paste0('Socs', c('1', '2', '3', '4' ,'5', '6', '7')), 
        group.by = 'time_celltype', scale = FALSE)+
  theme(panel.grid.major = element_line(colour = "lightgrey"))

#Socs in IPS
chimeric_mock_IPS_infected <- subset(chimeric_mock_ips, Treatment == 'rChLGTV')
chimeric_mock_IPS_infected_resident <- subset(chimeric_mock_IPS_infected, manualAnnotation %in% resident_celltypes)
chimeric_mock_IPS_infected_infiltrating<- subset(chimeric_mock_IPS_infected, manualAnnotation %in% infiltrating_celltypes)
DotPlot(chimeric_mock_IPS_infected_resident, features = paste0('Socs', c('1', '2', '3', '4' ,'5', '6', '7')), 
        group.by = 'celltype_time', scale = FALSE)+
  theme(panel.grid.major = element_line(colour = "lightgrey"),
        axis.text.x = element_text(angle = 45, vjust = 0.6))


