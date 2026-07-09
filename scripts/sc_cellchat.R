#Packages and functions
library(Seurat)
library(RColorBrewer)
library(CellChat)
library(ggpubr)
library(patchwork)
library(ComplexHeatmap)
library(tidyverse)
library(ggalluvial)
library(NMF) #Needed to make cellchat selectK work
source('~/Documents/ÖverbyLab/scripts/langatFunctions.R')
source('~/Documents/ÖverbyLab/scripts/cellChatBug.R')

#Load data
ParseSeuratObj_int <- LoadSeuratRds("~/Documents/ÖverbyLab/data/FilteredRpcaIntegratedDatNoDoublets.rds") 

#Check data
newCols <-  c(brewer.pal(12, 'Paired'), '#99FFE6', '#CE99FF', '#18662E','#737272',  '#FF8AEF')
newCols[11] =  '#FF8AEF'
DimPlot(ParseSeuratObj_int, label = FALSE, group.by = 'manualAnnotation', reduction = 'umap.integrated',
        cols = newCols)

ParseSeuratObj_int <- subset(ParseSeuratObj_int, manualAnnotation != 'unknown')

#Create sample column, split genotypes, and make cellchat object
ParseSeuratObj_int$samples = factor(ParseSeuratObj_int$orig.ident)
wt_cells <- subset(ParseSeuratObj_int, Genotype == 'WT' & Treatment == 'rChLGTV')
ips_cells <- subset(ParseSeuratObj_int, Genotype == 'IPS1' & Treatment == 'rChLGTV')
mock_cells <- subset(ParseSeuratObj_int, Treatment == 'PBS')
mock_wt_cells <- subset(ParseSeuratObj_int, Treatment == 'PBS' & Genotype == 'WT')
mock_ips_cells <- subset(ParseSeuratObj_int, Treatment == 'PBS' & Genotype == 'IPS1')

#Remove original data object to save space
rm(ParseSeuratObj_int)

#Split into timepoints as well
wt_cells_three <- subset(wt_cells, Timepoint == 'Day 3')
wt_cells_four <- subset(wt_cells, Timepoint == 'Day 4')
wt_cells_five <- subset(wt_cells, Timepoint == 'Day 5')

ips_cells_three <- subset(ips_cells, Timepoint == 'Day 3')
ips_cells_four <- subset(ips_cells, Timepoint == 'Day 4')
ips_cells_five <- subset(ips_cells, Timepoint == 'Day 5')

#Create ligand receptor database
CellChatDB <- CellChatDB.mouse
showDatabaseCategory(CellChatDB)

prep_cellchat_obj <- function(seu_obj){
  cellchat_obj <- createCellChat(object = seu_obj, group.by = "manualAnnotation", assay = "RNA")
  cellchat_obj <- setIdent(cellchat_obj, ident.use = "manualAnnotation")
  
  #set the used database in the object
  cellchat_obj@DB <- CellChatDB
  
  #Preprocess expression data
  cellchat_obj <- subsetData(cellchat_obj) # This step is necessary even if using the whole database
  
  future::plan("multisession", workers = 1) # doing in parellel leads to memory issues a few lines down. so 1 for now
  cellchat_obj <- identifyOverExpressedGenes(cellchat_obj)
  
  #Need to allow greater ram usage to run pca integration
  options(future.globals.maxSize = 3000 * 1024^2)
  cellchat_obj <- identifyOverExpressedInteractions(cellchat_obj)
  
  #reset max
  options(future.globals.maxSize = 500 * 1024^2)
  
  #Compute the communication probability - pretty slow step, takes around 10 mins on laptop
  cellchat_obj <- computeCommunProb(cellchat_obj, type = "triMean")
  
  #Infer pathway level
  cellchat_obj <- computeCommunProbPathway(cellchat_obj)
  
  cellchat_obj <- aggregateNet(cellchat_obj)
  
  #centrality scores
  cellchat_obj <- netAnalysis_computeCentrality(cellchat_obj, slot.name = "netP")
  return(cellchat_obj)
}

#Create cellchat objects for mock and infected cells and look at overall senders/receivers
custom_net_signal_scatter <- function(cc_obj, main = NULL, xlimit, ylimit){
  p1 <- netAnalysis_signalingRole_scatter(cc_obj)+
    xlim(c(0,xlimit))+
    ylim(0,ylimit)+
    ggtitle(main)+
    theme(plot.title = element_text(size = 25))
  plot(p1)
}



mock_wt_cells_cc <- prep_cellchat_obj(mock_wt_cells)
mock_ips_cells_cc <- prep_cellchat_obj(mock_ips_cells)

pdf('~/Documents/ÖverbyLab/scPlots/cellchat_plots/wt_mock_signal_scatter.pdf', height = 5, width = 5)
custom_net_signal_scatter(mock_wt_cells_cc, main = 'Mock wt', xlimit = 33, ylimit = 38)
dev.off()

pdf('~/Documents/ÖverbyLab/scPlots/cellchat_plots/ips_mock_signal_scatter.pdf', height = 5, width = 5)
custom_net_signal_scatter(mock_ips_cells_cc, main = 'Mock ips', xlimit = 33, ylimit = 38)
dev.off()

netAnalysis_signalingRole_heatmap(mock_wt_cells_cc, pattern = "outgoing", cluster.cols = TRUE)
netAnalysis_signalingRole_heatmap(mock_ips_cells_cc, pattern = "outgoing", cluster.cols = TRUE)

if(FALSE){
  #Plots split by genotype
  wt_cells_cc <- prep_cellchat_obj(wt_cells)
  custom_net_signal_scatter(wt_cells_cc, main = 'WT chLGTV', xlimit = 33, ylimit = 33)
  groupSize_wt <- as.numeric(table(wt_cells_cc@idents))
  
  ips_cells_cc <- prep_cellchat_obj(ips_cells)
  custom_net_signal_scatter(ips_cells_cc, main = 'IPS1 chLGTV', xlimit = 33, ylimit = 33)
  groupSize_ips <- as.numeric(table(ips_cells_cc@idents))
}

###### Create cellchat objects for each genotype at each timepoint ######
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

#Plots split by timepoint
wt_cells_three_cc <- prep_cellchat_obj(wt_cells_three)
wt_p3 <- custom_net_signal_scatter(wt_cells_three_cc, main = 'WT chLGTV Day 3', xlimit = 33, ylimit = 38)

wt_cells_four_cc <- prep_cellchat_obj(wt_cells_four)
pdf('~/Documents/ÖverbyLab/scPlots/cellchat_plots/wt_day_4_signal_scatter.pdf', height = 5, width = 5)
wt_p4 <- custom_net_signal_scatter(wt_cells_four_cc, main = 'WT chLGTV Day 4', xlimit = 33, ylimit = 38)
dev.off()

wt_cells_five_cc <- prep_cellchat_obj(wt_cells_five)
wt_p5 <- custom_net_signal_scatter(wt_cells_five_cc, main = 'WT chLGTV Day 5', xlimit = 33, ylimit = 38)

ggpubr::ggarrange(wt_p3, wt_p4, wt_p5, nrow = 1, ncol = 3)

ips_cells_three_cc <- prep_cellchat_obj(ips_cells_three)
ips_p3 <- custom_net_signal_scatter(ips_cells_three_cc, main = 'IPS chLGTV Day 3', xlimit = 33, ylimit = 38)

ips_cells_four_cc <- prep_cellchat_obj(ips_cells_four)

pdf('~/Documents/ÖverbyLab/scPlots/cellchat_plots/ips_day_4_signal_scatter.pdf', height = 5, width = 5)
ips_p4 <- custom_net_signal_scatter(ips_cells_four_cc, main = 'IPS chLGTV Day 4', xlimit = 33, ylimit = 38)
dev.off()

ips_cells_five_cc <- prep_cellchat_obj(ips_cells_five)
ips_p5 <- custom_net_signal_scatter(wt_cells_five_cc, main = 'IPS chLGTV Day 5', xlimit = 33, ylimit = 38)

ggpubr::ggarrange(ips_p3, ips_p4, ips_p5, nrow = 1, ncol = 3)
ggpubr::ggarrange(wt_p3, ips_p3, nrow = 1, ncol = 2)

#Plot general cell-cell trends
if(FALSE){
  netVisual_circle(wt_cells_cc@net$count, vertex.weight = groupSize_wt, weight.scale = T, label.edge= F, title.name = "Number of interactions")
  netVisual_circle(wt_cells_cc@net$weight, vertex.weight = groupSize_wt, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
  
  netVisual_circle(ips_cells_cc@net$count, vertex.weight = groupSize_ips, weight.scale = T, label.edge= F, title.name = "Number of interactions")
  netVisual_circle(ips_cells_cc@net$weight, vertex.weight = groupSize_ips, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
}

#Look at specific cell type interactions
#Had to recreate aggregateNet function in cellChatBug.R in order to use sources and targets. Bug reported on github
plot_source_targets <- function(cc_obj, target_cell = NULL, source_cell = NULL){
  cellchat_obj <- aggregateNet_bug_removed(cc_obj, source = source_cell, target = target_cell)
  groupSize <- as.numeric(table(cellchat_obj@idents))
  netVisual_circle(cellchat_obj@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
}

#Plot day 4 ips astros
pdf(file ="~/Documents/ÖverbyLab/scPlots/cellchat_plots/ips_day_4_astro_source.pdf", width = 6, height =6)
plot_source_targets(ips_cells_four_cc, target_cell = 'Astrocytes')
dev.off()

#Ips day 4 pathways
ips_4_paths <- netVisual_bubble(ips_cells_four_cc, targets.use = 'Astrocytes', remove.isolate = FALSE, 
                 sources.use = c('Choroid Plexus', 'Endothelial', 'Ependymal',
                                 'Granulocytes', 'Macrophage/Monocytes', 'Microglia', 'Nk cells',
                                 'Oligodendrocytes', 'T cells'), thresh = 0.01, return.data = TRUE)
unique(ips_4_paths$communication$pathway_name) %>% sort()

pdf(file ="~/Documents/ÖverbyLab/scPlots/cellchat_plots/infil_to_astro_paths.pdf", width = 6, height =9)
netVisual_bubble(ips_cells_four_cc, targets.use = 'Astrocytes', remove.isolate = FALSE, 
                 sources.use = c('Granulocytes', 'Macrophage/Monocytes', 'Nk cells','T cells'), thresh = 0.01, color.heatmap = 'viridis')+
  theme(text = element_text(size = 16))
dev.off()
#Astrocyte interactions
pdf(file ="~/Documents/ÖverbyLab/scPlots/cellchat_plots/mock_astro_source.pdf", width = 8, height =6)
plot_source_targets(mock_wt_cells_cc, source_cell = 'Astrocytes')
dev.off()

pdf(file ="~/Documents/ÖverbyLab/scPlots/cellchat_plots/wt_astro_source.pdf", width = 8, height =6)
plot_source_targets(wt_cells_three_cc, source_cell = 'Astrocytes')
dev.off()

pdf(file ="~/Documents/ÖverbyLab/scPlots/cellchat_plots/ips_astro_source.pdf", width = 8, height =6)
plot_source_targets(ips_cells_three_cc, source_cell = 'Astrocytes')
dev.off()

pdf(file ="~/Documents/ÖverbyLab/scPlots/cellchat_plots/mock_astro_target.pdf", width = 8, height =6)
plot_source_targets(mock_wt_cells_cc, target_cell = 'Astrocytes')
dev.off()

pdf(file ="~/Documents/ÖverbyLab/scPlots/cellchat_plots/wt_astro_target.pdf", width = 8, height =6)
plot_source_targets(wt_cells_three_cc, target_cell = 'Astrocytes')
dev.off()

pdf(file ="~/Documents/ÖverbyLab/scPlots/cellchat_plots/ips_astro_target.pdf", width = 8, height =6)
plot_source_targets(ips_cells_three_cc, target_cell = 'Astrocytes')
dev.off()

#Endothelial interactions
pdf(file ="~/Documents/ÖverbyLab/scPlots/cellchat_plots/mock_endo_source.pdf", width = 8, height = 6)
plot_source_targets(mock_wt_cells_cc, source_cell = 'Endothelial')
dev.off()

pdf(file ="~/Documents/ÖverbyLab/scPlots/cellchat_plots/wt_endo_source.pdf", width = 8, height = 6)
plot_source_targets(wt_cells_three_cc, source_cell = 'Endothelial')
dev.off()

pdf(file ="~/Documents/ÖverbyLab/scPlots/cellchat_plots/wt_endo_target.pdf", width = 8, height = 6)
plot_source_targets(wt_cells_three_cc, target_cell = 'Endothelial')
dev.off()

#Microglia interactions
pdf(file ="~/Documents/ÖverbyLab/scPlots/cellchat_plots/wt_micro_source.pdf", width = 8, height = 6)
plot_source_targets(wt_cells_three_cc, source_cell = 'Microglia')
dev.off()

pdf(file ="~/Documents/ÖverbyLab/scPlots/cellchat_plots/wt_micro_target.pdf", width = 8, height = 6)
plot_source_targets(wt_cells_three_cc, target_cell ='Microglia')
dev.off()

#Can view weights 
mock_wt_cells_cc@net

#Look at specific pathways
CellChatDB$interaction$pathway_name
mock_wt_cells_cc@netP$pathways %>% sort()

#Look at ligand-receptor pairs in select pathway
CellChatDB$interaction[CellChatDB$interaction$pathway_name == 'CCL',]$ligand %>% unique()

#Plot pathways overtime
path_over_time <- function(pathway, genotype){
  pathways.show <- c(pathway) 
  if(genotype == 'WT'){
    #Vertex weight being null makes dot sizes proportional to number of cells
    wt_p3 <- netVisual_aggregate(wt_cells_three_cc, signaling = pathways.show, vertex.weight = NULL)
    wt_p4 <- netVisual_aggregate(wt_cells_four_cc, signaling = pathways.show, vertex.weight = NULL)
    wt_p5 <- netVisual_aggregate(wt_cells_five_cc, signaling = pathways.show, vertex.weight = NULL)
    return(list(wt_p3, wt_p4, wt_p5))
  }
  if(genotype == 'IPS1'){
    ips_p3 <- netVisual_aggregate(ips_cells_three_cc, signaling = pathways.show, vertex.weight = NULL)
    ips_p4 <- netVisual_aggregate(ips_cells_four_cc, signaling = pathways.show, vertex.weight = NULL)
    ips_p5 <- netVisual_aggregate(ips_cells_five_cc, signaling = pathways.show, vertex.weight = NULL)
    return(list(ips_p3, ips_p4, ips_p5))
  }
}

if(FALSE){
  ccl_wt_plots <- path_over_time('TNF', genotype = 'WT')
  ccl_wt_plots[[1]]
  
  ccl_ips_plots <- path_over_time('TNF', genotype = 'IPS1')
  ccl_ips_plots[[1]]
  
  netVisual_aggregate(mock_wt_cells_cc, signaling = 'CCL', vertex.weight = NULL)
  
  #Aggregate pathway chord diagram
  pdf(file ="~/Documents/ÖverbyLab/scPlots/cellchat_plots/cellchat_agg_chord.pdf", width = 20, height =16)
  netVisual_aggregate(wt_cells_three_cc, signaling = pathways.show, layout = "chord")
  dev.off()
}


#Side by side pathway heatmaps for two datasets
#Have to use patchwork to get both heatmaps to show legends, therefore have to do this draw grid grab stuff
path_heatmap_comp <- function(cc_dat1, cc_dat2, pathway, plot_path){
  p1 <- netVisual_heatmap(cc_dat1, signaling = pathway, color.heatmap = "Reds") %>% 
    draw() %>% 
    grid.grabExpr()
  p2 <- netVisual_heatmap(cc_dat2, signaling = pathway, color.heatmap = "Reds") %>% 
    draw() %>% 
    grid.grabExpr()
  wrap_plots(list(p1, p2), ncol = 2)
}

#Cannot get legends to stay seperate if plot two together (looks like an issue in cellchat docs to me as well)
#So just plotting seperate for now
pdf('~/Documents/ÖverbyLab/scPlots/cellchat_plots/comp_heatmap_test.pdf', width = 10, height = 6)
path_heatmap_comp(wt_cells_three_cc, ips_cells_three_cc, pathway = 'CCL')
dev.off()
#Break pathway into specific lr pairs
netAnalysis_contribution(cellchat, signaling = pathways.show)

# show all the interactions sending from astrocytes to neurons
pdf(file ="~/Documents/ÖverbyLab/scPlots/cellchat_plots/mock_astro_to_neuron.pdf", width = 20, height =16)
netVisual_chord_gene(mock_cells_cc, sources.use = 1, targets.use = c(11), lab.cex = 1,legend.pos.y = 30)
dev.off()

pdf(file ="~/Documents/ÖverbyLab/scPlots/cellchat_plots/wt_astro_to_neuron.pdf", width = 20, height =16)
netVisual_chord_gene(wt_cells_cc, sources.use = 1, targets.use = c(11), lab.cex = 1,legend.pos.y = 30)
dev.off()

#Look at pathway centrality
netAnalysis_signalingRole_network(wt_cells_cc, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(wt_cells_cc, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(wt_cells_cc, pattern = "incoming")
ht1

#Look at top pathways between astros and neurons
int_path_heatmap <- function(cc_dat, source_cells = NULL, target_cells = NULL, main = "", prob_limit = 0.1){
  comm_subset <-  subsetCommunication(cc_dat, slot.name = "net",
                                       sources.use = source_cells, targets.use = target_cells,
                                       signaling = NULL,
                                       pairLR.use = NULL,
                                       thresh = 0.05)
  comm_subset$cell_int <- paste(comm_subset$source, comm_subset$target, sep = '_')
  
  #Order interaction name by probability for plotting
  comm_subset$interaction_name <- factor(comm_subset$interaction_name, levels = comm_subset$interaction_name[order(comm_subset$prob)])
  
  #Plot
  final_plot <- comm_subset %>% dplyr::arrange(prob) %>% dplyr::filter(prob > prob_limit) %>% 
    ggplot(aes(x = cell_int, y = interaction_name, fill = prob))+
    geom_tile(color = 'black')+
    ylab('')+
    xlab('')+
    ggtitle(main)+
    scale_fill_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white"), 
                         values = c(1.0,0.7,0.4,0),
                         limits = c(0,0.35))+
    theme_bw()
  
  print(final_plot)
}

int_path_heatmap(wt_cells_three_cc, source_cells = 'Astrocytes', target_cells = 'Neurons', main = 'mock')
int_path_heatmap(wt_cells_cc, source_cells = 'Astrocytes', target_cells = 'Neurons', main = 'wt')
int_path_heatmap(ips_cells_cc, source_cells = 'Astrocytes', target_cells = 'Neurons', main = 'ips1')

#Loot at specific gene expression
DotPlot(mock_cells, features = 'Sv2b', group.by = 'manualAnnotation')
DotPlot(wt_cells, features = 'Negr1', group.by = 'manualAnnotation')

#Try chord diagram
pdf(file ="~/Documents/ÖverbyLab/scPlots/cellchat_plots/astro_to_neuron_chord.pdf", width = 20, height =16)
netVisual_chord_gene(wt_cells_three_cc, sources.use = 'Astrocytes', targets.use = 'Neurons', lab.cex = 1.1,legend.pos.y = 30)
dev.off()

pdf(file ="~/Documents/ÖverbyLab/scPlots/cellchat_plots/neuron_to_astro_chord.pdf", width = 20, height =16)
netVisual_chord_gene(wt_cells_three_cc, sources.use = 'Neurons', targets.use = 'Astrocytes', lab.cex = 1.1,legend.pos.y = 30)
dev.off()

#Can plot select signalling pathways
plotGeneExpression(wt_cells_three_cc, signaling = "TNF", type = "violin")
plotGeneExpression(ips_cells_three_cc, signaling = "TNF", type = "violin")


########### Comparative analysis between wt and ips cellchat ########### 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

#Merge datasets
cellchat_three_merged <- mergeCellChat(list(wt = wt_cells_three_cc, ips = ips_cells_three_cc), add.names = c('WT', 'IPS1'))
compareInteractions(cellchat_three_merged, show.legend = F, group = c(1,2))
netVisual_diffInteraction(cellchat_three_merged, weight.scale = T, measure = 'count')
netVisual_heatmap(cellchat_three_merged)
netVisual_heatmap(cellchat_three_merged, comparison = c(1,2))

cellchat_four_merged <- mergeCellChat(list(wt = wt_cells_four_cc, ips = ips_cells_four_cc), add.names = c('WT', 'IPS1'))
compareInteractions(cellchat_four_merged, show.legend = F, group = c(1,2))
netVisual_diffInteraction(cellchat_four_merged, weight.scale = T)

cellchat_five_merged <- mergeCellChat(list(wt = wt_cells_five_cc, ips = ips_cells_five_cc), add.names = c('WT', 'IPS1'))
compareInteractions(cellchat_five_merged, show.legend = F, group = c(1,2))
netVisual_diffInteraction(cellchat_five_merged, weight.scale = T)
netVisual_heatmap(cellchat_five_merged)

#Identifying signalling changes of specific celltypes between datasets
netAnalysis_signalingChanges_scatter(cellchat_five_merged, idents.use = "Astrocytes")

###### Differential look at signalling networks ###### 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

reticulate::use_python("/Users/matthewbradley/miniconda3/bin/python", required=T) 
cellchat_three_merged <- computeNetSimilarityPairwise(cellchat_three_merged, type = "functional")
#Crashes R for some reason
#cellchat_three_merged <- netEmbedding(cellchat_three_merged, type = "functional")
#cellchat_three_merged <- netClustering(cellchat_three_merged, type = "functional")

#Pathway barplot comparisons for each timepoint
pdf(file ="~/Documents/ÖverbyLab/scPlots/cellchat_plots/day3_path_comp_bar.pdf", width = 8, height =12)
rankNet(cellchat_three_merged, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE)
dev.off()

rankNet(cellchat_four_merged, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE)

pdf(file ="~/Documents/ÖverbyLab/scPlots/cellchat_plots/day5_path_comp_bar.pdf", width = 8, height =12)
rankNet(cellchat_five_merged, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE)
dev.off()

#Up and downregulated ligand-receptor pairs by communication probabilties
netVisual_bubble(cellchat_three_merged, sources.use = 'Astrocytes', targets.use = c('Microglia', 'Neurons'),  comparison = c(1, 2), angle.x = 45)

#Compare pathways
pdf(file ="~/Documents/ÖverbyLab/scPlots/cellchat_plots/wt_cells_three_out_paths.pdf", width = 11, height = 16)
netAnalysis_signalingRole_heatmap(wt_cells_three_cc, pattern = "outgoing", width = 10, height = 25)
dev.off()


pdf(file ="~/Documents/ÖverbyLab/scPlots/cellchat_plots/ips_cells_three_out_paths.pdf", width = 11, height = 16)
netAnalysis_signalingRole_heatmap(ips_cells_three_cc, pattern = "outgoing", width = 10, height = 25)
dev.off()

object <- ips_cells_three_cc
pattern = 'outgoing'
slot.name = "netP"

########### Comparative analysis between day 3 (or 4) and mock by genotype ########### 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

cellchat_wt_3_merged <- mergeCellChat(list(mock = mock_wt_cells_cc, wt_inf = wt_cells_three_cc), add.names = c('PBS', 'WT_inf'))
cellchat_wt_4_merged <- mergeCellChat(list(mock = mock_wt_cells_cc, wt_inf = wt_cells_four_cc), add.names = c('PBS', 'WT_inf'))
compareInteractions(cellchat_wt_merged, show.legend = F, group = c(1,2))
netVisual_diffInteraction(cellchat_wt_merged, weight.scale = T, measure = 'count')

pdf('~/Documents/ÖverbyLab/scPlots/cellchat_plots/wt_4_vs_mock_astro_target.pdf', width = 5, height = 5)
netVisual_diffInteraction(cellchat_wt_3_merged, weight.scale = T, measure = 'count', targets.use = 'Astrocytes')
dev.off()

pdf('~/Documents/ÖverbyLab/scPlots/cellchat_plots/wt_4_vs_mock_astro_source.pdf', width = 5, height = 5)
netVisual_diffInteraction(cellchat_wt_4_merged, weight.scale = T, measure = 'count', sources.use = 'Astrocytes')
dev.off()

pdf('~/Documents/ÖverbyLab/scPlots/cellchat_plots/wt_4_vs_mock_signal_changes.pdf', width = 10, height = 8)
netAnalysis_signalingChanges_scatter(cellchat_wt_4_merged, idents.use = "Astrocytes", label.size = 9)+
  theme(text = element_text(size = 24))+
  xlim(c(-0.4, 3))+
  ylim(c(-0.4, 2.5))
dev.off()

netVisual_heatmap(cellchat_wt_merged, measure = 'count')

#Can subtract total interactions between mock and infected:
wt_celltype_count_diff <- cellchat_wt_merged@net[[2]][['count']] - cellchat_wt_merged@net[[1]][['count']]
rowSums(wt_celltype_count_diff)

#Same analysis with ips
cellchat_ips_3_merged<- mergeCellChat(list(mock = mock_ips_cells_cc, ips_inf = ips_cells_three_cc), add.names = c('PBS', 'IPS_inf'))
cellchat_ips_4_merged<- mergeCellChat(list(mock = mock_ips_cells_cc, ips_inf = ips_cells_four_cc), add.names = c('PBS', 'IPS_inf'))
cellchat_ips_5_merged<- mergeCellChat(list(mock = mock_ips_cells_cc, ips_inf = ips_cells_five_cc), add.names = c('PBS', 'IPS_inf'))

compareInteractions(cellchat_ips_3_merged, show.legend = F, group = c(1,2))
netVisual_diffInteraction(cellchat_ips_3_merged, weight.scale = T, measure = 'count')

pdf('~/Documents/ÖverbyLab/scPlots/cellchat_plots/ips_4_vs_mock_astro_target.pdf', width = 5, height = 5)
netVisual_diffInteraction(cellchat_ips_3_merged, weight.scale = T, measure = 'count', targets.use = 'Astrocytes')
dev.off()

pdf('~/Documents/ÖverbyLab/scPlots/cellchat_plots/ips_4_vs_mock_astro_source.pdf', width = 5, height = 5)
netVisual_diffInteraction(cellchat_ips_3_merged, weight.scale = T, measure = 'count', sources.use = 'Astrocytes')
dev.off()

pdf('~/Documents/ÖverbyLab/scPlots/cellchat_plots/ips_4_vs_mock_signal_changes.pdf', width = 10, height = 8)
netAnalysis_signalingChanges_scatter(cellchat_ips_4_merged, idents.use = "Astrocytes", label.size = 9)+
  theme(text = element_text(size = 24))+
  xlim(c(-0.4, 3))+
  ylim(c(-0.4, 2.5))
dev.off()

#Compare mocks
cellchat_mock_merged <- mergeCellChat(list(mock = mock_wt_cells_cc, wt_inf = mock_ips_cells_cc), add.names = c('wt_mock', 'ips_mock'))
compareInteractions(cellchat_mock_merged, show.legend = F, group = c(1,2))

netVisual_diffInteraction(cellchat_mock_merged, weight.scale = T, measure = 'count', targets.use = 'Astrocytes')
netVisual_diffInteraction(cellchat_mock_merged, weight.scale = T, measure = 'count', sources.use = 'Astrocytes')

pdf('~/Documents/ÖverbyLab/scPlots/cellchat_plots/ips_mock_vs_wt_mock_signal_changes.pdf', width = 10, height = 8)
netAnalysis_signalingChanges_scatter(cellchat_mock_merged, idents.use = "Astrocytes", label.size = 6)+
  theme(text = element_text(size = 24))+
  xlim(c(-0.4, 3.2))+
  ylim(c(-0.4, 2.5))
dev.off()
#Look at raw data for different pathways
mock_ips_cells_cc@netP$centr$Glutamate

#Do not use celltypes that have only a few cells at one point or another
netVisual_heatmap(cellchat_ips_merged, measure = 'count', targets.use = 'Astrocytes')

#total difference between inf and mock astrocytes (based on comapreInteractions function)
#Diff in outgoing
rowSums(cellchat_wt_merged@net$WT_inf$count)['Astrocytes'] - rowSums(cellchat_wt_merged@net$PBS$count)['Astrocytes']
rowSums(cellchat_ips_merged@net$IPS_inf$count)['Astrocytes'] - rowSums(cellchat_ips_merged@net$PBS$count)['Astrocytes']

#Diff in incoming
colSums(cellchat_wt_merged@net$WT_inf$count)['Astrocytes'] - colSums(cellchat_wt_merged@net$PBS$count)['Astrocytes']
colSums(cellchat_ips_merged@net$IPS_inf$count)['Astrocytes'] - colSums(cellchat_ips_merged@net$PBS$count)['Astrocytes']

#Identifying signalling changes of specific celltypes between datasets
netAnalysis_signalingChanges_scatter(cellchat_wt_merged, idents.use = "Astrocytes", xlims = c(-0.5,3.5))
netAnalysis_signalingChanges_scatter(cellchat_ips_merged, idents.use = "Astrocytes", xlims = c(-0.5,3.5))

#Plot dotplots of overall interaction together
merged_dot_comp <- function(cc_obj){
  num.link <- sapply(cc_obj, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
  weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
  gg <- list()
  for (i in 1:length(cc_obj)) {
    gg[[i]] <- netAnalysis_signalingRole_scatter(cc_obj[[i]], title = names(cc_obj)[i], weight.MinMax = weight.MinMax)+
      xlim(c(0, 40))+
      ylim(c(0,40))
  }
  #> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
  #> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
  patchwork::wrap_plots(plots = gg)
}

wt_cc_objects <- list(mock_wt = mock_wt_cells_cc, inf_wt = wt_cells_three_cc)
ips_cc_objects <- list(mock_ips= mock_ips_cells_cc, inf_ips = ips_cells_three_cc)
merged_dot_comp(wt_cc_objects)
merged_dot_comp(ips_cc_objects)
patchwork::wrap_plots(plots = list(merged_dot_comp(wt_cc_objects), merged_dot_comp(ips_cc_objects)),
                       nrow = 2)


########### Comparative analysis between mock wt and mock ips ########### 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
cellchat_mock_merged<- mergeCellChat(list(mock_wt = mock_wt_cells_cc, mock_ips = mock_ips_cells_cc), 
                                   add.names = c('mock_wt', 'mock_ips'))

#Mock ips has more interactions

pdf('~/Documents/ÖverbyLab/scPlots/cellchat_plots/mock_wt_vs_ips_counts.pdf', height = 5, width = 5)
compareInteractions(cellchat_mock_merged, show.legend = F, group = c(1,2))+
  theme(text = element_text(size = 20))
dev.off()
compareInteractions(cellchat_mock_merged, show.legend = F, group = c(1,2), measure = "weight")

#Red = increased in ips
netVisual_heatmap(cellchat_mock_merged)
netVisual_heatmap(cellchat_mock_merged, measure = "weight")

#Total strength by celltype
num.link <- sapply(list(mock_wt = mock_wt_cells_cc, mock_ips = mock_ips_cells_cc), 
                   function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link))

pdf('~/Documents/ÖverbyLab/scPlots/cellchat_plots/mock_signal_wt_dot.pdf', height = 5, width = 5)
netAnalysis_signalingRole_scatter(mock_wt_cells_cc, title = 'Mock WT', weight.MinMax = weight.MinMax)+
  ylim(c(0, 26))+
  theme(text = element_text(size = 20))
dev.off()

pdf('~/Documents/ÖverbyLab/scPlots/cellchat_plots/mock_signal_ips_dot.pdf', height = 5, width = 5)
netAnalysis_signalingRole_scatter(mock_ips_cells_cc, title = 'Mock IPS', weight.MinMax = weight.MinMax)+
  ylim(c(0, 26))+
  theme(text = element_text(size = 20))
dev.off()

#What changes between astrocytes
pdf('~/Documents/ÖverbyLab/scPlots/cellchat_plots/mock_signalling_comp_scatter.pdf', height = 5, width = 7)
netAnalysis_signalingChanges_scatter(cellchat_mock_merged, idents.use = "Astrocytes")+
  theme(text = element_text(size = 18))
dev.off()
#The glutamate difference to me seems to be due to the gls gene,
#since this is mentioned in all the pathways in the ips netAnalysis_contribution and seems to be
#lowly expressed in mock
netVisual_heatmap(mock_wt_cells_cc, signaling = 'Glutamate', color.heatmap = "Reds")
netAnalysis_contribution(mock_wt_cells_cc, signaling = 'Glutamate', sources.use = 'Astrocytes', targets.use = 'Neurons')
plotGeneExpression(mock_wt_cells_cc, signaling = "Glutamate", enriched.only = TRUE, type = "violin")

pdf('~/Documents/ÖverbyLab/scPlots/cellchat_plots/mock_wt_glut_genes.pdf', width = 6, height = 6)
plotGeneExpression(mock_wt_cells_cc, features = c('Slc1a2', 'Gls', 'Grik2', 'Gria4', 'Grm1'), enriched.only = TRUE, type = "violin")
dev.off()

netVisual_heatmap(mock_ips_cells_cc, signaling = 'Glutamate', color.heatmap = "Reds")
netAnalysis_contribution(mock_ips_cells_cc, signaling = 'Glutamate', sources.use = 'Astrocytes', targets.use = 'Neurons')
plotGeneExpression(mock_ips_cells_cc, signaling = "Glutamate", enriched.only = TRUE, type = "violin")

pdf('~/Documents/ÖverbyLab/scPlots/cellchat_plots/mock_ips_glut_genes.pdf', width = 6, height = 6)
plotGeneExpression(mock_ips_cells_cc, features = c('Slc1a2', 'Gls', 'Grik2', 'Gria4', 'Grm1'), enriched.only = TRUE, type = "violin")
dev.off()

#Look directly at path probs in data
mock_wt_cells_cc@netP$prob[,,'Glutamate']
mock_ips_cells_cc@netP$prob[,,'Glutamate']

#Look at genes in data
subsetCommunication(mock_wt_cells_cc) %>% dplyr::filter(source == 'Astrocytes' | target == 'Astrocytes') %>% 
  dplyr::filter(pathway_name == 'Glutamate') %>% 
  dplyr::arrange(source)

#dotplot genes
mock_astros <- subset(mock_cells, manualAnnotation == 'Astrocytes')
DotPlot(mock_astros, features = 'Gls', group.by = 'Genotype', scale = FALSE)

#River plots to compare astros between datasets
selectK(mock_wt_cells_cc, pattern = "incoming")
n_patterns = 5
mock_wt_cells_cc <- identifyCommunicationPatterns(mock_wt_cells_cc, pattern = "incoming", k = n_patterns)
netAnalysis_river(mock_wt_cells_cc, pattern = "incoming")

#SEMA4 pathway across groups
netVisual_aggregate(wt_cells_three_cc, signaling = 'SEMA4', layout = "circle", targets.use = 'Astrocytes')
netVisual_aggregate(ips_cells_three_cc, signaling = 'SEMA4', layout = "circle", targets.use = 'Astrocytes')

