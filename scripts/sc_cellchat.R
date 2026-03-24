#Packages and functions
library(Seurat)
library(RColorBrewer)
library(CellChat)
source('~/Documents/ÖverbyLab//scripts/langatFunctions.R')
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
mock_cells_cc <- prep_cellchat_obj(mock_cells)
netAnalysis_signalingRole_scatter(mock_cells_cc)

wt_cells_cc <- prep_cellchat_obj(wt_cells)
netAnalysis_signalingRole_scatter(wt_cells_cc)
groupSize_wt <- as.numeric(table(wt_cells_cc@idents))

ips_cells_cc <- prep_cellchat_obj(ips_cells)
netAnalysis_signalingRole_scatter(ips_cells_cc)
groupSize_ips <- as.numeric(table(ips_cells_cc@idents))

#Plot general cell-cell trends
netVisual_circle(wt_cells_cc@net$count, vertex.weight = groupSize_wt, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(wt_cells_cc@net$weight, vertex.weight = groupSize_wt, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

netVisual_circle(ips_cells_cc@net$count, vertex.weight = groupSize_ips, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(ips_cells_cc@net$weight, vertex.weight = groupSize_ips, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

#Look at specific cell type interactions
#Had to recreate aggregateNet function in cellChatBug.R in order to use sources and targets. Bug reported on github
plot_source_targets <- function(cc_obj, target_cell = NULL, source_cell = NULL){
  cellchat_obj <- aggregateNet_bug_removed(cc_obj, source = source_cell, target = target_cell)
  groupSize <- as.numeric(table(cellchat_obj@idents))
  netVisual_circle(cellchat_obj@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
}

#Astrocyte interactions
pdf(file ="~/Documents/ÖverbyLab/scPlots/cellchat_plots/mock_astro_source.pdf", width = 8, height =6)
plot_source_targets(mock_cells_cc, source_cell = 'Astrocytes')
dev.off()

pdf(file ="~/Documents/ÖverbyLab/scPlots/cellchat_plots/wt_astro_source.pdf", width = 8, height =6)
plot_source_targets(wt_cells_cc, source_cell = 'Astrocytes')
dev.off()

pdf(file ="~/Documents/ÖverbyLab/scPlots/cellchat_plots/ips_astro_source.pdf", width = 8, height =6)
plot_source_targets(ips_cells_cc, source_cell = 'Astrocytes')
dev.off()

pdf(file ="~/Documents/ÖverbyLab/scPlots/cellchat_plots/mock_astro_target.pdf", width = 8, height =6)
plot_source_targets(mock_cells_cc, target_cell = 'Astrocytes')
dev.off()

pdf(file ="~/Documents/ÖverbyLab/scPlots/cellchat_plots/wt_astro_target.pdf", width = 8, height =6)
plot_source_targets(wt_cells_cc, target_cell = 'Astrocytes')
dev.off()

pdf(file ="~/Documents/ÖverbyLab/scPlots/cellchat_plots/ips_astro_target.pdf", width = 8, height =6)
plot_source_targets(ips_cells_cc, target_cell = 'Astrocytes')
dev.off()

#Endothelial interactions
pdf(file ="~/Documents/ÖverbyLab/scPlots/cellchat_plots/mock_endo_source.pdf", width = 8, height = 6)
plot_source_targets(mock_cells_cc, source_cell = 'Endothelial')
dev.off()

pdf(file ="~/Documents/ÖverbyLab/scPlots/cellchat_plots/wt_endo_source.pdf", width = 8, height = 6)
plot_source_targets(wt_cells_cc, source_cell = 'Endothelial')
dev.off()

pdf(file ="~/Documents/ÖverbyLab/scPlots/cellchat_plots/wt_endo_target.pdf", width = 8, height = 6)
plot_source_targets(wt_cells_cc, target_cell = 'Endothelial')
dev.off()

#Microglia interactions
pdf(file ="~/Documents/ÖverbyLab/scPlots/cellchat_plots/wt_micro_source.pdf", width = 8, height = 6)
plot_source_targets(wt_cells_cc, source_cell = 'Microglia')
dev.off()

pdf(file ="~/Documents/ÖverbyLab/scPlots/cellchat_plots/wt_micro_target.pdf", width = 8, height = 6)
plot_source_targets(wt_cells_cc, target_cell ='Microglia')
dev.off()

#Can view weights 
wt_cells_cc@net

#Look at specific pathways
CellChatDB$interaction$pathway_name
wt_cells_cc@netP$pathways
#Look at ligand-receptor pairs in select pathwayy
CellChatDB$interaction[CellChatDB$interaction$pathway_name == 'Glutamate',]$ligand %>% unique()
#Plot pathway
pathways.show <- c("IL1") 
netVisual_aggregate(wt_cells_cc, signaling = pathways.show)

#Aggregate pathway chord diagram
pdf(file ="~/Documents/ÖverbyLab/scPlots/cellchat_plots/cellchat_agg_chord.pdf", width = 20, height =16)
netVisual_aggregate(wt_cells_cc, signaling = pathways.show, layout = "chord")
dev.off()

#Heatmap of il1 pathway
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")

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
int_path_heatmap <- function(cc_dat, source_cells = NULL, target_cells = NULL, main = ""){
  comm_subset <-  subsetCommunication(cc_dat, slot.name = "net",
                                       sources.use = source_cells, targets.use = target_cells,
                                       signaling = NULL,
                                       pairLR.use = NULL,
                                       thresh = 0.05)
  comm_subset$cell_int <- paste(comm_subset$source, comm_subset$target, sep = '_')
  
  #Order interaction name by probability for plotting
  comm_subset$interaction_name <- factor(comm_subset$interaction_name, levels = comm_subset$interaction_name[order(comm_subset$prob)])
  
  #Plot
  final_plot <- comm_subset %>% dplyr::arrange(prob) %>% dplyr::filter(prob > 0.1) %>% 
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

int_path_heatmap(mock_cells_cc, source_cells = 'Astrocytes', target_cells = 'Neurons', main = 'mock')
int_path_heatmap(wt_cells_cc, source_cells = 'Astrocytes', target_cells = 'Neurons', main = 'wt')
int_path_heatmap(ips_cells_cc, source_cells = 'Astrocytes', target_cells = 'Neurons', main = 'ips1')

#Loot at specific gene expression
DotPlot(mock_cells, features = 'Sv2b', group.by = 'manualAnnotation')
DotPlot(wt_cells, features = 'Negr1', group.by = 'manualAnnotation')

#Try chord diagram
pdf(file ="~/Documents/ÖverbyLab/scPlots/cellchat_plots/astro_to_neuron_chord.pdf", width = 20, height =16)
netVisual_chord_gene(wt_cells_cc, sources.use = 'Astrocytes', targets.use = 'Neurons', lab.cex = 1.1,legend.pos.y = 30)
dev.off()

pdf(file ="~/Documents/ÖverbyLab/scPlots/cellchat_plots/neuron_to_astro_chord.pdf", width = 20, height =16)
netVisual_chord_gene(wt_cells_cc, sources.use = 'Neurons', targets.use = 'Astrocytes', lab.cex = 1.1,legend.pos.y = 30)
dev.off()

#Can plot select signalling pathways
plotGeneExpression(wt_cells_cc, signaling = "NRXN", type = "violin")
plotGeneExpression(cellchat, signaling = "TNF", type = "violin")
