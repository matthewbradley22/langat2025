#Packages and functions
library(Seurat)
library(RColorBrewer)
library(CellChat)
source('~/Documents/ÖverbyLab//scripts/langatFunctions.R')

#Load data
ParseSeuratObj_int <- LoadSeuratRds("~/Documents/ÖverbyLab/data/FilteredRpcaIntegratedDatNoDoublets.rds") 

#Check data
newCols <-  c(brewer.pal(12, 'Paired'), '#99FFE6', '#CE99FF', '#18662E','#737272',  '#FF8AEF')
newCols[11] =  '#FF8AEF'
DimPlot(ParseSeuratObj_int, label = FALSE, group.by = 'manualAnnotation', reduction = 'umap.integrated',
        cols = newCols)

ParseSeuratObj_int <- subset(ParseSeuratObj_int, manualAnnotation != 'unknown')

#Create sample column and cellchat object
ParseSeuratObj_int$samples = factor(ParseSeuratObj_int$orig.ident)

#Test random matrix to see if error persists
cellchat <- createCellChat(object = ParseSeuratObj_int, group.by = "manualAnnotation", assay = "RNA")
cellchat <- setIdent(cellchat, ident.use = "manualAnnotation")
levels(cellchat@idents)

#Create ligand receptor database
CellChatDB <- CellChatDB.mouse
showDatabaseCategory(CellChatDB)

# set the used database in the object
cellchat@DB <- CellChatDB

#Preprocess expression data
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 1) # doing in parellel leads to memory issues a few lines down. so 1 for now
cellchat <- identifyOverExpressedGenes(cellchat)

#Need to allow greater ram usage to run pca integration
options(future.globals.maxSize = 3000 * 1024^2)
cellchat <- identifyOverExpressedInteractions(cellchat)

#reset max
options(future.globals.maxSize = 500 * 1024^2)

#Compute the communication probability - pretty slow step, takes around 10 mins on laptop
cellchat <- computeCommunProb(cellchat, type = "triMean")

#Infer pathway level
cellchat <- computeCommunProbPathway(cellchat)

#Aggregate cell-communication
groupSize <- as.numeric(table(cellchat@idents))
cellchat <- aggregateNet(cellchat)

#Plot general cell-cell trends
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

#Look at specific cell type interactions
#Had to recreate aggregateNet function in cellChatBug.R in order to use sources and targets. Bug reported on github
cellchat_astro <- aggregateNet_bug_removed(cellchat, source = "Astrocytes")

pdf(file ="~/Documents/ÖverbyLab/scPlots/cellchat_plots/astro_source.pdf", width = 8, height =6)
netVisual_circle(cellchat_astro@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()
#Can view weights 
cellchat_astro@net

cellchat_astro_target <- aggregateNet_bug_removed(cellchat, target = "Astrocytes")
groupSize <- as.numeric(table(cellchat_astro@idents))
pdf(file ="~/Documents/ÖverbyLab/scPlots/cellchat_plots/astro_target.pdf", width = 8, height =6)
netVisual_circle(cellchat_astro_target@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()
cellchat_astro_target@net

#Endothelial interactions
cellchat_endo <- aggregateNet_bug_removed(cellchat, source = "Endothelial")

pdf(file ="~/Documents/ÖverbyLab/scPlots/cellchat_plots/endo_source.pdf", width = 8, height =6)
netVisual_circle(cellchat_endo@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()

cellchat_endo_target <- aggregateNet_bug_removed(cellchat, target = "Endothelial")

pdf(file ="~/Documents/ÖverbyLab/scPlots/cellchat_plots/endo_target.pdf", width = 8, height =6)
netVisual_circle(cellchat_endo_target@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()

#Look at specific pathways
CellChatDB$interaction$pathway_name
cellchat@netP$pathways
#Look at ligand-receptor pairs in select pathwayy
CellChatDB$interaction[CellChatDB$interaction$pathway_name == 'Glutamate',]$ligand %>% unique()
#Plot pathway
pathways.show <- c("IL1") 
netVisual_aggregate(cellchat, signaling = pathways.show)

#Aggregate pathway chord diagram
pdf(file ="~/Documents/ÖverbyLab/scPlots/cellchat_plots/cellchat_agg_chord.pdf", width = 20, height =16)
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
dev.off()

#Heatmap of il1 pathway
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")

#Break pathway into specific lr pairs
netAnalysis_contribution(cellchat, signaling = pathways.show)

# show all the interactions sending from astrocytes to neurons
pdf(file ="~/Documents/ÖverbyLab/scPlots/cellchat_plots/astro_to_neuron.pdf", width = 20, height =16)
netVisual_chord_gene(cellchat, sources.use = 1, targets.use = c(11), lab.cex = 0.5,legend.pos.y = 30)
dev.off()

#Network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
#Look at pathway centrality
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

#Overall senders and recievers 
netAnalysis_signalingRole_scatter(cellchat)

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
netAnalysis_signalingRole_scatter(cellchat)

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1

#Look at top pathways between astros and neurons
astro_neuron <-  subsetCommunication(cellchat, slot.name = "net",
                                              sources.use = 'Astrocytes', targets.use = 'Neurons',
                                              signaling = NULL,
                                              pairLR.use = NULL,
                                              thresh = 0.05)
astro_neuron$cell_int <- paste(astro_neuron$source, astro_neuron$target, sep = '_')
astro_neuron %>% dplyr::filter(prob > 0.1) %>% dplyr::arrange(prob) %>% 
  ggplot(aes(x = cell_int, y = interaction_name, color = prob))+
  geom_point()+
  ylab('')+
  xlab('')+
  ggtitle('Top Astrocyte - Neuron interaction paths')+
  scale_color_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white"), 
                       values = c(1.0,0.7,0.4,0),
                       limits = c(0,0.32))+
  theme_bw()

#Try chord diagram
pdf(file ="~/Documents/ÖverbyLab/scPlots/cellchat_plots/astro_to_neuron_chord.pdf", width = 20, height =16)
netVisual_chord_gene(cellchat, sources.use = 'Astrocytes', targets.use = 'Neurons', lab.cex = 1.1,legend.pos.y = 30)
dev.off()

pdf(file ="~/Documents/ÖverbyLab/scPlots/cellchat_plots/neuron_to_astro_chord.pdf", width = 20, height =16)
netVisual_chord_gene(cellchat, sources.use = 'Neurons', targets.use = 'Astrocytes', lab.cex = 1.1,legend.pos.y = 30)
dev.off()

#Can plot select signalling pathways
plotGeneExpression(cellchat, signaling = "NRXN", type = "violin")
plotGeneExpression(cellchat, signaling = "TNF", type = "violin")
