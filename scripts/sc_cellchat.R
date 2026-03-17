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

#Perform cellchat analysis
cellchat_input <- ParseSeuratObj_int[["RNA"]]$data
labels <- ParseSeuratObj_int$manualAnnotation
meta <- data.frame(labels = labels, row.names = names(labels))

#Create sample column and cellchat object
ParseSeuratObj_int$samples = factor(ParseSeuratObj_int$orig.ident)
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
cellchat <- aggregateNet(cellchat)

#Plot general cell-cell trends
groupSize <- as.numeric(table(cellchat@idents))

netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

