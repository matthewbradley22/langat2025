#Testing bug
#Packages and functions
library(Seurat)
library(CellChat)


#Found the bug, using data from sc_cellchat.R script

aggregateNet_bug_removed <- function(cellchat_dat, target = NULL, source = NULL){
  net <- cellchat_dat@net
  
  df.net <- subsetCommunication(cellchat_dat, slot.name = "net",
                                sources.use = source, targets.use = target,
                                signaling = NULL,
                                pairLR.use = NULL,
                                thresh = 0.05)
  df.net$source_target <- paste(df.net$source, df.net$target, sep = "|")
  df.net2 <- df.net %>% group_by(source_target) %>% summarize(count = n(), .groups = 'drop')
  
  df.net3 <- df.net %>% group_by(source_target) %>% summarize(prob = sum(prob), .groups = 'drop')
  df.net2$prob <- df.net3$prob
  
  #NEW LINE HERE FOR DEFINING A, CELLCHAT LINE HAD A BUG
  a <- stringr::str_split_fixed(df.net2$source_target, "\\|", n = 2)
  df.net2$source <- as.character(a[, 1])
  df.net2$target <- as.character(a[, 2])
  cells.level <- levels(cellchat_dat@idents)
  
  df.net2$source <- factor(df.net2$source, levels = cells.level)
  df.net2$target <- factor(df.net2$target, levels = cells.level)
  
  count <- tapply(df.net2[["count"]], list(df.net2[["source"]], df.net2[["target"]]), sum)
  prob <- tapply(df.net2[["prob"]], list(df.net2[["source"]], df.net2[["target"]]), sum)
  net$count <- count
  net$weight <- prob
  net$weight[is.na(net$weight)] <- 0
  net$count[is.na(net$count)] <- 0
  
  cellchat_dat@net <- net
  return(cellchat_dat)
}


########## For a reprex using seurat pbmc data, uncomment below ########## 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
# 
# 
# set.seed(151)
# 
# #Load in 10x pbmc data
# pbmc.data <- Read10X(data.dir = "~/Downloads/filtered_gene_bc_matrices/hg19/")
# pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
# 
# #Standard seurat preprocessing (skipping some steps)
# pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)
# pbmc <- NormalizeData(pbmc)
# pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
# pbmc <- ScaleData(pbmc, features = rownames(pbmc))
# pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
# pbmc <- FindNeighbors(pbmc, dims = 1:10)
# pbmc <- FindClusters(pbmc, resolution = 0.5)
# pbmc <- RunUMAP(pbmc, dims = 1:10)
# DimPlot(pbmc, reduction = "umap", label = TRUE)
# 
# #Label random cell names for example
# pbmc[[]] <- pbmc[[]] %>% mutate(fake_celltypes = case_when(seurat_clusters %in% c(4,6) ~ 'NK Cells',
#                                                            seurat_clusters %in% c(2,0) ~ 'T Cells' ,
#                                                            seurat_clusters %in% c(3) ~ 'Macropahges',
#                                                            seurat_clusters %in% c(1,5,7) ~ 'Endothelial',
#                                                            seurat_clusters %in% c(8) ~ 'Neutrophils'))
# 
# #Create cellchat object from seurat object
# cellchat <- createCellChat(object = pbmc, group.by = "fake_celltypes", assay = 'RNA')
# cellchat <- setIdent(cellchat, ident.use = "fake_celltypes")
# levels(cellchat@idents)
# 
# # Initialize database
# CellChatDB <- CellChatDB.human
# showDatabaseCategory(CellChatDB)
# cellchat@DB <- CellChatDB
# 
# #Preprocess expression data
# cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
# future::plan("multisession", workers = 1) # doing in parellel leads to memory issues a few lines down. so 1 for now
# cellchat <- identifyOverExpressedGenes(cellchat)
# 
# #Need to allow greater ram usage to run pca integration
# cellchat <- identifyOverExpressedInteractions(cellchat)
# 
# #Compute the communication probability - pretty slow step, takes around 10 mins on laptop
# cellchat <- computeCommunProb(cellchat, type = "triMean")
# 
# #Infer pathway level
# cellchat <- computeCommunProbPathway(cellchat)
# 
# #Look at specific cell type interactions
# cellchat_NK <- aggregateNet_bug_removed(cellchat, target = "NK Cells")
# groupSize <- as.numeric(table(cellchat_NK@idents))
# 
# #ERROR HERE
# netVisual_circle(cellchat_NK@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
# 
