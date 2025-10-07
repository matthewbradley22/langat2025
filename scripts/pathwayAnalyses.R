#Load libraries
library(ReactomeGSA)
library(msigdbr)
library(fgsea)
library(dplyr)
library(Seurat)
library(RColorBrewer)
library(MAST)
library(ggplot2)
library(readr)
library(ReactomePA)
library(clusterProfiler)

#Loading to try new mousser analysis
library(mousseR)
library(stringr)

#Load data
ParseSeuratObj_int <- LoadSeuratRds("./data/FilteredRpcaIntegratedDatNoDoublets.rds") 
ParseSeuratObj_int$hasVirus = ifelse(ParseSeuratObj_int$virusCountPAdj >= 10, 1, 0)

#Check data
newCols <-  c(brewer.pal(12, 'Paired'), '#99FFE6', '#CE99FF', '#18662E',  '#FF8AEF','#737272')
DimPlot(ParseSeuratObj_int, label = FALSE, group.by = 'manualAnnotation', reduction = 'umap.integrated',
        cols = newCols)


#Molecular signatures database
all_gene_sets <- msigdbr(species = "Mus musculus")

gene_sets <- all_gene_sets %>%
  split(x = .$gene_symbol, f = .$gs_name)


#Maybe also try with subset to just significant genes and see if there's a bif difference
ParseSeuratObj_int$genotypeCellType <- paste0(ParseSeuratObj_int$Genotype, ParseSeuratObj_int$manualAnnotation)
uniqueGenoCell <- unique(ParseSeuratObj_int$genotypeCellType)

#Initiate with same columns as returned by fgsea
interferonDat <- data.frame(matrix(ncol = 9, nrow = 0))
colnames(interferonDat) = c("pathway", "pval", "padj", "log2err",  "ES" ,  "NES" ,"size" ,     
                            "leadingEdge", 'groupLab')

#Also should do this a different way to check
#Compare infected vs uninfected in both groups and see differences between groups

#Several gene sets to try with:
mouse_gene_sets[grep('INTERFERON' , names(mouse_gene_sets))]
#Want to try at least REACTOME_ANTIVIRAL_MECHANISM_BY_IFN_STIMULATED_GENES and HALLMARK_INTERFERON_ALPHA_RESPONSE ad
#REACTOME_INTERFERON_ALPHA_BETA_SIGNALING

for(i in 1:length(uniqueGenoCell)){
  ident1_label = uniqueGenoCell[i]
  
  #Need a ranked gene list for this, so create DEG list then rank by
  #most upregulated to most downregulated
  
  # Just using this for average log fold change, so don't think test used should matter, same below
  genoCellComp <- FindMarkers(ParseSeuratObj_int, 
              group.by = 'genotypeCellType', 
              ident.1 = ident1_label,
              min.pct = 0.05)
  genoCellComp <- genoCellComp[genoCellComp$pct.1 > 0 | genoCellComp$pct.2 > 0,]
  # genoCellComp <- FindMarkers(ParseSeuratObj_int, 
  #                             group.by = 'genotypeCellType', 
  #                             ident.1 = ident1_label,
  #                             min.pct = 0.1,
  #                             test.use = 'MAST')
  
  #Trying with only sig genes
  #Can also try with different ranking
  #genoCellComp = genoCellComp[genoCellComp$p_val_adj < 0.05,]
  genoCellComp <- genoCellComp[order(dplyr::desc(genoCellComp$avg_log2FC)),]
  rankings <- genoCellComp$avg_log2FC
  names(rankings) <- rownames(genoCellComp)
  
  #Run gsea
  GSEAres <- fgsea(pathways = gene_sets, # List of gene sets to check
                   stats = rankings,
                   scoreType = 'pos', # in this case we have both pos and neg rankings. if only pos or neg, set to 'pos', 'neg'
                   minSize = 15,
                   maxSize = 500) 
  interferonRow <- GSEAres[grep('REACTOME_ANTIVIRAL_MECHANISM_BY_IFN_STIMULATED_GENES', GSEAres$pathway),]
  interferonRow$groupLab <- ident1_label
  print(paste("Done with", ident1_label))
  interferonDat <- rbind(interferonDat, interferonRow)
  
} 

ggplot(interferonDat, aes(x = groupLab, y = pathway, size = -log10(padj),
                          color = NES))+
  geom_point()+
  theme(axis.text.x = element_text(angle = 90))

#Compare infected vs uninfected in each cell type split by genotype

cellTypes <- unique(ParseSeuratObj_int$manualAnnotation)
interferonDatInfected <- data.frame(matrix(ncol = 9, nrow = 0))
colnames(interferonDatInfected) = c("pathway", "pval", "padj", "log2err",  "ES" ,  "NES" ,"size" ,     
                            "leadingEdge", 'groupLab')

for(i in 1:length(uniqueGenoCell)){
  cellTypeObj <- subset(ParseSeuratObj_int, manualAnnotation ==  cellTypes[i])
  wt <- subset(cellTypeObj, Genotype == 'WT')
  IPS <- subset(cellTypeObj, Genotype == 'IPS1')
  
  #Need a ranked gene list for this, so create DEG list then rank by
  #most upregulated to most downregulated
  if(length(table(wt$hasVirus)) == 2 & all(table(wt$hasVirus) > 30)){
    # wtInfectComp <- Seurat::FoldChange(wt, 
    #                                    group.by = 'hasVirus', 
    #                                    assay = 'RNA',
    #                                    ident.1 = 1)
    
    #Using findMarkers because I don't think min.pct is working for foldChange
    wtInfectComp <- FindMarkers(wt, 
                                group.by = 'hasVirus', 
                                ident.1 = 1,
                                 min.pct = 0.05,
                                 test.use = 'MAST')
  }
  
  if(length(table(IPS$hasVirus)) == 2 & all(table(IPS$hasVirus) > 30)){
    IPSInfectComp <- FindMarkers(IPS, 
                                 group.by = 'hasVirus', 
                                 ident.1 = 1,
                                 min.pct = 0.05,
                                 test.use = 'MAST')
  }
 
  if(exists('wtInfectComp')){
    wtInfectComp <- wtInfectComp[wtInfectComp$pct.1 > 0.05 | wtInfectComp$pct.2 > 0.05,]
    wtInfectComp <- wtInfectComp[order(dplyr::desc(wtInfectComp$avg_log2FC)),]
    wtRankings <- wtInfectComp$avg_log2FC
    names(wtRankings) <- rownames(wtInfectComp)
    
    WtGSEAres <- fgsea(pathways = gene_sets, # List of gene sets to check
                       stats = wtRankings,
                       scoreType = 'pos', # in this case we have both pos and neg rankings. if only pos or neg, set to 'pos', 'neg'
                       minSize = 15,
                       maxSize = 500) 
    interferonRowWt <- WtGSEAres[grep('REACTOME_ANTIVIRAL_MECHANISM_BY_IFN_STIMULATED_GENES', WtGSEAres$pathway),]
    interferonRowWt$groupLab <- unique(paste0(wt$manualAnnotation, wt$Genotype))
    interferonDatInfected <- rbind(interferonDatInfected, interferonRowWt)
  }
  
  if(exists('IPSInfectComp')){
    IPSInfectComp <- IPSInfectComp[IPSInfectComp$pct.1 > 0 | IPSInfectComp$pct.2 > 0,]
    IPSInfectComp <- IPSInfectComp[order(dplyr::desc(IPSInfectComp$avg_log2FC)),]
    IpsRankings <- IPSInfectComp$avg_log2FC
    names(IpsRankings) <- rownames(IPSInfectComp)
    
    IpGSEAres <- fgsea(pathways = gene_sets, # List of gene sets to check
                       stats = IpsRankings,
                       scoreType = 'pos', # in this case we have both pos and neg rankings. if only pos or neg, set to 'pos', 'neg'
                       minSize = 15,
                       maxSize = 500) 
    interferonRowIps <- IpGSEAres[grep('REACTOME_ANTIVIRAL_MECHANISM_BY_IFN_STIMULATED_GENES', IpGSEAres$pathway),]
    interferonRowIps$groupLab <- unique(paste0(IPS$manualAnnotation, IPS$Genotype))
    interferonDatInfected <- rbind(interferonDatInfected, interferonRowIps)
  }
  
  
  # genoCellComp <- FindMarkers(ParseSeuratObj_int, 
  #                             group.by = 'genotypeCellType', 
  #                             ident.1 = ident1_label,
  #                             min.pct = 0.1,
  #                             test.use = 'MAST')
  
  #Trying with only sig genes
  #Can also try with different ranking
  #genoCellComp = genoCellComp[genoCellComp$p_val_adj < 0.05,]
  #Run gsea
  
  
  print(paste("Done with", cellTypes[i]))
  
} 

dplyr::arrange(interferonDatInfected, NES)

wtSamples <- interferonDatInfected[grep('WT',interferonDatInfected$groupLab),]
ipsSamples <-  interferonDatInfected[grep('IPS',interferonDatInfected$groupLab),]
plotOrder <- c(wtSamples$groupLab, ipsSamples$groupLab)
ggplot(interferonDatInfected, aes(x = factor(groupLab, levels = plotOrder), y = pathway, size = -log10(padj),
                          color = NES))+
  geom_point()+
  theme(axis.text.x = element_text(angle = 90))+
  xlab('cell genotype group')+
  ggtitle('Infected vs Uninfected Interferon Expression')


#Check simple plot
#messy import of genes from https://www.gsea-msigdb.org/gsea/msigdb/mouse/geneset/REACTOME_INTERFERON_SIGNALING.html
#gmt file
Reactome_interferon_genes <- read_table("REACTOME_INTERFERON_SIGNALING.v2025.1.Mm.gmt", 
                                                       skip = 1)
Reactome_interferon_genes <- colnames(Reactome_interferon_genes)
ParseSeuratObj_int <- AddModuleScore(ParseSeuratObj_int, features = list(Reactome_interferon_genes),
                                     assay = 'RNA', name = 'interferonReactomeScore')
infected <- subset(ParseSeuratObj_int, hasVirus == 1 & Treatment %in% c('rChLGTV', 'rLGTV'))
VlnPlot(ParseSeuratObj_int, features = 'interferonReactomeScore1', group.by = 'Treatment', split.by = 'Genotype',
        pt.size = 0) +
  ggtitle('Interferon reactome score relative expression')

VlnPlot(ParseSeuratObj_int, features = 'interferonReactomeScore1', group.by = 'manualAnnotation', split.by = 'Genotype',
        pt.size = 0) 

FeaturePlot(ParseSeuratObj_int, features = 'interferonReactomeScore1', reduction = 'umap.integrated')
VlnPlot(infected, features = 'interferonReactomeScore1', group.by = 'Treatment', split.by = 'Genotype',
        pt.size = 0) +
  ggtitle('Interferon reactome score relative expression')

#Look by celltype
table(ParseSeuratObj_int$Genotype, ParseSeuratObj_int$hasVirus, ParseSeuratObj_int$manualAnnotation)
plotGenotypeReactome <- function(dat, celltype){
  celltype <- subset(dat, manualAnnotation == celltype)
  print(table(celltype$Genotype, celltype$hasVirus))
  print(VlnPlot(celltype, features = 'interferonReactomeScore1', group.by = 'Genotype', split.by = 'hasVirus', 
                alpha = 0.2)) +
    ggtitle(paste('Interferon reactome score', unique(celltype$manualAnnotation)))
}
plotGenotypeReactome(ParseSeuratObj_int, cellTypes[1])


#Look for IFN genes and explore
IfnGenes <- rownames(ParseSeuratObj_int[['RNA']]$data)[grep('Ifn', rownames(ParseSeuratObj_int[['RNA']]$data))]
VlnPlot(ParseSeuratObj_int, features = IfnGenes, pt.size = 0, group.by = 'Genotype',
        assay = 'RNA')
DotPlot(ParseSeuratObj_int, features = IfnGenes,  group.by = 'Genotype', scale = FALSE)
ParseSeuratObj_int <- AddModuleScore(ParseSeuratObj_int, features = list(IfnGenes), name = 'ifnGenes')
VlnPlot(ParseSeuratObj_int, features = 'ifnGenes1', pt.size = 0, group.by = 'genotypeCellType')+
  theme(legend.position = 'None')

ifnExp <- ParseSeuratObj_int[['RNA']]$data[grep('Ifn', rownames(ParseSeuratObj_int[['RNA']]$data)),] 
rowSums(ifnExp >0) / ncol(ifnExp)
  
#ReactomeGSA try
pseudo_bulk_data <- ReactomeGSA::generate_pseudo_bulk_data(ParseSeuratObj_int, group_by = "Genotype")
pseudo_metadata <- generate_metadata(pseudo_bulk_data)

#create new request
my_request <- ReactomeAnalysisRequest(method = "Camera")
my_request <- set_parameters(request = my_request, max_missing_values = 0.5)

# add the pseudo-bulk data as a dataset
my_request <- add_dataset(request = my_request,
                          expression_values = pseudo_bulk_data,
                          sample_data = pseudo_metadata,
                          name = "Pseudo-bulk",
                          type = "rnaseq_counts",
                          comparison_factor = "Group",
                          comparison_group_1 = "IPS1",
                          comparison_group_2 = "WT")
my_request

#Submit analysis
quant_result <- perform_reactome_analysis(my_request, compress = FALSE)

quant_result

# get the pathway-level results
quant_pathways <- pathways(quant_result)
head(quant_pathways)
quant_pathways['R-HSA-913531',]
quant_pathways['R-HSA-1606341',]

#ReactomePA attempt
#Probably will turn this into loop and see common themes.
astrocytes <- subset(ParseSeuratObj_int, manualAnnotation == 'Astrocytes')
wt <- subset(astrocytes, Genotype == 'WT')
IPS <- subset(astrocytes, Genotype == 'IPS1')

infectionMarkers <- Seurat::FindMarkers(wt, group.by = 'hasVirus', ident.1 = 1, test.use = 'MAST', only.pos = TRUE)
infectionMarkersIPS <- Seurat::FindMarkers(IPS, group.by = 'hasVirus', ident.1 = 1, test.use = 'MAST', only.pos = TRUE)

infectionMarkersSig <- infectionMarkersIPS %>% filter(p_val_adj < 0.05)
eID <- bitr(rownames(infectionMarkersSig), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
upPaths <- enrichPathway(gene=eID$ENTREZID, pvalueCutoff = 0.05, readable=TRUE, 
                   organism = 'mouse')
upPaths@result

#Overall comparison of infected wt and infected ips
infected <- subset(ParseSeuratObj_int, hasVirus == 1)
table(infected$Genotype)
genotypeMarkers <- FindMarkers(infected, group.by = 'Genotype', ident.1 = 'WT',only.pos = TRUE, test.use = 'MAST')
genotypeMarkersIPS <- FindMarkers(infected, group.by = 'Genotype', ident.1 = 'IPS1',only.pos = TRUE, test.use = 'MAST')

#plot specific ifn genes
VlnPlot(ParseSeuratObj_int, features = 'Ifnar1', group.by = 'Genotype', split.by = 'hasVirus',
        pt.size = 0, assay = 'RNA')

#try new mousser library
counts.matrix <- as.data.frame(t(as.matrix(ParseSeuratObj_int@assays$RNA@layers$data)))
rownames(counts.matrix) <- colnames(ParseSeuratObj_int)
colnames(counts.matrix) <- str_to_title(rownames(ParseSeuratObj_int))
mousse.scores <- mousse(counts.matrix = counts.matrix, numGenes = 60)

#prepare for plotting
ParseSeuratObj_int[["Mousse"]] <- CreateAssayObject(data = (t(mousse.scores)))
ParseSeuratObj_int <- ScaleData(ParseSeuratObj_int, assay = 'Mousse')

mark <- FindAllMarkers(ParseSeuratObj_int, min.pct = 0.05, test.use = "roc", assay = 'Mousse', group.by = 'hasVirus')
GOI_niche <- mark %>% 
  group_by(cluster) %>% 
  top_n(5, myAUC)
DoHeatmap(ParseSeuratObj_int, features = unique(GOI_niche$gene), size = 2.5, group.by = 'Treatment', assay = 'Mousse')

#LRP8
FeaturePlot(ParseSeuratObj_int, reduction = 'umap.integrated', features = 'Lrp8')
FeaturePlot(ParseSeuratObj_int, reduction = 'umap.integrated', features = 'Rsad2')

VlnPlot(ParseSeuratObj_int, features = 'Rsad2', group.by = 'hasVirus', pt.size = 0)
VlnPlot(ParseSeuratObj_int, features = 'Lrp8', group.by = 'hasVirus', pt.size = 0)

