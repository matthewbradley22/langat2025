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

#How many wt infected cells from each cell type
ParseSeuratObj_int[[]] %>% dplyr::group_by(Genotype, manualAnnotation, hasVirus) %>% dplyr::summarise(total = n()) %>% 
  dplyr::filter(Genotype == 'WT') %>% ggplot(aes(x = manualAnnotation, y = total, fill = factor(hasVirus)))+
  geom_bar(stat = 'identity', position = 'dodge')+
  geom_text(aes(label=total), vjust=0, position = position_dodge(1)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  xlab('')+
  labs(fill = 'has virus')


#Molecular signatures database
all_gene_sets <- msigdbr(species = "Mus musculus")

mouse_gene_sets <- all_gene_sets %>%
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
    #Don't think findMarkers test matters because we are only using Log fold change for rankigns
    wtInfectComp <- FindMarkers(wt, 
                                group.by = 'hasVirus', 
                                ident.1 = 1)
                                # min.pct = 0.05)
                                 #test.use = 'MAST')
  }
  
  if(length(table(IPS$hasVirus)) == 2 & all(table(IPS$hasVirus) > 30)){
    IPSInfectComp <- FindMarkers(IPS, 
                                 group.by = 'hasVirus', 
                                 ident.1 = 1)
                                 #min.pct = 0.05)
                                 #test.use = 'MAST')
  }
 
  if(exists('wtInfectComp')){
    wtInfectComp <- wtInfectComp[wtInfectComp$pct.1 > 0 | wtInfectComp$pct.2 > 0,]
    wtInfectComp <- wtInfectComp[order(dplyr::desc(wtInfectComp$avg_log2FC)),]
    wtRankings <- wtInfectComp$avg_log2FC
    names(wtRankings) <- rownames(wtInfectComp)
    
    WtGSEAres <- fgsea(pathways = gene_sets, # List of gene sets to check
                       stats = wtRankings,
                       scoreType = 'pos', # in this case we have both pos and neg rankings. if only pos or neg, set to 'pos', 'neg'
                       minSize = 15,
                       maxSize = 500) 
    interferonRowWt <- WtGSEAres[grep('HALLMARK_INTERFERON_ALPHA_RESPONSE', WtGSEAres$pathway),]
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
    interferonRowIps <- IpGSEAres[grep('HALLMARK_INTERFERON_ALPHA_RESPONSE', IpGSEAres$pathway),]
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


#Do above, but not split by celltype
wt <- subset(ParseSeuratObj_int, Genotype == 'WT')
ips1 <- subset(ParseSeuratObj_int, Genotype == 'IPS1')
genoCellComp <- FindMarkers(ParseSeuratObj_int, 
                            group.by = 'hasVirus', 
                            ident.1 = ident1_label,
                            min.pct = 0.05)
genoCellComp <- genoCellComp[genoCellComp$pct.1 > 0 | genoCellComp$pct.2 > 0,]
DotPlot(ips1Astro, group.by = 'hasVirus', features = c('Rsad2', 'Ifitm3', 'Oas2', 'Isg20'), scale = FALSE)

#Maybe, for a few cell types (split by wt and ips1), just look at degs between infected and uninfected,
#subset to positive, and map how many are in ISG list

wtImNeu = subset(wt, manualAnnotation == 'Immature Neurons')
ips1ImNeu = subset(ips1, manualAnnotation == 'Immature Neurons')
wtEpen = subset(wt, manualAnnotation == 'Ependymal')
ips1Epen = subset(ips1, manualAnnotation == 'Ependymal')

#Get an idea if wt infections are increasing isg's
manualInspectionIfnGenes <- function(celltype){
  if(!exists("mouse_gene_sets")){
    return("No gene set loaded")
  }
  wtCell = subset(wt, manualAnnotation == celltype)
  ips1Cell = subset(ips1, manualAnnotation == celltype)
  
  wtCellComp <- FindMarkers(wtCell, 
                             group.by = 'hasVirus', 
                             ident.1 = 1,
                             test.use = 'MAST')
  ips1CellComp <- FindMarkers(ips1Cell, 
                               group.by = 'hasVirus', 
                               ident.1 = 1,
                               test.use = 'MAST')
  #Get total number of positive DEGs between infected and uninfected for both genotypes
  wtCellComp_pos <- wtCellComp[wtCellComp$p_val_adj < 0.05 & wtCellComp$avg_log2FC > 0,]
  ips1CellComp_pos <- ips1CellComp[ips1CellComp$p_val_adj < 0.05 & ips1CellComp$avg_log2FC > 0,]
  
  newRow = data.frame(cellType = celltype, wtPos = nrow(wtCellComp_pos), ipsPos = nrow(ips1CellComp_pos))

  
  #Get total number of positive DEGs as above, but subset to ISGs
  ifnA_response <- mouse_gene_sets$HALLMARK_INTERFERON_ALPHA_RESPONSE
  ifnA_GOBP_response <- mouse_gene_sets$GOBP_RESPONSE_TO_INTERFERON_ALPHA
  type1_response <- mouse_gene_sets$GOBP_RESPONSE_TO_TYPE_I_INTERFERON
  reactome_ifn_antiviral <- mouse_gene_sets$REACTOME_ANTIVIRAL_MECHANISM_BY_IFN_STIMULATED_GENES
  all_ISGs = unique(c(ifnA_response, ifnA_GOBP_response, type1_response, reactome_ifn_antiviral))
  
  wt_upreg_isgs <- sum(rownames(wtCellComp_pos) %in% all_ISGs)
  ips_upreg_isgs <- sum(rownames(ips1CellComp_pos) %in% all_ISGs)
  newRow_isg = data.frame(cellType = celltype, wtPos = wt_upreg_isgs, ipsPos = ips_upreg_isgs)
  
  return(list('allDEGS' = newRow, 'ISG_DEGs' = newRow_isg))
  print(paste('Done with', celltype))
  #Just running beginning for now
  if(FALSE){
   
    print(paste("Number of wt positive DEGs", nrow(wtCellComp_pos)))
    print(paste("Number of IPS positive DEGs", nrow(ips1CellComp_pos)))
    
    
    print(paste("Number of wt DEGs in ifn response list", sum(rownames(wtCellComp_pos) %in% ifnA_response)))
    print(paste("Number of IPS DEGs in ifn response list", sum(rownames(ips1CellComp_pos) %in% ifnA_response)))
    print(paste("Number of wt DEGs in ifn GOBP response list", sum(rownames(wtCellComp_pos) %in% ifnA_GOBP_response)))
    print(paste("Number of IPS DEGs in ifn GOBP response list", sum(rownames(ips1CellComp_pos) %in% ifnA_GOBP_response)))
    print(paste("Number of wt DEGs in ifn ISG Viral response list", sum(rownames(wtCellComp_pos) %in% reactome_ifn_antiviral)))
    print(paste("Number of IPS DEGs in ifn ISG Viral response list", sum(rownames(ips1CellComp_pos) %in% reactome_ifn_antiviral)))
    
    geneLocations <- which(rownames(dplyr::arrange(wtAstroComp, desc(avg_log2FC))) %in% ifnA_response)
    geneLocaitonsIps <- which(rownames(dplyr::arrange(ips1CellComp_pos, desc(avg_log2FC))) %in% ifnA_response)
    
    p1 = ggplot(data.frame(geneLocations, y = 1), aes(x = geneLocations, y = y))+
      geom_point()+
      xlim(0, nrow(wtAstroComp))+
      xlab('IFN gene LFC Ranking')+
      ylab('')
    p2 = ggplot(data.frame(geneLocations, y = 1), aes(x = geneLocations, y = y))+
      geom_point()+
      xlim(0, nrow(ips1AstroComp))+
      xlab('IFN gene LFC Ranking')+
      ylab('')
    
    #fgsea
    WtRankings <- wtCellComp$avg_log2FC
    names(WtRankings) <- rownames(wtCellComp)
    WtRankings <- WtRankings[!is.na(WtRankings)]
    WtRankings <- WtRankings[names(WtRankings) != '']
    WtRankings <- WtRankings[!duplicated(names(WtRankings))]
    
    GSEAres <- fgsea(pathways = mouse_gene_sets, # List of gene sets to check
                     stats = WtRankings,
                     scoreType = 'pos', # in this case we have both pos and neg rankings. if only pos or neg, set to 'pos', 'neg'
                     minSize = 10,
                     maxSize = 500)
    
    IpsRankings <- ips1CellComp$avg_log2FC
    names(IpsRankings) <- rownames(ips1CellComp)
    IpsRankings <- IpsRankings[!is.na(IpsRankings)]
    IpsRankings <- IpsRankings[names(IpsRankings) != '']
    IpsRankings <- IpsRankings[!duplicated(names(IpsRankings))]
    
    Ips_GSEAres <- fgsea(pathways = mouse_gene_sets, # List of gene sets to check
                         stats = IpsRankings,
                         scoreType = 'pos', # in this case we have both pos and neg rankings. if only pos or neg, set to 'pos', 'neg'
                         minSize = 10,
                         maxSize = 500)
    
    View(GSEAres)
    
  }
  
}

celltypes <- unique(ParseSeuratObj_int$manualAnnotation)
#Run one at at time - for loop crashed r for some reason
EndoUp <- manualInspectionIfnGenes(celltypes[1])
MicroUp <- manualInspectionIfnGenes(celltypes[2])
oligoUp <- manualInspectionIfnGenes(celltypes[3])
astroUp <- manualInspectionIfnGenes(celltypes[4])
manualInspectionIfnGenes(celltypes[5])
manualInspectionIfnGenes(celltypes[6])
epenUp <- manualInspectionIfnGenes(celltypes[7])
manualInspectionIfnGenes(celltypes[8])
manualInspectionIfnGenes(celltypes[9])
manualInspectionIfnGenes(celltypes[10])
manualInspectionIfnGenes(celltypes[11])
manualInspectionIfnGenes(celltypes[12])
manualInspectionIfnGenes(celltypes[13])
manualInspectionIfnGenes(celltypes[14])
manualInspectionIfnGenes(celltypes[15])
manualInspectionIfnGenes(celltypes[16])

#Compare wt isg directly w ips isgs, all infected
testing <- subset(ParseSeuratObj_int, manualAnnotation == 'Astrocytes' & hasVirus == 1)
infectedISGComp <- FindMarkers(testing, 
                            group.by = 'Genotype', 
                            ident.1 = 'WT',
                            test.use = 'MAST')
infectedISGComp_pos <- infectedISGComp[infectedISGComp$p_val_adj < 0.05 & infectedISGComp$avg_log2FC > 0,]
infectedISGComp_neg <- infectedISGComp[infectedISGComp$p_val_adj < 0.05 & infectedISGComp$avg_log2FC < 0,]
sum(rownames(infectedISGComp_pos) %in% all_ISGs)
sum(rownames(infectedISGComp_neg) %in% all_ISGs)
#Check simple plot
ifnStimGenes <- mouse_gene_sets$REACTOME_ANTIVIRAL_MECHANISM_BY_IFN_STIMULATED_GENES
ParseSeuratObj_int <- AddModuleScore(ParseSeuratObj_int, features = list(ifnStimGenes),
                                     assay = 'RNA', name = 'ifnStimGenes')
infected <- subset(ParseSeuratObj_int, hasVirus == 1 & Treatment %in% c('rChLGTV', 'rLGTV'))
VlnPlot(ParseSeuratObj_int, features = 'ifnStimGenes1', group.by = 'Treatment', split.by = 'Genotype',
        pt.size = 0) +
  ggtitle('Interferon reactome score relative expression')
ParseSeuratObj_int[[]] %>% filter(Genotype == 'WT', hasVirus == 1, Treatment == 'rLGTV') %>% count(manualAnnotation)

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

#Look at select isg in cells
ifnA_response <- mouse_gene_sets$HALLMARK_INTERFERON_ALPHA_RESPONSE
ifnA_GOBP_response <- mouse_gene_sets$GOBP_RESPONSE_TO_INTERFERON_ALPHA
type1_response <- mouse_gene_sets$GOBP_RESPONSE_TO_TYPE_I_INTERFERON
reactome_ifn_antiviral <- mouse_gene_sets$REACTOME_ANTIVIRAL_MECHANISM_BY_IFN_STIMULATED_GENES
all_ISGs = unique(c(ifnA_response, ifnA_GOBP_response, type1_response, reactome_ifn_antiviral))

ParseSeuratObj_int = AddModuleScore(ParseSeuratObj_int, features = list(all_ISGs), name = 'ISG_score') 
VlnPlot(subset(ParseSeuratObj_int, Timepoint == 'Day 5'), features = 'ISG_score1', group.by = 'Genotype', pt.size = 0,
        split.by = 'hasVirus')+
  ggtitle('ISG Score day 5')
  
parseSub = subset(ParseSeuratObj_int, Timepoint == 'Day 5')
table(parseSub$Genotype, parseSub$hasVirus)

FeaturePlot(ParseSeuratObj_int, features = 'ISG_score1', reduction = 'umap.integrated') + ggtitle('ISG Expression')

