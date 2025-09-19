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

cat(paste("Loaded", length(gene_sets), "gene sets\n"))

#Maybe also try with subset to just significant genes and see if there's a bif difference
ParseSeuratObj_int$genotypeCellType <- paste0(ParseSeuratObj_int$Genotype, ParseSeuratObj_int$manualAnnotation)
uniqueGenoCell <- unique(ParseSeuratObj_int$genotypeCellType)

#Initiate with same columns as returned by fgsea
interferonDat <- data.frame(matrix(ncol = 9, nrow = 0))
colnames(interferonDat) = c("pathway", "pval", "padj", "log2err",  "ES" ,  "NES" ,"size" ,     
                            "leadingEdge", 'groupLab')

for(i in 1:length(uniqueGenoCell)){
  ident1_label = uniqueGenoCell[i]
  
  #Need a ranked gene list for this, so create DEG list then rank by
  #most upregulated to most downregulated
  genoCellComp <- Seurat::FoldChange(ParseSeuratObj_int, ident.1 = ident1_label, 
                                     group.by = 'genotypeCellType')
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
  interferonRow <- GSEAres[grep('REACTOME_INTERFERON_SIGNALING', GSEAres$pathway),]
  interferonRow$groupLab <- ident1_label
  print(paste("Done with", ident1_label))
  interferonDat <- rbind(interferonDat, interferonRow)
  
} 

interferonDat %>% head()
ggplot(interferonDat, aes(x = groupLab, y = pathway, size = -log10(padj),
                          color = NES))+
  geom_point()+
  theme(axis.text.x = element_text(angle = 90))


(GSEAres[order(pval), ]) %>% as.data.frame()
GSEAres[grep('INTERFERON', GSEAres$pathway),] %>% dplyr::arrange(padj)
GSEAres[grep('ANTIVIRAL', GSEAres$pathway),] %>% dplyr::arrange(padj)

plotEnrichment(gene_sets[["REACTOME_INTERFERON_SIGNALING"]],
               rankings) + labs(title="REACTOME_INTERFERON_SIGNALING")

#Check simple plot
#messy import of genes from https://www.gsea-msigdb.org/gsea/msigdb/mouse/geneset/REACTOME_INTERFERON_SIGNALING.html
#gmt file
Reactome_interferon_genes <- read_table("REACTOME_INTERFERON_SIGNALING.v2025.1.Mm.gmt", 
                                                       skip = 1)
Reactome_interferon_genes <- colnames(Reactome_interferon_genes)
ParseSeuratObj_int <- AddModuleScore(ParseSeuratObj_int, features = list(Reactome_interferon_genes),
                                     assay = 'RNA', name = 'interferonReactomeScore')
infected <- subset(ParseSeuratObj_int, hasVirus == 1)
VlnPlot(ParseSeuratObj_int, features = 'interferonReactomeScore1', group.by = 'genotypeCellType',
        pt.size = 0) +
  theme(legend.position = 'None')
VlnPlot(infected, features = 'interferonReactomeScore1', group.by = 'Genotype',
        pt.size = 0) +
  theme(legend.position = 'None')

#Look for IFN genes and explore
IfnGenes <- rownames(ParseSeuratObj_int[['RNA']]$data)[grep('Ifn', rownames(ParseSeuratObj_int[['RNA']]$data))]
VlnPlot(ParseSeuratObj_int, features = IfnGenes, pt.size = 0, group.by = 'Genotype',
        assay = 'RNA')
DotPlot(ParseSeuratObj_int, features = IfnGenes,  group.by = 'Genotype', scale = FALSE)
ParseSeuratObj_int <- AddModuleScore(ParseSeuratObj_int, features = list(IfnGenes), name = 'ifnGenes')
VlnPlot(ParseSeuratObj_int, features = 'ifnGenes1', pt.size = 0, group.by = 'genotypeCellType')+
  theme(legend.position = 'None')


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

 