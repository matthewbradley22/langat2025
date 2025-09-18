#Load libraries
library(ReactomeGSA)
library(msigdbr)
library(fgsea)
library(dplyr)

#Load data
ParseSeuratObj_int <- LoadSeuratRds("./data/FilteredRpcaIntegratedDatNoDoublets.rds") 
ParseSeuratObj_int$hasVirus = ifelse(ParseSeuratObj_int$virusCountPAdj >= 10, 1, 0)

#Check data
newCols <-  c(brewer.pal(12, 'Paired'), '#99FFE6', '#CE99FF', '#18662E',  '#FF8AEF','#737272')
DimPlot(ParseSeuratObj_int, label = FALSE, group.by = 'manualAnnotation', reduction = 'umap.integrated',
        cols = newCols)

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


#Molecular signatures database
all_gene_sets <- msigdbr(species = "Mus musculus")

gene_sets <- all_gene_sets %>%
  split(x = .$gene_symbol, f = .$gs_name)

cat(paste("Loaded", length(gene_sets), "gene sets\n"))

#Need a ranked gene list for this, so create DEG list then rank by
#most upregulated to most downregulated
#Think I need to make a new column that combines genotype + manual annotation cell type
#Then for each group run Findmarkers against all others, run to fgsea, and get NES value
#Will also check against simple addModuleScore for reactome pathway 
ParseSeuratObj_int$genotypeCellType <- paste0(ParseSeuratObj_int$Genotype, ParseSeuratObj_int$manualAnnotation)
genotypeComp <- FindMarkers(ParseSeuratObj_int, group.by = 'Genotype', ident.1 = 'IPS1', ident.2 = 'WT')
genotypeComp <- genotypeComp[order(desc(genotypeComp$avg_log2FC)),]
rankings <- genotypeComp$avg_log2FC
names(rankings) <- rownames(genotypeComp)
head(rankings)

#Run gsea
GSEAres <- fgsea(pathways = gene_sets, # List of gene sets to check
                 stats = rankings,
                 scoreType = 'std', # in this case we have both pos and neg rankings. if only pos or neg, set to 'pos', 'neg'
                 minSize = 15,
                 maxSize = 500) 
(GSEAres[order(pval), ]) %>% as.data.frame()
GSEAres[grep('INTERFERON', GSEAres$pathway),] %>% dplyr::arrange(padj)

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
VlnPlot(ParseSeuratObj_int, features = 'interferonReactomeScore1', group.by = 'Genotype',
        pt.size = 0) +
  theme(legend.position = 'None')

