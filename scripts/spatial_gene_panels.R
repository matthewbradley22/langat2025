#Read in xenium mouse 5k panel
xenium_5k_genes <- read.csv("~/Documents/ÖverbyLab/spatial_seq_genes/XeniumPrimeMouse5Kpan_tissue_pathways_metadata.csv")

#Compare panel to upregulated genes in data, and known gene lists (relevant KEGG / GO / msigdbr pathways)

#First, msigdbr
all_gene_sets <- msigdbr(species = "Mus musculus", db_species = "MM")

mouse_gene_sets <- all_gene_sets %>%
  split(x = .$gene_symbol, f = .$gs_name)

#Get total number of positive DEGs as above, but subset to ISGs
ifnA_response <- mouse_gene_sets$HALLMARK_INTERFERON_ALPHA_RESPONSE
ifnA_GOBP_response <- mouse_gene_sets$GOBP_RESPONSE_TO_INTERFERON_ALPHA
type1_response <- mouse_gene_sets$GOBP_RESPONSE_TO_TYPE_I_INTERFERON
#reactome_ifn_antiviral <- mouse_gene_sets$REACTOME_ANTIVIRAL_MECHANISM_BY_IFN_STIMULATED_GENES
all_ISGs_type1 = unique(c(ifnA_response, ifnA_GOBP_response, type1_response))

sum(all_ISGs_type1 %in% xenium_5k_genes$gene_name)/length(all_ISGs_type1)
all_ISGs_type1[all_ISGs_type1 %in% xenium_5k_genes$gene_name]
all_ISGs_type1[!all_ISGs_type1 %in% xenium_5k_genes$gene_name]

#Missing - Cd74, Ifit and Ifitm genes, Rsad2, oas genes, Ifna genes (only has receptors)

#Check other genes 
#ccl present (missing ccl4?)
xenium_5k_genes[grep('ccl', xenium_5k_genes$gene_name, ignore.case = TRUE),]

#Tnf present
xenium_5k_genes[grep('tnf', xenium_5k_genes$gene_name, ignore.case = TRUE),]

#Has cxcl and cxcr (check if all)
xenium_5k_genes[grep('Cxc', xenium_5k_genes$gene_name, ignore.case = TRUE),]

#Has Cx3cr1 and l1
xenium_5k_genes[grep('Cx3', xenium_5k_genes$gene_name, ignore.case = TRUE),]

#Has some mmp genes
xenium_5k_genes[grep('Mmp', xenium_5k_genes$gene_name, ignore.case = TRUE),]

#Has ripk genes
xenium_5k_genes[grep('Ripk', xenium_5k_genes$gene_name, ignore.case = TRUE),]

#Cell annotation markers
xenium_5k_genes[grep('Snap25', xenium_5k_genes$gene_name, ignore.case = TRUE),]
xenium_5k_genes[grep('Flt1', xenium_5k_genes$gene_name, ignore.case = TRUE),]
xenium_5k_genes[grep('Pecam1', xenium_5k_genes$gene_name, ignore.case = TRUE),]
xenium_5k_genes[grep('Gfap', xenium_5k_genes$gene_name, ignore.case = TRUE),]
xenium_5k_genes[grep('Aqp4', xenium_5k_genes$gene_name, ignore.case = TRUE),]
xenium_5k_genes[grep('Gli3', xenium_5k_genes$gene_name, ignore.case = TRUE),]
xenium_5k_genes[grep('Ctss', xenium_5k_genes$gene_name, ignore.case = TRUE),] #No Ctss
xenium_5k_genes[grep('P2ry12', xenium_5k_genes$gene_name, ignore.case = TRUE),]#No P2ry12
xenium_5k_genes[grep('Cx3cr1', xenium_5k_genes$gene_name, ignore.case = TRUE),]
xenium_5k_genes[grep('Cd74', xenium_5k_genes$gene_name, ignore.case = TRUE),] #No Cd74
xenium_5k_genes[grep('Klra2', xenium_5k_genes$gene_name, ignore.case = TRUE),]#No Klra2




