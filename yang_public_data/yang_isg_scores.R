library(Seurat)
library(dplyr)
library(ggplot2)

#Load in data
#Load in data
yang_data <- LoadSeuratRds("~/Documents/ÖverbyLab/yang_public_data/yang_rpca_integrated_obj.RDS")
newCols <-  c(brewer.pal(12, 'Paired'), '#99FFE6', '#CE99FF', '#18662E','#737272',  '#FF8AEF')
newCols[11] =  '#FF8AEF'

#Plot data
DimPlot(yang_data, reduction = "umap.rpca", label = FALSE, group.by = 'manualAnnotation', cols = newCols)

#Load in ISGs
#Molecular signatures database
#Should compare using db_species vs not using it
all_gene_sets <- msigdbr(species = "Mus musculus", db_species = "MM")

mouse_gene_sets <- all_gene_sets %>%
  split(x = .$gene_symbol, f = .$gs_name)

#Get total number of positive DEGs as above, but subset to ISGs
ifnA_response <- mouse_gene_sets$HALLMARK_INTERFERON_ALPHA_RESPONSE
ifnA_GOBP_response <- mouse_gene_sets$GOBP_RESPONSE_TO_INTERFERON_ALPHA
type1_response <- mouse_gene_sets$GOBP_RESPONSE_TO_TYPE_I_INTERFERON
#reactome_ifn_antiviral <- mouse_gene_sets$REACTOME_ANTIVIRAL_MECHANISM_BY_IFN_STIMULATED_GENES
all_ISGs_type1 = unique(c(ifnA_response, ifnA_GOBP_response, type1_response))

#Should also try ucell module scores
yang_data <- AddModuleScore(yang_data, features = list(all_ISGs_type1), name = 'ifna_response')

#Add column to split by treatment and celltype
yang_data
DotPlot(yang_data, features = 'ifna_response1', group.by = 'manualAnnotation')$data



