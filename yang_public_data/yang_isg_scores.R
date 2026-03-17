library(Seurat)
library(dplyr)
library(ggplot2)
library(UCell)

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
yang_data <- AddModuleScore_UCell(yang_data, features = list(all_ISGs_type1), name = 'ifna_response')

#Add column to split by treatment and celltype
yang_data$treatment_celltype <- paste(yang_data$treatment, yang_data$manualAnnotation, sep = '_')
yang_isg_dat <- DotPlot(yang_data, features = 'signature_1ifna_response', group.by = 'treatment_celltype', scale = FALSE)$data

#Split id column into two
isg_meta <- str_split_fixed(yang_isg_dat$id, "_", 2)
colnames(isg_meta) <- c('treatment', 'celltype')
yang_isg_dat <- cbind(yang_isg_dat, isg_meta)

#Plot isg scores
ggplot(yang_isg_dat, aes(x = treatment, y = celltype, fill = avg.exp.scaled))+
  geom_tile()+
  scale_fill_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white"), 
                       values = c(1.0,0.7,0.4,0),
                       limits = c(0,0.25))+
  geom_text(aes(label= round(avg.exp.scaled, digits = 2)))+
  ggtitle('ISG scores Yang data')

table(yang_data$treatment, yang_data$manualAnnotation)

#Which genes are most correlated with gene score?
yang_data <- JoinLayers(yang_data)
genes_for_cor <- t(yang_data[['RNA']]$data)
gene_isg_cors <- cor(yang_data$signature_1ifna_response, genes_for_cor, method = 'pearson')
gene_isg_cors_df <- t(gene_isg_cors) %>% as.data.frame() %>% dplyr::arrange(desc(V1))
head(gene_isg_cors_df, n = 50)


########### Look at specific genes ########### 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## #


yang_isg_gene_dat <- DotPlot(yang_data, features = c('Isg15'), group.by = 'treatment_celltype', scale = FALSE)$data

#Split id column into two
isg_meta <- str_split_fixed(yang_isg_gene_dat$id, "_", 2)
colnames(isg_meta) <- c('treatment', 'celltype')
yang_isg_gene_dat <- cbind(yang_isg_gene_dat, isg_meta)

#Plot isg scores
ggplot(yang_isg_gene_dat, aes(x = treatment, y = celltype, fill = avg.exp.scaled))+
  geom_tile()+
  scale_fill_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white"), 
                       values = c(1.0,0.7,0.4,0),
                       limits = c(0,5))+
  geom_text(aes(label= round(avg.exp.scaled, digits = 2)))+
  ggtitle('Isg15 Yang data')
