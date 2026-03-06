#Packages and functions
library(ggpubr)
library(Seurat)
library(UCell)
library(RColorBrewer)
source('~/Documents/ÖverbyLab//scripts/langatFunctions.R')

#Load data
ParseSeuratObj_int <- LoadSeuratRds("~/Documents/ÖverbyLab/data/FilteredRpcaIntegratedDatNoDoublets.rds") 

#Check data
newCols <-  c(brewer.pal(12, 'Paired'), '#99FFE6', '#CE99FF', '#18662E','#737272',  '#FF8AEF')
newCols[11] =  '#FF8AEF'
DimPlot(ParseSeuratObj_int, label = FALSE, group.by = 'manualAnnotation', reduction = 'umap.integrated',
        cols = newCols)

#or load in from csv
#Labelled these macs in another script, need to find which one
false_macs_to_remove <- read.csv("~/Documents/ÖverbyLab/false_macs_to_remove.csv")[[2]]
#should be 409 cells
length(false_macs_to_remove)

ParseSeuratObj_int <- subset(ParseSeuratObj_int, cells = false_macs_to_remove, invert = TRUE)

#Wt cerebrum celltypes across times
wt_cerebrum <-  subset(ParseSeuratObj_int, Treatment %in% c('PBS', 'rLGTV') & Organ == 'Cerebrum' & Genotype == 'WT')

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
wt_cerebrum <- AddModuleScore_UCell(wt_cerebrum, features = list(all_ISGs_type1), name = 'ifna_response')

#Add column to split by treatment and celltype
wt_cerebrum$timepoint_celltype_treatment <- paste(wt_cerebrum$Timepoint, wt_cerebrum$manualAnnotation, wt_cerebrum$Treatment, sep = '_')
wt_cerebrum_isg_dat <- DotPlot(wt_cerebrum, features = 'signature_1ifna_response', group.by = 'timepoint_celltype_treatment', scale = FALSE)$data

#Split id column into two
isg_meta <- str_split_fixed(wt_cerebrum_isg_dat$id, "_", 3)
colnames(isg_meta) <- c('timepoint', 'celltype', 'treatment')
wt_cerebrum_isg_dat <- cbind(wt_cerebrum_isg_dat, isg_meta)

#Plot isg scores
ggplot(wt_cerebrum_isg_dat, aes(x = timepoint, y = celltype, fill = avg.exp.scaled))+
  geom_tile()+
  facet_grid(~treatment, scales = 'free')+
  scale_fill_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white"), 
                       values = c(1.0,0.7,0.4,0),
                       limits = c(0,0.25))+
  geom_text(aes(label= round(avg.exp.scaled, digits = 2)))+
  ggtitle('ISG scores sc wt data')


