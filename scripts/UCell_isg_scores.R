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
wt <-  subset(ParseSeuratObj_int, Treatment %in% c('PBS', 'rChLGTV') & Genotype == 'WT')

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
wt <- AddModuleScore_UCell(wt, features = list(all_ISGs_type1), name = 'ifna_response')

#Add column to split by treatment and celltype
wt$timepoint_celltype_treatment <- paste(wt$Timepoint, wt$manualAnnotation, wt$Treatment, sep = '_')
wt_isg_dat <- DotPlot(wt, features = 'signature_1ifna_response', group.by = 'timepoint_celltype_treatment', scale = FALSE)$data

#Split id column into three
isg_meta <- str_split_fixed(wt_isg_dat$id, "_", 3)
colnames(isg_meta) <- c('timepoint', 'celltype', 'treatment')
wt_isg_dat <- cbind(wt_isg_dat, isg_meta)

#Reorder celltypes
wt_isg_dat <- dplyr::filter(wt_isg_dat, celltype != 'unknown')

wt_isg_dat$celltype <- factor(wt_isg_dat$celltype, levels = rev(c('Astrocytes', 'Choroid Plexus', 'Endothelial', 'Ependymal',
                                                                  'Immature Neurons', 'Microglia', 'Muscle cells', 'Neurons',
                                                                  'Oligodendrocytes', 'Pericytes', 'B Cells', 'Granulocytes',
                                                                  'Macrophage/Monocytes', 'Nk cells', 'T cells')))

#Plot isg scores
ggplot(wt_isg_dat, aes(x = timepoint, y = celltype, fill = avg.exp.scaled))+
  geom_tile()+
  facet_grid(~treatment, scales = 'free')+
  scale_fill_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white"), 
                       values = c(1.0,0.7,0.4,0),
                       limits = c(0,0.25))+
  geom_text(aes(label= round(avg.exp.scaled, digits = 2)))+
  ggtitle('ISG scores sc wt data')+
  theme(panel.background = element_rect(fill = 'white', colour = 'white'))

#Compare wt mock and ips mock across celltypes
mock <- subset(ParseSeuratObj_int, Treatment %in% c('PBS'))
mock <- AddModuleScore_UCell(mock, features = list(all_ISGs_type1), name = 'ifna_response')

#Add column to split by treatment and celltype
mock$genotype_celltype_treatment <- paste(mock$Genotype, mock$manualAnnotation, mock$Treatment, sep = '_')
mock_isg_dat <- DotPlot(mock, features = 'signature_1ifna_response', group.by = 'genotype_celltype_treatment', scale = FALSE)$data

#Split id column into three
isg_mock_meta <- str_split_fixed(mock_isg_dat$id, "_", 3)
colnames(isg_mock_meta) <- c('genotype', 'celltype', 'treatment')
mock_isg_dat <- cbind(mock_isg_dat, isg_mock_meta)

#Reorder celltypes
mock_isg_dat <- dplyr::filter(mock_isg_dat, celltype != 'unknown')

mock_isg_dat$celltype <- factor(mock_isg_dat$celltype, levels = rev(c('Astrocytes', 'Choroid Plexus', 'Endothelial', 'Ependymal',
                                                                  'Immature Neurons', 'Microglia', 'Muscle cells', 'Neurons',
                                                                  'Oligodendrocytes', 'Pericytes', 'B Cells', 'Granulocytes',
                                                                  'Macrophage/Monocytes', 'Nk cells', 'T cells')))

#Plot isg scores
ggplot(mock_isg_dat, aes(x = genotype, y = celltype, fill = avg.exp.scaled))+
  geom_tile()+
  facet_grid(~treatment, scales = 'free')+
  scale_fill_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white"), 
                       values = c(1.0,0.7,0.4,0),
                       limits = c(0,0.25))+
  geom_text(aes(label= round(avg.exp.scaled, digits = 2)))+
  ggtitle('ISG scores sc wt mock data')+
  theme(panel.background = element_rect(fill = 'white', colour = 'white'))


