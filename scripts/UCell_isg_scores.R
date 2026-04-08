#Packages and functions
library(ggpubr)
library(Seurat)
library(UCell)
library(msigdbr)
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
ips <-  subset(ParseSeuratObj_int, Treatment %in% c('PBS', 'rChLGTV') & Genotype == 'IPS1')

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
ips <- AddModuleScore_UCell(ips, features = list(all_ISGs_type1), name = 'ifna_response')

#Add column to split by treatment and celltype
wt$timepoint_celltype_treatment <- paste(wt$Timepoint, wt$manualAnnotation, wt$Treatment, sep = '_')
ips$timepoint_celltype_treatment <- paste(ips$Timepoint, ips$manualAnnotation, ips$Treatment, sep = '_')

create_isg_heatmap <- function(dat, group_by_arg, title = NULL){
  isg_dat <- DotPlot(dat, features = 'signature_1ifna_response', group.by = group_by_arg, scale = FALSE)$data
  #Split id column into three
  isg_meta <- str_split_fixed(isg_dat$id, "_", 3)
  colnames(isg_meta) <- c('timepoint', 'celltype', 'treatment')
  isg_dat <- cbind(isg_dat, isg_meta)
  
  #Reorder celltypes
  isg_dat <- dplyr::filter(isg_dat, celltype != 'unknown')
  
  isg_dat$celltype <- factor(isg_dat$celltype, levels = rev(c('Astrocytes', 'Choroid Plexus', 'Endothelial', 'Ependymal',
                                                                    'Immature Neurons', 'Microglia', 'Muscle cells', 'Neurons',
                                                                    'Oligodendrocytes', 'Pericytes', 'B Cells', 'Granulocytes',
                                                                    'Macrophage/Monocytes', 'Nk cells', 'T cells')))
  #Plot isg scores
  isg_heatmap <- ggplot(isg_dat, aes(x = timepoint, y = celltype, fill = avg.exp.scaled))+
    geom_tile()+
    facet_grid(~treatment, scales = 'free')+
    scale_fill_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white"), 
                         values = c(1.0,0.7,0.4,0),
                         limits = c(0,0.25))+
    geom_text(aes(label= round(avg.exp.scaled, digits = 2)))+
    ggtitle(title)+
    theme(panel.background = element_rect(fill = 'white', colour = 'white'))
  print(isg_heatmap)
}

create_isg_heatmap(wt, group_by_arg = 'timepoint_celltype_treatment')

wt_cerebrum <- subset(wt, Organ == 'Cerebrum')
wt_cerebellum <- subset(wt, Organ == 'Cerebellum')

pdf("~/Documents/ÖverbyLab/single_cell_ISG_figures/isg_fig_plots/wt_cerebrum_isg.pdf", width = 7, height = 6)
create_isg_heatmap(wt_cerebrum, group_by_arg = 'timepoint_celltype_treatment', title = 'wt cerebrum ISG scores')
dev.off()

pdf("~/Documents/ÖverbyLab/single_cell_ISG_figures/isg_fig_plots/wt_cerebellum_isg.pdf", width = 7, height = 6)
create_isg_heatmap(wt_cerebellum, group_by_arg = 'timepoint_celltype_treatment', title = 'wt cerebellum ISG scores')
dev.off()

pdf("~/Documents/ÖverbyLab/single_cell_ISG_figures/isg_fig_plots/ips_isg.pdf", width = 7, height = 6)
create_isg_heatmap(ips, group_by_arg = 'timepoint_celltype_treatment', title = 'IPS1 ISG scores')
dev.off()

#Compare wt mock and ips mock across celltypes
mock <- subset(ParseSeuratObj_int, Treatment %in% c('PBS'))
mock <- AddModuleScore_UCell(mock, features = list(all_ISGs_type1), name = 'ifna_response')

#Add column to split by treatment and celltype
mock$genotype_celltype_treatment <- paste(mock$Genotype, mock$manualAnnotation, mock$Treatment, sep = '_')

pdf("~/Documents/ÖverbyLab/single_cell_ISG_figures/isg_fig_plots/mock_isg.pdf", width = 7, height = 6)
create_isg_heatmap(mock, group_by_arg = 'genotype_celltype_treatment', title = 'Mock ISG scores')
dev.off()

#Which isgs are up in wt, which in ips, which both. Also split by time
#And ideally by celltype

#Data objects created in scRNA_celltype_analysis.R script
degs_by_celltype_wt <- readRDS('~/Documents/ÖverbyLab/cell_upsetPlots/wt_celltype_degs.rds')
degs_by_celltype_ips <- readRDS('~/Documents/ÖverbyLab/cell_upsetPlots/ips_celltype_degs.rds')

#Subset deg lists to just isgs
isg_degs_wt <- lapply(degs_by_celltype_wt, FUN = function(x){
  dat <- as.data.frame(x)
  isgs_only <- x[all_ISGs_type1,]
  isgs_only[!is.na(isgs_only$p_val),]
})

isg_degs_ips <- lapply(degs_by_celltype_ips, FUN = function(x){
  dat <- as.data.frame(x)
  isgs_only <- x[all_ISGs_type1,]
  isgs_only[!is.na(isgs_only$p_val),]
})

#For each celltype, check which genes are upregulated
#in wt or ips infected
celltypes <- names(isg_degs_wt)

for(i in 1:length(isg_degs_wt)){
  celltype_wt_genes <- isg_degs_wt[[celltypes[[i]]]]
  celltype_ips_genes <- isg_degs_ips[[celltypes[[i]]]]
  
  celltype_wt_genes_sig <- subset(celltype_wt_genes, avg_log2FC > 1 & p_val_adj < 0.01)
  celltype_ips_genes_sig <- subset(celltype_ips_genes, avg_log2FC > 1 & p_val_adj < 0.01)
  
  
}


