#Packages and functions
library(Seurat)
library(dplyr)
library(Augur)
source('~/Documents/ÖverbyLab//scripts/langatFunctions.R')

#Load data
ParseSeuratObj_int <- LoadSeuratRds("~/Documents/ÖverbyLab/data/FilteredRpcaIntegratedDatNoDoublets.rds") 
ParseSeuratObj_int$manualAnnotation <- factor(ParseSeuratObj_int$manualAnnotation, 
                                              levels = c('Astrocytes', 'Choroid Plexus', 'Endothelial',
                                                         'Ependymal', 'Immature Neurons', 'Microglia', 'Muscle cells',
                                                         'Neurons', 'Oligodendrocytes', 'Pericytes', 'B Cells',
                                                         'Granulocytes', 'Macrophage/Monocytes', 'Nk cells', 
                                                         'T cells', 'unknown'))

#Check data
umap_color_list <- c( "#7047A1", "#B370AE","#292270",  "#166DF0","#6D92F8",  "#6DC3F8", "#8a0000","#F76363", "#FF96A2", 
                      "#D6644B", "#F08C3A", "#fdc087","#074F00", "#208d1f","#7bcd79", 
                      "gray")

DimPlot(ParseSeuratObj_int, label = FALSE, group.by = 'manualAnnotation', reduction = 'umap.integrated',
        cols = umap_color_list)

#Split by organ + other groups and do augur analysis 
#to see which celltypes are most affected by infection

#Can skip loop with 
augur_result_list <- readRDS(file = '~/Documents/ÖverbyLab/augur_results/augur_result_list.rds')

organs <- unique(ParseSeuratObj_int$Organ)
treatment <- c('rLGTV', 'rChLGTV')
genotypes = unique(ParseSeuratObj_int$Genotype)
all_combos <- do.call(expand.grid, list(organs, treatment, genotypes))
colnames(all_combos) = c('organ', 'treatment', 'genotype')

#There are no cerebellum langat samples, so remove these
all_combos <- dplyr::filter(all_combos, organ == 'Cerebrum' | treatment != 'rLGTV')
augur_result_list <- list()

for(i in 1:nrow(all_combos)){
  cur_treatment = as.character(all_combos[i,]$treatment)
  cur_genotype = as.character(all_combos[i,]$genotype)
  cur_organ = all_combos[i,]$organ
  cur_dat <- subset(ParseSeuratObj_int, Treatment %in% c(cur_treatment, 'PBS') & Organ == cur_organ & Genotype == cur_genotype)
  
  #Add columns so that augur package knows what to look for
  cur_dat$cell_type <- factor(cur_dat$manualAnnotation)
  cur_dat$label <- factor(cur_dat$Treatment)
  
  cur_augur = calculate_auc(cur_dat)
  
  augur_result_list[[length(augur_result_list) + 1]] = cur_augur
}

#saveRDS(augur_result_list, file = '~/Documents/ÖverbyLab/augur_results/augur_result_list.rds')
all_combos$label <- paste(all_combos$organ, all_combos$treatment, all_combos$genotype, sep = '_')

celltype_auc <- lapply(augur_result_list, FUN = function(x){
  dplyr::filter(x$AUC, !cell_type %in% c('Immature Neurons', 'Muscle cells', 'Neurons', 'Pericytes',
                                         'B Cells', 'unknown'))
  })

names(celltype_auc) = all_combos$label

augur_plot_list <- lapply(celltype_auc, FUN = function(x){
  ggplot(x, aes(x = reorder(cell_type,auc), y = auc))+
    geom_bar(stat = 'identity', color = 'black', fill = 'lightgrey')+
    theme_classic()+
    theme(axis.text.x = element_text(angle =90))+
    xlab('')
})

pdf('~/Documents/ÖverbyLab/augur_results/augur_plot.pdf', width = 14, height = 8)
ggarrange(plotlist=augur_plot_list, labels = names(celltype_auc))
dev.off()
