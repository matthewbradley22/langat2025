library(Seurat)
library(ggplot2)
library(dplyr)
library(networkD3)
library(ggpubr)

#Load data
ParseSeuratObj_int <- LoadSeuratRds("~/Documents/ÖverbyLab/data/FilteredRpcaIntegratedDatNoDoublets.rds") 
ParseSeuratObj_int$hasVirus = ifelse(ParseSeuratObj_int$virusCountPAdj >= 10, 1, 0)

#Check data
umap_color_list <- c( "#7047A1", "#B370AE","#292270",  "#166DF0","#6D92F8",  "#6DC3F8", "#8a0000","#F76363", "#FF96A2", "#D6644B", 
                      "#F08C3A", "#fdc087","#074F00", "#208d1f","#7bcd79", 
                      "gray")

DimPlot(ParseSeuratObj_int, label = FALSE, group.by = 'manualAnnotation', reduction = 'umap.integrated',
        cols = umap_color_list)+
  theme(axis.ticks = element_blank(),
        axis.text=element_blank(),
        legend.text=element_text(size=17))+
  xlab('Umap1')+
  ylab('Umap2')+
  guides(color=guide_legend(override.aes=list(size=8)))+
  ggtitle('')

#Remove lgtv samples and unlabelled celltypes
chimeric_mock <- subset(ParseSeuratObj_int, Treatment != 'rLGTV' & manualAnnotation != 'unknown')

#Create column for future plotting purposes
chimeric_mock$treatment_timepoint_celltype <- paste(chimeric_mock$Treatment, chimeric_mock$Timepoint, chimeric_mock$manualAnnotation, sep = '_')

#Subset by genotype for plotting
wt_chimeric_mock <- subset(chimeric_mock, Genotype == 'WT' | Treatment == 'PBS')
ips_chimeric_mock <- subset(chimeric_mock, Genotype == 'IPS1' | Treatment == 'PBS')

#Look at markers that may downregulate neuroinflammation
single_gene_dotplot <- function(dat, gene, size_lim = NULL, fill_lims = NULL){
  wt_ccl_dat <- DotPlot(dat, features = c(gene), group.by = 'treatment_timepoint_celltype', scale = FALSE)$data
  gene_plot = wt_ccl_dat %>% 
    tidyr::separate(col = id, into = c('treatment', 'time', 'celltype'), sep = '_') %>% 
    ggplot(aes(x = time, y = celltype, size = pct.exp, fill = avg.exp.scaled))+
    facet_wrap(~treatment, scales = 'free_x')+
    geom_point(pch = 21)+
    scale_fill_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white"),
                         values = c(1.0,0.7,0.4,0), 
                         limits = fill_lims)+
    theme_classic()+
    ggtitle(gene)+
    scale_size(limits = size_lim)
  print(gene_plot)
}

#Possible means of lowering inflammation
#Decreased ccl2, ccr2, cxcl10 - not seen in wt
single_gene_dotplot(wt_chimeric_mock, 'Ccl2')
single_gene_dotplot(ips_chimeric_mock, 'Ccl2')
single_gene_dotplot(wt_chimeric_mock, 'Ccr2')
single_gene_dotplot(wt_chimeric_mock, 'Cxcl10')

#Higher il10, il1rn, tgfb1, il27.
#Il10 goes up in both genotypes, indicating that both begin to regulate the immune response. only slightly more in wt
single_gene_dotplot(wt_chimeric_mock, 'Il10', size_lim = c(0, 24), fill_lims = c(0, 0.75))
single_gene_dotplot(ips_chimeric_mock, 'Il10', size_lim = c(0, 24), fill_lims = c(0, 0.75))
single_gene_dotplot(wt_chimeric_mock, 'Il1rn')
single_gene_dotplot(wt_chimeric_mock, 'Tgfb1')
single_gene_dotplot(wt_chimeric_mock, 'Il27')

#Socs genes. Increases in both
single_gene_dotplot(wt_chimeric_mock, 'Socs3')
single_gene_dotplot(ips_chimeric_mock, 'Socs3')

#Can do gene ontology analysis of day 5 vs day 3 cells in each genotype and see if
#any recovery pathways are unique to one genotype or the other
#maybe not late enough in sampling though
