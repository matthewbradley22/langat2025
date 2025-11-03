library(Seurat)
library(RColorBrewer)
library(msigdbr)
library(dplyr)
library(matrixStats)
library(ggplot2)

#Load data
ParseSeuratObj_int <- LoadSeuratRds("./data/FilteredRpcaIntegratedDatNoDoublets.rds") 
ParseSeuratObj_int$hasVirus = ifelse(ParseSeuratObj_int$virusCountPAdj >= 10, 1, 0)

#Create new column for plotting later 
ParseSeuratObj_int$time_treatment <- paste(ParseSeuratObj_int$Treatment, ParseSeuratObj_int$Timepoint, sep = '_')

#Check data
newCols <-  c(brewer.pal(12, 'Paired'), '#99FFE6', '#CE99FF', '#18662E',  '#FF8AEF','#737272')
DimPlot(ParseSeuratObj_int, label = FALSE, group.by = 'manualAnnotation', reduction = 'umap.integrated',
        cols = newCols)

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

#Look at select isg in cells
ParseSeuratObj_int = AddModuleScore(ParseSeuratObj_int, features = list(all_ISGs_type1), name = 'ISG_score',
                                    slot = 'data', assay = 'RNA') 

#Subset by genotype for plotting later
ips <- subset(ParseSeuratObj_int, Genotype == 'IPS1')
wt <- subset(ParseSeuratObj_int, Genotype == 'WT')

#Plot ISGs
pdf('~/Documents/ÖverbyLab/scPlots/ISG_vln_ips.pdf', width = 7, height = 4)
VlnPlot(ips, features = 'ISG_score1', group.by = 'Treatment', split.by = 'Timepoint',
        pt.size = 0)+
  ggtitle('ISG Module Scores IPS1')
dev.off()

isg_data <- ParseSeuratObj_int[['RNA']]$data[which(rownames(ParseSeuratObj_int[['RNA']]$data) %in% all_ISGs_type1),]
isg_variances <- rowVars(as.matrix(isg_data))
high_var_ISGs <- rev(names(tail(sort(isg_variances), n = 20)))

pdf('~/Documents/ÖverbyLab/scPlots/ISG_dotplot.pdf', width = 9, height = 4)
DotPlot(ParseSeuratObj_int, group.by = 'time_treatment', features = high_var_ISGs, assay = 'RNA',
        scale = FALSE)+
  theme(axis.text.x = element_text(angle = 90))+
  ggtitle('High variance ISGs')
dev.off()

pdf('~/Documents/ÖverbyLab/scPlots/ISG_dotplot_ips.pdf', width = 9, height = 4)
DotPlot(ips, group.by = 'time_treatment', features = high_var_ISGs, assay = 'RNA',
        scale = FALSE)+
  theme(axis.text.x = element_text(angle = 90))+
  ggtitle('High variance ISGs - IPS1')
dev.off()

pdf('~/Documents/ÖverbyLab/scPlots/ISG_dotplot_wt.pdf', width = 9, height = 4)
DotPlot(wt, group.by = 'time_treatment', features = high_var_ISGs, assay = 'RNA',
        scale = FALSE)+
  theme(axis.text.x = element_text(angle = 90))+
  ggtitle('High variance ISGs - WT')
dev.off()

pdf('~/Documents/ÖverbyLab/scPlots/ISG_Umap.pdf', width = 7, height = 5)
FeaturePlot(ParseSeuratObj_int, features = 'ISG_score1', reduction = 'umap.integrated')
dev.off()

#Compare wt and IPS
ParseSeuratObj_int$time_treatment_genotype <- paste(ParseSeuratObj_int$time_treatment, ParseSeuratObj_int$Genotype, sep = '_')
DotPlot(subset(ParseSeuratObj_int, Treatment != 'PBS'), group.by = 'time_treatment_genotype', features = high_var_ISGs, assay = 'RNA',
        scale = FALSE)+
  theme(axis.text.x = element_text(angle = 90))+
  ggtitle('High variance ISGs')


###############################################################################################
#Other gene markers

#IFN genes
ifn_genes <- rownames(ParseSeuratObj_int[['RNA']]$data)[grep('Ifn', rownames(ParseSeuratObj_int[['RNA']]$data))]
ifn_genes <- ifn_genes[grep('Ifna', ifn_genes)]

pdf('~/Documents/ÖverbyLab/scPlots/IfnDotPlot.pdf', width = 8, height = 4)
DotPlot(ParseSeuratObj_int, features = ifn_genes, scale = FALSE, group.by = 'Timepoint')
dev.off()

#Reln
dplot <- DotPlot(chimeric_and_mock, features = 'Reln', group.by = 'cellType_treatment') 
dplot_dat_Reln <- dplot$data
dplot_Reln_meta <- str_split_fixed(dplot_dat_Reln$id, "_", 2)
colnames(dplot_Reln_meta) = c('cellType', 'treatment')
dplot_dat_Reln <- cbind(dplot_dat_Reln, dplot_Reln_meta)
ggplot(dplot_dat_Reln, aes(x = treatment, y = cellType, color = avg.exp.scaled, size = pct.exp))+
  geom_point()+
  ggtitle('Reln expression single-cell')
VlnPlot(chimeric_and_mock, 'Reln', group.by = 'cellType_treatment') + theme(legend.position = 'None')

chi_mock_pseudo <- createPseudoBulk(chimeric_and_mock, c('Genotype', 'Timepoint', 'Treatment'))
chi_mock_pseudo <- DESeq(chi_mock_pseudo)
resultsNames(chi_mock_pseudo)
treatment_res <- results(chi_mock_pseudo, name = 'Treatment_rChLGTV_vs_PBS')
treatment_res['Reln',]

#LRP8
FeaturePlot(ParseSeuratObj_int, reduction = 'umap.integrated', features = 'Lrp8')
FeaturePlot(ParseSeuratObj_int, reduction = 'umap.integrated', features = 'Rsad2')

VlnPlot(ParseSeuratObj_int, features = 'Rsad2', group.by = 'hasVirus', pt.size = 0)
DotPlot(subset(ParseSeuratObj_int, Treatment == 'rChLGTV' | Treatment == 'PBS'), features = 'Lrp8', 
        group.by = 'manualAnnotation',split.by = 'Treatment', scale = FALSE)

chimeric_and_mock <- subset(ParseSeuratObj_int, Treatment == 'rChLGTV' | Treatment == 'PBS')
chimeric_and_mock$cellType_treatment <- paste(chimeric_and_mock$manualAnnotation, chimeric_and_mock$Treatment, sep = '_')
chimeric_and_mock_ips <- subset(chimeric_and_mock, Genotype == 'IPS1')
chimeric_and_mock_wt <- subset(chimeric_and_mock, Genotype == 'WT')
dplot <- DotPlot(chimeric_and_mock_wt, features = 'rna_Lrp8', 
                 group.by = 'cellType_treatment', scale = TRUE)

dplot_dat <- dplot$data
dplot_meta <- str_split_fixed(dplot_dat$id, "_", 2)
colnames(dplot_meta) = c('cellType', 'treatment')
dplot_dat <- cbind(dplot_dat, dplot_meta)
ggplot(dplot_dat, aes(x = treatment, y = cellType, color = avg.exp.scaled, size = pct.exp))+
  geom_point()+
  ggtitle('Lrp8 expression WT')+
  scale_color_gradient2(low = '#3DB9FF', mid = 'white', high = 'red', midpoint = 0)+
  scale_size(range = c(2, 10))+
  theme_bw()+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 


FeaturePlot(ParseSeuratObj_int, features = 'Lrp8', reduction = 'umap.integrated')