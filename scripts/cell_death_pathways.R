#Look at necroptosis/pyroptosis genes from email
source('~/Documents/ÖverbyLab//scripts/langatFunctions.R')

#Load data
ParseSeuratObj_int <- LoadSeuratRds("~/Documents/ÖverbyLab/data/FilteredRpcaIntegratedDatNoDoublets.rds") 
ParseSeuratObj_int$hasVirus = ifelse(ParseSeuratObj_int$virusCountPAdj >= 10, 1, 0)

#Check data
newCols <-  c(brewer.pal(12, 'Paired'), '#99FFE6', '#CE99FF', '#18662E','#737272',  '#FF8AEF')
newCols[11] =  '#FF8AEF'
DimPlot(ParseSeuratObj_int, label = FALSE, group.by = 'manualAnnotation', reduction = 'umap.integrated',
        cols = newCols)
ParseSeuratObj_int$treatment_celltype <- paste(ParseSeuratObj_int$Treatment, ParseSeuratObj_int$manualAnnotation, sep = '_')

#Gene lists from email
anti_apoptosis <- c("Bcl2", "Bcl2l1", 'Bcl2l2', 'Bag1', 'Bag2', 'Bag3', 'Bag4')
pro_apoptosis <- c('Apaf1', 'Casp9', 'Casp8', 'Bax', 'Bak',
                   'Bid', 'Bad', 'Bim', 'Bcl10', 'Bik',
                   'Blk', 'Fas', 'Fasl', 'Tnfrsf1a', 'Tnf', 'Tyro3',
                   'Axl', 'Mertk', 'Tnfsf10', 'Tnfrsf10b', 'Casp3', 'Casp6', 'Casp7')

#Also lists Ferritin, Fth1? Ftl1? Think the map genes are right, should double check though
ferroptosis <- c('Gpx4', 'Acsl4', 'Ptgs2', 'Slc39a14', 'Prnp', 'Steap3', 'Vdac2', 'Vdac3', 'Alox15', 'Atf3')
autophagy <- c('Atg3', 'Atg5', 'Atg7', 'Atg10', 'Atg12', 'Atg13', 'Atg14', 'Ulk1', 'Becn1',
               'Ambra1', 'Map1lc3a', 'Map1lc3a')
cuproptosis <- c('Fdx1', 'Lias', 'Lipt1', 'Dld', 'Dlat', 'Pdha1', 'Pdhb', 'Mtf1',
                 'Gls', 'Cdkn2a', 'Atp7b', 'Slc31a1', 'Atp7a', 'Dlst', 'Dbt', 'Gcsh')

############Analyze gene lists in scRNA data#############
ParseSeuratObj_int <- AddModuleScore(ParseSeuratObj_int, features = list(pro_apoptosis), name = 'pro_apoptosis_score', assay = 'RNA')
ParseSeuratObj_int <- AddModuleScore(ParseSeuratObj_int, features = list(anti_apoptosis), name = 'anti_apoptosis_score', assay = 'RNA')
ParseSeuratObj_int <- AddModuleScore(ParseSeuratObj_int, features = list(ferroptosis), name = 'ferroptosis_score', assay = 'RNA')
ParseSeuratObj_int <- AddModuleScore(ParseSeuratObj_int, features = list(autophagy), name = 'autophagy_score', assay = 'RNA')
ParseSeuratObj_int <- AddModuleScore(ParseSeuratObj_int, features = list(cuproptosis), name = 'cuproptosis_score', assay = 'RNA')

DotPlot(ParseSeuratObj_int, features = 'pro_apoptosis_score1', group.by = 'manualAnnotation', scale = FALSE)
DotPlot(ParseSeuratObj_int, features = pro_apoptosis, group.by = 'manualAnnotation', scale = FALSE)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

#Data subsets
wt_cerebrum <- subset(ParseSeuratObj_int, Genotype == 'WT' & Organ == 'Cerebrum')
wt_cerebellum <- subset(ParseSeuratObj_int, Genotype == 'WT' & Organ == 'Cerebellum')
ips_cerebrum <- subset(ParseSeuratObj_int, Genotype == 'IPS1' & Organ == 'Cerebrum')
ips_cerebellum <- subset(ParseSeuratObj_int, Genotype == 'IPS1' & Organ == 'Cerebellum')

#Ccl12 analysis - expression across celltypes after infection
DotPlot(ParseSeuratObj_int, features = 'Ccl12', group.by = 'manualAnnotation', scale = FALSE)

draw_gene_dotplot <- function(dat, gene){
  dotplot_dat <- DotPlot(dat, features = gene, group.by = 'treatment_celltype', scale = FALSE)$data
  dotplot_dat %>% separate_wider_delim(id, delim = "_", names = c("Treatment", "celltype")) %>% 
    ggplot(aes(x = Treatment, y = celltype, size = pct.exp, color = avg.exp))+
    geom_point()+
    scale_colour_gradient(low = "gray", high = "red")+
    theme(text = element_text(size = 14))
}

draw_gene_dotplot(wt_cerebrum, 'Ccl12') + ggtitle('Ccl12 wt cerebrum')
draw_gene_dotplot(wt_cerebellum, 'Ccl12')+ ggtitle('Ccl12 wt cerebellum')
draw_gene_dotplot(ips_cerebrum, 'Ccl12')+ ggtitle('Ccl12 IPS cerebrum')
draw_gene_dotplot(ips_cerebellum, 'Ccl12')+ ggtitle('Ccl12 IPS cerebrum')

#Gene list dotplots 
#Pro apoptosis
draw_gene_dotplot(wt_cerebrum, 'pro_apoptosis_score1') + ggtitle('pro apoptosis score wt cerebrum')
draw_gene_dotplot(wt_cerebellum, 'pro_apoptosis_score1') + ggtitle('pro apoptosis score wt cerebellum')
draw_gene_dotplot(ips_cerebrum, 'pro_apoptosis_score1') + ggtitle('pro apoptosis score IPS cerebrum')
draw_gene_dotplot(ips_cerebellum, 'pro_apoptosis_score1') + ggtitle('pro apoptosis score IPS cerebellum')

#anti apoptosis
draw_gene_dotplot(wt_cerebrum, 'anti_apoptosis_score1') + ggtitle('anti apoptosis score wt cerebrum')
draw_gene_dotplot(wt_cerebellum, 'anti_apoptosis_score1') + ggtitle('anti apoptosis score wt cerebellum')
draw_gene_dotplot(ips_cerebrum, 'anti_apoptosis_score1') + ggtitle('anti apoptosis score IPS cerebrum')
draw_gene_dotplot(ips_cerebellum, 'anti_apoptosis_score1') + ggtitle('anti apoptosis score IPS cerebellum')

#ferroptosis_score
draw_gene_dotplot(wt_cerebrum, 'ferroptosis_score1') + ggtitle('ferroptosis score wt cerebrum')
draw_gene_dotplot(wt_cerebellum, 'ferroptosis_score1') + ggtitle('ferroptosis apoptosis score wt cerebellum')
draw_gene_dotplot(ips_cerebrum, 'ferroptosis_score1') + ggtitle('ferroptosis apoptosis score IPS cerebrum')
draw_gene_dotplot(ips_cerebellum, 'ferroptosis_score1') + ggtitle('ferroptosis apoptosis score IPS cerebellum')


