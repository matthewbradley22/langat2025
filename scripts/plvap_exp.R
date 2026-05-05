#Plvap expression in single cell data. WT vs IPS

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
newCols <-  c( "#7047A1", "#B370AE","#292270",  "#166DF0","#6D92F8",  "#6DC3F8", "#8a0000","#F76363", "#FF96A2", "#D6644B", 
               "#F08C3A", "#fdc087","#074F00", "#208d1f","#7bcd79", 
               "gray")
DimPlot(ParseSeuratObj_int, label = FALSE, group.by = 'manualAnnotation', reduction = 'umap.integrated',
        cols = newCols)

#Look at just langat infection
langat_mock <- subset(ParseSeuratObj_int, Treatment != 'rChLGTV')
table(langat_mock$Treatment)

# Prep to compare between treatments, timepoints, and genotypes
langat_mock$treatment_time_geno <- paste(langat_mock$Treatment, langat_mock$Timepoint, langat_mock$Genotype, sep = '_')

#Look at plvap expression across data
#Split by organ
langat_mock_cerebrum <- subset(langat_mock, Organ == 'Cerebrum')
langat_mock_cerebrum_cp <- subset(langat_mock_cerebrum, manualAnnotation == 'Choroid Plexus')
langat_mock_cerebellum <- subset(langat_mock, Organ == 'Cerebellum')

#Prep data for plotting
cerebrum_plvap <- DotPlot(langat_mock_cerebrum, features = 'Plvap', group.by = 'treatment_time_geno', scale = FALSE)$data

pdf('~/Documents/ÖverbyLab/chp_bulk_rna_amanda/plvap_plots/single_cell_all_cerebrum.pdf', width = 8, height = 5)
cerebrum_plvap %>% tidyr::separate(col = 'id', into = c('treatment', 'timepoint', 'genotype'), sep = '_') %>% 
  ggplot(aes(x = timepoint, y = genotype, fill = avg.exp.scaled, size = pct.exp))+
  facet_wrap(~treatment)+
  geom_point(pch = 21)+
  scale_fill_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white"), 
                       values = c(1.0,0.7,0.4,0))+
  ggtitle('Cerebrum Plvap Expression')+
  theme_classic()+
  theme(text = element_text(size = 20))
dev.off()

cerebrum_langat_markers <- FindMarkers(langat_mock_cerebrum, group.by = 'Treatment', ident.1 = 'rLGTV', test.use = 'MAST')
cerebrum_langat_markers['Plvap',] #returns NA, no sig difference

#Choroid plexus only
cerebrum_plvap_cp <- DotPlot(langat_mock_cerebrum_cp, features = 'Plvap', group.by = 'treatment_time_geno', scale = FALSE)$data

pdf('~/Documents/ÖverbyLab/chp_bulk_rna_amanda/plvap_plots/single_cell_cp_cerebrum.pdf', width = 8, height = 5)
cerebrum_plvap_cp %>% tidyr::separate(col = 'id', into = c('treatment', 'timepoint', 'genotype'), sep = '_') %>% 
  ggplot(aes(x = timepoint, y = genotype, fill = avg.exp.scaled, size = pct.exp))+
  facet_wrap(~treatment)+
  geom_point(pch = 21)+
  scale_fill_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white"), 
                       values = c(1.0,0.7,0.4,0))+
  ggtitle('Cerebrum CP Plvap Expression')+
  theme_classic()+
  theme(text = element_text(size = 20))
dev.off()

cp_langat_markers <- FindMarkers(langat_mock_cerebrum_cp, group.by = 'Treatment', ident.1 = 'rLGTV', test.use = 'MAST')
cp_langat_markers['Plvap',] 
#Do not have infected cerebellum samples with langat infection
cerebellum_plvap <- DotPlot(langat_mock_cerebellum, features = 'Plvap', group.by = 'treatment_time_geno', scale = FALSE)$data
cerebellum_plvap %>% tidyr::separate(col = 'id', into = c('treatment', 'timepoint', 'genotype'), sep = '_') %>% 
  ggplot(aes(x = timepoint, y = genotype, fill = avg.exp.scaled, size = pct.exp))+
  facet_wrap(~treatment)+
  geom_point(pch = 21)+
  scale_fill_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white"), 
                       values = c(1.0,0.7,0.4,0))+
  ggtitle('Cerebellum Plvap Expression')



