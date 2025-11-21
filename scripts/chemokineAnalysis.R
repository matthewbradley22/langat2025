#Analyse changes in chemokin expression over time

#Load data
ParseSeuratObj_int <- LoadSeuratRds("~/Documents/ÖverbyLab/data/FilteredRpcaIntegratedDatNoDoublets.rds") 
ParseSeuratObj_int$hasVirus = ifelse(ParseSeuratObj_int$virusCountPAdj >= 10, 1, 0)

#Create new column for plotting later 
ParseSeuratObj_int$time_treatment <- paste(ParseSeuratObj_int$Treatment, ParseSeuratObj_int$Timepoint, sep = '_')
ParseSeuratObj_int$time_celltype <-  paste(ParseSeuratObj_int$Timepoint, ParseSeuratObj_int$manualAnnotation, sep = '_')

#Check data
newCols <-  c(brewer.pal(12, 'Paired'), '#99FFE6', '#CE99FF', '#18662E','#737272',  '#FF8AEF')
newCols[11] =  '#FF8AEF'
DimPlot(ParseSeuratObj_int, label = FALSE, group.by = 'manualAnnotation', reduction = 'umap.integrated',
        cols = newCols)

#Chemokines
ccl_chemokines <- rownames(ParseSeuratObj_int@assays$RNA$data)[grep('Ccl', rownames(ParseSeuratObj_int@assays$RNA$data))]
cxcl_chemokines <- rownames(ParseSeuratObj_int@assays$RNA$data)[grep('Cxcl', rownames(ParseSeuratObj_int@assays$RNA$data))]

#subset out lgtv
chimeric_mock <- subset(ParseSeuratObj_int, Treatment != 'rLGTV')
chimeric_mock_wt_infected <- subset(chimeric_mock, Genotype == 'WT' & Treatment == 'rChLGTV')
chimeric_mock_wt_mock <- subset(chimeric_mock, Genotype == 'WT' & Treatment == 'PBS')
chimeric_mock_ips_infected<- subset(chimeric_mock, Genotype == 'IPS1' & Treatment == 'rChLGTV')
chimeric_mock_ips_mock<- subset(chimeric_mock, Genotype == 'IPS1' & Treatment == 'PBS')

#Cytokine dotplots
#WT infected
DotPlot(chimeric_mock_wt_infected, ccl_chemokines, group.by = 'Timepoint', scale = FALSE)+
  coord_flip()+
  scale_color_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white","lightblue"),
                       values = c(1.0,0.7,0.4,0.2,0),
                       limits = c(0,1.6))+
  scale_size_continuous(range = c(1,6))+
  ggtitle("WT rChLGTV chemokines")

DotPlot(chimeric_mock_wt_infected, cxcl_chemokines, group.by = 'Timepoint', scale = FALSE)+
  coord_flip()+
  scale_color_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white","lightblue"),
                        values = c(1.0,0.7,0.4,0.2,0),
                        limits = c(0,2.8))+
  scale_size_continuous(range = c(1,6))+
  ggtitle("WT rChLGTV chemokines")

#IPS
DotPlot(chimeric_mock_ips_infected, ccl_chemokines, group.by = 'Timepoint', scale = FALSE)+
  coord_flip()+
  scale_color_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white","lightblue"),
                        values = c(1.0,0.7,0.4,0.2,0),
                        limits = c(0,1.6))+
  scale_size_continuous(range = c(1,6))+
  ggtitle("IPS rChLGTV chemokines")

DotPlot(chimeric_mock_ips_infected, cxcl_chemokines, group.by = 'Timepoint', scale = FALSE)+
  coord_flip()+
  scale_color_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white","lightblue"),
                        values = c(1.0,0.7,0.4,0.2,0),
                        limits = c(0,2.8))+
  scale_size_continuous(range = c(1,6))+
  ggtitle("IPS rChLGTV chemokines")


