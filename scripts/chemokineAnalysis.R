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
chimeric_mock_wt <- subset(ParseSeuratObj_int, Treatment != 'rLGTV' & Genotype == 'WT')
chimeric_wt_infected <- subset(chimeric_mock_wt, Treatment == 'rChLGTV')
chimeric_wt_mock <- subset(chimeric_mock_wt, Treatment == 'PBS')
chimeric_ips_infected<- subset(ParseSeuratObj_int, Genotype == 'IPS1' & Treatment == 'rChLGTV')
chimeric_ips_mock<- subset(ParseSeuratObj_int, Genotype == 'IPS1' & Treatment == 'PBS')

#Cytokine dotplots
#WT infected

#Considering subtracting pbs expression from wt expression with this dataframe,
#but not sure what to do with percent expressing?
wt_ccl <- DotPlot(chimeric_mock_wt, ccl_chemokines, group.by = 'time_treatment', scale = FALSE)$data
wt_ccl_meta <- str_split_fixed(wt_ccl$id, "_", 2)
colnames(wt_ccl_meta) <- c('treatment', 'timepoint')
wt_ccl <- cbind(wt_ccl, wt_ccl_meta)
wt_ccl <- rownames_to_column(wt_ccl, var = 'gene')

#Ccl chemokines, wt infected
DotPlot(chimeric_wt_infected, ccl_chemokines, group.by = 'Timepoint', scale = FALSE)+
  coord_flip()+
  scale_color_gradientn(colours = c("#F03C0C","#F57456","#FFB975","lightgrey"),
                       values = c(1.0,0.7,0.4,0),
                       limits = c(0,1.6))+
  scale_size_continuous(range = c(1,6))+
  ggtitle("WT rChLGTV chemokines")

#Why does ccl12 decrease?
ccl12_dat <- DotPlot(chimeric_wt_infected, features = 'Ccl12', group.by = 'time_celltype', scale = FALSE)$data
ccl12_dat_meta <- stringr::str_split_fixed(ccl12_dat$id, "_", 2)
colnames(ccl12_dat_meta) = c('timepoint', 'celltype')
ccl12_dat <- cbind(ccl12_dat, ccl12_dat_meta)
ggplot(ccl12_dat, aes(x = timepoint, y = celltype, color = avg.exp, size = pct.exp))+
  geom_point()+
  scale_color_gradientn(colours = c("#F03C0C","#F57456","#FFB975","lightgrey"),
                        values = c(1.0,0.7,0.4,0),
                        limits = c(0,3))+
  scale_size_continuous(range = c(1,6))+
  theme_bw()+
  ggtitle('WT rChLGTV Ccl12 Expression')

table(chimeric_wt_infected$Timepoint, chimeric_wt_infected$manualAnnotation) %>% 
  as.data.frame() %>% dplyr::group_by(Var1) %>% dplyr::mutate(freq_props = Freq/sum(Freq))%>% 
  ggplot(aes(x = Var1, y = freq_props, fill = Var2))+
  geom_bar(stat = 'identity', position = 'stack', width = 0.6)+
  scale_fill_manual(values = newCols)+
  theme_classic()+
  theme(axis.text.x = element_text(size = 22),
        axis.text.y = element_text(size = 22),
        axis.title=element_text(size=22),
        plot.title = element_text(size = 22))+
  xlab('')+
  ggtitle('WT rChLGTV cell proportions')+
  ylab('Proportion of cells')

#Ccl chemokines, wt mock#Ccl chemokines, wt mockFALSE
DotPlot(chimeric_wt_mock, ccl_chemokines, group.by = 'Timepoint', scale = FALSE)+
  coord_flip()+
  scale_color_gradientn(colours = c("#F03C0C","#F57456","#FFB975","lightgrey"),
                        values = c(1.0,0.7,0.4,0),
                        limits = c(0,1.6))+
  scale_size_continuous(range = c(1,6))+
  ggtitle("WT mock chemokines")

#Cxcl wt infected
DotPlot(chimeric_wt_infected, cxcl_chemokines, group.by = 'Timepoint', scale = FALSE)+
  coord_flip()+
  scale_color_gradientn(colours = c("#F03C0C","#F57456","#FFB975","lightgrey"),
                        values = c(1.0,0.7,0.4,0),
                        limits = c(0,1.6))+
  scale_size_continuous(range = c(1,6))+
  ggtitle("WT rChLGTV chemokines")

#Cxcl wt mock
DotPlot(chimeric_wt_mock, cxcl_chemokines, group.by = 'Timepoint', scale = FALSE)+
  coord_flip()+
  scale_color_gradientn(colours = c("#F03C0C","#F57456","#FFB975","lightgrey"),
                        values = c(1.0,0.7,0.4,0),
                        limits = c(0,1.6))+
  scale_size_continuous(range = c(1,6))+
  ggtitle("WT mock chemokines")

#Ccl chemokines, ips infected
DotPlot(chimeric_ips_infected, ccl_chemokines, group.by = 'Timepoint', scale = FALSE)+
  coord_flip()+
  scale_color_gradientn(colours = c("#F03C0C","#F57456","#FFB975","lightgrey"),
                        values = c(1.0,0.7,0.4,0),
                        limits = c(0,1.6))+
  scale_size_continuous(range = c(1,6))+
  ggtitle("IPS rChLGTV chemokines")

#Ccl chemokines, ips mock
DotPlot(chimeric_ips_mock, ccl_chemokines, group.by = 'Timepoint', scale = FALSE)+
  coord_flip()+
  scale_color_gradientn(colours = c("#F03C0C","#F57456","#FFB975","lightgrey"),
                        values = c(1.0,0.7,0.4,0),
                        limits = c(0,1.6))+
  scale_size_continuous(range = c(1,6))+
  ggtitle("IPS mock chemokines")

#Cxcl chemokines, ips infected
DotPlot(chimeric_ips_infected, cxcl_chemokines, group.by = 'Timepoint', scale = FALSE)+
  coord_flip()+
  scale_color_gradientn(colours = c("#F03C0C","#F57456","#FFB975","lightgrey"),
                        values = c(1.0,0.7,0.4,0),
                        limits = c(0,1.6))+
  scale_size_continuous(range = c(1,6))+
  ggtitle("IPS rChLGTV chemokines")

#Cxcl chemokines, ips mcok
DotPlot(chimeric_ips_mock, cxcl_chemokines, group.by = 'Timepoint', scale = FALSE)+
  coord_flip()+
  scale_color_gradientn(colours = c("#F03C0C","#F57456","#FFB975","lightgrey"),
                        values = c(1.0,0.7,0.4,0),
                        limits = c(0,1.6))+
  scale_size_continuous(range = c(1,6))+
  ggtitle("IPS mock chemokines")


