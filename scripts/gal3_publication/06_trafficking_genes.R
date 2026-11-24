library(Seurat)
source('./langatFunctions.R')

#Load data
ParseSeuratObj_int <- LoadSeuratRds("~/Documents/ÖverbyLab/data/FilteredRpcaIntegratedDatNoDoublets.rds") 

#Add new column for future plotting
ParseSeuratObj_int$treatment_time <- paste(ParseSeuratObj_int$Treatment, ParseSeuratObj_int$Timepoint, sep = '_')

#Combine PBS timepoints for plotting
ParseSeuratObj_int$treatment_time = ifelse(ParseSeuratObj_int$Treatment == 'PBS', 
                                          yes = 'PBS_3+5',
                                          no = ParseSeuratObj_int$treatment_time)

# - - - - - - - - - - - - - - - 
#### Endothelial signalling ####
# - - - - - - - - - - - - - - - 


#Genes involved in endothelial signalling
wt_endothelial <-  subset(ParseSeuratObj_int, Treatment %in% c('PBS', 'rLGTV') & Organ == 'Cerebrum' & 
                         Genotype == 'WT' & manualAnnotation == 'Endothelial')

endo_trafficking_dot_dat <- DotPlot(wt_endothelial, 
                                    features = c('Vcam1', 'Icam1', 'Icam2', 'Sele', 'Selp', 
                                                 'Mcam', 'F11r', 'Jam2', 'Pecam1', 'Pvr', 'Cd99l2',
                                                 'Cdh5', 'Cxcl12'), group.by = 'treatment_time', 
                                    scale = FALSE)$data

#Split id column into treatment and timepoint for plotting
endo_trafficking_dot_dat <- endo_trafficking_dot_dat %>% 
  tidyr::separate(id, into = c('treatment', 'time'), sep = '_')


ggplot(endo_trafficking_dot_dat, aes(x = time, y = features.plot))+
  geom_point(aes(size = pct.exp, fill = avg.exp.scaled), pch = 21)+
  facet_grid(cols = vars(treatment), scales = "free")+
  #coord_flip()+
  scale_size_continuous(range = c(0.5,6), limits = c(0,100))+
  scale_fill_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white"), 
                       values = c(1.0,0.7,0.4,0),
                       limits = c(0,2.5))+
  ggtitle('Endothelial trafficking genes')+
  theme_bw()+
  theme(panel.grid = element_blank(),
        #panel.border = element_rect(colour = "white", fill = NA),
        panel.spacing = unit(0, "line"),
        strip.background = element_blank())+
  scale_size(range = c(0,9), limits = c(0, 100))+
  ylab('')+
  xlab('')


# - - - - - - - - - - - - - - - 
#### Monocyte migration ####
# - - - - - - - - - - - - - - - 

#Genes involved in monocyte migration
wt_mono <-  subset(ParseSeuratObj_int, Treatment %in% c('PBS', 'rLGTV') & Organ == 'Cerebrum' & 
                            Genotype == 'WT' & manualAnnotation == 'Macrophage/Monocytes')

mono_trafficking_dot_dat <- DotPlot(wt_mono, features = c('Pecam1', 'Pvr', 'Cd99l2', 'Epha1', 'Ephb1', 
                                                              'Itgal', 'Itga4', 'Itgb1', 'Ccr1', 'Ccr2', 'Ccr5',
                                                              'Cxcr4', 'Ackr3'), 
                                    group.by = 'treatment_time', scale = FALSE)$data

#Split id column into treatment and timepoint for plotting
mono_trafficking_dot_dat <- mono_trafficking_dot_dat %>% 
  tidyr::separate(id, into = c('treatment', 'time'), sep = '_')



ggplot(mono_trafficking_dot_dat, aes(x = time, y = features.plot))+
  geom_point(aes(size = pct.exp, fill = avg.exp.scaled), pch = 21)+
  facet_grid(cols = vars(treatment), scales = "free")+
  #coord_flip()+
  scale_size_continuous(range = c(0.5,6), limits = c(0,100))+
  scale_fill_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white"), 
                       values = c(1.0,0.7,0.4,0),
                       limits = c(0,2.5))+
  ggtitle('Monocyte trafficking genes')+
  theme_bw()+
  theme(panel.grid = element_blank(),
        #panel.border = element_rect(colour = "white", fill = NA),
        panel.spacing = unit(0, "line"),
        strip.background = element_blank())+
  scale_size(range = c(0,9), limits = c(0, 100))+
  ylab('')+
  xlab('')
