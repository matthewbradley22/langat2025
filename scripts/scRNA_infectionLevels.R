#Load data
ParseSeuratObj_int <- LoadSeuratRds("./data/FilteredRpcaIntegratedDatNoDoublets.rds") 
ParseSeuratObj_int$hasVirus = ifelse(ParseSeuratObj_int$virusCountPAdj >= 10, 1, 0)

#Check data
newCols <-  c(brewer.pal(12, 'Paired'), '#99FFE6', '#CE99FF', '#18662E',  '#FF8AEF','#737272')

pdf('~/Documents/ÖverbyLab/scPlots/scUMAP.pdf', width = 8, height = 6)
DimPlot(ParseSeuratObj_int, label = FALSE, group.by = 'manualAnnotation', reduction = 'umap.integrated',
        cols = newCols)
dev.off()

#Infection levels over time
lgtv_ips <- subset(ParseSeuratObj_int, Treatment == 'rLGTV' & Genotype == 'IPS1')
lgtv_wt <- subset(ParseSeuratObj_int, Treatment == 'rLGTV' & Genotype == 'WT')
chLgtv_ips <- subset(ParseSeuratObj_int, Treatment == 'rChLGTV' & Genotype == 'IPS1')
chLgtv_wt <- subset(ParseSeuratObj_int, Treatment == 'rChLGTV' & Genotype == 'WT')

lgtv_ips[[]] %>% dplyr::group_by(Timepoint, Genotype, hasVirus) %>% 
  dplyr::summarise(total = n()) %>% 
  ggplot(aes(x = Timepoint, y = total, fill = factor(hasVirus)))+
  geom_bar(stat = 'identity', position = 'dodge')+
  ggtitle('LGTV IPS1')+
  geom_text(aes(label=total), vjust=0, position = position_dodge(width = .9))

lgtv_wt[[]] %>% dplyr::group_by(Timepoint, Genotype, hasVirus) %>% 
  dplyr::summarise(total = n()) %>% 
  ggplot(aes(x = Timepoint, y = total, fill = factor(hasVirus)))+
  geom_bar(stat = 'identity', position = 'dodge')+
  ggtitle('LGTV WT')+
  geom_text(aes(label=total), vjust=0, position = position_dodge(width = .9))

chLgtv_ips[[]] %>% dplyr::group_by(Timepoint, Genotype, hasVirus) %>% 
  dplyr::summarise(total = n()) %>% 
  ggplot(aes(x = Timepoint, y = total, fill = factor(hasVirus)))+
  geom_bar(stat = 'identity', position = 'dodge')+
  ggtitle('ChLGTV IPS1')+
  geom_text(aes(label=total), vjust=0, position = position_dodge(width = .9))

chLgtv_wt[[]] %>% dplyr::group_by(Timepoint, Genotype, hasVirus) %>% 
  dplyr::summarise(total = n()) %>% 
  ggplot(aes(x = Timepoint, y = total, fill = factor(hasVirus)))+
  geom_bar(stat = 'identity', position = 'dodge')+
  ggtitle('ChLGTV WT')+
  geom_text(aes(label=total), vjust=0, position = position_dodge(width = .9))

#How many wt infected cells from each cell type
ParseSeuratObj_int[[]] %>% dplyr::group_by(Genotype, manualAnnotation, hasVirus) %>% dplyr::summarise(total = n()) %>% 
  dplyr::filter(Genotype == 'WT') %>% ggplot(aes(x = manualAnnotation, y = total, fill = factor(hasVirus)))+
  geom_bar(stat = 'identity', position = 'dodge')+
  geom_text(aes(label=total), vjust=0, position = position_dodge(1)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  xlab('')+
  labs(fill = 'has virus')

ParseSeuratObj_int[[]] %>% dplyr::group_by(Genotype, manualAnnotation, hasVirus) %>% dplyr::summarise(total = n()) %>% 
  dplyr::filter(Genotype == 'IPS1') %>% ggplot(aes(x = manualAnnotation, y = total, fill = factor(hasVirus)))+
  geom_bar(stat = 'identity', position = 'dodge')+
  geom_text(aes(label=total), vjust=0, position = position_dodge(1)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  xlab('')+
  labs(fill = 'has virus')
