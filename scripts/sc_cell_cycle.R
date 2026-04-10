library(gprofiler2)
#Load data
ParseSeuratObj_int <- LoadSeuratRds("~/Documents/ÖverbyLab/data/FilteredRpcaIntegratedDatNoDoublets.rds") 

#Check data
newCols <-  c(brewer.pal(12, 'Paired'), '#99FFE6', '#CE99FF', '#18662E','#737272',  '#FF8AEF')
newCols[11] =  '#FF8AEF'

DimPlot(ParseSeuratObj_int, label = FALSE, group.by = 'manualAnnotation', reduction = 'umap.integrated',
        cols = newCols)+
  theme(axis.ticks = element_blank(),
        axis.text=element_blank(),
        legend.text=element_text(size=17))+
  xlab('Umap1')+
  ylab('Umap2')+
  guides(color=guide_legend(override.aes=list(size=8)))+
  ggtitle('')

#Compare infected ips vs infected wt over all celltypes
chimeric_mock <- subset(ParseSeuratObj_int, Treatment != 'rLGTV')

#Load cell cycle genes
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

#Convert to mouse genes
chimeric_mock$time_geno_treatment <- paste(chimeric_mock$Timepoint, chimeric_mock$Genotype, chimeric_mock$Treatment, sep = '_')
mmus_s = gorth(s.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name
mmus_g2m = gorth(g2m.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name

chimeric_mock <- CellCycleScoring(chimeric_mock, s.features = mmus_s, g2m.features = mmus_g2m, set.ident = TRUE)
DimPlot(chimeric_mock, group.by = 'Phase', reduction = 'umap.integrated')

DotPlot(chimeric_mock, features = 'S.Score', group.by = 'time_geno_treatment')
table(chimeric_mock$Genotype, chimeric_mock$Treatment, chimeric_mock$Phase)

chimeric_mock[[]] %>% dplyr::group_by(Genotype, Treatment, Timepoint) %>% 
  dplyr::summarise(s_score = mean(S.Score), g2m_score = mean(G2M.Score))

#Is this effect only because of standout microglia
chimeric_mock_no_microglia <- subset(chimeric_mock, manualAnnotation != 'Microglia')
chimeric_mock_no_microglia[[]] %>% dplyr::group_by(Genotype, Treatment, Timepoint) %>% 
  dplyr::summarise(s_score = mean(S.Score), g2m_score = mean(G2M.Score))

