library(tidyr)
library(Seurat)
library(dplyr)
library(ggplot2)

#Check that viral reads cover whole virus genome
 
coverageList <- list()

#List all files with genome locus coverages
covFiles <- list.files("~/Documents/ÖverbyLab/bamLGTVCoverages/viralLocusDepths/")
for(i in 1:length(covFiles)){
  coverageList[[i]] = read.delim(paste0("~/Documents/ÖverbyLab/bamLGTVCoverages/viralLocusDepths/",
                                    covFiles[i]))
}

for(i in 1:length(coverageList)){
  #Get total count of positions with at least one read mapped
  totalPosWithReads <- sum(coverageList[[2]][3] > 0)
  print(totalPosWithReads)
}

plot(coverageList[[4]][,2],coverageList[[4]][,3],
     xlab="Chromosome position",
     ylab = 'Number of reads',
     main = 'Sublibrary 3')

#Plot # bases with each coverage
baseCovList <- list()
for(j in 1:length(coverageList)){
  curList <- coverageList[[j]]
  numBasesWithCov <- c()
  for(i in 1:max(curList[,3])){
    totalBases <- sum(curList[,3] >= i)
    numBasesWithCov = c(numBasesWithCov, totalBases)
  }
  baseCovList[[j]] <- numBasesWithCov
}

mostReadsSeen <- lapply(baseCovList, length) %>% unlist() %>% max()
baseCovList <- lapply(baseCovList, FUN = function(x){
  dif = mostReadsSeen - length(x) 
  x = c(x, rep(0, dif))
})

depthListDf <- data.frame(do.call(cbind, baseCovList))
depthListDf$readDepth <- seq(1:nrow(depthListDf))
depthListDf <- depthListDf %>% pivot_longer(!readDepth,
                             values_to = 'value',
                             names_to = 'name')
depthListDf$fraction = depthListDf$value / 10943
ggplot(depthListDf, aes(x = readDepth, y = fraction, color = name))+
  geom_smooth(se = FALSE)+
  xlab('Read depth')+
  ylab('Fraction of bases >= depth')+
  ggtitle('LGTV Virus Coverage')+
  guides(color=guide_legend(title="SubLibrary"))

#Does sublib 8 stand out in total features/counts? No.
VlnPlot(ParseSeuratObj_int, features = c('nFeature_RNA', 'nCount_RNA'), 
        group.by = 'subLib', assay = 'RNA', pt.size = 0)

#Infected PBS samples across sublibs
infectedPBS <- subset(ParseSeuratObj_int, Treatment == 'PBS' & virusCountPAdj > 0)
infectedPBSHigh <- subset(infectedPBS, virusCountPAdj > 10)
table(infectedPBSHigh$subLib)

###### Create plots of viral counts using fully adjusted values ###### 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

#Read in fully corrected viral levels
adjusted_viral_reads <- list.files('~/Documents/ÖverbyLab/viralReadCountsFullyAdjusted/', full.names = TRUE)
adjusted_viral_reads_df <- lapply(adjusted_viral_reads, FUN = function(x){
  dat <- readr::read_table(x, col_names = FALSE)
  colnames(dat) = c('virusCountFinalAdj', 'barcode')
  dat$barcode = substring(dat$barcode, 6)
  dat
})

for(i in 1:length(adjusted_viral_reads_df)){
  adjusted_viral_reads_df[[i]]$sublib = paste0('seuObj', i)
}

adjusted_viral_reads_binded <- do.call(rbind, adjusted_viral_reads_df)

#Load data
ParseSeuratObj_int <- LoadSeuratRds("~/Documents/ÖverbyLab/data/FilteredRpcaIntegratedDatNoDoublets.rds") 

ParseSeuratObj_int[[]] <- left_join(ParseSeuratObj_int[[]], adjusted_viral_reads_binded, by = c('cell' = 'barcode', 'subLib' = 'sublib'))
ParseSeuratObj_int$virusCountFinalAdj <- ifelse(is.na(ParseSeuratObj_int$virusCountFinalAdj), 0, ParseSeuratObj_int$virusCountFinalAdj)
ParseSeuratObj_int[[]][which(ParseSeuratObj_int$virusCountFinalAdj != ParseSeuratObj_int$virusCountPAdj),][c('virusCount', 'virusCountPAdj', 'virusCountFinalAdj')]

#Normalize viral count by cell
parse_counts <- ParseSeuratObj_int[['RNA']]$counts
viral_read_matrix <- matrix(data = ParseSeuratObj_int$virusCountFinalAdj, ncol = length(ParseSeuratObj_int$virusCountFinalAdj))
rownames(viral_read_matrix) = 'lgtv_final'
parse_counts_with_virus <- rbind(parse_counts, viral_read_matrix)

parse_counts_with_virus_norm <- NormalizeData(parse_counts_with_virus)

norm_vs_unnorm <- data.frame(unnorm = parse_counts_with_virus['lgtv_final',], norm = parse_counts_with_virus_norm['lgtv_final',] )
ParseSeuratObj_int$virus_count_normalized <- parse_counts_with_virus_norm['lgtv_final',]
ParseSeuratObj_int$Genotype = factor(ParseSeuratObj_int$Genotype, levels = c('WT', 'IPS1'))

#Look at viral expression in pbs samples
pdf('~/Documents/ÖverbyLab/single_cell_ISG_figures/fig_1_plots/uncorrected_viral_counts.pdf', width = 8, height = 5)
ParseSeuratObj_int[[]] %>% 
  dplyr::filter(Treatment != 'rLGTV') %>% 
  ggplot(aes(virus_count_normalized, fill = Treatment))+
  geom_histogram(position="identity",, alpha=0.5)+
  theme_classic()+
  facet_wrap(~Genotype)+
  scale_fill_manual(values=c("black", "#56B4E9"))+
  theme(text = element_text(size = 16))+
  ggtitle('uncorrected viral expression')
dev.off()

chlgtv_viral_hist <- hist(subset(ParseSeuratObj_int, Treatment == 'rChLGTV')$virus_count_normalized)
pbs_viral_hist <- hist(subset(ParseSeuratObj_int, Treatment == 'PBS')$virus_count_normalized)

plot(chlgtv_viral_hist, col=rgb(0,0,1,1/4), xlim=c(0,10))  # first histogram
plot(pbs_viral_hist, col=rgb(1,0,0,1/4), xlim=c(0,10), add=T)  # second
quantile(subset(ParseSeuratObj_int, Treatment == 'PBS')$virus_count_normalized, 0.98)
viral_level_to_remove <- sd(subset(ParseSeuratObj_int, Treatment == 'PBS')$virus_count_normalized) * 3
viral_level_to_remove <- 2

#Remove amount of virus that seems to account for most contamination
ParseSeuratObj_int$virusCountFinalAdj_corrected = ParseSeuratObj_int$virus_count_normalized - viral_level_to_remove
ParseSeuratObj_int$virusCountFinalAdj_corrected[ParseSeuratObj_int$virusCountFinalAdj_corrected < 0] = 0

#Prepare celltype order for plotting
ParseSeuratObj_int$manualAnnotation <- factor(ParseSeuratObj_int$manualAnnotation, levels = rev(c('Astrocytes', 'Choroid Plexus', 'Endothelial', 'Ependymal',
                                        'Immature Neurons', 'Microglia', 'Muscle cells', 'Neurons',
                                        'Oligodendrocytes', 'Pericytes', 'B Cells', 'Granulocytes',
                                        'Macrophage/Monocytes', 'Nk cells', 'T cells')))

#Do not keep cell pops with very few cells
cells_with_not_enough <- ParseSeuratObj_int[[]] %>% dplyr::filter(Treatment == 'PBS') %>% 
  dplyr::group_by(Timepoint, Genotype, manualAnnotation) %>% 
  dplyr::summarise(count = n()) %>% 
  dplyr::filter(count < 50)
groups_to_not_keep <- paste(cells_with_not_enough$Timepoint, cells_with_not_enough$Genotype, cells_with_not_enough$manualAnnotation, sep = '_')
#Should filter by cell count so we don't have some groups that are 10 cells total
pdf('~/Documents/ÖverbyLab/single_cell_ISG_figures/fig_1_plots/corrected_viral_counts_mock.pdf', width = 7, height = 5)
ParseSeuratObj_int[[]] %>% dplyr::filter(Treatment == 'PBS') %>% 
  dplyr::filter(manualAnnotation != 'unknown') %>% 
  dplyr::mutate(time_geno_celltype = paste(Timepoint, Genotype, manualAnnotation, sep = '_')) %>% 
  dplyr::filter(!time_geno_celltype %in% groups_to_not_keep) %>% 
  dplyr::mutate(viral_presence <- ifelse(virusCountFinalAdj_corrected > 0, 1, 0)) %>% 
  dplyr::group_by(Genotype, Timepoint, manualAnnotation) %>% 
  dplyr::summarise(pct.exp = sum(virusCountFinalAdj_corrected > 0)*100/length(virusCountFinalAdj_corrected),
                   corrected.viral.counts = mean(virusCountFinalAdj_corrected)) %>% 
  ggplot(aes(x = Timepoint, y = manualAnnotation))+
  facet_wrap(~Genotype)+
  geom_point(aes(size = pct.exp, fill = corrected.viral.counts), pch = 21)+
  scale_fill_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white"), 
                       values = c(1.0,0.7,0.4,0),
                       limits = c(0, 1.6))+
  scale_size_continuous(range = c(1, 2), breaks = c(0, 6))+
  ggtitle('Mock viral count')+
  theme_classic()+
  theme(text = element_text(size = 16))+
  ylab('')+
  xlab('Day')
dev.off()

#Do same for infected samples
cells_with_not_enough_inf <- ParseSeuratObj_int[[]] %>% dplyr::filter(Treatment == 'rChLGTV') %>% 
  dplyr::group_by(Timepoint, Genotype, manualAnnotation) %>% 
  dplyr::summarise(count = n()) %>% 
  dplyr::filter(count < 50)

groups_to_not_keep_inf <- paste(cells_with_not_enough_inf$Timepoint, cells_with_not_enough_inf$Genotype, cells_with_not_enough_inf$manualAnnotation, sep = '_')

pdf('~/Documents/ÖverbyLab/single_cell_ISG_figures/fig_1_plots/corrected_viral_counts.pdf', width = 8, height = 5)
ParseSeuratObj_int[[]] %>% dplyr::filter(Treatment == 'rChLGTV') %>% 
  dplyr::filter(manualAnnotation != 'unknown') %>% 
  dplyr::group_by(Genotype, Timepoint, manualAnnotation) %>%
  dplyr::mutate(time_geno_celltype = paste(Timepoint, Genotype, manualAnnotation, sep = '_')) %>% 
  dplyr::filter(!time_geno_celltype %in% groups_to_not_keep_inf) %>% 
  dplyr::summarise(pct.exp = sum(virusCountFinalAdj_corrected > 0)*100/length(virusCountFinalAdj_corrected),
                   corrected.viral.counts = mean(virusCountFinalAdj_corrected)) %>% 
  ggplot(aes(x = Timepoint, y = manualAnnotation))+
  facet_wrap(~Genotype)+
  geom_point(aes(size = pct.exp, fill = corrected.viral.counts), pch = 21)+
  scale_fill_gradientn(colours = c("#F03C0C","#F57456","#FFB975","white"), 
                       values = c(1.0,0.7,0.4,0),
                       limits = c(0, 1.6))+
  scale_size_continuous(range = c(1,6))+
  ggtitle('rChLGTV viral count')+
  theme_classic()+
  theme(text = element_text(size = 16))+
  ylab('')+
  xlab('Day')+
  scale_x_discrete(labels = c('3', '4' , '5'))
dev.off()

#Check a group as confirmation, day 4 ips astros should have 30-50% infection
infected_astros_ips_day4 <- subset(ParseSeuratObj_int, Treatment == 'rChLGTV' & Timepoint == 'Day 4' & Genotype == 'IPS1' &
                                     manualAnnotation == 'Astrocytes')
table(infected_astros_ips_day4$virusCountFinalAdj_corrected > 0) / nrow(infected_astros_ips_day4[[]])
