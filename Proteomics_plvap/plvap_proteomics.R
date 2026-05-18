library(readr)
library(stringr)
library(QFeatures)
library(ggplot2)
library(dplyr)
library(ComplexHeatmap)
library(limma)

third_vent = TRUE
fourth_vent = FALSE

if(fourth_vent){
  #read in proteomics data - reads in 4th ventricle samples
  prot_dat <- read_csv("Proteomics_plvap/proteomics_csv.csv")
}

if(third_vent){
  #Read in 3rd ventricle instead 
  prot_dat <- read_csv("Proteomics_plvap/S21-S40_Proteins.csv")
  #Filter out outlier from pca later
  prot_dat <- prot_dat[names(prot_dat) != 'Abundances (Normalized): S40, LGTV, D5, F']
}

grep('Sym', colnames(prot_dat))
colnames(prot_dat)[13] = 'symbol'

#Quick look at plvap
#pdf('~/Documents/ÖverbyLab/for_anna_plots/Plvap_proteomics_time.pdf', width = 6, height = 5)
dplyr::filter(prot_dat, symbol == 'Plvap') %>%
  dplyr::select(starts_with('Abundances')) %>% 
  tidyr::pivot_longer(everything(),names_to = 'sample', values_to = 'exp_levels') %>% 
  dplyr::arrange(sample) %>% 
  dplyr::mutate(treatment = factor(stringr::str_extract(string = sample, pattern = '(PBS)|(LGTV)|(ChLGTV)'), levels = c('PBS', 'LGTV'))) %>% 
  dplyr::mutate(timepoint = stringr::str_extract(string = sample, pattern = 'D.')) %>% 
  ggplot(aes(x = treatment, y = exp_levels, fill = timepoint))+
  geom_boxplot()+
  geom_point(aes(group = timepoint), position = position_dodge(width = .75))+
  theme_classic()+
  ggtitle('Plvap')+
  ylim(c(0, 170000000))
#dev.off()

plvap_exp <- dplyr::filter(prot_dat, symbol == 'Plvap') %>%
  dplyr::select(starts_with('Abundances')) %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'sample') 
colnames(plvap_exp) = c('sample', 'plvap abundances normalized')
plvap_exp$sample = stringr::str_sub(plvap_exp$sample, start = 25, end = 38)
plvap_exp %>% tidyr::separate(col = sample, into = c('sex', 'timepoint', 'treatment'), sep = ',') %>% 
  dplyr::mutate(timepoint = ifelse(grepl('D', timepoint), yes = timepoint, no = NA)) %>% 
  write.csv(file = '~/Documents/ÖverbyLab/Proteomics_plvap/plvap_exp.csv', quote = FALSE,
          row.names = FALSE)
  
dplyr::filter(prot_dat, symbol == 'Stat1') %>%
  dplyr::select(starts_with('Abundances')) %>% 
  tidyr::pivot_longer(everything(),names_to = 'sample', values_to = 'exp_levels') %>% 
  dplyr::arrange(sample) %>% 
  dplyr::mutate(treatment = factor(stringr::str_extract(string = sample, pattern = '(PBS)|(LGTV)|(ChLGTV)'), levels = c('PBS', 'LGTV'))) %>% 
  ggplot(aes(x = treatment, y = exp_levels, color = treatment))+
  geom_point()+
  geom_boxplot()+
  theme_classic()+
  ggtitle('Stat1')

#Try with q features
quant_cols <-  grep('Normalized', names(prot_dat))

prot_q <- readQFeatures(prot_dat, quantCols = quant_cols, name = 'symbol')

sexes = unlist(lapply(names(prot_dat)[quant_cols], FUN = function(x){
  split_meta <- stringr::str_split(x, pattern = ' ')
  sex <- stringr::str_sub(split_meta[[1]][3], start = 1, end = 1) 
  sex
}))

treatments = unlist(lapply(names(prot_dat)[quant_cols], FUN = function(x){
  treatment = factor(stringr::str_extract(string = x, pattern = '(PBS)|(LGTV)|(ChLGTV)'), levels = c('PBS', 'LGTV'))
  if(grepl('wt', treatment)){
    treatment = stringr::str_sub(treatment, start = 1, end = 6)
  }
  treatment
}))

timepoints = unlist(lapply(names(prot_dat)[quant_cols], FUN = function(x){
 stringr::str_extract(string = x, pattern = 'D.')
}))
timepoints[is.na(timepoints)] = 'Mo'

prot_q$sex = sexes
prot_q$treatment = treatments
colData(prot_q)


#Which variables must we keep
rowDataNames(prot_q)

rowvars <- c("symbol")
prot_q <- selectRowData(prot_q, rowvars)

#Check missing values
prot_q <- zeroIsNA(prot_q, i = seq_along(prot_q))
na_vals <- nNA(prot_q, i = seq_along(prot_q))
table(na_vals$nNArows$nNA == 0)

#Filter out any peptides with missing values, can come back to this and change later
prot_q <- filterNA(prot_q, i = seq_along(prot_q), pNA = 0)

#Count unique features
prot_q <- countUniqueFeatures(prot_q,
                             i = "symbol",
                             colDataName = "prot_counts")
colData(prot_q)


#Check if data needs a log transformation
addAssay(prot_q,
         logTransform(prot_q[[1]]),
         name = "peptides_log")
prot_q <- logTransform(prot_q,
                      i = "symbol",
                      name = "prot_log")

par(mfrow = c(1, 2))
limma::plotDensities(assay(prot_q[[1]]), legend = FALSE)
limma::plotDensities(assay(prot_q[[2]]), legend = FALSE)


head(assay(prot_q[['prot_log']]))
(rowData(prot_q[["prot_log"]])$symbol)

#Check outliers
pc <- prcomp(assay(prot_q[['prot_log']]),
             center = TRUE,
             scale. = TRUE)
samp_names <- stringr::str_sub(rownames(pc$rotation), start = 26, end = -1)
meta_data_clean <- data.frame(day = stringr::str_extract(samp_names, pattern = 'D.'),
                              treatment = stringr::str_extract(samp_names, pattern = '(PBS)|(LGTV)'),
                              sex = stringr::str_extract(samp_names, pattern = 'M|F'))
meta_data_clean$day[is.na(meta_data_clean$day)] = 'Mock'

pc$rotation %>% as.data.frame() %>% 
  cbind(meta_data_clean) %>% 
  ggplot(aes(x = PC1, y = PC2, color = treatment, shape = sex))+
  geom_point(size = 3)+
  theme_classic()

#VIsualrownames_to_column()#VIsualize proteins of interest
plot(prot_q)

assay(prot_q, "prot_log")
Heatmap(matrix = assay(prot_q, "prot_log"),
        show_row_names = FALSE)

#try to use limma for degs
norm_prot <- prot_q[['prot_log']]
abundance_matrix <- assay(norm_prot)
sample_meta <- colData(prot_q)
protein_names <- rowData(prot_q[['prot_log']])

rownames(abundance_matrix) <- protein_names$symbol

treatments <- factor(treatments)
timepoints <- factor(timepoints)

treatment_comp = TRUE
timepoint_comp = FALSE

if(timepoint_comp){
  design <- model.matrix(~ 0 + timepoints)
  colnames(design) <- levels(timepoints)
  
  contrast_matrix <- makeContrasts(
    d2_vs_mo =  D2 - Mo,    # adjust group names to match your levels
    d3_vs_mo =  D3 - Mo, 
    d4_vs_mo =  D4 - Mo, 
    d5_vs_mo =  D5 - Mo, 
    d4_vs_d2 =  D4 - D2, 
    levels = design
  )
}

if(treatment_comp){
  design <- model.matrix(~ 0 + treatments)
  colnames(design) <- levels(treatments)
  
  contrast_matrix <- makeContrasts(
    lgtv_vs_pbs =  LGTV - PBS,
    levels = design)
}


# Fit
fit <- lmFit(abundance_matrix, design)
fit <- contrasts.fit(fit, contrast_matrix)
fit <- eBayes(fit, trend = TRUE)

# Quick results without DEqMS
results_limma <- topTable(fit, coef = "lgtv_vs_pbs", 
                          number = Inf, 
                          adjust.method = "BH")

results_limma = dplyr::mutate(results_limma, 
                              sig = ifelse(adj.P.Val < 0.05 & logFC >1,  yes = 'sig', no ='not_sig'))%>% 
  mutate(my_label = ifelse(sig == 'sig', ID, ''))
ggplot(results_limma, aes(x = logFC, y = -log10(adj.P.Val), color = sig))+
  geom_point()+
  scale_color_manual(values = c('black', 'red'))+
  geom_text_repel(aes(label = my_label))+
  theme_classic()+
  theme(legend.position = 'none')

dplyr::filter(results_limma, ID == 'Plvap')
