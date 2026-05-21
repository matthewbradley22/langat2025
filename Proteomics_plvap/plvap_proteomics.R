library(readr)
library(stringr)
library(QFeatures)
library(ggplot2)
library(dplyr)
library(ComplexHeatmap)
library(limma)

# Load data --------------------------------------------------------------------

#Check true or false depending on which data you want to be read in 
third_vent = T
fourth_vent = F

#Read in 3rd ventricle data if third_vent is true (change to whatever path data is in on your computer)
if(third_vent){
  prot_dat <- read_csv("Proteomics_plvap/S21-S40_Proteins.csv")
}

#read in 4th ventricle data if fourth_vent is true (change to whatever path data is in on your computer)
if(fourth_vent){
  prot_dat <- read_csv("Proteomics_plvap/proteomics_csv.csv") 
}

#Change column name "Gene Symbol" to just "symbol" for easier access
gene_symbol_col <- grep('Sym', colnames(prot_dat))
colnames(prot_dat)[gene_symbol_col] = 'symbol'

# Construct QFeatures object for analyis ---------------------------------------

#Using normalized abundance columns for analysis
quant_cols <-  grep('Normalized', names(prot_dat))
prot_dat[quant_cols]

#Read data in as a QFeatures object for analysis
prot_q <- readQFeatures(prot_dat, quantCols = quant_cols, name = 'symbol')

#The metadata in written in the abundances column names (names(prot_dat)[quant_cols])
#Pull out relevant metadata for labelling
sexes = unlist(lapply(names(prot_dat)[quant_cols], FUN = function(x){
  factor(stringr::str_extract(string = x, pattern = '(M$)|(F$)'))
}))

treatments = unlist(lapply(names(prot_dat)[quant_cols], FUN = function(x){
  treatment = factor(stringr::str_extract(string = x, pattern = '(PBS)|(LGTV)'), levels = c('PBS', 'LGTV'))
  if(grepl('wt', treatment)){
    treatment = stringr::str_sub(treatment, start = 1, end = 6)
  }
  treatment
}))

timepoints = unlist(lapply(names(prot_dat)[quant_cols], FUN = function(x){
  stringr::str_extract(string = x, pattern = 'D.')
}))

#Will label mock timepoint as 'Mo', but could make it anything
timepoints[is.na(timepoints)] = 'Mo'

#Add metadata to QFeature object
prot_q$sex = sexes
prot_q$treatment = treatments
prot_q$timepoints = timepoints

#Can check that metadata matches samples
colData(prot_q) %>% as.data.frame()

#Using the 'symbol' column with protein names as our row data in qfeature object
rowvars <- c("symbol")
prot_q <- selectRowData(prot_q, rowvars)
rowData(prot_q)[["symbol"]]

#Quality control steps --------------------------------------------------------
#Not filtering out contaminants or reverse hits, do we need to?
#Following this guide somewhat https://www.bioconductor.org/packages//release/bioc/vignettes/QFeatures/inst/doc/Processing.html

#Check missing values in data
prot_q <- zeroIsNA(prot_q, i = seq_along(prot_q))
na_vals <- nNA(prot_q, i = seq_along(prot_q))

#If I understand right, about 22% of reads have missing values for third ventricle
#26.5% 4th ventricle
na_vals$nNA

#Most peptides have 0 missing value, some have 18/20 ie all samples are missing?
table(na_vals$nNArows$nNA)

#Filter out any peptides that have over 30% values missing. 
#A somewhat high threshold because Bst2 seems to missing values and I do not want to filter it out for now at least
prot_q <- filterNA(prot_q, i = seq_along(prot_q), pNA = 0.3)

#Count unique features across samples to see if it is similar
prot_q <- countUniqueFeatures(prot_q,
                             i = "symbol",
                             colDataName = "prot_counts")
hist(colData(prot_q)$prot_counts, main = 'unique protein count')

#Create a new assay of log transformed protein counts so that data is closer to normal
prot_q <- logTransform(prot_q,
                      i = "symbol",
                      name = "prot_log")

#Plot raw data and log transformed data to see change
par(mfrow = c(1, 2))
limma::plotDensities(assay(prot_q[[1]]), legend = FALSE)
limma::plotDensities(assay(prot_q[[2]]), legend = FALSE)

#Check fpr outliers. pca cannot handle missing values so creating a new object with no missing values
prot_q_for_pca <- filterNA(prot_q, i = seq_along(prot_q), pNA = 0)

#Run PCA
pc <- prcomp(assay(prot_q_for_pca[['prot_log']]),
             center = TRUE,
             scale. = TRUE)

#Create metadata for pca
samp_names <- stringr::str_sub(rownames(pc$rotation), start = 26, end = -1)

#Third ventricle appears to have s40 as a big outlier
if(third_vent){
  meta_data_clean <- data.frame(day = stringr::str_extract(samp_names, pattern = 'D.'),
                                treatment = stringr::str_extract(samp_names, pattern = '(PBS)|(LGTV)'),
                                sex = stringr::str_extract(samp_names, pattern = '(M$)|(F$)'),
                                sample = stringr::str_extract(samp_names, pattern = 'S..'))
  meta_data_clean$day[is.na(meta_data_clean$day)] = 'Mock'
  
 
}

#Fourth ventricle appears to have a sample 9 as big outlier
if(fourth_vent){
  meta_data_clean <- data.frame(day = stringr::str_extract(samp_names, pattern = 'D.'),
                                treatment = stringr::str_extract(samp_names, pattern = '(PBS)|(LGTV)'),
                                sex = stringr::str_extract(samp_names, pattern = '(^M)|(^F)'))
  meta_data_clean$day[is.na(meta_data_clean$day)] = 'Mock'
  meta_data_clean$sample = seq(1:nrow(meta_data_clean))
}

#Plot pca
pc$rotation %>% as.data.frame() %>% 
  cbind(meta_data_clean) %>% 
  ggplot(aes(x = PC1, y = PC2, color = treatment, shape = sex))+
  geom_point(size = 3)+
  geom_text(aes(label = sample), color = 'black')+
  theme_classic()

#For third ventricle: remove S40 as it appears to be a big outlier 
if(third_vent){
  prot_q$sample = seq(21, 40)
  prot_q <- prot_q[,!prot_q$sample == 40]
}

if(fourth_vent){
  prot_q$sample = seq(1, 18)
  prot_q <- prot_q[,!prot_q$sample == 9]
}

#rerun pca with outliers removed
#5 appears a possible outlier
if(fourth_vent){
  prot_q_for_pca <- filterNA(prot_q, i = seq_along(prot_q), pNA = 0)
  pc <- prcomp(assay(prot_q_for_pca[['prot_log']]),
               center = TRUE,
               scale. = TRUE)
  samp_names <- stringr::str_sub(rownames(pc$rotation), start = 26, end = -1)
  meta_data_clean <- data.frame(day = stringr::str_extract(samp_names, pattern = 'D.'),
                                treatment = stringr::str_extract(samp_names, pattern = '(PBS)|(LGTV)'),
                                sex = stringr::str_extract(samp_names, pattern = '(^M)|(^F)'))
  meta_data_clean$day[is.na(meta_data_clean$day)] = 'Mock'
  meta_data_clean$sample = seq(1:nrow(meta_data_clean))
  
}

pc$rotation %>% as.data.frame() %>% 
  cbind(meta_data_clean) %>% 
  ggplot(aes(x = PC1, y = PC2, color = treatment, shape = sex))+
  geom_point(size = 3)+
  geom_text(aes(label = sample), color = 'black')+
  theme_classic()

#Sample to sample heatmap with normalized values
cormat <- round(cor(assay(prot_q_for_pca, "prot_log")), 2) 
cormat %>% as.data.frame() %>% rownames_to_column(var = 'sample_1') %>% 
  tidyr::pivot_longer(cols = starts_with('Abundances'), names_to = 'sample_2', values_to = 'cor') %>% 
  dplyr::mutate(sample_1 = stringr::str_replace(sample_1, pattern = 'Abundances \\(.+\\)\\:', replacement = ''),
                sample_2 = stringr::str_replace(sample_2, pattern = 'Abundances \\(.+\\)\\:', replacement = '')) %>% 
  ggplot(aes(x = sample_1, y = sample_2, fill = cor))+
  geom_tile()+
  theme(axis.text.x = element_text(angle = 90))

#Use limma model to test for degs ----------------------------------------------

#Prepare matrix with log2 normalized values so highly expressed proteins don't dominate model
norm_prot <- prot_q[['prot_log']]

#Should we be imputing missing values? Only about 4% missing after filtering for third vent, 2% for fourth vent
#Imputation increases power for finding degs
nNA(prot_q, i = seq_along(prot_q))$nNAcols$pNA %>% mean()
norm_prot_imp <- impute(norm_prot, method = "min")

abundance_matrix <- assay(norm_prot_imp) #Can use non imputed data here too
sample_meta <- colData(prot_q)
protein_names <- rowData(prot_q[['prot_log']])
rownames(abundance_matrix) <- protein_names$symbol

#We know data is normalized, but also check that distributions are similar across samples
#Looks good
abundance_matrix %>% as.data.frame() %>% 
  rownames_to_column(var = 'protein_symbol') %>% 
  dplyr::filter(protein_symbol != 'NA.') %>% 
  tidyr::pivot_longer(cols = starts_with('Abundances'), names_to = 'sample', values_to = 'value') %>% 
  ggplot(aes(x = sample, y = value))+
  geom_boxplot()+
  theme(axis.text.x = element_blank())

#check how many missing values are left in data
mean(is.na(abundance_matrix))*100 #Only 4% for third ventricle data

#Make variables of interest factors for model
treatments <- factor(prot_q$treatment)
timepoints <- factor(prot_q$timepoints)

#Can either compare treatments or timepoints
treatment_comp = FALSE
timepoint_comp = TRUE

#Treatment comparison does not account for time
if(treatment_comp){
  design <- model.matrix(~ 0 + treatments)
  colnames(design) <- levels(treatments)
  
  contrast_matrix <- makeContrasts(
    lgtv_vs_pbs =  LGTV - PBS,
    levels = design)
}

#Time comparison can compare infected timepoints vs eachother or vs mock
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

#Fit a linear model for each protein
fit <- lmFit(abundance_matrix, design)
#estimate coefficients and standard errors for the linear model
fit <- contrasts.fit(fit, contrast_matrix)
#moderation of the standard errors towards a common value
fit <- eBayes(fit, trend = TRUE)

#Now get results, can change the "coef" to any comparison from the contrast_matrix
#Doing default p value correction and choosing number = inf to return all genes
results_limma <- topTable(fit, coef = "d2_vs_mo", 
                          number = Inf,
                          adjust.method = 'BH')

#Prep data for plotting
results_limma = dplyr::mutate(results_limma, 
                              sig = ifelse(adj.P.Val < 0.05 & abs(logFC) >1,  yes = 'sig', no ='not_sig'))%>% 
  mutate(my_label = ifelse(sig == 'sig', ID, ''))

ggplot(results_limma, aes(x = logFC, y = -log10(adj.P.Val), color = sig))+
  geom_point()+
  scale_color_manual(values = c('black', 'red'))+
  geom_text_repel(aes(label = my_label))+
  theme_classic()+
  theme(legend.position = 'none')

dplyr::filter(results_limma, ID == 'Plvap')

results_limma %>% dplyr::arrange(desc(logFC))


#Quick look at plvap
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
  ggtitle('Plvap')


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

dplyr::filter(prot_dat, symbol == 'Npc1') %>%
  dplyr::select(starts_with('Abundances')) %>% 
  tidyr::pivot_longer(everything(),names_to = 'sample', values_to = 'exp_levels') %>% 
  dplyr::arrange(sample) %>% 
  dplyr::mutate(treatment = factor(stringr::str_extract(string = sample, pattern = '(PBS)|(LGTV)|(ChLGTV)'), levels = c('PBS', 'LGTV'))) %>% 
  ggplot(aes(x = treatment, y = exp_levels, color = treatment))+
  geom_point()+
  geom_boxplot()+
  theme_classic()+
  ggtitle('Npc1')

