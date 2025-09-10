library(tidyr)

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
  ggtitle('LGTV Virus Coverage')

