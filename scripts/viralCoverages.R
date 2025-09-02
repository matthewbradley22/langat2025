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

plot(coverageList[[4]][,2],coverageList[[4]][,3])
