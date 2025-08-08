#SoupX attempt, the problem here is that the virus count is not actually in the 
#gene matrix, so would need to add to cells and filtered out droplets

#Load packages
library(SoupX)

#Read in already processed data to provide clustering info
ParseSeuratObj_int <- LoadSeuratRds("./FilteredRpcaIntegratedDat.rds")

#Read in data
parseOutput <- list.files("./FilteredParseOutput/")

ParseMatricies <- lapply(parseOutput, FUN = function(x){
  ReadParseBio(paste0("./FilteredParseOutput/", x))
})

ParseMatricies_unfiltered <- lapply(parseOutput, FUN = function(x){
  ReadParseBio(paste0("./unfilteredParseOutput/", x))
})


#Run soupx on data
#Load in filtered and unfilitered data
soupDat <- SoupChannel(tod = ParseMatricies_unfiltered[[1]], toc = ParseMatricies[[1]])

#Assign clustering metadata
soupDat = setClusters(soupDat, as.character(subset(ParseSeuratObj_int, subLib == 'seuObj1')$seurat_clusters))
sc = autoEstCont(soupDat)



