#Load in libraries
library(Seurat)

#List out our data files, and read them into R
parseOutput <- list.files("./unfilteredParseOutput")

ParseMatricies <- lapply(parseOutput, FUN = function(x){
  ReadParseBio(paste0("unfilteredParseOutput/", x))
})

#Convert data to Seurat objects  
ParseSeuratObj <- lapply(ParseMatricies, FUN = function(x){
  CreateSeuratObject(counts = x, project = "Langat", 
                     min.cells = 3, min.features = 200)
})

#Run through data processing before integration
ParseSeuratObj <- lapply(ParseSeuratObj, FUN = function(x){
  dat <- NormalizeData(x)
  dat <- FindVariableFeatures(dat)
  dat <- ScaleData(dat)
  dat <- RunPCA(dat)
  dat
})

pbmc.big <- merge(pbmc3k, y = c(pbmc4k, pbmc8k), add.cell.ids = c("3K", "4K", "8K"), project = "PBMC15K")

