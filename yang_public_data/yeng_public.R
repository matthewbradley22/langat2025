library(dplyr)
library(Seurat)
library(ggplot2)

#Load in public yang data from https://link.springer.com/article/10.1186/s12974-024-03071-1#Sec2
yang_data_files <- list.dirs(path = "~/Documents/ÖverbyLab/yang_public_data/yang_encephalitis_public_data", full.names = TRUE, recursive = FALSE)

#Error is in second file (healthy_2) + warnings in a couple other
#skipping healthy_2 for now due to data not being complete and healthy_1 as it does not seem to have all entries either
yang_data <- lapply(yang_data_files[c(3:8)], FUN = function(x){
  Read10X(data.dir = x)
})

#Add name prefixes so there are no cell merge issues
for(i in 1:length(yang_data)){
  colnames(yang_data[[i]]) = paste(i, colnames(yang_data[[i]]), sep = '_')
}

#Create seurat object
yang_data <- CreateSeuratObject(counts = yang_data, project = "yang_public", min.cells = 3, min.features = 200)
yang_data[["percent.mt"]] <- PercentageFeatureSet(yang_data, pattern = "^mt-")

#In paper they say they kept umis within 2 standard devs, checking what cutoff that would be
feature_dev <- sd(yang_data$nFeature_RNA)
feature_mean <- mean(yang_data$nFeature_RNA)
feature_mean + 2*feature_dev
feature_mean - 2*feature_dev

#Find reasonable cutoffs for initial filtering. Can be liberal here and then remove later any weird clusters
VlnPlot(yang_data, features = c("nFeature_RNA"),  pt.size = 0)+
  geom_hline(yintercept = 5000)

#Will make mito cutoff at 30 since that's what the paper did
VlnPlot(yang_data, features = c("percent.mt"),  pt.size = 0)+
  geom_hline(yintercept = 30)

#Filter data
yang_data <- subset(yang_data, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 30)

#Standard preprocessing
yang_data <- NormalizeData(yang_data)
yang_data <- FindVariableFeatures(yang_data)
yang_data <- ScaleData(yang_data)
yang_data <- RunPCA(yang_data)

#Need to allow greater ram usage to run pca integration
options(future.globals.maxSize = 14500 * 1024^2)

#Integrate data
yang_data <- IntegrateLayers(object = yang_data, method = RPCAIntegration, orig.reduction = "pca", new.reduction = "integrated.rpca",
                        verbose = FALSE)

#Reset max ram to not accidentally use a bunch later without realizing
options(future.globals.maxSize = 500 * 1024^2)

#Prepare integrated data to plot
ElbowPlot(yang_data, ndims = 40)
yang_data <- FindNeighbors(yang_data, reduction = "integrated.rpca", dims = 1:30)
yang_data <- FindClusters(yang_data, resolution = 2, cluster.name = "rpca_clusters")
yang_data <- RunUMAP(yang_data, reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap.rpca")

#Plot cells
DimPlot(yang_data,reduction = "umap.rpca")
SaveSeuratRds(yang_data, file = "~/Documents/ÖverbyLab/yang_public_data/yang_rpca_integrated_obj.RDS")



######## Debugging ############
## ## ## ## ## ## ## ## ## ## ##

#Bug found, last row of matrix data in healthy_2 sample only contains one value
#rather than all 3 columns. Need to figure out how to fix
#Can recreate with: 
file_here = file("/Users/matthewbradley/Documents/ÖverbyLab/yang_public_data/yang_encephalitis_public_data/healthy_2/matrix.mtx.gz")
Matrix::readMM(file_here)

#Have to redefine for it to continue working
file_here = file("/Users/matthewbradley/Documents/ÖverbyLab/yang_public_data/yang_encephalitis_public_data/healthy_2/matrix.mtx.gz")
open(file_here)

scan1 <- function(what, ...){
  scan(file_here, nmax = 1, what = what, quiet = TRUE, ...)
}

if (scan1(character()) != "%%MatrixMarket")# hdr
  stop("file is not a MatrixMarket file")
if (!(typ <- tolower(scan1(character()))) %in% "matrix")
  stop(gettextf("type '%s' not recognized", typ), domain = NA)
if (!(repr <- tolower(scan1(character()))) %in% c("coordinate", "array"))
  stop(gettextf("representation '%s' not recognized", repr), domain = NA)
elt <- tolower(scan1(character()))
if (!elt %in% c("real", "complex", "integer", "pattern"))
  stop(gettextf("element type '%s' not recognized", elt), domain = NA)

sym <- tolower(scan1(character()))

nr <- scan1(integer(), comment.char = "%")
nc <- scan1(integer())
nz <- scan1(integer())


checkIJ <- function(els) {
  if((nz. <- length(els$i)) < nz)
    warning(gettextf("readMM(): expected %d entries but found only %d",
                     nz, nz.), call. = FALSE, domain = NA)
  if(any(is.na(els$i) | els$i < 1L | els$i > nr))
    stop(gettextf("readMM(): row indices 'i' are not in 1:nrow[=%d]",
                  nr), call. = FALSE, domain = NA)
  if(any(is.na(els$j) | els$j < 1L | els$j > nc))
    stop(gettextf("readMM(): column indices 'j' are not in 1:ncol[=%d]",
                  nc), call. = FALSE, domain = NA)
}

if (repr == "coordinate") {
  switch(elt,
         "real" = ,
         "integer" = {
           ## TODO: the "integer" element type should be returned as
           ##       an object of an "iMatrix" subclass--once there are
           els <- scan(file_here, nmax = nz, quiet = TRUE,
                       what= list(i= integer(), j= integer(), x= numeric()))
           checkIJ(els)
           switch(sym,
                  "general" = {
                    new("dgTMatrix", Dim = c(nr, nc), i = els$i - 1L,
                        j = els$j - 1L, x = els$x)
                  },
                  "symmetric" = {
                    new("dsTMatrix", uplo = "L", Dim = c(nr, nc),
                        i = els$i - 1L, j = els$j - 1L, x = els$x)
                  },
                  "skew-symmetric" = {
                    stop("symmetry form 'skew-symmetric' not yet implemented for reading")
                    ## FIXME: use dgT... but must expand the (i,j,x) slots!
                    new("dgTMatrix", uplo = "L", Dim = c(nr, nc),
                        i = els$i - 1L, j = els$j - 1L, x = els$x)
                    
                  },
                  "hermitian" = {
                    stop("symmetry form 'hermitian' not yet implemented for reading")
                  },
                  ## otherwise (not possible; just defensive programming):
                  stop(gettextf("symmetry form '%s' is not yet implemented",
                                sym), domain = NA)
           )
         },
         "pattern" = {
           els <- scan(file, nmax = nz, quiet = TRUE,
                       what = list(i = integer(), j = integer()))
           checkIJ(els)
           switch(sym,
                  "general" = {
                    new("ngTMatrix", Dim = c(nr, nc),
                        i = els$i - 1L, j = els$j - 1L)
                  },
                  "symmetric" = {
                    new("nsTMatrix", uplo = "L", Dim = c(nr, nc),
                        i = els$i - 1L, j = els$j - 1L)
                  },
                  "skew-symmetric" = {
                    stop("symmetry form 'skew-symmetric' not yet implemented for reading")
                    ## FIXME: use dgT... but must expand the (i,j,x) slots!
                    new("ngTMatrix", uplo = "L", Dim = c(nr, nc),
                        i = els$i - 1L, j = els$j - 1L)
                    
                  },
                  "hermitian" = {
                    stop("symmetry form 'hermitian' not yet implemented for reading")
                  },
                  ## otherwise (not possible; just defensive programming):
                  stop(gettextf("symmetry form '%s' is not yet implemented",
                                sym), domain = NA)
           )
         },
         "complex" = {
           stop("element type 'complex' not yet implemented")
         },
         ## otherwise (not possible currently):
         stop(gettextf("'%s()' is not yet implemented for element type '%s'",
                       "readMM", elt), domain = NA))
}

close(file_here)
