# This script creates the preprocessig raw object from raw idat files in a directory and performs background correction
# datadir is the directory where raw idat files are stored
# return an MethylSet object
library(minfi)
library(wateRmelon)
background.correction<-function(datadir){
  file.names <-
    list.files(path =  datadir,
               pattern = 'Grn.idat',
               full.names = T)
  if (length(file.names) == 0) {
    stop(paste0('no idat file found'))
  }
  targets <- data.frame('Basename' = sub('_Grn.idat.*', "", file.names))
  RGset <- read.metharray.exp(targets = targets)
  GRset.noob <-
    preprocessNoob(
      RGset,
      offset = 15,
      dyeCorr = TRUE,
      verbose = FALSE,
      dyeMethod = c("single", "reference")
    )
  return(GRset.noob)
}
