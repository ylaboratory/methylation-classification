# This script creates the preprocessig raw object from raw idat files in a directory
# datadir is the directory where raw idat files are stored
# return an RGset raw object
create.preprocessing.raw<-function(datadir){
  file.names <-
    list.files(path =  datadir,
               pattern = 'Grn.idat',
               full.names = T)
  if (length(file.names) == 0) {
    stop(paste0('no idat file found for ', accession.num))
  }
  targets <- data.frame('Basename' = sub('_Grn.idat*', "", file.names))
  RGset <- read.metharray.exp(targets = targets)
  return(RGset)
}
