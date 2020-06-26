# This script is an example of preprocessing the geo and ENCODEmicroarray data given an accession number
source('./src/download-data-geo-microarray.R')
source('./src/get-metadata-geo-microarray.R')
source('./src/background-correction.R')
source('./src/beta-value-normalization.R')
source('./src/liftover.R')
source('./src/write-to-file-microarray.R')
source('./src/download-encode-microarray.R')
source('./src/get-metadata-encode-microarray.R')


# The example dataset in GEO used is GSE146179

accession.num <- 'GSE146179'
download.data.geo.microarray(accession.num, ignore.exisit = T)
download.geo.metadata(accession.num)
corrected.data <-
  background.correction(paste0('./raw/GEO/', accession.num))
normalized.data <- normalization(corrected.data)
lifted.data <-
  liftover('./annotation/hg19ToHg38.over.chain', normalized.data)
write.to.file(lifted.data,
              paste0('./data/GEO/', accession.num, '_sample_metadata.txt'))


#For ENCODE data we use ENCSR420WUN as an example
 
accession.name<-'ENCSR420WUN'
download.encode(accession.name)
convert2target(accession.name)
get.metadata.encode(accession.name)
corrected.data <-
  background.correction(paste0('./raw/ENCODE/', accession.name))
normalized.data <- normalization(corrected.data)
lifted.data <-
  liftover('./annotation/hg19ToHg38.over.chain', normalized.data)
write.to.file(lifted.data,
              paste0('./data/ENCODE/', accession.name, '_metadata.txt'))
