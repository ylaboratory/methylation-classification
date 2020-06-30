# This script is an example of preprocessing the geo and ENCODEmicroarray data given an accession number
source('./src/download-data-microarray.R')
source('./src/preprocess-microarray.R')
source('./src/liftover.R')
source('./src/output-data-microarray.R')


# The example dataset in GEO used is GSE146179

GEO_accession <- 'GSE146179'
download_data_geo_microarray(GEO_accession, ignore_exist = T)
download_geo_metadata(GEO_accession)
corrected_data <-
  background_correction(paste0('./raw/GEO/', GEO_accession))
normalized_data <- normalization(corrected_data, paste0('./data/GEO/', GEO_accession, '_sample_metadata.txt'))
lifted_data <-
  liftover('./annotation/hg19ToHg38.over.chain', normalized_data)
write_to_file(lifted_data,
              paste0('./data/GEO/', GEO_accession, '_sample_metadata.txt'))


#For ENCODE data we use ENCSR420WUN as an example
 
Encode_accesion<-'ENCSR420WUN'
download_encode(Encode_accesion)
convert2target(Encode_accesion)
get_metadata_encode(Encode_accesion)
corrected_data <-
  background_correction(paste0('./raw/ENCODE/', Encode_accesion))
normalized_data <- normalization(corrected_data,paste0('./data/ENCODE/', Encode_accesion, '_metadata.txt'))
lifted_data <-
  liftover('./annotation/hg19ToHg38.over.chain', normalized_data)
write_to_file(lifted_data,
              paste0('./data/ENCODE/', Encode_accesion, '_metadata.txt'))
