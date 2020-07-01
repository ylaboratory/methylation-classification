# Take in accession and database from the command line
args <- commandArgs(trailingOnly = TRUE)
type_arg = args[1]
if (type_arg == 'GEO') {
  GEO_accession <- args[2]
} else if (type_arg == 'ENCODE') {
  Encode_accesion <- args[2]
}

# This script is an example of preprocessing the geo and ENCODEmicroarray data given an accession number
# r = getOption("repos")
# r["CRAN"] = "http://cran.rstudio.com"
# options(repos = r)


if (!require('RPMM')) {
  install.packages("RPMM")
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (!require('minfi')) {
  install.packages("minfi")
}
if (!require('GEOquery')) {
  install.packages("GEOquery")
}
if (!require('wateRmelon')) {
  install.packages("wateRmelon")
}
if (!require("IlluminaHumanMethylationEPICmanifest")) {
  install.packages("IlluminaHumanMethylationEPICmanifest")
}
if (!require("IlluminaHumanMethylation450kmanifest")) {
  install.packages("IlluminaHumanMethylation450kmanifest")
}
if (!require("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")) {
  install.packages("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
}
if (!require("IlluminaHumanMethylation450kanno.ilmn12.hg19")) {
  install.packages("IlluminaHumanMethylation450kanno.ilmn12.hg19")
}
if (!require("ENCODExplorer")) {
  install.packages("ENCODExplorer")
}
if (!require("rtracklayer")) {
  install.packages("rtracklayer")
}
if (!require('data.table')) {
  install.packages('data.table')
}
source('./src/download-data-microarray.R')
source('./src/process-microarray.R')
source('./src/liftover.R')
source('./src/output-data-microarray.R')


# The example dataset in GEO used is GSE146179

if (type_arg == 'GEO') {
  download_data_geo_microarray(GEO_accession, ignore_exist = T)
  download_geo_metadata(GEO_accession)
  corrected_data <-
    background_correction(paste0('./raw/GEO/', GEO_accession))
  normalized_data <-
    normalization(corrected_data,
                  paste0('./data/GEO/', GEO_accession, '_sample_metadata.txt'))
  lifted_data <-
    liftover('./annotation/hg19ToHg38.over.chain', normalized_data)
  write_to_file(lifted_data,
                paste0('./data/GEO/', GEO_accession, '_sample_metadata.txt'))
}



#For ENCODE data we use ENCSR420WUN as an example

if (type_arg == 'ENCODE') {
  download_encode(Encode_accesion)
  convert2target(Encode_accesion)
  get_metadata_encode(Encode_accesion)
  corrected_data <-
    background_correction(paste0('raw/ENCODE/', Encode_accesion))
  normalized_data <-
    normalization(corrected_data,
                  paste0('data/ENCODE/', Encode_accesion, '_metadata.txt'))
  lifted_data <-
    liftover('annotation/hg19ToHg38.over.chain', normalized_data)
  write_to_file(lifted_data,
                paste0('data/ENCODE/', Encode_accesion, '_metadata.txt'))
  
}

# This script merges all sample metadata from the corresponding dataset
all_metadata_name <-
  list.files(path = paste0('data/', type_arg, '/'), pattern = '_sample_metadata.txt')
metadatadf <-
  lapply(paste0('data/', type_arg, '/', all_metadata_name), function(x) {
    read.table(file = x,
               header = T,
               sep = "\t")
  })
metadata_total <- do.call('rbind', lapply(metadatadf, as.data.frame))
write.table(
  metadata_total,
  file = paste0('data/', type_arg, '/all_samplemetadata.txt'),
  quote = F,
  row.names = F,
  sep = "\t"
)
all_betavalue_name <-
  list.files(path = paste0('data/', type_arg, '/'), pattern = '_beta_values.txt')
betavaluedf <-
  lapply(paste0('data/', type_arg, '/', all_betavalue_name), function(x) {
    read.table(file = x,
               header = T,
               sep = "\t")
  })
betavalue_total <-
  do.call('merge', lapply(betavaluedf, as.data.table))
write.table(
  betavalue_total,
  file = paste0('data/', type_arg, '/all_betavalues.txt'),
  quote = F,
  row.names = F,
  sep = "\t"
)
