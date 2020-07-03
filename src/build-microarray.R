# Take in accession and database from the command line
args <- commandArgs(trailingOnly = TRUE)
database_type<- args[1]
accession <- args[2]
if (args[3]=="F") {
  ignore_exist_state<-F
} else if (args[3]=="T"){
  ignore_exist_state<-T
}

# This script is an example of preprocessing the geo and ENCODEmicroarray data given an accession number
if (!require('RPMM')) {
  install.packages("RPMM")
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (!require('minfi')) {
  BiocManager::install("minfi")
}
if (!require('GEOquery')) {
  BiocManager::install("GEOquery")
}
if (!require('wateRmelon')) {
  BiocManager::install("wateRmelon")
}
if (!require("IlluminaHumanMethylationEPICmanifest")) {
  BiocManager::install("IlluminaHumanMethylationEPICmanifest")
}
if (!require("IlluminaHumanMethylation450kmanifest")) {
  BiocManager::install("IlluminaHumanMethylation450kmanifest")
}
if (!require("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")) {
  BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
}
if (!require("IlluminaHumanMethylation450kanno.ilmn12.hg19")) {
  BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")
}
if (!require("ENCODExplorer")) {
  BiocManager::install("ENCODExplorer")
}
if (!require("rtracklayer")) {
  BiocManager::install("rtracklayer")
}
if (!require('data.table')) {
  install.packages('data.table')
}
source('./src/download-data-microarray.R')
source('./src/process-microarray.R')
source('./src/liftover.R')

# The example dataset in GEO used is GSE146179
# For ENCODE data we use ENCSR420WUN as an example
if (database_type == 'GEO') {
  download_data_geo_microarray(accession, ignore_exist = ignore_exist_state)
  download_geo_metadata(accession)
} else if (database_type=='ENCODE'){
  download_encode(accession, ignore_exist = ignore_exist_state)
  convert2target(accession)
  get_metadata_encode(accession)
}
corrected_data <-
    background_correction(paste0('raw/', database_type, '/', accession))
normalized_data <-
    normalization(corrected_data,
                  paste0('data/', database_type,'/', accession, '_sample_metadata.txt'))
lifted_data <-
    liftover('annotation/hg19ToHg38.over.chain', normalized_data)
write.table(lifted_data, paste0('data/', database_type,'/', accession, '_beta_values.txt'),
            quote = F,
            sep = '\t', row.names = F, col.names = T)

# This script merges all sample metadata from the corresponding dataset
all_metadata_name <-
  list.files(path = paste0('data/', database_type), pattern = '_sample_metadata.txt', full.names = T)
metadata_total <-
  rbindlist(lapply(all_metadata_name, function(x) {
    fread(x)
  }), use.names = T, fill = T)
write.table(
  metadata_total,
  file = paste0('data/', database_type, '/all_samplemetadata.txt'),
  quote = F,
  row.names = F,
  sep = "\t"
)
all_betavalue_name <-
  list.files(path = paste0('data/', database_type), pattern = '_beta_values.txt', full.names = T)
betavaluedf <-
  lapply(all_betavalue_name, function(x) {
    fread(x)
  })
betavalue_total<-Reduce(function(x, y) merge(x, y, by=c('chr','loci')), betavaluedf)
write.table(
  betavalue_total,
  file = paste0('data/', database_type, '/all_betavalues.txt'),
  quote = F,
  row.names = F,
  sep = "\t"
)
