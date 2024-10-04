#last updated: Sept 2024

# Take in accession and database from the command line to output raw files to raw/database/
# The working directory needs to be changed accordingly

# if updates don't work: options(repos="https://CRAN.R-project.org")
# if installs don't work: lib="/home/mk98/R/x86_64-redhat-linux-gnu-library/4.0"


args <- commandArgs(trailingOnly = TRUE)

# manual arguments for now start
args[1]<-'GEO'
args[2]<-'F'
args[3]<-'GEO_microarray_accession_450K.txt'
# args[3] <- 'gdc_manifest.2024-04-19.beta.value.txt'
# manual arguments for now end

# The example dataset in GEO used is GEO_microarray_accession.txt
# The example dataset in ENCODE used is ENCODE_microarray_accession.txt
# The example dataset in TCGA used is TCGA_microarray_manifest_test.txt
# The full dataset in GEO used is GEO_microarray_accession_full.txt
database_type<- args[1]
manifest_file<- args[3]

if (args[2]=="F") {
  ignore_exist_state<-F
} else if (args[2]=="T"){
  ignore_exist_state<-T
}

setwd('/grain/mk98/methyl/methylation-classification') #needs to change accordingly
print(.libPaths())
.libPaths( c( "/home/mk98/R/x86_64-redhat-linux-gnu-library/4.0", .libPaths() ) ) #needs to change accordingly

if (!require('readr')) {
  install.packages('readr')
}
if (!require('RPMM')) {
  install.packages("RPMM", repos = "http://cran.us.r-project.org")
}
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "http://cran.us.r-project.org", version="3.14")
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
# if (!require("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")) {
#   BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
# }
# if (!require("IlluminaHumanMethylation450kanno.ilmn12.hg19")) {
#   BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")
# }
# if (!require("ENCODExplorer")) {
#   BiocManager::install("ENCODExplorer")
# }
if (!require("rtracklayer")) {
  BiocManager::install("rtracklayer")
}
if (!require("TCGAutils")) {
  BiocManager::install("TCGAutils")
}
if (!require("TCGAbiolinks")) {
  BiocManager::install("TCGAbiolinks")
}
if (!require('data.table')) {
  install.packages('data.table', repos = "http://cran.us.r-project.org")
}

source('./src/download-data-microarray.R')
source('./src/process-microarray.R')
source('./src/liftover.R')
# manifest<-read_tsv('./annotation/HM450.hg38.manifest.gencode.v36.tsv')
# manifest$CpG<-manifest$CpG_beg+1

if (database_type == 'GEO') {
  accession_list<-read.table(paste0('annotation/',manifest_file), header = F)
  for (accession in as.character(accession_list$V1)) {
    print(accession)
    tryCatch({
      download_data_geo_microarray(accession, ignore_exist = ignore_exist_state)
      },
      error=function(e) {
        print(paste("BETAVALUE ERROR: ",accession))
        message(e)
        writeLines('\n')
      }
    )
    tryCatch({
      download_geo_metadata(accession, ignore_exist= ignore_exist_state)
      if (file.exists(paste0('./data/GEO/',accession, "_beta_values_probe.txt",sep = "")) == FALSE){
        corrected_data <-
          background_correction(paste0('raw/', database_type, '/', accession))
        normalized_data <-
          normalization(corrected_data,
                        paste0('data/', database_type,'/', accession, '_sample_metadata.txt'))
        # lifted_data <-
        #   liftover('annotation/hg19ToHg38.over.chain', normalized_data)
        write.table(normalized_data, paste0('data/', database_type,'/', accession, '_beta_values_probe.txt'),
                    quote = F,
                    sep = '\t', row.names = F, col.names = T)
        # unlink(paste0('raw/GEO/', accession,"/*"), force=TRUE)
      }
      else {print(paste(accession, " probe file already exists"))}
      },
      error=function(e) {
        print(paste("METADATA ERROR: ",accession))
        message(e)
        writeLines('\n')
      }
    )
  }
} else if (database_type == 'ENCODE') {
  accession_list<-read.table(paste0('annotation/',manifest_file), header = F)
  accession_list=accession_list$V1
  other_accession=list()
  for (i in 1:length(accession_list)) {
    accession=as.character(accession_list[i])
    print(accession)
    download_encode(accession, ignore_exist = ignore_exist_state)
    signal=convert2target(accession)
    if (signal!='No_red' & accession!='ENCSR719GFJ'){
      platform=get_metadata_encode(accession)
      if (platform=='Illumina Infinium Omni5Exome-4 Kit'){
        other_accession=c(other_accession,accession)
        print(accession)
      }
      else{
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
      }
    }
    }
    
  accession_list=accession_list[! accession_list %in% other_accession]
} else if (database_type == 'TCGA'){
  TCGA_manifest<-read.table(paste0('annotation/',manifest_file), header = T, sep = '\t')
  print(TCGA_manifest$id)
  # file2case<-UUIDtoUUID(as.character(TCGA_manifest$id), to_type="case_id")
  # print(file2case)
  file2bar<-UUIDtoBarcode(TCGA_manifest$id, from_type = "file_id", legacy = F)
  # download_TCGA(manifest_file)
  reorganize_TCGA(file2bar)
  get_metadata_TCGA(file2bar)
  accession_list<-unique(file2bar[,'associated_entities.entity_submitter_id'])
  for (accession in as.character(accession_list)){
      print(accession)
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
    }
}
