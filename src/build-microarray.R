# Take in accession and database from the command line
# The working directory needs to be changed accordingly
args <- commandArgs(trailingOnly = TRUE)

#fornow_start
args[1]<-'GEO'
args[2]<-'F'
args[3]<-'GEO_microarray_accession_full.txt'
#fornow_end

database_type<- args[1]
manifest_file<- args[3]

if (args[2]=="F") {
  ignore_exist_state<-F
} else if (args[2]=="T"){
  ignore_exist_state<-T
}
setwd('/grain/mk98/methyl/methylation-classification')
# This script is an example of preprocessing the geo and ENCODEmicroarray data given an accession number
if (!require('RPMM')) {
  install.packages("RPMM", repos = "http://cran.us.r-project.org")
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "http://cran.us.r-project.org")
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

# The example dataset in GEO used is GEO_microarray_accession.txt
# The example dataset in ENCODE used is ENCODE_microarray_accession.txt
# The example dataset in TCGA used is TCGA_microarray_manifest_test.txt
if (database_type == 'GEO') {
  accession_list<-read.table(paste0('annotation/',manifest_file), header = F)
  for (accession in as.character(accession_list$V1)) {
    print('Here')
    print(accession)
    tryCatch(
      {
        download_data_geo_microarray(accession, ignore_exist = ignore_exist_state)
        download_geo_metadata(accession)
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
        unlink(paste0('raw/GEO/', accession,"/*"), force=TRUE)
      },
      error=function(e) print(paste("ERROR: ",accession)))
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
  print(as.character(TCGA_manifest$id))
  file2bar<-UUIDtoBarcode(as.character(TCGA_manifest$id), from_type = "file_id", legacy = T)
  download_TCGA(manifest_file)
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
database_type='GEO'
# This script merges all sample metadata from the corresponding dataset
all_metadata_name <-
  list.files(path = paste0('data/', database_type), pattern = '_sample_metadata.txt', full.names = T)
metadata_total <-
  rbindlist(lapply(all_metadata_name, function(x) {
    fread(x)
  }), use.names = T, fill = T)
write.table(
  metadata_total,
  file = paste0('data/', database_type, '/all_metadata.txt'),
  quote = F,
  row.names = F,
  sep = "\t"
)
library(R.utils)
library(data.table)
gzip(filename=paste0('data/', database_type, '/all_metadata.txt'), 
     destname=paste0('data/', database_type, '/all_metadata.txt.gz'), overwrite=TRUE, remove=TRUE)

# This script merges all beta_value files together and only keeps the loci that are present in all the files
# The script only merges 450k data(485344 loci) and EpicBeadchip data (865613 loci)
all_betavalue_name <-
  list.files(path = paste0('data/', database_type), pattern = '_beta_values.txt', full.names = T, recursive = F)
betavaluedf <-
  lapply(all_betavalue_name, function(x) {
    fread(x, quote = "")
  })
total_matrix<-list()

# This function removes the duplicated loci in 450k and Beadchip data
for (i in 1:length(betavaluedf)){
  dataset_sample<-betavaluedf[[i]]
  datatype1=betavaluedf[[1]]#Here put in a dataset in the matrix list that is Beadchip data
  loci1=paste(datatype1$chr,datatype1$loci, sep = ' ')
  duplicate1<-which(loci1%in%loci1[duplicated(loci1)])
  datatype2=betavaluedf[[3]]#Here put in a dataset in the matrix list that is 450k data
  loci2=paste(datatype2$chr,datatype2$loci, sep = ' ')
  duplicate2<-which(loci2%in%loci2[duplicated(loci2)])
  if (nrow(dataset_sample)>800000){
    dataset_sample_removed<-dataset_sample[-duplicate1,]
    total_matrix[[i]]<-dataset_sample_removed
  } else if(nrow(dataset_sample)<500000){
    dataset_sample_removed<-dataset_sample[-duplicate2,]
    total_matrix[[i]]<-dataset_sample_removed
  }
}
#Here put in a dataset(like the second in total_matrix) in the deduplicated matrix list that is Beadchip data
loci1_rm=paste(total_matrix[[1]]$chr,total_matrix[[1]]$loci, sep = ' ')
#Here put in a dataset in the deduplicated matrix list that is 450k data
loci2_rm=paste(total_matrix[[3]]$chr,total_matrix[[3]]$loci, sep = ' ')
loci_location1<-which(loci1_rm%in%intersect(loci1_rm,loci2_rm))
loci_location2<-which(loci2_rm%in%intersect(loci1_rm,loci2_rm))
intersect_matrix<-cbind(datatype1$chr[loci_location1],datatype1$loci[loci_location1])
#This function merges 450k and Beadchip data into one big matrix
for (i in 1:length(total_matrix)){
  dataset_sample<-total_matrix[[i]]
  if (nrow(dataset_sample)==865581){
    dataset_sample_intersect<-dataset_sample[loci_location1,-c(1,2)]
    intersect_matrix<-cbind(intersect_matrix, dataset_sample_intersect)
  } else if(nrow(dataset_sample)==485306){
    dataset_sample_intersect<-dataset_sample[loci_location2,-c(1,2)]
    intersect_matrix<-cbind(intersect_matrix, dataset_sample_intersect)
  }

}
write.table(
  intersect_matrix,
  file = paste0('data/', database_type, '/all_betavalues.txt'),
  quote = F,
  row.names = F,
  sep = "\t"
)
gzip(filename=paste0('data/', database_type, '/all_betavalues.txt'), 
     destname=paste0('data/', database_type, '/all_betavalues.txt.gz'), overwrite=TRUE, remove=TRUE, BFR.SIZE=1e+07)
