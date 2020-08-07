# Take in accession and database from the command line
args <- commandArgs(trailingOnly = TRUE)
accession_num<- args[1]
# Install packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (!require('GEOquery')) {
  BiocManager::install("GEOquery")
}
if (!require('data.table')) {
  install.packages('data.table')
}
library(data.table)
library(GEOquery)
# This file extracts the metadata for a GSE series in GEO database for microarray data
# assume the working directory is the master folder Methylation-classfication
out_directory='data/GEO/'
gse_series <- getGEO(accession_num, GSEMatrix = F)
if (length(gse_series) > 1) {
  stop(paste0(
    'This is a super series, please use the series list: ',
    names(gse_series)
  ))
}
gsm_names <- names(GSMList(gse_series))
platform_all <-
  rep(gse_series@header$platform_id, length(gsm_names))
series_all <- rep(accession_num, length(gsm_names))
database_all <- rep('GEO', length(gsm_names))
if (dir.exists(out_directory)
    == FALSE) {
  dir.create(out_directory)
}
gsm_source_all <- vector()
gsm_status_all <- vector()
gsm_datatype_all <- vector()
sra_acc_all<- vector()
for (i in 1:length(gsm_names)) {
  gsm_sample <- getGEO(gsm_names[i])
  gsm_source <- gsm_sample@header$source_name_ch1
  gsm_status <- gsm_sample@header$title
  gsm_datatype <- gsm_sample@header$library_selection
  sra_string <- gsm_sample@header$relation[which(grep("SRA:", gsm_sample@header$relation, fixed = T)==T)]
  sra_acc<- sub("SRA.*term=", "", sra_string)
  sra_acc_all<- c(sra_acc_all, sra_acc)
  gsm_source_all <- c(gsm_source_all, gsm_source)
  gsm_status_all <- c(gsm_status_all, gsm_status)
  gsm_datatype_all<- c(gsm_datatype_all, gsm_datatype)
}
metadata <-
  data.table(
    'Samples' = gsm_names,
    'Source' = gsm_source_all,
    'Title' = gsm_status_all,
    'Series' = series_all,
    'Platform' = platform_all,
    'Database' = database_all,
    'Datatype' = gsm_datatype_all,
    "Sra_accession" = sra_acc_all
  )
write.table(
  metadata,
  paste0(out_directory, accession_num, "/",
         accession_num,
         '_sample_metadata.txt',
         sep = "")
  ,
  sep = "\t",
  row.names = FALSE,
  quote = F
)
series_relation <- gse_series@header$relation
series_design <- gse_series@header$overall_design
series_name <- gse_series@header$geo_accession
series_supp <- gse_series@header$supplementary_file
series_title <- gse_series@header$title
series_info <-
  data.table(
    'name' = series_name,
    'design' = series_design,
    'relation' = series_relation,
    'supplement' = series_supp,
    'title' = series_title,
    key = c('name', 'design', 'relation', 'supplement', 'title')
  )
write.table(
  series_info,
  paste0(out_directory,accession_num, "/",
         accession_num, 
         '_series_metadata.txt',
         sep = "")
  ,
  sep = "\t",
  row.names = FALSE,
  quote = F
)

