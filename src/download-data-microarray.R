library(data.table)
library(GEOquery)
library(ENCODExplorer)
library(TCGAutils)
library(TCGAbiolinks)
# This script is used for downloading and extracting raw microarray files from the GEO database
# To donwload data even if directory exists, use download.data.geo.microarray(accession.num, ignore.exisiting = T)
download_data_geo_microarray <- function(accession_num, ignore_exist=F, download_directory='raw/GEO/') {
  if (dir.exists(
    paste(
      download_directory,
      accession_num,
      sep = ""
    )
  ) == FALSE | ignore_exist) {
    result <-
      getGEOSuppFiles(
        accession_num,
        makeDirectory = T,
        baseDir = download_directory,
        fetch_files = T,
        filter_regex = '.tar'
      )
    if (is.null(result)) {
      stop(paste0('no .tar supplement file found in ', accession_num))
    }
    data_dir<-paste0(
      download_directory,
      accession_num
    )
    file_name <-
      list.files(
        path = data_dir,
        pattern = '.tar'
      )
    print('uncompressing the supplement file')
    untar(
      paste0(
        data_dir ,
        '/',
        file_name
      ),
      exdir = data_dir
    )
    
  }
  print(paste0('Executed downloading command for ', accession_num))
  
}

# This file extracts the metadata for a GSE series in GEO database for microarray data
# assume the working directory is the master folder Methylation-classfication
download_geo_metadata <- function(accession_num, out_directory='data/GEO/') {
  gse_series <- getGEO(accession_num, GSEMatrix = F)
  if (length(gse_series)>1){
    stop(paste0('This is a super series, please use the series list: ', names(gse_series)))
  }
  gsm_names <- names(GSMList(gse_series))
  platform_all<-rep(gse_series@header$platform_id, length(gsm_names))
  series_all<-rep(accession_num,length(gsm_names))
  database_all<-rep('GEO',length(gsm_names))
  if (dir.exists(
    out_directory
  )
  == FALSE) {
    dir.create(
      out_directory
    )
  }
  gsm_source_all <- vector()
  gsm_status_all <- vector()
  for (i in 1:length(gsm_names)) {
    gsm_sample <- getGEO(gsm_names[i])
    gsm_source <- gsm_sample@header$source_name_ch1
    gsm_status <- gsm_sample@header$title
    gsm_source_all <- c(gsm_source_all, gsm_source)
    gsm_status_all <- c(gsm_status_all, gsm_status)
  }
  metadata <-
    data.table(
      'Samples' = gsm_names,
      'Source' = gsm_source_all,
      'Title' = gsm_status_all,
      'Series'= series_all,
      'Platform' = platform_all,
      'Database' = database_all
    )
  write.table(
    metadata,
    paste0(out_directory,
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
    paste0(out_directory,
           accession_num,
           '_series_metadata.txt',
           sep = "")
    ,
    sep = "\t",
    row.names = FALSE,
    quote = F
  )
  
}
# This function downloads the microarray data in ENCODE
# Input is the accession name of the experiment staring with ENCS
download_encode <- function(accession_name, download_directory='raw/ENCODE/', ignore_exist=F) {
  if (dir.exists(paste0(download_directory, accession_name))==F | ignore_exist){
    if (dir.exists(paste0(download_directory, accession_name))==F){
      dir.create(paste0(download_directory, accession_name))
    }
    downloadEncode(file_acc = accession_name,
                   dir = paste0(download_directory, accession_name))
  }
  
  
}

# This function renames the downloaded file 
# Input is the accession name starting with ENCS and the files will be renamed as ENCLB_Red/ENCS_Grn
convert2target <- function(accession_name,download_directory='raw/ENCODE/') {
  encode_df<-get_encode_df()
  files_all_rep <- encode_df[accession == accession_name]
  sample_names<-unique(files_all_rep[,replicate_libraries])
  for (i in 1: length(sample_names)){
    files_all<-files_all_rep[which(files_all_rep[,replicate_libraries]==sample_names[i]),]
    pattern_name_red <-
      files_all[which(files_all[, output_type] == 'idat red channel'), file_accession]
    replace_name_red <-
      paste0(files_all[which(files_all[, output_type] == 'idat red channel'), replicate_libraries], '_Red')
    pattern_name_green <-
      files_all[which(files_all[, output_type] == 'idat green channel'), file_accession]
    replace_name_green <-
      paste0(files_all[which(files_all[, output_type] == 'idat green channel'), replicate_libraries], '_Grn')
    file.rename(
      from =  paste0(download_directory,accession_name,'/',pattern_name_red, '.idat'),
      to = paste0(download_directory,accession_name,'/', replace_name_red, '.idat')
    )
    file.rename(
      from =  paste0(download_directory,accession_name,'/', pattern_name_green, '.idat'),
      to = paste0(download_directory,accession_name,'/', replace_name_green, '.idat')
    )
    
  }
  
}

# This file extracts the metadata for all microarray data in the encode database
# Input is the accession name (experiment)
# Output is the encode metadata of the experiment (might contain biological replicates)
get_metadata_encode<-function(accession_name,out_directory='data/ENCODE/') {
  encode_df<-get_encode_df()
  metadata <- encode_df[accession == accession_name]
  metadata <- metadata[which(metadata[,output_type == 'idat green channel']), ]
  sample <- metadata[, replicate_libraries]
  experiment <- metadata[, accession]
  database <- rep('ENCODE', nrow(metadata))
  source_type <- metadata[, biosample_type]
  source <- metadata[, biosample_name]
  platform <- metadata[, platform]
  ENCODE_metadata <-
    data.frame(
      'Samples' = sample,
      'Assay_type' = platform,
      'Source' = source,
      'Source_type' = source_type,
      'Series' = experiment,
      'Database' = database
    )
  
  write.table(
    ENCODE_metadata,
    paste0(out_directory, accession_name, '_sample_metadata.txt'),
    quote = F,
    sep = '\t',
    row.names = F
  )
}

# This script downloads the TCGA microarray data
download_TCGA<-function(manifest_file_name){
  system(paste0('src/download-data-TCGA-microarray.sh -m ', manifest_file_name))
  
}
# This function reorganizes the TCGA microarray samples into folders of barcodes
# Input: manifest file containing all file names
# Output: new folders containing the idat files corresponding to each sample barcode
reorganize_TCGA<-function(file2bar){
  for (i in unique(file2bar[,'associated_entities.entity_submitter_id'])){
    file_name<-file2bar[which(file2bar[,'associated_entities.entity_submitter_id']==i), 'file_id'] 
    for (j in as.character(file_name)){
      target<-list.files(path = paste0('raw/TCGA/',j), pattern = '.idat', full.names = T)
      target_short<-list.files(path = paste0('raw/TCGA/',j), pattern = '.idat')
      if (dir.exists(paste0('raw/TCGA/',i))==F){
        dir.create(paste0('raw/TCGA/',i))
      }
      file.copy(target, paste0('raw/TCGA/',i, '/',target_short))
    }
  }
  
}
# This script gets the metadata of a TCGA microarray sample identified by the sample barcode
# Input: manifest file containing all microarray file names ; directory to write the metadata to
# Output: metadata of each sample
get_metadata_TCGA<-function(file2bar, output_directory='data/TCGA'){
  for (i in unique(file2bar[,'associated_entities.entity_submitter_id'])){
    sample<-colDataPrepare(i)
    extracted_info<-sample[c('barcode','patient','definition','sample_id','sample_type_id','tissue_or_organ_of_origin','name')]
    file_name<-as.character(file2bar[which(file2bar[,'associated_entities.entity_submitter_id']==i), 'file_id'])
    colnames(extracted_info)<-c('Samples','patient','definition','sample_id','sample_type_id','tissue_of_origin','disease_name')
    metadata_TCGA<-data.table(extracted_info, 'file1' = file_name[1], 'file2' = file_name[2])
    write.table(metadata_TCGA, paste0(output_directory,'/', i, '_sample_metadata.txt'), quote = F, row.names = F, sep = '\t')
  }
  
}
