# This function downloads the microarray data in ENCODE
# Input is the accession name of the experiment staring with ENCS
library(ENCODExplorer)
download.encode <- function(accession.name) {
  if (dir.exists(paste0("./raw/ENCODE/", accession.name))==F){
    dir.create(paste0("./raw/ENCODE/", accession.name))
  }
  downloadEncode(file_acc = accession.name,
                 dir = paste0("./raw/ENCODE/", accession.name))
}

# This function renames the downloaded file 
# Input is the accession name starting with ENCS and the files will be renamed as ENCLB_Red/ENCS_Grn
library(ENCODExplorer)
convert2target <- function(accession.name) {
  encode_df<-get_encode_df()
  files.all.rep <- encode_df[accession == accession.name]
  sample.names<-unique(files.all.rep[,replicate_libraries])
  for (i in 1: length(sample.names)){
    files.all<-files.all.rep[which(files.all.rep[,replicate_libraries]==sample.names[i]),]
    pattern.name.red <-
      files.all[which(files.all[, output_type] == 'idat red channel'), file_accession]
    replace.name.red <-
      paste0(files.all[which(files.all[, output_type] == 'idat red channel'), replicate_libraries], '_Red')
    pattern.name.green <-
      files.all[which(files.all[, output_type] == 'idat green channel'), file_accession]
    replace.name.green <-
      paste0(files.all[which(files.all[, output_type] == 'idat green channel'), replicate_libraries], '_Grn')
    file.rename(
      from =  paste0('./raw/ENCODE/',accession.name,'/',pattern.name.red, '.idat'),
      to = paste0('./raw/ENCODE/',accession.name,'/', replace.name.red, '.idat')
    )
    file.rename(
      from =  paste0('./raw/ENCODE/',accession.name,'/', pattern.name.green, '.idat'),
      to = paste0('./raw/ENCODE/',accession.name,'/', replace.name.green, '.idat')
    )
    
  }
  
}