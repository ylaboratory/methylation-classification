# This script renames the downloaded file 
# Input is the accession name starting with ENCS and the files will be renamed as ENCLB_Red/ENCS_Grn
convert2target <- function(accession.name) {
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
      from =  paste0('./raw/ENCODE', pattern.name.red, '.idat'),
      to = paste0('./raw/ENCODE', replace.name.red, '.idat')
    )
    file.rename(
      from =  paste0('./raw/ENCODE', pattern.name.green, '.idat'),
      to = paste0('./raw/ENCODE', replace.name.green, '.idat')
    )
    
  }
  
}