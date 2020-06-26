# This script is used for downloading and extracting raw microarray files from the GEO database
# To donwload data even if directory exists, use download.data.geo.microarray(accession.num, ignore.exisiting = T)
download.data.geo.microarray <- function(accession.num, ignore.exisit=F) {
  if (dir.exists(
    paste(
      './raw/GEO/',
      accession.num,
      sep = ""
    )
  ) == FALSE | ignore.exisit) {
    result <-
      getGEOSuppFiles(
        accession.num,
        makeDirectory = T,
        baseDir = paste(
          './raw/GEO',
          sep = ""
        ),
        fetch_files = T,
        filter_regex = '.tar'
      )
    if (is.null(result)) {
      stop(paste0('no .tar supplement file found in ', accession.num))
    }
    data.dir<-paste0(
      './raw/GEO/',
      accession.num
    )
    file.name <-
      list.files(
        path = data.dir,
        pattern = '.tar'
      )
    print('uncompressing the supplement file')
    untar(
      paste0(
        data.dir ,
        '/',
        file.name
      ),
      exdir = data.dir
    )
    
  }
  print(paste0('Executed downloading command for ', accession.num))
  
}
