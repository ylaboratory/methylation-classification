# This script is used for downloading raw microarray files from the GEO database
download.data.geo.microarray <- function(accession.num) {
  if (dir.exists(
    paste(
      './raw/GEO/',
      accession.num,
      sep = ""
    )
  ) == FALSE) {
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
  }
  print(paste0('Executed downloading command for ', accession.num))
}
