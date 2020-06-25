# This file extracts the .tar microarray files downloaded from the GEO database

file.extract <- function(accession.num) {
  data.dir<-paste0(
    './raw/GEO/',
    accession.num
  )
  file.name <-
    list.files(
      path = data.dir,
      pattern = '.tar'
    )
  untar(
    paste0(
      data.dir ,
      '/',
      file.name
    ),
    exdir = data.dir
  )
}
