# This file extracts the metadata for a GSE series in GEO database for microarray data
# assume the working directory is the master folder Methylation-classfication
download.geo.series <- function(accession.num) {
  gse.series <- getGEO(accession.num, GSEMatrix = F)
  gsm.names <- names(GSMList(gse.series))
  platform.table <- c(850, 850, 450, 450)
  names(platform.table) <-
    c('GPL21145', 'GPL23976', 'GPL16304', 'GPL13534')
  assay.type <- platform.table[gse.series@header$platform_id]
  assay.type.all <- rep(assay.type, length(gsm.names))
  platform.all<-rep(gse.series@header$platform_id, length(gsm.names))
  series.all<-rep(accession.num,length(gsm.names))
  database.all<-rep('GEO',length(gsm.names))
  if (dir.exists(
      './data/GEO/'
    )
   == FALSE) {
    dir.create(
        './data/GEO/'
      )
  }
  gsm.source.all <- vector()
  gsm.status.all <- vector()
  for (i in 1:length(gsm.names)) {
    gsm.sample <- getGEO(gsm.names[i])
    gsm.source <- gsm.sample@header$source_name_ch1
    gsm.status <- gsm.sample@header$title
    gsm.source.all <- c(gsm.source.all, gsm.source)
    gsm.status.all <- c(gsm.status.all, gsm.status)
  }
  metadata <-
    data.table(
      'Samples' = gsm.names,
      'Assay_type' = assay.type.all,
      'Source' = gsm.source.all,
      'Title' = gsm.status.all,
      'Series'= series.all,
      'Platform' = platform.all,
      'Database' = database.all
    )
  write.table(
    metadata,
      paste0('./data/GEO/',
      accession.num,
      '_sample_metadata.txt',
      sep = "")
    ,
    sep = "\t",
    row.names = FALSE,
    quote = F
  )
  series.relation <- gse.series@header$relation
  series.design <- gse.series@header$overall_design
  series.name <- gse.series@header$geo_accession
  series.supp <- gse.series@header$supplementary_file
  series.title <- gse.series@header$title
  series.info <-
    data.table(
      'name' = series.name,
      'design' = series.design,
      'relation' = series.relation,
      'supplement' = series.supp,
      'title' = series.title,
      key = c('name', 'design', 'relation', 'supplement', 'title')
    )
  write.table(
    series.info,
    paste0('./data/GEO/',
           accession.num,
           '_series_metadata.txt',
           sep = "")
    ,
    sep = "\t",
    row.names = FALSE,
    quote = F
  )
  
}