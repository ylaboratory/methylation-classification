# create a directory for illumina microarray data on GEO
createdir <- function (target.directory) {
  target.dir <<- target.directory
  if (dir.exists(paste0(
    target.directory,
    '/Methylation-classification-preprocessing'
  )) == FALSE) {
    dir.create(paste0(
      target.directory,
      '/Methylation-classification-preprocessing'
    ))
  }
  if (dir.exists(paste0(
    target.directory,
    '/Methylation-classification-preprocessing/raw'
  )) == FALSE) {
    dir.create(paste0(
      target.directory,
      '/Methylation-classification-preprocessing/raw'
    ))
  }
  if (dir.exists(paste0(
    target.directory,
    '/Methylation-classification-preprocessing/data'
  )) == FALSE) {
    dir.create(paste0(
      target.directory,
      '/Methylation-classification-preprocessing/data'
    ))
  }
  if (dir.exists(
    paste0(
      target.directory,
      '/Methylation-classification-preprocessing/annotation'
    )
  ) == FALSE) {
    dir.create(
      paste0(
        target.directory,
        '/Methylation-classification-preprocessing/annotation'
      )
    )
  }
  dir.create(
    paste0(
      target.directory,
      '/Methylation-classification-preprocessing/data/GEO'
    )
  )
  dir.create(paste0(
    target.directory,
    '/Methylation-classification-preprocessing/raw/GEO'
  ))
}


# main function for preprocessing microarray data on GEO
# Input: GEO series accession number (ex. GSE151355)
# Output: metadata and methylation beta values of each sample (saved in local directory)
library(GEOquery)
library(minfi)
library(wateRmelon)
library(data.table)
target.directory <- 'local'
accession.num <- 'GSE146179'
main.microarray.geo(accession.num, target.directory)
main.microarray.geo <- function(accession.num, target.directory) {
  createdir(target.directory)
  download.geo.series(accession.num)
  file.extract(accession.num)
  preprocess.data(accession.num)
}


# Download all raw idat file from GEO
# Input: GSE accession number
# Output: downloaded .idat file in a compressed format under local/usr/GEO/microarray/accession number, series and sample metadata in RDS object form saved under local/usr/GEO/microarray/accession number
# Assumptions: the .idat files are compressed in .tar format and the .tar file is provided in series supplement record
download.geo.series <- function(accession.num) {
  gse.series <- getGEO(accession.num, GSEMatrix = F)
  gsm.names <- names(GSMList(gse.series))
  platform.table <- c(850, 850, 450, 450)
  names(platform.table) <-
    c('GPL21145', 'GPL23976', 'GPL16304', 'GPL13534')
  assay.type <- platform.table[gse.series@header$platform_id]
  if (dir.exists(
    paste(
      target.dir,
      '/Methylation-classification-preprocessing/raw/GEO/',
      accession.num,
      sep = ""
    )
  ) == FALSE) {
    result <-
      getGEOSuppFiles(
        accession.num,
        makeDirectory = T,
        baseDir = paste(
          target.dir,
          '/Methylation-classification-preprocessing/raw/GEO',
          sep = ""
        ),
        fetch_files = T,
        filter_regex = '.tar'
      )
    if (is.null(result)) {
      stop(paste0('no .tar supplement file found in ', accession.num))
    }
  }
  assay.type.all <- rep(assay.type, length(gsm.names))
  if (dir.exists(
    paste0(
      target.dir,
      '/Methylation-classification-preprocessing/data/GEO/',
      accession.num
    )
  ) == FALSE) {
    dir.create(
      paste0(
        target.dir,
        '/Methylation-classification-preprocessing/data/GEO/',
        accession.num
      )
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
      'source' = gsm.source.all,
      'title' = gsm.status.all
    )
  write.table(
    metadata,
    paste(
      target.dir,
      '/Methylation-classification-preprocessing/data/GEO/',
      accession.num,
      '/sample_metadata.txt',
      sep = ""
    ),
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
    paste(
      target.dir,
      '/Methylation-classification-preprocessing/data/GEO/',
      accession.num,
      '/series_metadata.txt',
      sep = ""
    ),
    sep = "\t",
    row.names = FALSE,
    quote = F
  )
  
}




# extract the files
# Input: accession number
# Output: idat files for each sample, saved in local/usr/GEO/microarray/[accession.num]
file.extract <- function(accession.num) {
  file.name <-
    list.files(
      path = paste(
        target.dir,
        '/Methylation-classification-preprocessing/raw/GEO/',
        accession.num ,
        sep = ""
      ),
      pattern = '.tar'
    )
  untar(
    paste(
      target.dir,
      '/Methylation-classification-preprocessing/raw/GEO/',
      accession.num ,
      '/',
      file.name,
      sep = ""
    ),
    exdir = paste(
      target.dir,
      '/Methylation-classification-preprocessing/raw/GEO/',
      accession.num ,
      sep = ""
    )
  )
}


# Preprocess the data
# Input: accession number
# Output: methylation beta values saved in local/usr/GEO/microarray/[accession.num]/[accession.num]_processed
# Assumptions: the .idat files are named as XX_Grn(Red).idatXX
# Major functions: preprocessNoob()--use Noob method for background correction BMIQ()--use BMIQ to normalize the beta values
preprocess.data <- function(accession.num) {
  series.metadata <-
    read.table(
      paste(
        target.dir,
        '/Methylation-classification-preprocessing/data/GEO/',
        accession.num,
        '/sample_metadata.txt',
        sep = ""
      ),
      header = T,
      sep = '\t'
    )
  gsm.names <- series.metadata[, 1]
  datadir <-
    paste(
      target.dir,
      '/Methylation-classification-preprocessing/raw/GEO/',
      accession.num,
      sep = ""
    )
  file.names <-
    list.files(path =  datadir,
               pattern = 'Grn.idat',
               full.names = T)
  if (length(file.names) == 0) {
    stop(paste0('no idat file found for ', accession.num))
  }
  targets <- data.frame('Basename' = sub('_Grn.idat.*', "", file.names))
  RGset <- read.metharray.exp(targets = targets)
  GRset.noob <-
    preprocessNoob(
      RGset,
      offset = 15,
      dyeCorr = TRUE,
      verbose = FALSE,
      dyeMethod = c("single", "reference")
    )
  ratioSet <- ratioConvert(GRset.noob, what = "both", keepCN = TRUE)
  gset <- mapToGenome(ratioSet)
  # Prepare liftover parameters
  if (!file.exists(
    paste0(
      target.dir,
      '/Methylation-classification-preprocessing/annotation/',
      'hg19_',
      series.metadata$Assay_type[1],
      '_coordinate.txt'
    )
  )) {
    liftover1(gset)
    stop(
      'hg19 to hg38 liftover files are needed. please use the extracted hg19 coordinate file saved as local/usr/hg19_[450/850]_coordinate.txt'
    )
  } else {
    if (series.metadata$Assay_type[1] == 450) {
      if (!file.exists(
        paste0(
          target.dir,
          '/Methylation-classification-preprocessing/annotation/',
          'hg19_450_deleted.txt'
        )
      ) ||
      !file.exists(
        paste0(
          target.dir,
          '/Methylation-classification-preprocessing/annotation/',
          'hg38_450_converted_coordinate.txt'
        )
      )) {
        stop('No hg19 transfer failed loci or converted hg38 loci found for 450k platform')
        
      } else {
        hg19.delete <-
          read.table(
            paste0(
              target.dir,
              '/Methylation-classification-preprocessing/annotation/',
              'hg19_450_deleted.txt'
            )
          )
        hg19.lift <-
          read.table(
            paste0(
              target.dir,
              '/Methylation-classification-preprocessing/annotation/',
              'hg38_450_converted_coordinate.txt'
            ),
            header = T
          )
      }
      
    } else {
      if (!file.exists(
        paste0(
          target.dir,
          '/Methylation-classification-preprocessing/annotation/',
          'hg19_850_deleted.txt'
        )
      ) ||
      !file.exists(
        paste0(
          target.dir,
          '/Methylation-classification-preprocessing/annotation/',
          'hg38_850_converted_coordinate.txt'
        )
      )) {
        stop('No hg19 transfer failed loci or converted hg38 loci found for 850k platform')
      } else {
        hg19.delete <-
          read.table(
            paste0(
              target.dir,
              '/Methylation-classification-preprocessing/annotation/',
              'hg19_850_deleted.txt'
            )
          )
        hg19.lift <-
          read.table(
            paste0(
              target.dir,
              '/Methylation-classification-preprocessing/annotation/',
              'hg38_850_converted_coordinate.txt'
            ),
            header = T
          )
      }
      
    }
    chr <- as.character(gset@rowRanges@seqnames)
    loci.start <- gset@rowRanges@ranges@start
    loci.end <- gset@rowRanges@ranges@start
    hg19.coordinate <- paste0(chr, ':', loci.start, '-', loci.end)
    coordinate.delete <- which(hg19.coordinate %in% hg19.delete$V1)
    coordinate.original <- 1:length(hg19.coordinate)
    # BMIQ normalization of methylation beta values
    GRset.BMIQ <- BMIQ(GRset.noob)
    genome_loci <-
      which(row.names(GRset.BMIQ) %in% as.character(gset@rowRanges@ranges@NAMES))
    GRset.BMIQ.genome_loci <- GRset.BMIQ[genome_loci, ]
    # liftover from hg19 to hg38
    GRset.BMIQ.genome_loci_lift <-
      GRset.BMIQ.genome_loci[setdiff(coordinate.original, coordinate.delete), ]
    for (i in 1:ncol(GRset.BMIQ.genome_loci_lift)) {
      Output <-
        cbind(
          'chr' = hg19.lift$chr,
          'loci' = hg19.lift$start,
          GRset.BMIQ.genome_loci_lift[, i]
        )
      colnames(Output)[3] <- colnames(GRset.BMIQ.genome_loci_lift)[i]
      write.table(
        Output,
        paste0(
          target.dir,
          '/Methylation-classification-preprocessing/data/GEO/',
          accession.num,
          '/',
          gsm.names[i],
          '_betavalues.txt'
        ),
        quote = F,
        sep = '\t'
      )
    }
  }
  
  
}

# liftover of 850/450k beadchip CpG sites from hg19 to hg38
# gset is the GRange object converted from minfi package processed Methylset object
# Input: Genomic ranges associated with 450/850k platform
# Output: text file containing CpG genome loci mapped to hg19 (saved in /local/usr)
liftover1 <- function (gset) {
  chr <- as.character(gset@rowRanges@seqnames)
  loci.start <- gset@rowRanges@ranges@start
  loci.end <- gset@rowRanges@ranges@start
  hg19.coordinate <- paste0(chr, ':', loci.start, '-', loci.end)
  platform <- get.platform(length(loci.start))
  write.table(
    paste(chr, ':', loci.start, '-', loci.end, sep = ""),
    file = paste0(
      target.dir,
      '/Methylation-classification-preprocessing/annotation/',
      'hg19_',
      platform,
      '_coordinate.txt'
    ),
    quote = F,
    col.names = F,
    row.names = F
  )
}

get.platform <- function(len) {
  if (len > 500000) {
    return(850)
  } else {
    return(450)
  }
}
