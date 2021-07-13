# This script creates the preprocessig raw object from raw idat files in a accession directory (experiment for ENCODE and series for GEO) and performs background correction
# datadir is the directory where raw idat files are stored
# return an MethylSet object
library(minfi)
library(wateRmelon)
background_correction<-function(accession_datadir){
  file_names <-
    list.files(path =  accession_datadir,
               pattern = 'Grn.idat',
               full.names = T)
  if (length(file_names) == 0) {
    stop(paste0('no idat file found'))
  }
  targets <- data.frame('Basename' = sub('_Grn.idat.*', "", file_names))
  print(targets)
  RGset <- read.metharray.exp(targets = targets, force = T)
  GRset_noob <-
    preprocessNoob(
      RGset,
      offset = 15,
      dyeCorr = TRUE,
      verbose = FALSE,
      dyeMethod = c("single", "reference")
    )
  return(GRset_noob)
}

# This script does the BMIQ normalization of data
# Input is the Methylset object after background correction
# returns a matrix that contains the normalized beta values with hg19 coordinate
normalization<-function(GRset_noob,dir2metadata) {
  print(dir2metadata);
  metadata_table <- read.table(dir2metadata, header = T, sep = '\t', fill = T)
  sample_names <- metadata_table[, 'Samples']
  ratioSet <- ratioConvert(GRset_noob, what = "both", keepCN = TRUE)
  gset <- mapToGenome(ratioSet)
  GRset_BMIQ <- BMIQ(GRset_noob)
  genome_loci <-
    which(row.names(GRset_BMIQ) %in% as.character(gset@rowRanges@ranges@NAMES))
  GRset_BMIQ_genome_loci <- as.matrix(GRset_BMIQ[genome_loci,])
  colnames(GRset_BMIQ_genome_loci) <-colnames(GRset_BMIQ)
  GRset_BMIQ_genome_loci <-
    data.table('chr' = as.character(gset@rowRanges@seqnames),
               'loci' = gset@rowRanges@ranges@start,GRset_BMIQ_genome_loci)
  setnames(GRset_BMIQ_genome_loci, colnames(GRset_BMIQ), as.character(sample_names), skip_absent = T)
  return(GRset_BMIQ_genome_loci)
}