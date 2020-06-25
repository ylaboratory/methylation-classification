# This script does liftover from hg19 to hg38
# Input is the downloaded chain file from hg19 to hg38 on UCSC website and the beta values with hg19 coordinate
# Output is the beta values with hg38 coordinate
library(rtracklayer)
liftover <- function(dir2chain, GRset.BMIQ.genome_loci) {
  chain <- import.chain(dir2chain)
  seqname <- GRset.BMIQ.genome_loci[, 'chr']
  rownames(seqname) <- NULL
  start <- GRset.BMIQ.genome_loci[, 'loci']
  rownames(start) <- NULL
  Grange19 <-
    GRanges(seqnames = seqname,
            ranges = IRanges(start = as.numeric(start), width = 1))
  Grange38 <- liftOver(Grange19, chain)
  GRset.BMIQ.genome_loci[, 'chr'] <- as.character(seqnames(Grange38))
  GRset.BMIQ.genome_loci[, 'loci'] <- as.numeric(start(Grange38))
  start.new <- GRset.BMIQ.genome_loci[, 'loci']
  rownames(start.new) <- NULL
  GRset.BMIQ.genome_loci_lift <-
    GRset.BMIQ.genome_loci[which(!is.na(start.new)), ]
  return(GRset.BMIQ.genome_loci_lift)
}
