# This script does liftover from hg19 to hg38
# Input is the downloaded chain file from hg19 to hg38 on UCSC website and the beta values with hg19 coordinate
# Output is the beta values with hg38 coordinate
library(rtracklayer)
liftover <- function(dir2chain, GRset_BMIQ_genome_loci) {
  chain <- import.chain(dir2chain)
  seqname <- GRset_BMIQ_genome_loci$chr
  rownames(seqname) <- NULL
  start <- GRset_BMIQ_genome_loci$loci
  rownames(start) <- NULL
  Grange19 <-
    GRanges(seqnames = seqname,
            ranges = IRanges(start = as.numeric(start), width = 1))
  Grange38 <- liftOver(Grange19, chain)
  GRset_BMIQ_genome_loci$chr <- as.character(seqnames(Grange38))
  GRset_BMIQ_genome_loci$loci <- as.numeric(start(Grange38))
  start_new <- GRset_BMIQ_genome_loci$loci
  rownames(start_new) <- NULL
  GRset_BMIQ_genome_loci_lift <-
    GRset_BMIQ_genome_loci[which(!is.na(start_new)), ]
  return(GRset_BMIQ_genome_loci_lift)
}
