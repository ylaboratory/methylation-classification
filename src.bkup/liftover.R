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
  Grange38_elementlen<-elementNROWS(Grange38)
  rownames(Grange38_elementlen)=NULL
  unmatched<-which(Grange38_elementlen==0)
  GRset_BMIQ_genome_loci_lift=GRset_BMIQ_genome_loci[-unmatched,]
  GRset_BMIQ_genome_loci_lift$chr=as.character(data.frame(Grange38@unlistData)$seqnames)
  GRset_BMIQ_genome_loci_lift$loci=as.character(data.frame(Grange38@unlistData)$start)
  return(GRset_BMIQ_genome_loci_lift)
}
