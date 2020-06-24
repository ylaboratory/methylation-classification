# This script does the BMIQ normalization of data
# Input is the Methylset object after background correction
# returns a matrix that contains the normalized beta values with hg19 coordinate
function(GRset.noob) {
  ratioSet <- ratioConvert(GRset.noob, what = "both", keepCN = TRUE)
  gset <- mapToGenome(ratioSet)
  GRset.BMIQ <- BMIQ(GRset.noob)
  genome_loci <-
    which(row.names(GRset.BMIQ) %in% as.character(gset@rowRanges@ranges@NAMES))
  GRset.BMIQ.genome_loci <- GRset.BMIQ[genome_loci,]
  GRset.BMIQ.genome_loci <-
    cbind(
      'chr' = as.character(gset@rowRanges@seqnames),
      'loci' = gset@rowRanges@ranges@start,
      GRset.BMIQ.genome_loci
    )
  return(GRset.BMIQ.genome_loci)
}