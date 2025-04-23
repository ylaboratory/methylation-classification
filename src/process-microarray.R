# This script contains functios for background correction and normalization for raw idat files
# background correction: preprocessNoob from minfi
# normalization: BMIQ from wateRmelon

suppressMessages({
  library(minfi)
  library(wateRmelon)
  library(readr)
})

background_correction<-function(RGset){
  Mset_noob <- preprocessNoob(
    RGset,
    offset = 15,
    dyeCorr = TRUE,
    verbose = FALSE,
    dyeMethod = c("single", "reference")
  )
  return(Mset_noob)
}

normalization <- function(Mset_noob) {
  ratioSet <- ratioConvert(Mset_noob, what = "both", keepCN = TRUE)
  beta <- BMIQ(Mset_noob)
  beta <- beta[rownames(ratioSet), colnames(ratioSet)]
  assay(ratioSet, 'Beta') <- beta
  return(ratioSet)
}
