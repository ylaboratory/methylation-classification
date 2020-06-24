# This file does background correction of microarray data
background.correction<-function(RGset){
  GRset.noob <-
    preprocessNoob(
      RGset,
      offset = 15,
      dyeCorr = TRUE,
      verbose = FALSE,
      dyeMethod = c("single", "reference")
    )
  return(GRset.noob)
}
