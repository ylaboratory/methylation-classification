# This script downloads the microarray data in ENCODE
# Input is the accession name of the experiment staring with ENCS
library(ENCODExplorer)
download.encode <- function(accession.name) {
  if (dir.exists(paste0("./raw/ENCODE/"), accession.name)==F){
    dir.create(paste0("./raw/ENCODE/"), accession.name)
  }
  downloadEncode(file_acc = accession.name,
                 df = get_encode_df(),
                 dir = paste0("./raw/ENCODE/"), accession.name)
}