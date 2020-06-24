# This script downloads the microarray data in ENCODE
# Input is the accession name of the experiment staring with ENCS
library(ENCODExplorer)
download.encode <- function(accession.name) {
  downloadEncode(file_acc = accession.name,
                 df = get_encode_df(),
                 dir = "./raw/ENCODE")
}