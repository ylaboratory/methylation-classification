# This script saves the beta values for each sample (GSM) in directory ./data/GEO
# Input is the lifted beta values to hg38 for a series and the directory to the metadata file
# Output is methylation beta value files for each GSM sample
function(GRset.BMIQ.genome_loci_lift,
         dir2metadata) {
  metadata.table <- read.table(dir2metadata, header = T, sep = '\t')
  sample.names <- metadata.table[, 'Sample']
  for (i in 3:ncol(GRset.BMIQ.genome_loci_lift)) {
    write.table(
      cbind(
        GRset.BMIQ.genome_loci_lift[, 1:3],
        GRset.BMIQ.genome_loci_lift[, i]
      ),
      paste0('./data/', metadata.table[1,'Database'],'/', gsm.names[i - 2], '_beta_values.txt'),
      quote = F,
      sep = '\t'
    )
  }
  
}