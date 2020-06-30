# This script saves the beta values for each sample (GSM) in directory ./data/GEO
# Input is the lifted beta values to hg38 for a series and the directory to the metadata file
# Output is methylation beta value files for each GSM sample
write_to_file<-function(GRset_BMIQ_genome_loci_lift, dir2metadata) {
  metadata_table <- read.table(dir2metadata, header = T, sep = '\t')
  col_names<-colnames(GRset_BMIQ_genome_loci_lift)
  # for (i in 3:ncol(GRset_BMIQ_genome_loci_lift)) {
  #   write.table(
  #     data.table(
  #       GRset_BMIQ_genome_loci_lift[, 1:2],
  #       GRset_BMIQ_genome_loci_lift$col_names[i],keep.rownames = T
  #     ),
  #     paste0('./data/', metadata.table[1,'Database'],'/', col_names[i], '_beta_values.txt'),
  #     quote = F,
  #     sep = '\t'
  #   )
  write.table(GRset_BMIQ_genome_loci_lift, paste0('./data/', metadata_table[1,'Database'],'/', metadata_table[1,'Series'], '_beta_values.txt'),
              quote = F,
              sep = '\t', row.names = F)
  }
  

