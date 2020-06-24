# This file extracts the metadata for all microarray data in the encode database
# Input is the accession name (experiment)
# Output is the encode metadata of the experiment (might contain biological replicates)
function(accession.name) {
  metadata <- encode_df[accession == accession.name]
  metadata <- metadata[which(metadata[,output_type == 'idat green channel']), ]
  sample <- metadata[, replicate_libraries]
  experiment <- metadata[, accession]
  database <- rep('ENCODE', nrow(metadata))
  source.type <- metadata[, biosample_type]
  source <- metadata[, biosample_name]
  platform <- metadata[, platform]
  ENCODE.metadata <-
    data.frame(
      'Samples' = sample,
      'Assay_type' = platform,
      'Source' = source,
      'Source_type' = source.type,
      'Experiment' = experiment,
      'Database' = database
    )
  
  write.table(
    ENCODE.metadata,
    paste0('./data/ENCODE/', accession.name, '_metadata.txt'),
    quote = F,
    sep = '\t',
    row.names = F
  )
}