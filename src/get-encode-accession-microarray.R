# This script queries the data_set accessions of all ENCODE microarray data
library(ENCODExplorer)
query_results <- queryEncode(organism = "Homo sapiens", 
                             file_format = "idat",
                             fixed = TRUE, fuzzy = F)
write.table(unique(query_results$accession), 'annotation/ENCODE_microarray_accession_full_raw.txt', quote = F, sep = '\n', row.names = F, col.names = F)

