# This script renames the Blueprint file and makes the metadata file
args <- commandArgs( trailingOnly = TRUE )
manifest <- args[1]
sample_manifest <- args[2]
base_directory <- args[3]
database <- args[4]
if (!require('data.table')) {
  install.packages('data.table', repos = "http://cran.us.r-project.org")
}
library(data.table)
url_list <- read.table(paste0(base_directory,'/../annotation/', manifest), sep = '\n', header = F )
url_list <- as.vector(url_list[,1])
sample_csv<- read.table(paste0(base_directory,'/../annotation/',sample_manifest), sep =  ',', header = T)
for (i in 1:nrow(sample_csv)){
  filenames <- url_list[grep(sample_csv[i, 'id'], url_list)]
  if (file.exists(paste0(base_directory, '/../data/',database,'/', i, '_beta_values.bigWig'))==FALSE){
    if (length(filenames) == 2){
      beta_value_filename<- filenames[grep('methylation_profile', filenames)]
      original_name<- sub(paste0(".*", database, "/hg38/"), "", beta_value_filename)
      replace_name<- paste0(sample_csv[i, 'id'], '_beta_values.bigWig')
      file.rename(paste0(base_directory,'/../data/',database, '/', original_name), paste0(base_directory,'/../data/',database,'/', replace_name))
    }
    else if (length(filenames) == 4){
      signal_filename<- filenames[grep('signal_unstranded.bigWig', filenames)]
      identify_names <- sub(paste0(".*", database, "/hg38/"), "", signal_filename)
      identify_number <- sub(paste0(".", database,".*"),"", identify_names)
      original_name <- identify_names[which(identify_number==min(as.numeric(identify_number)))] 
      replace_name <- paste0(sample_csv[i, 'id'], '_beta_values.bigWig')
      file.rename(paste0(base_directory,'/../data/',database, '/', original_name), paste0(base_directory,'/../data/',database,'/', replace_name))
    }
  }
  
}
sample_name <- sample_csv$id
source <- sample_csv$cell_type
source_class <- sample_csv$biomaterial_type
disease<- sample_csv$disease_ontology_uri
donor_sex <- sample_csv$donor_sex
database_all <- rep(database, nrow(sample_csv))
datatype_all <- rep('WGS', nrow(sample_csv))
metadata<-data.table(
  'Samples' = sample_name,
  'Source' = source,
  'Database' = database_all,
  'Datatype' = datatype_all,
  'Sample_class' = source_class,
  'Disease' = disease,
  'Donor_sex' = donor_sex
)
write.table(
  metadata,
  paste0(base_directory,"/../data/",database,'/',database,
         '_sample_metadata.txt')
  ,
  sep = "\t",
  row.names = FALSE,
  quote = F
)
