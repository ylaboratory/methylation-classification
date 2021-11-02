database_type='GEO'
# This script merges all sample metadata from the corresponding dataset
all_metadata_name <-
  list.files(path = paste0('data/', database_type), pattern = '_sample_metadata.txt', full.names = T)
metadata_total <-
  rbindlist(lapply(all_metadata_name, function(x) {
    fread(x)
  }), use.names = T, fill = T)
write.table(
  metadata_total,
  file = paste0('data/', database_type, '/all_metadata.txt'),
  quote = F,
  row.names = F,
  sep = "\t"
)
library(R.utils)
library(data.table)
gzip(filename=paste0('data/', database_type, '/all_metadata.txt'), 
     destname=paste0('data/', database_type, '/all_metadata.txt.gz'), overwrite=TRUE, remove=TRUE)
# This script merges all beta_value files together and only keeps the loci that are present in all the files
# The script only merges 450k data(485344 loci) and EpicBeadchip data (865613 loci)
all_betavalue_name <-
  list.files(path = paste0('data/', database_type), pattern = '_beta_values.txt', full.names = T, recursive = F)
betavaluedf <-
  lapply(all_betavalue_name, function(x) {
    fread(x, quote = "")
  })
total_matrix<-list()

# This function removes the duplicated loci in 450k and Beadchip data
for (i in 1:length(betavaluedf)){
  dataset_sample<-betavaluedf[[i]]
  datatype1=betavaluedf[[1]]#Here put in a dataset in the matrix list that is Beadchip data
  loci1=paste(datatype1$chr,datatype1$loci, sep = ' ')
  duplicate1<-which(loci1%in%loci1[duplicated(loci1)])
  datatype2=betavaluedf[[3]]#Here put in a dataset in the matrix list that is 450k data
  loci2=paste(datatype2$chr,datatype2$loci, sep = ' ')
  duplicate2<-which(loci2%in%loci2[duplicated(loci2)])
  if (nrow(dataset_sample)>800000){
    dataset_sample_removed<-dataset_sample[-duplicate1,]
    total_matrix[[i]]<-dataset_sample_removed
  } else if(nrow(dataset_sample)<500000){
    dataset_sample_removed<-dataset_sample[-duplicate2,]
    total_matrix[[i]]<-dataset_sample_removed
  }
}
#Here put in a dataset(like the second in total_matrix) in the deduplicated matrix list that is Beadchip data
loci1_rm=paste(total_matrix[[1]]$chr,total_matrix[[1]]$loci, sep = ' ')
#Here put in a dataset in the deduplicated matrix list that is 450k data
loci2_rm=paste(total_matrix[[3]]$chr,total_matrix[[3]]$loci, sep = ' ')
# This script merges all beta_value files together and only keeps the loci that are present in all the files
# The script only merges 450k data(485344 loci) and EpicBeadchip data (865613 loci)
all_betavalue_name <-
  list.files(path = paste0('data/', database_type), pattern = '_beta_values.txt', full.names = T, recursive = F)
betavaluedf <-
  lapply(all_betavalue_name, function(x) {
    fread(x, quote = "")
  })
total_matrix<-list()

# This function removes the duplicated loci in 450k and Beadchip data
for (i in 1:length(betavaluedf)){
  dataset_sample<-betavaluedf[[i]]
  datatype1=betavaluedf[[1]]#Here put in a dataset in the matrix list that is Beadchip data
  loci1=paste(datatype1$chr,datatype1$loci, sep = ' ')
  duplicate1<-which(loci1%in%loci1[duplicated(loci1)])
  datatype2=betavaluedf[[3]]#Here put in a dataset in the matrix list that is 450k data
  loci2=paste(datatype2$chr,datatype2$loci, sep = ' ')
  duplicate2<-which(loci2%in%loci2[duplicated(loci2)])
  if (nrow(dataset_sample)>800000){
    dataset_sample_removed<-dataset_sample[-duplicate1,]
    total_matrix[[i]]<-dataset_sample_removed
  } else if(nrow(dataset_sample)<500000){
    dataset_sample_removed<-dataset_sample[-duplicate2,]
    total_matrix[[i]]<-dataset_sample_removed
  }
}
#Here put in a dataset(like the second in total_matrix) in the deduplicated matrix list that is Beadchip data
loci1_rm=paste(total_matrix[[1]]$chr,total_matrix[[1]]$loci, sep = ' ')
#Here put in a dataset in the deduplicated matrix list that is 450k data
loci2_rm=paste(total_matrix[[3]]$chr,total_matrix[[3]]$loci, sep = ' ')
loci_location1<-which(loci1_rm%in%intersect(loci1_rm,loci2_rm))
loci_location2<-which(loci2_rm%in%intersect(loci1_rm,loci2_rm))
intersect_matrix<-cbind(datatype1$chr[loci_location1],datatype1$loci[loci_location1])
#This function merges 450k and Beadchip data into one big matrix
for (i in 1:length(total_matrix)){
  dataset_sample<-total_matrix[[i]]
  if (nrow(dataset_sample)==865581){
    dataset_sample_intersect<-dataset_sample[loci_location1,-c(1,2)]
    intersect_matrix<-cbind(intersect_matrix, dataset_sample_intersect)
  } else if(nrow(dataset_sample)==485306){
    dataset_sample_intersect<-dataset_sample[loci_location2,-c(1,2)]
    intersect_matrix<-cbind(intersect_matrix, dataset_sample_intersect)
  }

}
write.table(
  intersect_matrix,
  file = paste0('data/', database_type, '/all_betavalues.txt'),
  quote = F,
  row.names = F,
  sep = "\t"
)
gzip(filename=paste0('data/', database_type, '/all_betavalues.txt'), 
     destname=paste0('data/', database_type, '/all_betavalues.txt.gz'), overwrite=TRUE, remove=TRUE, BFR.SIZE=1e+07)

