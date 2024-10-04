#date updated: Jun 2022

setwd('/grain/mk98/methyl/methylation-classification')
print(.libPaths())
.libPaths( c( "/home/mk98/R/x86_64-redhat-linux-gnu-library/4.0", .libPaths() ) )

if (!require('R.utils')) {
  install.packages('R.utils') }
if (!require('data.table')) {
  install.packages('data.table') }
library(R.utils)
library(data.table)

database_type='GEO'

# This script merges all beta_value files together and only keeps the loci that are present in all the files
# The script only merges 450k data(485344 loci) and EpicBeadchip data (865613 loci)
all_betavalue_name <-
  # list.files(path = paste0('data/', database_type), pattern = '_beta_values.txt', full.names = T, recursive = F)
  list.files(path = paste0('data/', database_type), pattern = '_beta_values_probe.txt', full.names = T, recursive = F)
print(paste0("betavalue names acquired: length of ",length(all_betavalue_name)))
betavaluedf <-
  lapply(all_betavalue_name, function(x) {
    fread(x, quote = "")
  })
total_matrix<-list()

#Find a representative sample for each Beadchip and 450K
for (i in 1:length(all_betavalue_name)){
  if (dim(fread(all_betavalue_name[i]))[1]>500000){
    # print(i)
    beta850<-i
    for (j in 1:length(all_betavalue_name)){
      if (dim(fread(all_betavalue_name[j]))[1]<500000){
        # print(j)
        beta450<-j
        break
      }
    }
    break
  }
  
}

# This function removes the duplicated loci in 450k and Beadchip data
# datatype1=betavaluedf[[beta850]]#Here put in a dataset in the matrix list that is Beadchip data
# loci1=paste(datatype1$chr,datatype1$loci, sep = ' ')
# duplicate1<-which(loci1%in%loci1[duplicated(loci1)])
# datatype2=betavaluedf[[beta450]]#Here put in a dataset in the matrix list that is 450k data
# loci2=paste(datatype2$chr,datatype2$loci, sep = ' ')
# duplicate2<-which(loci2%in%loci2[duplicated(loci2)])

# datatype1=betavaluedf[[beta850]]#Here put in a dataset in the matrix list that is Beadchip data
# loci1=datatype1$probe
# duplicate1<-which(loci1%in%loci1[duplicated(loci1)])
datatype2=betavaluedf[[beta450]]#Here put in a dataset in the matrix list that is 450k data
loci2=datatype2$probe
duplicate2<-which(loci2%in%loci2[duplicated(loci2)])

# for (i in 1:length(betavaluedf)){
#   dataset_sample<-betavaluedf[[i]]
#   if (nrow(dataset_sample)>800000){
#     if (length(duplicate1)>0){
#       dataset_sample_removed<-dataset_sample[-duplicate1,]
#       total_matrix[[i]]<-dataset_sample_removed
#     }
#     else{
#       total_matrix[[i]]<-dataset_sample
#     }
#   } 
#   else if(nrow(dataset_sample)<500000){
#     if (length(duplicate2)>0){
#       dataset_sample_removed<-dataset_sample[-duplicate2,]
#       total_matrix[[i]]<-dataset_sample_removed
#     }
#     else{
#     total_matrix[[i]]<-dataset_sample
#     }
#   }
# }
to_remove<-list()
for (i in 1:length(betavaluedf)){
  dataset_sample<-betavaluedf[[i]]
  if (nrow(dataset_sample)>800000){
    to_remove<-append(to_remove, i)
  }
}
betavaluedf<-betavaluedf[-strtoi(to_remove)]

for (i in 1:length(betavaluedf)){
  dataset_sample<-betavaluedf[[i]]
  if(nrow(dataset_sample)<500000){
    if (length(duplicate2)>0){
      dataset_sample_removed<-dataset_sample[-duplicate2,]
      total_matrix[[i]]<-dataset_sample_removed
    }
    else{
    total_matrix[[i]]<-dataset_sample
    }
  }
}



# #Here put in a dataset in the deduplicated matrix list that is Beadchip data
# loci1_rm=paste(total_matrix[[beta850]]$chr,total_matrix[[beta850]]$loci, sep = ' ')
# #Here put in a dataset in the deduplicated matrix list that is 450k data
# loci2_rm=paste(total_matrix[[beta450]]$chr,total_matrix[[beta450]]$loci, sep = ' ')

#Here put in a dataset in the deduplicated matrix list that is Beadchip data
# loci1_rm=total_matrix[[beta850]]$probe
#Here put in a dataset in the deduplicated matrix list that is 450k data
# loci2_rm=total_matrix[[beta450]]$probe

#Intersection loci of deduplicated matrices from both Beadchip and 450K
# loci_location1<-which(loci1_rm%in%intersect(loci1_rm,loci2_rm))
# loci_location2<-which(loci2_rm%in%intersect(loci1_rm,loci2_rm))
# intersect_matrix<-cbind(datatype1$chr[loci_location1],datatype1$loci[loci_location1])

#This function merges 450k and Beadchip data into one big matrix
# for (i in 1:length(total_matrix)){
#   dataset_sample<-total_matrix[[i]]
#   if (nrow(dataset_sample)==865581){
#     dataset_sample_intersect<-dataset_sample[loci_location1,-c(1,2)]
#     intersect_matrix<-cbind(intersect_matrix, dataset_sample_intersect)
#   } else if(nrow(dataset_sample)==485306){
#     dataset_sample_intersect<-dataset_sample[loci_location2,-c(1,2)]
#     intersect_matrix<-cbind(intersect_matrix, dataset_sample_intersect)
#   }
# 
# }

# write.table(
#   intersect_matrix,
#   file = paste0('data/', database_type, '/all_betavalues_450K_probe.txt'),
#   quote = F,
#   row.names = F,
#   sep = "\t"
# )

write.table(
  total_matrix,
  file = paste0('data/', database_type, '/450K_betavalues_probe.txt'),
  quote = F,
  row.names = F,
  sep = "\t"
)

gzip(filename=paste0('data/', database_type, '/450K_betavalues_probe.txt'), 
     destname=paste0('data/', database_type, '/450K_betavalues_probe.txt.gz'), 
     overwrite=TRUE, remove=TRUE, BFR.SIZE=1e+07)
