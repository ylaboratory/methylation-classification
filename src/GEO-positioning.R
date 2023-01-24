#date updated: Jun 2022

setwd('/grain/mk98/methyl/methylation-classification')
print(.libPaths())
.libPaths( c( "/home/mk98/R/x86_64-redhat-linux-gnu-library/4.0", .libPaths() ) )

if (!require('R.utils')) {
  install.packages('R.utils') }
if (!require('readr')) {
  install.packages('readr') }
if (!require('data.table')) {
  install.packages('data.table') }
# library(R.utils)
# library(data.table)
# library(readr)

database_type='GEO'
manifest<-read_tsv('./annotation/HM450.hg38.manifest.gencode.v36.tsv')
manifest$CpG<-manifest$CpG_beg+1
sample_beta_name <-
  list.files(path = paste0('data/', database_type), pattern = '_beta_values_probe.txt', full.names = T, recursive = F)[1]
betavaluedf <- fread(sample_beta_name, quote = "")

chr=list()
pos=list()
isl_pos=list()
isl=list()
iter=1
range=seq(1,dim(betavaluedf)[1], by=1)
for (p in betavaluedf$probe[range]){
  probe_row<-manifest[manifest$probeID==p,]
  chr<-append(chr,probe_row$CpG_chrm)
  pos<-append(pos,probe_row$CpG)
  isl_pos<-append(isl_pos, probe_row$CGIposition)
  isl<-append(isl, probe_row$CGI)
  if (!iter%%10000){print(iter)}
  iter=iter+1
}

characteristics<-data.table(
  'probe'<-betavaluedf$probe[range],
  'chr'<-chr,
  'position'<-pos,
  'CGI'<-isl,
  'CGI type'<-isl_pos
)
print(characteristics[dim(betavaluedf)[1]])

characteristics_df<-data.frame(lapply(characteristics, as.character), stringsAsFactors=FALSE)
print(characteristics_df[dim(betavaluedf)[1],])

write.table(
  characteristics_df,
  file = paste0('data/', database_type, '/450K_cpg.txt'),
  quote = F,
  row.names = F,
  sep = "\t"
)

gzip(filename=paste0('data/', database_type, '/450K_cpg.txt'), 
     destname=paste0('data/', database_type, '/450K_cpg.txt.gz'), 
     overwrite=TRUE, remove=TRUE, BFR.SIZE=1e+07)
