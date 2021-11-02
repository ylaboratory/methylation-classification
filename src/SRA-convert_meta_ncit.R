setwd('/grain/mk98/methyl/methylation-classification')

if (!require('R.uilts')) {
  install.packages("R.utils")
}
library(R.utils)

meta<-read.csv('data/SRA/SRA-RRBS-metadata_manual.txt', header=TRUE, sep='\t')
unique(meta$tissue)
meta[meta=='Nerves']<-'Nerves'
meta[meta=='Thyroid']<-'Thyroid Gland'
meta[meta=='Stem Cells']<-'Stem Cell'
meta[meta=='Prostate']<-'Prostate Gland'
write.table(meta, 'data/SRA/total_metadata.txt',
          row.names=F,
          sep="\t",
          quote=F)
gzip(filename='data/SRA/total_metadata.txt', 
     destname='data/SRA/total_metadata.txt.gz',
     overwrite=TRUE, remove=TRUE)
