setwd('/grain/mk98/methyl/methylation-classification')

geo_features=read.csv("annotation/features.txt", header=FALSE)
counter=1
loci_start=c()
loci_end=c()

#union of all loci starts and loci ends from all SRA samples (sept 2021)
for (folder in list.files('data/SRA/RRBS/', pattern='SRR')){
  print(counter)
  print(folder)
  tryCatch(
    {filename=paste0('data/SRA/RRBS/',folder,"/",folder,"_trimmed_bismark_bt2.bedGraph.gz")
  file<-read.csv(filename, sep="\t", header=FALSE)[-c(1),] #first row omitted because it just says bedGraph
  names(file)<-c('chr', 'loci_start','loci_end','beta')
  loci_start=union(loci_start, paste(paste('chr',file$chr, sep=''),file$loci_start, sep = ' '))
  loci_end=union(loci_end, paste(paste('chr',file$chr, sep=''),file$loci_end, sep = ' '))
  },
  error=function(e) print(paste("ERROR: ",folder)))
  counter=counter+1
}
intersect_start=intersect(geo_features$V1, loci_start)
intersect_end=intersect(geo_features$V1, loci_end)
print(length(intersect_start))
print(length(intersect_end))

write.table(intersect_start, 'annotation/SRA/intersect_start.csv', row.names=FALSE, col.names=FALSE)
write.table(intersect_end, 'annotation/SRA/intersect_end.csv', row.names=FALSE, col.names=FALSE)
write.table(loci_start, 'annotation/SRA/union_start.csv', row.names=FALSE, col.names=FALSE)
write.table(loci_end, 'annotation/SRA/union_end.csv', row.names=FALSE, col.names=FALSE)
