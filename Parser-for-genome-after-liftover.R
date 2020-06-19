# The hg19_850/450_coordinate text file is lifted over to hg38 using the UCSC liftOver tool
# https://genome.ucsc.edu/cgi-bin/hgLiftOver
# Input: coordinate of CpG sites on hg19
# Output: converted coordinate of CpG sites on hg38 and coordinate that failed conversion
# Failed coordinates saved as 'hg19_[450/850]_deleted.txt', converted coordinate saved as 'hg38_[450/850]_coordinate.txt'
# All Outputs saved in /local/usr

# Parse the converted coordinate of hg38
# Input: 'local/hg38_[850/450]_coordinate.txt'
# Output: 'hg38_[850/450]_converted_coordinate.txt'
# liftover2('local/Methylation-classification-preprocessing/annotation/hg38_450_coordinate.txt')
liftover2<-function (dir.to.lifted_genome){
  a<-read.table(dir.to.lifted_genome)
  b<-tidyr::separate(a,V1,c('chr','range'),sep=':',convert=T)
  c<-tidyr::separate(b,'range',c('start','end'),sep='-',convert=T)
  if (nrow(a)>500000){
    write.table(c,paste0(target.dir,'/Methylation-classification-preprocessing/annotation/','hg38_850_converted_coordinate.txt'),quote = F, col.names = T, row.names = F)
  } else{
    write.table(c,paste0(target.dir,'/Methylation-classification-preprocessing/annotation/','hg38_450_converted_coordinate.txt'),quote = F, col.names = T, row.names = F)
    
  }
}
