# liftover of 850/450k beadchip CpG sites from hg19 to hg38
# gset is the GRange object converted from minfi package processed Methylset object
# Input: Genomic ranges associated with 450/850k platform
# Output: text file containing CpG genome loci mapped to hg19 (saved in /local/usr)
chr<-as.character(gset@rowRanges@seqnames)
loci.start<-gset@rowRanges@ranges@start
loci.end<-gset@rowRanges@ranges@start
hg19.coordinate<-paste0(chr,':',loci.start,'-',loci.end)
write.table(paste(chr,':',loci.start,'-',loci.end,sep = ""),file = 'local/hg19_[850/450]_coordinate.txt',quote = F,col.names = F, row.names = F )

# The hg19_850/450_coordinate text file is lifted over to hg38 using the UCSC liftOver tool
# https://genome.ucsc.edu/cgi-bin/hgLiftOver
# Input: coordinate of CpG sites on hg19
# Output: converted coordinate of CpG sites on hg38 and coordinate that failed conversion
# Failed coordinates saved as 'hg19_[450/850]_deleted.txt', converted coordinate saved as 'hg38_[450/850]_coordinate.txt'
# All Outputs saved in /local/usr

# Parse the converted coordinate of hg38
# Input: 'local/hg38_[850/450]_coordinate.txt'
# Output: 'hg38_[850/450]_converted_coordinate.txt'
a<-read.table('local/hg38_[850/450]_coordinate.txt')
b<-tidyr::separate(a,V1,c('chr','range'),sep=':',convert=T)
c<-tidyr::separate(b,'range',c('start','end'),sep='-',convert=T)
write.table(c,'hg38_850_converted_coordinate.txt',quote = F, col.names = T, row.names = F)