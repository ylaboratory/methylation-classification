# create a directory for illumina microarray data on GEO 
if (dir.exists("local/GEO")==FALSE){
  dir.create("local/GEO")
}
if (dir.exists("local/GEO/microarray")==FALSE){
  dir.create("local/GEO/microarray")
}

# main function for preprocessing microarray data on GEO
# Input: GEO series accession number (ex. GSE151355)
# Output: metadata and methylation beta values of each sample (saved in local directory)
library(GEOquery)
library(minfi)
library(wateRmelon)
main.microarray.geo<-function(accession.num){
  download.geo.series(accession.num)
  file.extract(accession.num)
  preprocess.data(accession.num)
}


# Download all raw idat file from GEO
# Input: GSE accession number 
# Output: downloaded .idat file in a compressed format under local/usr/GEO/microarray/accession number, series and sample metadata in RDS object form saved under local/usr/GEO/microarray/accession number 
# Assumptions: the .idat files are compressed in .tar format and the .tar file is provided in series supplement record
download.geo.series<-function(accession.num){
  gse.series<-getGEO(accession.num, GSEMatrix=F)
  gsm.names<-names(GSMList(gse.series))
  if (gse.series@header$platform_id== 'GPL21145'|gse.series@header$platform_id=='GPL23976'){
    assay.type<-850
  }
  if (gse.series@header$platform_id== 'GPL16304'|gse.series@header$platform_id=='GPL13534'){
    assay.type<-450
  }
  if (dir.exists(paste("local/GEO/microarray/",accession.num,sep = ""))==FALSE){
    getGEOSuppFiles(accession.num, makeDirectory = T, baseDir = paste("local/GEO/microarray",sep = ""),
                    fetch_files = T, filter_regex = '.tar')
  }
  assay.type.all<-rep(assay.type,length(gsm.names))
  metadata<-data.frame('Samples'=gsm.names, 'Assay_type'=assay.type.all)
  write.table(metadata,paste("local/GEO/microarray/",accession.num,'/metadata.txt',sep = ""),sep="\t",row.names=FALSE)
  saveRDS(gse.series,file = paste("local/GEO/microarray/",accession.num,'/',accession.num, '_metadata.rds',sep = "") )
  for (i in 1:length(gsm.names)){
    gsm.info<-getGEO(gsm.names[i])
    saveRDS(gsm.info, file = paste("local/GEO/microarray/",accession.num,'/',gsm.names[i], '_metadata.rds',sep = ""))
  }
}


# extract the files 
# Input: ccession number
# Output: idat files for each sample, saved in local/usr/GEO/microarray/[accession.num]
file.extract<-function(accession.num){
  file.name<-list.files(path = paste("local/GEO/microarray/",accession.num ,sep = ""), pattern = '.tar')
  untar(paste("local/GEO/microarray/",accession.num ,'/',file.name,sep = ""),exdir = paste("local/GEO/microarray/",accession.num ,sep = ""))
}


# Preprocess the data 
# Input: accession number
# Output: methylation beta values saved in local/usr/GEO/microarray/[accession.num]/[accession.num]_processed
# Assumptions: the .idat files are named as XX_Grn(Red).idatXX
# Major functions: preprocessNoob()--use Noob method for background correction BMIQ()--use BMIQ to normalize the beta values
preprocess.data<-function(accession.num){
  series.metadata<-read.table(paste("local/GEO/microarray/",accession.num,'/metadata.txt',sep = ""),header = T, sep = '\t')
  gsm.names<-series.metadata[,1]
  datadir<-paste('local/GEO/microarray/',accession.num, sep = "")
  file.names<-list.files(path =  datadir, pattern = 'Grn.idat', full.names = T)
  targets<-data.frame('Basename'=sub('_Grn.idat.*',"",file.names))
  RGset <- read.metharray.exp(targets = targets)
  GRset.noob<-preprocessNoob(RGset, offset = 15, dyeCorr = TRUE, verbose = FALSE,
                             dyeMethod=c("single", "reference"))
  ratioSet <- ratioConvert(GRset.noob, what = "both", keepCN = TRUE)
  gset <- mapToGenome(ratioSet)
  if (!file.exists(paste0('local/hg19_',series.metadata$Assay_type[1],'_coordinate.txt'))){
    liftover1(gset)
    print('hg19 to hg38 liftover files are needed. please use the extracted hg19 coordinate file saved as local/usr/hg19_[450/850]_coordinate.txt')
    return(NULL)
    } else {
    if (series.metadata$Assay_type[1]==450){
      if (!file.exists('local/hg19_450_deleted.txt') || !file.exists('local/hg38_450_converted_coordinate.txt')){
        print('No hg19 transfer failed loci or converted hg38 loci found for 450k platform')
        return(NULL)
      } else {
        hg19.delete<-read.table('local/hg19_450_deleted.txt')
        hg19.lift<-read.table('local/hg38_450_converted_coordinate.txt',header = T)
      }
      
    } else {
      if (!file.exists('local/hg19_850_deleted.txt') || !file.exists('local/hg38_850_converted_coordinate.txt')){
        print('No hg19 transfer failed loci or converted hg38 loci found for 850k platform')
        return(NULL)
      } else {
        hg19.delete<-read.table('local/hg19_850_deleted.txt')
        hg19.lift<-read.table('local/hg38_850_converted_coordinate.txt',header = T)
      }
      
    }
    chr<-as.character(gset@rowRanges@seqnames)
    loci.start<-gset@rowRanges@ranges@start
    loci.end<-gset@rowRanges@ranges@start
    hg19.coordinate<-paste0(chr,':',loci.start,'-',loci.end)
    GRset.BMIQ<-BMIQ(GRset.noob)
    genome_loci<-which(row.names(GRset.BMIQ) %in% as.character(gset@rowRanges@ranges@NAMES))
    GRset.BMIQ.genome_loci<-GRset.BMIQ[genome_loci,]
    coordinate.delete<-which(hg19.coordinate %in% hg19.delete$V1)
    coordinate.original<-1:length(hg19.coordinate)
    GRset.BMIQ.genome_loci_lift<-GRset.BMIQ.genome_loci[setdiff(coordinate.original,coordinate.delete),]
    if (dir.exists(paste("local/GEO/microarray/",accession.num,'/', accession.num, '_processed/',sep = ""))==FALSE){
      dir.create(paste("local/GEO/microarray/",accession.num,'/', accession.num, '_processed/',sep = ""))
    }
    for (i in 1:ncol(GRset.BMIQ.genome_loci_lift)){
      Output<-cbind('chr'=hg19.lift$chr, 'loci'=hg19.lift$start, GRset.BMIQ.genome_loci_lift[,i])
      colnames(Output)[3]<-colnames(GRset.BMIQ.genome_loci_lift)[i]
      write.table(Output,paste0(datadir,'/', accession.num, '_processed/', gsm.names[i],'.txt'),quote = F, sep = '\t')
    }
  }
  
  
}

# liftover of 850/450k beadchip CpG sites from hg19 to hg38
# gset is the GRange object converted from minfi package processed Methylset object
# Input: Genomic ranges associated with 450/850k platform
# Output: text file containing CpG genome loci mapped to hg19 (saved in /local/usr)
liftover1<-function (gset){
  chr<-as.character(gset@rowRanges@seqnames)
  loci.start<-gset@rowRanges@ranges@start
  loci.end<-gset@rowRanges@ranges@start
  hg19.coordinate<-paste0(chr,':',loci.start,'-',loci.end)
  if (length(loci.start)>500000){
    write.table(paste(chr,':',loci.start,'-',loci.end,sep = ""),file = 'local/hg19_850_coordinate.txt',quote = F,col.names = F, row.names = F )
  } else{
    write.table(paste(chr,':',loci.start,'-',loci.end,sep = ""),file = 'local/hg19_450_coordinate.txt',quote = F,col.names = F, row.names = F )
    
  }
}

