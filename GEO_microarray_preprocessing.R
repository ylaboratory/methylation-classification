# create a directory for illumina microarray data on GEO 
if (dir.exists("local/GEO")==FALSE){
  dir.create("local/GEO")
}
if (dir.exists("local/GEO/microarray")==FALSE){
  dir.create("local/GEO/microarray")
}

# main function for preprocessing microarray data on GEO
library(GEOquery)
library(minfi)
library(wateRmelon)
main.microarray.geo<-function(accession.num){
  download.geo.series(accession.num)
  file.extract(accession.num)
  preprocess.data(accession.num)
}


# download all raw idat file from GEO
# Input: GSE accession number 
# Output: downloaded .idat file in a compressed format under local/GEO/microarray

# accession.num<-"GSE146179"
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
}

# extract the files 
# Input: the accession number
# Output: zipped files for each sample and a constructed csv file
file.extract<-function(accession.num){
  file.name<-list.files(path = paste("local/GEO/microarray/",accession.num ,sep = ""), pattern = '.tar')
  untar(paste("local/GEO/microarray/",accession.num ,'/',file.name,sep = ""),exdir = paste("local/GEO/microarray/",accession.num ,sep = ""))
}


# Preprocess the data 
# Input: directory of the dataset (csv and idat file should be available in the directory)
# Output: methylation beta values stored in the same directory
preprocess.data<-function(accession.num){
  series.metadata<-read.table(paste("local/GEO/microarray/",accession.num,'/metadata.txt',sep = ""),header = T, sep = '\t')
  gsm.names<-series.metadata[,1]
  datadir<-paste('local/GEO/microarray/',accession.num, sep = "")
  file.names<-list.files(path =  datadir, pattern = 'Grn.idat', full.names = T)
  # targets<-read.metharray.sheet(datadir, pattern ='csv$', recursive = T)
  # targets['Basename']=sub('_Grn.idat.*',"",file.names)
  targets<-data.frame('Basename'=sub('_Grn.idat.*',"",file.names))
  RGset <- read.metharray.exp(targets = targets)
  GRset.noob<-preprocessNoob(RGset, offset = 15, dyeCorr = TRUE, verbose = FALSE,
                             dyeMethod=c("single", "reference"))
  ratioSet <- ratioConvert(GRset.noob, what = "both", keepCN = TRUE)
  gset <- mapToGenome(ratioSet)
  GRset.BMIQ<-BMIQ(GRset.noob)
  genome_loci<-which(row.names(GRset.BMIQ) %in% as.character(gset@rowRanges@ranges@NAMES))
  GRset.BMIQ.genome_loci<-GRset.BMIQ[genome_loci,]
  if (dir.exists(paste("local/GEO/microarray/",accession.num,'/processed',sep = ""))==FALSE){
    dir.create(paste("local/GEO/microarray/",accession.num,'/processed',sep = ""))
  }
  for (i in 1:ncol(GRset.BMIQ.genome_loci)){
    Output<-cbind('chr'=as.character(gset@rowRanges@seqnames), 'loci'=gset@rowRanges@ranges@start, GRset.BMIQ.genome_loci[,i])
    colnames(Output)[3]<-colnames(GRset.BMIQ.genome_loci)[i]
    if (dir.exists(paste("local/GEO/microarray/",accession.num,'/processed/',gsm.names[i],sep = ""))==FALSE){
      dir.create(paste("local/GEO/microarray/",accession.num,'/processed/',gsm.names[i],sep = ""))
    }
    write.table(Output,paste0(datadir,'/processed/',gsm.names[i],'/',gsm.names[i],'.txt'),quote = F, sep = '\t')
  }
  
}
