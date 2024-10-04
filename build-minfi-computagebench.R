# Mirae Sunny Kim
# Aug 2024

### settings for vsc ###
# local(
#   {r <- getOption("repos")
#        r["CRAN"] <- "http://cran.r-project.org"
#        options(repos = r)
# }
# )
# library(vscDebugger)
######

#accession list:  annotation/GEO_microarray_accession_450K.txt
#database:        GEO

manifest_file <- "computagebench_hc.txt"
database_type <- "GEO"

#wd needs to change accordingly
setwd("/grain/mk98/methyl/methylation-classification")
print(.libPaths())
.libPaths(c("/home/mk98/R/x86_64-redhat-linux-gnu-library/4.2", .libPaths()))

if (!require("readr")) {
  install.packages("readr")
}
if (!require("RPMM")) {
  install.packages("RPMM", repos = "http://cran.us.r-project.org")
}
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "http://cran.us.r-project.org", version = "3.16")
}
if (!require("minfi")) {
  BiocManager::install("minfi")
}
if (!require("GEOquery")) {
  BiocManager::install("GEOquery")
}
if (!require("wateRmelon")) {
  BiocManager::install("wateRmelon")
}
if (!require("IlluminaHumanMethylationEPICmanifest")) {
  BiocManager::install("IlluminaHumanMethylationEPICmanifest")
}
if (!require("IlluminaHumanMethylation450kmanifest")) {
  BiocManager::install("IlluminaHumanMethylation450kmanifest")
}
if (!require("rtracklayer")) {
  BiocManager::install("rtracklayer")
}
# if (!require("TCGAutils")) {
#   BiocManager::install("TCGAutils")
# }
# if (!require("TCGAbiolinks")) {
#   BiocManager::install("TCGAbiolinks")
# }
if (!require("data.table")) {
  install.packages("data.table", repos = "http://cran.us.r-project.org")
}

if (!require("wateRmelon")) {
  install.packages("wateRmelon")
}

library("R.utils")

accession_list <- read.table("annotation/computagebench_hc.txt", header = F)
sample_list <- read.table("annotation/computagebench_samples.csv", header = F)
probe_list <- read.table("annotation/computagebench_probes.csv", header = F)
# metadata<-read.table('data/GEO/compiled/450K_meta_atleast2.csv.gz', sep=',', header=T)

targets <- data.frame()
sample_names_all <- list()
first = TRUE

for (accession in as.character(accession_list$V1)) {
  # print(accession)
  accession_datadir <- paste0("raw/", database_type, "/", accession)
  gr_file_names <- sub("_Grn.idat.*", "", list.files(path =  accession_datadir,
                                                     pattern = "Grn.idat",
                                                     full.names = T)
  )
  rd_file_names <- sub("_Red.idat.*", "", list.files(path =  accession_datadir,
                                                     pattern = "Red.idat",
                                                     full.names = T)
  )
  file_names<- intersect(gr_file_names, rd_file_names)
  sample_names<- sub("_.*", "", sub(".*GSM", "GSM", file_names))
  file_names = file_names[sample_names %in% sample_list$V1]
  sample_names = sample_names[sample_names %in% sample_list$V1]
  
  targets <- rbind(targets, data.frame('Basename' = file_names))
  sample_names_all<- c(sample_names_all,sample_names)
}

RGset <- read.metharray(basenames= targets$Basename, force=TRUE, verbose=FALSE)

GRset_noob <-preprocessNoob(RGset)
# GRset_funnorm <-preprocessFunnorm(RGset)

ratioSet <- ratioConvert(GRset_noob, what = "both", keepCN = TRUE)
GRset_BMIQ <- BMIQ(GRset_noob)

datatable_BMIQ<-data.table(GRset_BMIQ)
setnames(datatable_BMIQ, colnames(datatable_BMIQ), as.character(sample_names_all), skip_absent = T)

datatable_BMIQ$probe<-as.character(ratioSet@NAMES)
datatable_filtered<-datatable_BMIQ[datatable_BMIQ$probe %in% probe_list$V1]

##### disregard below - saving for comparison purposes #####

write.table(datatable_filtered, paste0('data/ComputAgeBench/compiled/raw_hc.csv'),
            quote = F,
            sep = '\t', row.names = T, col.names = T)

gzip(filename=paste0('data/ComputAgeBench/compiled/raw_hc.csv'), 
     destname=paste0('data/ComputAgeBench/compiled/raw_hc.csv.gz'), 
     overwrite=TRUE, remove=TRUE, BFR.SIZE=1e+07)
