
setwd('/grain/mk98/methyl/methylation-classification') #needs to change accordingly
print(.libPaths())
.libPaths( c( "/home/mk98/R/x86_64-redhat-linux-gnu-library/4.0", .libPaths() ) ) #needs to change accordingly
library("TCGAbiolinks")
# if (!require("TCGAbiolinks")) {
#   BiocManager::install("BioinformaticsFMRP/TCGAbiolinks")
# }
# BiocManager::install("BioinformaticsFMRP/TCGAbiolinks")

# source('./src/download-data-microarray.R')
# source('./src/process-microarray.R')
# source('./src/liftover.R')

# manifest_file<-'gdc_manifest.2024-04-19.beta.value.txt'
# TCGA_manifest<-read.table(paste0('annotation/', manifest_file), header = T, sep = '\t')
# print(TCGA_manifest$id)

# cancer_project = "TCGA-LIHC"

cat("Please enter a cancer (e.g. TCGA-BRCA):")
user_input <- readLines(con = "stdin", n = 1)
cancer_project = paste0("", user_input)

query <- GDCquery(project = cancer_project, data.category = "DNA Methylation",
                  data.type = "Methylation Beta Value",
                  platform = "Illumina Human Methylation 450",
                  # barcode = TCGA_manifest$id,
)

num_samples<-dim(query$results[[1]])[1] #how many samples

cases <- query$results[[1]]$cases[1:num_samples]
query <- GDCquery(project = cancer_project, data.category = "DNA Methylation",
                  data.type = "Methylation Beta Value",
                  platform = "Illumina Human Methylation 450",
                  barcode = cases,
)
GDCdownload(query, 
            directory='/grain/mk98/methyl/methylation-classification/data/TCGA/',
            files.per.chunk=20
            )
data <- GDCprepare(query,
                   directory='/grain/mk98/methyl/methylation-classification/data/TCGA/'
                   )

# beta values in ...data/TCGA/[project]/harmonized/[data category]/[data type]/[id]/[file_name]
# sample metadata in ...data/TCGA/[cases]_sample_metadata.txt

if (grepl("BRCA", cancer_project)){
  betavalues_folder <- paste0('/grain/mk98/methyl/methylation-classification/data/TCGA/',cancer_project,'/harmonized/DNA_Methylation/Methylation_Beta_Value/')
} else {
  betavalues_folder <- paste0('/grain/mk98/methyl/methylation-classification/data/TCGA/',cancer_project,'/DNA_Methylation/Methylation_Beta_Value/')
}
metadata_folder <- '/grain/mk98/methyl/methylation-classification/data/TCGA/'

# all_beta <- data.frame()
# all_meta <- data.frame()

# library(readr)
library(data.table)
barcodes<-c()
wanted_columns<-c("barcode", "sample_id", "primary_site", "name", "patient", "sample_type")

for (i in seq_along(cases)){
  if (i%%200==0){print(i)}
  case <- cases[i]
  
  case_query <- query$results[[1]][query$results[[1]]$cases==case, ]
  case_beta_file <- paste0(betavalues_folder, case_query$id,'/', case_query$file_name)
  
  # if (grepl("BRCA", cancer_project) && file.exists(paste0(metadata_folder, case_query$cases,'_sample_metadata.txt'))){
  #   case_meta <- read.table(paste0(metadata_folder, case_query$cases,'_sample_metadata.txt'), header = TRUE, sep='\t')
  # } else 
  if (case %in% data@colData$barcode){
    case_meta <- data@colData[data@colData$barcode==case,]
  } else {
    print(paste("The case ", case, " meta does not exist. Continuing to the next case_meta_file."))
    next
  }
  
  if (file.exists(case_beta_file)){
    case_beta <- read.table(case_beta_file, header = FALSE, row.names = 1)
  } else {
    print(paste("The case ", case, " beta does not exist. Continuing to the next case_meta_file."))
    next
  }
  if ("paper_Purity.ABSOLUTE.calls" %in% colnames(case_meta)){
    # print("purity information exists")
    case_meta$paper_Purity.ABSOLUTE.calls <- as.numeric(gsub(",", ".", case_meta$paper_Purity.ABSOLUTE.calls))
    wanted_columns <- wanted_columns + c("paper_Purity.ABSOLUTE.calls")
  } else {
    # print("purity information DOES NOT exist")
  }
  
  if (i==1){
    all_beta <- case_beta
    all_meta <- as.data.table(case_meta)
    barcodes <- c(case_meta$barcode)
  } else { 
    all_beta <- cbind(all_beta, case_beta)
    all_meta <- rbind(all_meta, as.data.table(case_meta))
    barcodes <- c(barcodes, case_meta$barcode)
  }
}

colnames(all_beta) <- barcodes
rownames(all_meta) <- all_meta$barcode
all_meta$primary_site <- sapply(all_meta$primary_site, function(x) x[[1]])
# all_meta$treatments <- sapply(all_meta$treatments, function(x) x[[1]])
all_meta <- all_meta[, !"treatments", with = FALSE]
all_meta$disease_type <- sapply(all_meta$disease_type, function(x) x[[1]])

write.table(all_beta, gzfile(paste0('data/TCGA/compiled/',cancer_project,'_beta_',num_samples,'.txt.gz')),
            quote = F,
            sep = '\t', row.names = T, col.names = T)
write.table(all_meta, gzfile(paste0('data/TCGA/compiled/',cancer_project,'_meta_',num_samples,'.txt.gz')),
            quote = F,
            sep = '\t', row.names = T, col.names = T)

# sample <- read.table(
#   paste0(betavalues_folder,'/03594190-9bbf-44c4-b519-e7dd1c2434eb/34b6d871-e421-4c0b-965e-97113d998c44.methylation_array.sesame.level3betas.txt'),
#   header = F, sep = '\t')
# 
# sample_metadata_folder <-  '/grain/mk98/methyl/methylation-classification/data/TCGA/'
