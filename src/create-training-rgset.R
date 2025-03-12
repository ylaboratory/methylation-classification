# Configuration
WORK_DIR <- '/grain/mk98/methyl/methylation-classification'
R_LIB_PATH <- "/home/mk98/R/x86_64-redhat-linux-gnu-library/4.0"
DATABASE_TYPE <- 'GEO'
# META_FILE <- 'data/GEO/preprocessed/training_meta.csv'
# META_FILE <- "annotation/filtered_training_meta.csv"
META_FILE <- "annotation/combined_metadata.csv"
RAW_DATA_DIR <- 'raw/GEO'
MAKE_RGSET <- FALSE
MAKE_MSET <- FALSE
QC <- TRUE

required_packages <- list(
  cran = c("R.utils", "data.table", "reticulate", "progress", "remotes"),
  bioc = c("minfi", "rhdf5")
)

# Set up logging function
log_file <- "logs/R_training_process_log.txt"
log <- file(log_file, open = "wt")  # Open log file in write mode

# Log function to handle console and file output
log_print <- function(message, level = "INFO") {
  timestamp <- Sys.time()  # Add timestamp
  log_message <- paste0("[", timestamp, "] [", level, "] ", message)
  cat(log_message, "\n") #console
  cat(log_message, "\n", file = log, append = TRUE) #file
}

# Initialize environment
setup_environment <- function() {
  setwd(WORK_DIR)
  .libPaths(c(R_LIB_PATH, .libPaths()))
}

install_dependencies <- function() {
  log_print("Installing dependencies...")
  options(warn = -1)  # Suppress warnings
  
  suppressMessages({
    if (packageVersion("matrixStats") != "1.1.0") {
      remove.packages("matrixStats")  # Remove the package
      remotes::install_version("matrixStats", version = "1.1.0", repos = "http://cran.us.r-project.org", force = TRUE, quiet = TRUE)
    } else {
      log_print("matrixStats is already at version 1.1.0.")
    }
    
    # Install CRAN packages if needed
    cran_packages <- required_packages$cran
    for (pkg in cran_packages) {
      if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
        install.packages(pkg, repos = "http://cran.us.r-project.org", quiet = TRUE)
      }
    }
    
    # Install BiocManager if needed
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager", repos = "http://cran.us.r-project.org", quiet = TRUE)
    }
    
    # Install Bioconductor packages
    bioc_packages <- required_packages$bioc
    for (pkg in bioc_packages) {
      if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
        BiocManager::install(pkg, quiet = TRUE)
      }
    }
  })
}

# Efficiently filter matching files (optimized with pre-allocated list)
filter_files <- function(files, tissue_samples) {
  grn_files <- vector("list", length(tissue_samples$Sample))
  red_files <- vector("list", length(tissue_samples$Sample))
  
  for (sample in tissue_samples$Sample) {
    grn_file_pattern <- paste0(sample, ".*_(Grn|Gn)?\\.idat(\\.gz)?$")
    red_file_pattern <- paste0(sample, ".*_(Red|Rd)?\\.idat(\\.gz)?$")
    
    grn_files[[sample]] <- grep(grn_file_pattern, files, value = TRUE)
    red_files[[sample]] <- grep(red_file_pattern, files, value = TRUE)
  }
  
  return(list(grn_files = unlist(grn_files), red_files = unlist(red_files)))
}

# Main execution
log_print("Setting up environment...")
setup_environment()
install_dependencies()
source("src/process-microarray.R")

# Read metadata
log_print("Reading metadata...")
meta_data <- fread(META_FILE)

# Read id_to_name
log_print("Reading id_to_name...")
id_to_name_raw <- read.csv("src/id_to_name.csv", header = FALSE)
id_to_name <- as.data.frame(t(id_to_name_raw), stringsAsFactors = FALSE)
colnames(id_to_name) <- c("id", "name")
id_to_name_dict <- setNames(id_to_name$name, id_to_name$id)
name_to_id_dict <- setNames(id_to_name$id, id_to_name$name)

# Process 450k data using minfi
unique_tissue_ids <- unique(meta_data$`merged.ID`)

if (MAKE_RGSET){
  # Read raw files
  # log_print("Reading all files...")
  # all_files <- list.files(RAW_DATA_DIR, pattern = "\\.idat(\\.gz)?$", recursive = TRUE, full.names = TRUE)
  # file_match <- filter_files(all_files, meta_data)
  # 
  # grn_files <- file_match$grn_files
  # red_files <- file_match$red_files
  
  
  # log_print("Length of green and red files are equal, proceeding to processing...")
  # sample_files <- data.frame(
  #   # Basename = gsub("_(Grn|Gn)?\\.idat(\\.gz)?$", "", basename(grn_files)),
  #   Basename = basename(grn_files),
  #   stringsAsFactors = FALSE
  # )
  # sample_files$Basename <- all_files[match(sample_files$Basename, basename(all_files))]
  # sample_files <- sample_files[!duplicated(sample_files$Basename), , drop = FALSE]
  # sample_files$Basename <- gsub("_(Grn|Gn)?\\.idat(\\.gz)?$", "", sample_files$Basename)
  # sample_files$Sample <- sub("(_.*)?$", "", basename(sample_files$Basename))
  # sample_files$Series <- sub(".*/(GSE\\d+).*/.*", "\\1", sample_files$Basename)
  
  sample_files <- data.frame(
    Basename = meta_data$File,
    Series = meta_data$Dataset,
    Sample = meta_data$Sample
  )
  
  # Create sample_series data frame
  # sample_series <- data.frame(Sample = meta_data$Sample, Series = meta_data$Dataset, stringsAsFactors = FALSE)
  # sample_series <- sample_series[match(sample_files$Sample, sample_series$Sample), ]
  # 
  # sample_files$Basename <- file.path(RAW_DATA_DIR, sample_series$Series, sample_files$Basename)
  # sample_files <- sample_files[sample_files$Basename %in% all_basenames, ]
  
  log_print(paste("n sample files:", nrow(sample_files)))
  # 
  # missing_samples <- setdiff(meta_data$Sample, sample_files$Sample)
  # if (length(missing_samples) > 0) {
  #   log_print(paste0("There are ", length(missing_samples), " missing samples"))
  # }
  
  # Create RGSet
  log_print("Creating RGSet...")
  RGSet <- read.metharray.exp(targets = sample_files)
  
  log_print(" Saving RGSet...")
  saveRDS(RGSet, file = "raw/GEO/RGSet_training.rds")
  log_print(paste("n sample files:", nrow(sample_files)))
} else {
  log_print("Reading RGSet...")
  # RGSet_add <- readRDS("raw/GEO/training_RGSet_additional.rds")
  # RGSet <- readRDS("raw/GEO/RGSet_training.rds")
  # RGSet <- combineArrays(RGSet, RGSet_add, "IlluminaHumanMethylation450k", FALSE)
  # RGSet <- RGSet[,RGSet$Sample %in% meta_data$Sample]
  # log_print(paste(" Dimensions of RGSet: ",dim(RGSet)))
  # log_print(" Saving RGSet_added...")
  # saveRDS(RGSet, file = "raw/GEO/training_RGSet_added.rds")
  
  RGSet <- readRDS("raw/GEO/training_RGSet_added.rds")
  log_print(paste(" Dimensions of RGSet: ",dim(RGSet)))
}

if (MAKE_MSET){
  log_print("Creating MSet...")
  MSet_noob <- background_correction(RGSet)
  log_print(" Saving MSet...")
  saveRDS(MSet_noob, file = "raw/GEO/training_MSet.rds")
  
  qc <- getQC(MSet_noob)
  plotQC(qc)
  meds <- (qc$mMed + qc$uMed)/2
  whichBad <- which((meds < 10.5))
  failed_samples <- MSet_noob$Sample[whichBad]
  log_print(paste(" QC Failed samples:", length(failed_samples)))
  
  log_print(paste(" QC Failed samples saved at: data/GEO/compiled/R_training_qc_failed.csv"))
  write.csv(failed_samples, 'data/GEO/compiled/R_training_qc_failed.csv')
  
  log_print(paste(" QC Passed samples saved at: data/GEO/compiled/R_training_qc_passed.csv"))
  write.csv(MSet_noob$Sample[which((meds >= 10.5))], 'data/GEO/compiled/R_training_qc_passed.csv')
  
  passed_samples <- if(length(whichBad) > 0) rownames(qc[-whichBad, ]) else rownames(qc)
  MSet_qc <- MSet_noob[, passed_samples]
  RGSet_qc <- RGSet[, passed_samples]
  
  log_print(paste(" Dimensions of MSet_qc:",dim(MSet_qc)))
  log_print(" Saving MSet_qc...")
  saveRDS(MSet_qc, file = "raw/GEO/training_MSet_qc.rds")
} else {
  # log_print(paste(" hpv_positive samples filtered: GSM2319954,GSM2319956,GSM2319961,GSM2319969"))
  # samples_to_remove <- c("GSM2319954", "GSM2319956", "GSM2319961", "GSM2319969")
  # passed_samples <- read.csv('data/GEO/compiled/R_training_qc_passed.csv')$x
  # passed_samples <- setdiff(passed_samples, samples_to_remove)
  
  log_print("Reading MSet_qc...")
  MSet_qc <- readRDS("raw/GEO/training_MSet_qc.rds")
  log_print(paste(" Dimensions of MSet_qc:",dim(MSet_qc)))
  
  log_print("Reading traininig_meta...")
  training_meta <- read.csv('data/GEO/preprocessed/training_meta.csv')
  log_print(paste(" Dimensions of training_meta:",dim(training_meta)))
  
  MSet_qc <- MSet_qc[, match(training_meta$Sample, MSet_qc$Sample)]
  RGSet_qc <- RGSet[, match(training_meta$Sample, RGSet$Sample)]
  
  log_print(paste(" Dimensions of MSet_qc:",dim(MSet_qc)))
  log_print(paste(" Dimensions of RGSet_qc:",dim(RGSet_qc)))
}

if (QC){
  # Detection P-value
  log_print("Calculating Detection P-values...")
  detP <- detectionP(RGSet_qc)
  log_print(paste(" Dimensions of detP (initial):", dim(detP)))
  
  # Detection P-value barplot
  log_print(" Generating Detection P-values Barplot...")
  pdf(file.path(WORK_DIR, "logs/Detection_P_Values_Barplot.pdf"))
  barplot(colMeans(detP), las=2, cex.names=0.8, ylab="Mean detection p-values")
  abline(h=0.05, col="red")
  dev.off()
  
  # Match the probes' order between RGSet and GRSet_bmiq
  detP <- detP[match(MSet_qc@NAMES, rownames(detP)),]
  log_print(paste(" Dimensions of detP after matching:", dim(detP)))
  
  # Remove any probes with p-value >= 0.01
  keep <- rowSums(detP < 0.01) == ncol(RGSet_qc)
  log_print(paste(" Number of probes to keep (p-value < 0.01):", sum(keep)))
  
  # Subset the RGSet based on the probes to keep
  # RGSet_qc <- RGSet_qc[keep,]
  MSet_qc <- MSet_qc[keep,]
  # log_print(paste(" Dimensions of RGSet_qc:", dim(RGSet_qc)))
  log_print(paste(" Dimensions of MSet_qc:", dim(MSet_qc))) 
  write.csv(rownames(MSet_qc), "annotation/detP_probes.csv", row.names = FALSE)
} else {
  detP_probes <- read.csv("annotation/detP_probes.csv")$x
  MSet_qc <- MSet_qc[match(detP_probes, rownames(MSet_qc)), ]
  log_print(paste(" Dimensions of MSet_qc:",dim(MSet_qc)))
}

# Perform normalization
log_print("Normalizing...")
RSet <- normalization(MSet_qc)
log_print(" Saving RSet...")
log_print(paste(" Dimensions of RSet:", paste(dim(RSet), collapse = " x ")))
saveRDS(RSet, file = "raw/GEO/training_RSet.rds")

log_print("Creating GRSet...")
GRSet_bmiq <- mapToGenome(RSet)
log_print(" Saving GRSet_bmiq...")
log_print(paste(" Dimensions of GRSet_bmiq:", paste(dim(GRSet_bmiq), collapse = " x ")))
saveRDS(GRSet_bmiq, file = "raw/GEO/training_GRSet_bmiq.rds")

# Prepare BMIQ-adjusted GRSet
# GRSet_bmiq <- copy(GRset)
# library(dplyr)
# GRSet_bmiq$Beta <- beta_values
# log_print(" Saving GRSet_bmiq...")
# log_print(paste(" Dimensions of GRSet_bmiq: ", paste(dim(GRSet_bmiq), collapse = " x ")))
# saveRDS(GRSet_bmiq, file = "raw/GEO/training_GRSet_bmiq.rds")
# GRSet_bmiq <- readRDS("raw/GEO/training_GRSet_bmiq.rds")

# Filtering and plot post-filtering Beta and M values
log_print("Starting Filtering Process...")

# Drop SNP loci
GRSet_flt <- dropLociWithSnps(GRSet_bmiq)
log_print(paste(" Dimensions of GRSet_flt after SNP loci removal:", dim(GRSet_flt)))

# Exclude cross-reactive probes
xReactiveProbes <- read.csv(file="_annotation/cross_reactive_probes.csv", stringsAsFactors=FALSE, col.names = c('TargetID'))
keep <- !(featureNames(GRSet_flt) %in% xReactiveProbes$TargetID)
log_print(paste(" Number of probes after excluding cross-reactive probes:", sum(keep)))
GRSet_flt <- GRSet_flt[keep,]
log_print(paste(" Dimensions of GRSet_flt after chromosome removal:", dim(GRSet_flt)))

# Remove X and Y chromosomes
ann450k = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
keep <- !(featureNames(GRSet_flt) %in% ann450k$Name[ann450k$chr %in% c("chrX", "chrY")])
log_print(paste(" Number of probes after X/Y removal:", sum(keep)))
GRSet_flt <- GRSet_flt[keep,]
log_print(paste(" Dimensions of GRSet_flt after chromosome removal:", dim(GRSet_flt)))

# Log2 transformation of Beta values
assay(GRSet_flt, 'M') <- log2(getBeta(GRSet_flt) / (1 - getBeta(GRSet_flt)))
saveRDS(GRSet_flt, file = "raw/GEO/training_GRSet_flt.rds")

# # Dynamically adjust xlim based on the range of the data
# range_beta <- range(GRSet_flt$Beta, na.rm = TRUE)
# range_M <- range(GRSet_flt$M, na.rm = TRUE)
# log_print(paste(" Beta range:", range_beta))
# log_print(paste(" M range:", range_M))

# # Plot post-filtering Beta and M values
# log_print("Generating Post-Filtering Beta and M Density Plots...")
# pdf(file.path(WORK_DIR, "logs/Density_Plot_Post_Filtering.pdf"))
# par(mfrow=c(1,2))
# densityPlot(GRSet_flt$Beta, sampGroups = RGSet$Series, main = "Beta", xlim = range_beta, legend = FALSE)
# densityPlot(GRSet_flt$M, sampGroups = RGSet$Series, main = "M", xlim = range_M, legend = FALSE)
# legend("right", legend = levels(factor(GRSet_flt$Series)),
#        text.col=brewer.pal(8,"Dark2"))
# dev.off()
# 
# # Final summary plot
# log_print("Generating Final Summary Plots...")
# pdf(file.path(WORK_DIR, "logs/Final_Summary_Plots.pdf"))
# par(mfrow=c(1,2))
# densityPlot(GRSet_flt$Beta, sampGroups=RGSet$Series, main="Beta", legend=FALSE)
# densityPlot(GRSet_flt$M, sampGroups=RGSet$Series, main="M", legend=FALSE)
# legend("top", legend = levels(factor(GRSet_flt$Series)),
#        text.col=brewer.pal(8,"Dark2"))
# dev.off()

log_print("All data are generated and saved.")

# Close log file
close(log)