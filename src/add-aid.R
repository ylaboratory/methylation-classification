# Configuration
WORK_DIR <- '/grain/mk98/methyl/methylation-classification'
R_LIB_PATH <- "/home/mk98/R/x86_64-redhat-linux-gnu-library/4.0"
DATABASE_TYPE <- 'GEO'
META_FILE <- 'annotation/aid_meta.csv'
RAW_DATA_DIR <- 'raw/GEO'
log_file <- "logs/R_add_aid.txt"

required_packages <- list(
  cran = c("R.utils", "data.table", "reticulate", "progress", "remotes"),
  bioc = c("minfi", "rhdf5")
)

# Set up logging function
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

# Main execution
log_print("Setting up environment...")
setup_environment()
install_dependencies()
source("src/process-microarray.R")

# Load both GRSet
log_print("Loading GRSets...")

GRSet_rest <- readRDS("raw/GEO/training_GRSet_flt.rds")
GRSet_aid <- readRDS("raw/GEO/aid_GRSet_flt.rds")

log_print(paste0("GRSet_rest:", dim(GRSet_rest)))
log_print(paste0("GRSet_aid:", dim(GRSet_aid)))

keep_probes <- intersect(featureNames(GRSet_rest), featureNames(GRSet_aid))
log_print(paste0("probes to keep:", length(keep_probes)))

GRSet_rest_flt <- GRSet_rest[keep_probes,]
GRSet_aid_flt <- GRSet_aid[keep_probes,]

log_print(paste0("GRSet_rest_flt:", dim(GRSet_rest_flt)))
log_print(paste0("GRSet_aid_flt:", dim(GRSet_aid_flt)))

GRSet_combined <- combineArrays(GRSet_rest_flt, GRSet_aid_flt, "IlluminaHumanMethylation450k", FALSE)
log_print(paste0("GRSet_combined:", dim(GRSet_combined)))
log_print(paste0("GRSet_combined beta:", dim(getBeta(GRSet_combined))))
log_print(paste0("GRSet_combined M:", dim(getM(GRSet_combined))))

log_print("Saving GRSet_combined...")
saveRDS(GRSet_combined, file = "raw/GEO/combined_GRSet.rds")

training_meta <- read.csv('annotation/training_meta.csv')
aid_meta <- read.csv('annotation/aid_meta.csv')
aid_meta$X <-NULL

system_aids <- read.csv("annotation/system_aids.csv", col.names=c('from.name','from.id','to.id','to.name'), header=FALSE)
aid_meta <- merge(aid_meta, system_aids[, c("from.id", "to.id")], by.x = "UBERON.ID", by.y = "from.id", all.x = TRUE)
colnames(aid_meta)[which(names(aid_meta) == "to.id")] <- "training.ID"
aid_meta$training.ID <- ifelse(is.na(aid_meta$training.ID), aid_meta$merged.ID, aid_meta$training.ID)
aid_meta

anno <- read.csv("annotation/gsm_450k_annotations.csv")
meta_anno <- anno[match(aid_meta$Sample, anno$sample), ]
aid_meta <- merge(aid_meta, meta_anno, by.x = "Sample", by.y = "sample", all.x = TRUE)
aid_meta

library(dplyr)
wanted_cols <- c('Sample','Dataset','Annotated.tissue','UBERON.ID','UBERON.Name','Display.Name','merged.ID','training.ID','File','FileSeries')
common_cols <- intersect(wanted_cols, names(aid_meta))
concatenated_data <- bind_rows(training_meta[wanted_cols], aid_meta[common_cols])

write.csv(concatenated_data, 'annotation/training_meta.csv')

close(log)
