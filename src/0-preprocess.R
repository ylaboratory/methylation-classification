# Configuration
WORK_DIR <- dirname(dirname(normalizePath(commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs())])))
META_FILE <- file.path(WORK_DIR, "annotation/metadata.csv") # Metadata from GEO listing filename ($File), GSE ($Dataset), and GSM ($Sample) of wanted samples
RAW_DATA_DIR <- file.path(WORK_DIR, "raw/GEO") # Path to raw idat files
OUTPUT_DIR <- file.path(WORK_DIR, "preprocessed/")
LOG_FILE <- file.path(WORK_DIR, "logs/preprocessing_log.txt")
REACTIVE_PROBE_FILE <- file.path(WORK_DIR, "annotation/cross_reactive_probes.csv")

required_packages <- list(
  cran = c("R.utils", "data.table", "reticulate", "progress", "remotes"),
  bioc = c("minfi", "rhdf5")
)

# Logging setup
log <- file(LOG_FILE, open = "wt")

log_print <- function(message, level = "INFO") {
  timestamp <- Sys.time()
  log_message <- paste0("[", timestamp, "] [", level, "] ", message)
  cat(log_message, "\n")
  cat(log_message, "\n", file = log, append = TRUE)
}

# Environment setup
setup_environment <- function() {
  setwd(WORK_DIR)
  if (!dir.exists(OUTPUT_DIR)) {
    dir.create(OUTPUT_DIR, recursive = TRUE)
  }
}

# Install dependencies
install_dependencies <- function() {
  log_print("Installing dependencies...")
  suppressMessages({
    if (packageVersion("matrixStats") != "1.1.0") {
      remove.packages("matrixStats")
      remotes::install_version("matrixStats", version = "1.1.0", repos = "http://cran.us.r-project.org", force = TRUE, quiet = TRUE)
    }
    for (pkg in required_packages$cran) {
      if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
        install.packages(pkg, repos = "http://cran.us.r-project.org", quiet = TRUE)
      }
    }
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager", repos = "http://cran.us.r-project.org", quiet = TRUE)
    }
    for (pkg in required_packages$bioc) {
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

log_print("Reading metadata...")
meta_data <- fread(META_FILE)

# Process RGSet
log_print("Creating RGSet...")
sample_files <- data.frame(
  Basename = meta_data$File,
  Series = meta_data$Dataset,
  Sample = meta_data$Sample
)
RGSet <- read.metharray.exp(targets = sample_files)
saveRDS(RGSet, file = file.path(OUTPUT_DIR, "RGSet.rds"))
log_print("RGSet saved.")

# Process MSet
log_print("Background correction...")
MSet_noob <- background_correction(RGSet)
saveRDS(MSet_noob, file = file.path(OUTPUT_DIR, "MSet.rds"))
log_print("MSet saved.")

# Quality control
log_print("Performing quality control...")
qc <- getQC(MSet_noob)
plotQC(qc)
meds <- (qc$mMed + qc$uMed)/2
whichBad <- which((meds < 10.5))
failed_samples <- MSet_noob$Sample[whichBad]
log_print(paste("QC Failed samples:", length(failed_samples)))

log_print(paste("QC Failed samples saved as qc_failed.csv"))
write.csv(failed_samples, paste0(OUTPUT_DIR, "qc_failed.csv"))

log_print(paste("QC Passed samples saved as qc_passed.csv"))
write.csv(MSet_noob$Sample[which((meds >= 10.5))], paste0(OUTPUT_DIR, "qc_passed.csv"))

passed_samples <- if(length(whichBad) > 0) rownames(qc[-whichBad, ]) else rownames(qc)
MSet_qc <- MSet_noob[, passed_samples]
RGSet_qc <- RGSet[, passed_samples]

detP <- detectionP(RGSet_qc)
keep <- rowSums(detP < 0.01) == ncol(RGSet_qc)
MSet_qc <- MSet_qc[keep, ]
saveRDS(MSet_qc, file = file.path(OUTPUT_DIR, "MSet_qc.rds"))
log_print("QC completed and MSet_qc saved.")

# Normalization
log_print("Normalizing data...")
RSet <- normalization(MSet_qc)
saveRDS(RSet, file = file.path(OUTPUT_DIR, "RSet.rds"))
log_print("RSet saved.")

# Genome mapping
log_print("Mapping to genome...")
GRSet_bmiq <- mapToGenome(RSet)
saveRDS(GRSet_bmiq, file = file.path(OUTPUT_DIR, "GRSet_bmiq.rds"))
log_print("GRSet_bmiq saved.")

# Filtering SNPs
log_print("Filtering probes...")
GRSet_flt <- dropLociWithSnps(GRSet_bmiq)

# Filtering cross reactive probes
xReactiveProbes <- read.csv(REACTIVE_PROBE_FILE, stringsAsFactors = FALSE, col.names = c('TargetID'))
keep <- !(featureNames(GRSet_flt) %in% xReactiveProbes$TargetID)
GRSet_flt <- GRSet_flt[keep, ]

# Filtering X, Y probes
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
keep <- !(featureNames(GRSet_flt) %in% ann450k$Name[ann450k$chr %in% c("chrX", "chrY")])
GRSet_flt <- GRSet_flt[keep, ]

# Convert to M-values
assay(GRSet_flt, 'M') <- log2(getBeta(GRSet_flt) / (1 - getBeta(GRSet_flt)))
saveRDS(GRSet_flt, file = file.path(OUTPUT_DIR, "GRSet_flt.rds"))
log_print("GRSet_flt saved.")

# Convert the matrix to a data frame and add row names and column names
beta_df <- as.data.frame(getBeta(GRSet_flt))
beta_df$probe <- rownames(beta_df)  # Add row names as a new column
print(rownames(beta_df)[1:5])

# Convert the matrix to a data frame and add row names and column names
m_df <- as.data.frame(getM(GRSet_flt))
m_df$probe <- rownames(m_df)  # Add row names as a new column
print(rownames(m_df)[1:5])

log_print("writing parquet")
write_parquet(beta_df, paste0(OUTPUT_DIR,"betavalues.parquet"))
write_parquet(m_df, paste0(OUTPUT_DIR,"Mvalues.parquet"))

# Finalize
log_print("All data preprocessing completed.")
close(log)