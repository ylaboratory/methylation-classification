# Configuration
WORK_DIR <- '/grain/mk98/methyl/methylation-classification'
R_LIB_PATH <- "/home/mk98/R/x86_64-redhat-linux-gnu-library/4.0"
DATABASE_TYPE <- 'GEO'
META_FILE <- "annotation/training_meta.csv"
RAW_DATA_DIR <- 'raw/GEO'

required_packages <- list(
  cran = c("R.utils", "data.table", "reticulate", "progress", "remotes", "arrow"),
  bioc = c("minfi", "rhdf5")
)

# Initialize environment
setup_environment <- function() {
  setwd(WORK_DIR)
  .libPaths(c(R_LIB_PATH, .libPaths()))
}

install_dependencies <- function() {
  suppressMessages({
    library(RColorBrewer)
    library(dplyr)
    if (packageVersion("matrixStats") != "1.1.0") {
      remove.packages("matrixStats")
      remotes::install_version("matrixStats", version = "1.1.0", repos = "http://cran.us.r-project.org", force = TRUE, quiet = TRUE)
    }
    cran_packages <- required_packages$cran
    for (pkg in cran_packages) {
      if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
        install.packages(pkg, repos = "http://cran.us.r-project.org", quiet = TRUE)
      }
    }
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager", repos = "http://cran.us.r-project.org", quiet = TRUE)
    }
    bioc_packages <- required_packages$bioc
    for (pkg in bioc_packages) {
      if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
        BiocManager::install(pkg, quiet = TRUE)
      }
    }
  })
}

# Main execution
setup_environment()
install_dependencies()

cat("loading GRSet\n")
GRSet <- readRDS(paste0("raw/GEO/combined_GRSet.rds"))
meta <- read.csv(META_FILE)

print(dim(getBeta(GRSet)))
print(dim(getM(GRSet)))
print(dim(meta))

# Convert the matrix to a data frame and add row names and column names
beta_df <- as.data.frame(getBeta(GRSet))
beta_df$probe <- rownames(beta_df)  # Add row names as a new column
print(rownames(beta_df)[1:5])

# Convert the matrix to a data frame and add row names and column names
m_df <- as.data.frame(getM(GRSet))
m_df$probe <- rownames(m_df)  # Add row names as a new column
print(rownames(m_df)[1:5])

# check if meta matches that of GRSet
print("Discrepancy between GRSet samples and metadata samples:")
print(length(setdiff(GRSet$Sample, meta$Sample)))

# cat("writing parquet")
write_parquet(beta_df, "data/GEO/preprocessed/training_betavalues.parquet")
write_parquet(m_df, "data/GEO/preprocessed/training_Mvalues.parquet")
