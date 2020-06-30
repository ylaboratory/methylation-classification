# Methylation-classification
This project uses whole-genome DNA methylation data to contruct a computational model for classifying human tissue/cell type or diseases.

## Setup directory
First, run directory-setup script in shell to create the data, raw, processed and annotation directory. The data, raw and processed directories all have subdirectories of different databases.
```
bash directory-setup.sh
```
## Installing R packages for preprocessing microarray data from GEO and ENCODE database (and potentially other databases)
Set-up in RStudio container (in shell): 
```
$ podman exec -it [container-name or container-id] bash 
$ sudo apt-get update 
$ sudo apt-get install libbz2-dev 
$ sudo apt-get install -y liblzma-dev 
$ sudo apt-get install curl
```
Set-up in RStudio (in R console):
```
>install.packages("RPMM")
>if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager") 
>BiocManager::install("minfi") 
>BiocManager::install("GEOquery") 
>BiocManager::install("wateRmelon") 
>BiocManager::install("IlluminaHumanMethylationEPICmanifest") 
>BiocManager::install("IlluminaHumanMethylation450kmanifest") 
>BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
>BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")
>BiocManager::install("ENCODExplorer")
>BiocManager::install("rtracklayer")
```
## Loading all functions for microarray data preprocessing (GEO and ENCODE)
```
source('./src/download-data-microarray.R')
source('./src/preprocess-microarray.R')
source('./src/liftover.R')
source('./src/output-data-microarray.R')
```
## Downloading data for microarray
All the functions for downloading raw data and metadata of GEO and ENCODE microarray datasets are housed under 'download-data-microarray.R'. 
### Download microarray data from GEO 
First, the package GEOquery is used to download raw microarray data saved in the supplemental material from each GEO series accession. Run the download_data_geo_microarray function and input the accession number. The .tar file containing all raw .idat files corresponding to the GSE series will be stored in ./raw/GEO/accession.number \
Then raw idat files will be automatically extracted in the same directory as the tar file. 
### Download microarray data from ENCODE
The R package ENCODEexplorer is used to download raw microarray data from each experiment (which often has one sample and possible replicates). Run the download_encode function and input the accession name (ENCS...). The .idat files will be downloaded to ./raw/ENCODE/accession.name by default.
## Getting metadata for microarray
### Get metadata from GEO
The package GEOquery is also used to get the metadata of each series. There will be two metadata files for each GSE series, one contains series information (name, design, title, relation and supplementary files), while the other one contains sample specific information (sample name, assay type, platform code, source, title and database). To get the metadata, run download_geo_metadat.R and input the accession number and the metadata will be saved as text files in ./data/GEO
### Get metadata from ENCODE
The package ENCODEexplorer is also used to download the metadata from each experiment. There will only be one metadata file for each experiment. Run get_metadata_encode function and input the accession name and the metadata will be saved as text file in ./data/ENCODE by default.
## Get ready for microarray data preprocessing
### GEO
No additional step is needed. Can proceed directly to next step
### ENCODE
The .idat files need to be renamed for passing into the preprocessing steps. Originally, ENCODE has a specific name (ENCF...) for each file (both red and green channel), which makes it hard to pair up the red and green channel data for each sample. Thus, we rename the files using the library name (ENCLB...), which is unique to one sample (similar to GSM in ENCODE). Run convert2target function and input the accession name to rename the files under that accession name. 
## Preprocessing microarray data
The minfi package is used to create the object for microarray preprocessing and the waterRmelon object is used for performing background correction and normalization. The following instruction is uniform for all microarray data regardless of database.\
First, to create the preprocessing raw object and perform background correction, run background_correction function and input the directory where the to-be-processed .idat files are stored. All .idat files will be read, and the output is a Methylset object after background correction with the reference probe.\
Next, run normalization function, which uses the BMIQ method to normalize the methylation beta values. The output is a matrix that contains the methylation beta values at each probe on hg19 coordinate.\
## Liftover
To convert the microarray data from hg19 assembly to hg38 assembly, package rtracklayer is used. Additionally, chain files of hg19 to hg38 from UCSC (https://hgdownload.soe.ucsc.edu/downloads.html#human) should be downloaded and saved (optimally in ./annotation). The script download-chain.sh can be used to automatically download the chain file to the target location.
```
bash download-chain-file.sh
```
Run liftover.R and input the matrix of beta values on hg19 coordinate and the directory to the chain file, and the output would be a matrix of beta values on hg38 coordinate (some CpG sites are lost during mapping). 
## Save beta values from microarray data
Finally, the beta values for each sample is written to its own text file and saved in the right directory (GEO/ENCODE..). Run write_to_file fucntion and input the beta value matrix and the path to the metadata file corresponding to this batch of samples. The bata values for each series will be saved as "series-ID_beta_values.txt" in side ./data/database



