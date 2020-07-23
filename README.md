# Methylation-classification
This project uses whole-genome DNA methylation data to contruct a computational model for classifying human tissue/cell type or diseases.
## Directory setup and preparation
First, run setup.sh script in shell to create the data, raw, processed and annotation directory. The data, raw and processed directories all have subdirectories of different databases. The script also downloads the annotation files needed in the /annotation directory (for now only hg19toHg38.chain file will be downloaded). The script then reads in the text files of datasets accession numbers from different databases that will be processed, and runs build-microarray.R file which does the processing of microarray data.

```
bash setup.sh
```
 
## Microarray data
### Installing R packages for preprocessing microarray data from GEO and ENCODE database (and potentially other databases)
The script build-microarray.R installs all R packages needed for microarray data preprocessing and handles data downloading, processing and outputting given the list of accession number of datasets. This script is called by build-dataset-microarray.R to run in command line. Three input arguments are needed (-i for whether existing datasets are ignored, -d for choosing which database to process data from and -m for the accessionlist/manifest file that contains datasets to be processed). Below are the setups needed to run the build-microarray.R script. The script also calls downloading, processing and outputting components, which are different scripts and will be explained in the following parts. After preprocessing all dataset inputs, the metadata and beta values for each database are integrated and written to a single file. 
Set-up in RStudio container (in shell): 
```
$ podman exec -it [container-name or container-id] bash 
$ sudo apt-get update 
$ sudo apt-get install libbz2-dev 
$ sudo apt-get install -y liblzma-dev 
$ sudo apt-get install curl
```
Component scripts for microarray data preprocessing (GEO and ENCODE)
```
/src/download-data-microarray.R
/src/preprocess-microarray.R
/src/liftover.R'
/src/output-data-microarray.R
```
### Microarray data downloading
All the functions for downloading raw data and metadata of GEO and ENCODE microarray datasets are housed under 'download-data-microarray.R'. 
### GEO 
First, the package GEOquery is used to download raw microarray data saved in the supplemental material from each GEO series accession. Run the download_data_geo_microarray function and input the accession number. The .tar file containing all raw .idat files corresponding to the GSE series will be stored in raw/GEO/accession.number \
Then raw idat files will be automatically extracted in the same directory as the tar file. \
The package GEOquery is also used to get the metadata of each series. There will be two metadata files for each GSE series, one contains series information (name, design, title, relation and supplementary files), while the other one contains sample specific information (sample name, assay type, platform code, source, title and database). To get the metadata, run download_geo_metadata function and input the accession number and the metadata will be saved as text files in data/GEO
### ENCODE
The R package ENCODEexplorer is used to download raw microarray data from each experiment (which often has one sample and possible replicates). Run the download_encode function and input the accession name (ENCS*). The .idat files will be downloaded to raw/ENCODE/accession.name by default. \
The package ENCODEexplorer is also used to download the metadata from each experiment. There will only be one metadata file for each experiment. Run get_metadata_encode function and input the accession name and the metadata will be saved as text file in data/ENCODE by default. \
For ENCODE data, the .idat files need to be renamed for passing into the preprocessing steps. Originally, ENCODE has a specific name (ENCF*) for each file (both red and green channel), which makes it hard to pair up the red and green channel data for each sample. Thus, we rename the files using the library name (ENCLB*), which is unique to one sample (similar to GSM in ENCODE). Run convert2target function and input the accession name to rename the files under that accession name.
### TCGA
The gdc-client is used to download TCGA microarray data using manifest file as input. The manifest files are donwloaded on the TCGA website. TCGAutils and TCGABiolinks packages are used to retrieve the metadata file given the barcode of each sample.
### Microarray data preprocessing
The minfi package is used to create the object for microarray preprocessing and the waterRmelon object is used for performing background correction and normalization. The following instruction is uniform for all microarray data regardless of database.\
First, to create the preprocessing raw object and perform background correction, run background_correction function and input the directory where the to-be-processed .idat files are stored. All .idat files will be read, and the output is a Methylset object after background correction with the reference probe.\
Next, run normalization function, which uses the BMIQ method to normalize the methylation beta values. The output is a matrix that contains the methylation beta values at each probe on hg19 coordinate.\
To convert the microarray data from hg19 assembly to hg38 assembly, package rtracklayer is used. Additionally, chain files of hg19 to hg38 from UCSC (https://hgdownload.soe.ucsc.edu/downloads.html#human) is used for reference. Run liftover.R and input the matrix of beta values on hg19 coordinate and the directory to the chain file, and the output would be a matrix of beta values on hg38 coordinate (some CpG sites are lost during mapping). 
### Microarray data output
Finally, the beta values for each sample is written to its own text file and saved in the right directory (GEO/ENCODE/..). Run write_to_file function and input the beta value matrix and the path to the metadata file corresponding to this batch of samples. The bata values for each series will be saved as "series-ID_beta_values.txt" in side data/database

## Sequencing data
### Install tools for sequnencing data processing 
The following tools need to be installed for running the RRBS and WGBS data processing scripts.
```
conda install -c bioconda entrez-direct
conda install -c bioconda bismark
conda install -c bioconda trim-galore
```
### Genome preparation
Before processing, configure-reference-genome.sh will download and prepare the reference genome for bisulfite sequencing alignment.
### Data downloading
#### GEO
For RRBS data processing, Entrez-direct is first used to find the SRR accessions that corresponds to a GEO sample and then download the raw sequence file. The donwload-files-GEO-sequencing.sh takes in command line inputs -i for whether to ignore existing datasets, -t for type of data to be processed (WGBS or RRBS), -d for database and -c for number of cores used for processing. The list of accessions in /annotation will be downloaded accordingly (filename: database-sequencing-accession-WGBS/RRBS.txt)
### Data processing
#### RRBS data
The process-geo-sequencing-RRBS.sh is used to process downloaded RRBS data. First, trim-galore is used to trim any adaptor sequence found and then bismark is used to align sequence to the reference genome and extract methylation percentage at each CpG sites. The outputs are written into /data folder and each sample results are grouped into a folder. 


