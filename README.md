# Methylation-classification
This project uses whole-genome DNA methylation data to contruct a computational model for classifying human tissue/cell type or diseases.

 
## Microarray data
### Directory setup and preparation
First, run setup.sh script in shell to create the data, raw, processed and annotation directory. The data, raw and processed directories all have subdirectories of different databases. The script also downloads the annotation files needed in the /annotation directory (for now only hg19toHg38.chain file will be downloaded). The script then reads in the text files of datasets accession numbers from different databases that will be processed, and runs build-microarray.R file which does the processing of microarray data.

```
bash setup.sh
```
### Installing R packages for preprocessing microarray data from GEO and ENCODE database (and potentially other databases)
The script build-microarray.R installs all R packages needed for microarray data preprocessing and handles data downloading, processing and outputting given the accession number of datasets. Below are the setups needed to run the build-microarray.R script. The script also calls downloading, processing and outputting components, which are different scripts and will be explained in the following parts. After preprocessing all dataset inputs, the metadata and beta values for each database are integrated and written to a single file. 
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
 
### Microarray data preprocessing
The minfi package is used to create the object for microarray preprocessing and the waterRmelon object is used for performing background correction and normalization. The following instruction is uniform for all microarray data regardless of database.\
First, to create the preprocessing raw object and perform background correction, run background_correction function and input the directory where the to-be-processed .idat files are stored. All .idat files will be read, and the output is a Methylset object after background correction with the reference probe.\
Next, run normalization function, which uses the BMIQ method to normalize the methylation beta values. The output is a matrix that contains the methylation beta values at each probe on hg19 coordinate.\
To convert the microarray data from hg19 assembly to hg38 assembly, package rtracklayer is used. Additionally, chain files of hg19 to hg38 from UCSC (https://hgdownload.soe.ucsc.edu/downloads.html#human) is used for reference. Run liftover.R and input the matrix of beta values on hg19 coordinate and the directory to the chain file, and the output would be a matrix of beta values on hg38 coordinate (some CpG sites are lost during mapping). 
### Microarray data output
Finally, the beta values for each sample is written to its own text file and saved in the right directory (GEO/ENCODE/..). Run write_to_file function and input the beta value matrix and the path to the metadata file corresponding to this batch of samples. The bata values for each series will be saved as "series-ID_beta_values.txt" in side data/database



