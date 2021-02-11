# Methylation-classification

This project uses whole-genome DNA methylation data to contruct a computational
model for classifying human tissue/cell type or diseases. Corresponding
elements of the pipeline download and process both microarray and sequencing
data for use.

## Directory setup and preparation

To setup the directory structure of data, annotation, raw and processed as well
as download chain file for data processing run the initial setup script:

```sh
bash src/initial_setup.sh
```

Several packages will also need to be installed to run the various command
line tools for data processing. All of the relevant packages are contained
in the `env.yml` file.

```sh
conda env create -f env.yml
conda activate methyl
```

## Download data

To get started downloading and preprocessing the relevant data run the build
all script: MAKE THIS SCRIPT!

```sh
bash src/build-all-sequencing.sh
```

THE FOLLOWING IS OLD AND WILL BE RESTRUCTURED...

## Microarray data

### Environment setup

To set-up the RStudio container (in shell) for running microarray
data processing scripts:

```sh
podman exec -it [container-name or container-id] bash
sudo apt-get update
sudo apt-get install libbz2-dev
sudo apt-get install -y liblzma-dev
sudo apt-get install curl
```

### Overall processing microarray data from GEO, ENCODE and TCGA database

To install R packages needed for processing microarray data, read in a
list of accession numbers, download and process the corresponding datasets
and output methylation beta values and metadata:

```sh
bash build-dataset-microarray.sh [-i] [-m manifest_file] [-d database]
```

* `-i`: Overide existing dataset
* `-d`: Database name (TCGA, GEO or ENCODE)
* `-m`: Text file containing the dataset accession number information

Component scripts for microarray data processing

```sh
download-data-microarray.R
preprocess-microarray.R
liftover.R
output-data-microarray.R
```

### Microarray data downloading

All the functions for downloading raw data and metadata of GEO and ENCODE
microarray datasets are housed under 'download-data-microarray.R'.

### GEO

The package GEOquery is used to download raw microarray data saved in the
supplemental material from each GEO series accession.

To download the .idat microarray data from GEO:

```sh
download_data_geo_microarray(accession_number, ignore_exist, download_directory)
```

* `accession number`: GEO dataset accession number (GSE*)
* `ignore_exist`: If ture, ignore existing dataset when downloading. Default: False
* `download_directory`: Specify directory to download the data. Default: `/raw/GEO`
There will be two metadata files for each GSE series, one contains series
information (name, design, title, relation and supplementary files), while the
other one contains sample specific information (sample name, assay type,
platform code, source, title and database). To get the metadata:

```sh
download_geo_metadata(accession_number, output_directory)
```

* `accession_number`: GEO dataset accession number (GSE*)
* `output_directory`: directory of output metadata. Default: /data/GEO

### ENCODE

The R package ENCODEexplorer is used to download raw microarray data from
each experiment (which often has one sample and possible replicates).
To download datasets from ENCODE:

```sh
download_encode(accession_name,  download_directory, ignore_exist)
```

* `accession_name`: ENCODE dataset accession name (ENCSR*)
* `ignore_exist`: If ture, ignore existing dataset when downloading. Default: False
* `download_directory`: Specify directory to download the data. Default: /raw/ENCODE

For ENCODE data, the .idat files need to be renamed for passing into the
processing steps. Originally, ENCODE has a specific name (ENCF*) for each file
(both red and green channel), which makes it hard to pair up the red and green
channel data for each sample. Thus, we rename the files using the library name (ENCLB*),
which is unique to one sample (similar to GSM in ENCODE). To rename the files:

```sh
convert2target(accession_name, download_directory)
```

* `accession_name`: ENCODE dataset accession name (ENCSR*) to be renamed
* `download_directory`: Specify directory that the data were downloaded to.
Default: /raw/ENCODE

To get the metadata of ENCODE dataset:

```sh
get_metadata_encode(accession_name, out_directory)
```

* `accession_name`: ENCODE dataset accession name (ENCSR*)
* `out_directory`: Directory to output metadata. Default: /data/ENCODE

### TCGA

The gdc-client is used to download TCGA microarray data using manifest file
as input. The manifest files are donwloaded from the TCGA website.
TCGAutils and TCGABiolinks packages are used to retrieve the metadata file given
the barcode of each sample. To downloaded the TCGA dataset:

```sh
download-data-TCGA-microarray.sh [-m manifest_file]
```

* `-m`: The manifest text file that has selected datasets
information from the TCGA website

### Microarray data preprocessing

The minfi package is used to create the object for microarray preprocessing and
the waterRmelon object is used for performing background correction and
normalization. All functions are in process-microarray.R.
To perform background correction:

```sh
GRset_noob<-background_correction(accession_datadir)
```

* `accession_datadir`: The directory to read .idat data to be processed from.
All files would be read.
* `GRset_noob`: Methylation beta values at each probe site after single
color background correction
To perform BMIQ normalization:

```sh
BMIQ_genome_loci<-normalization(GRset_noob,dir2metadata)
```

* `GRset_noob`: Methylation beta values at each probe site after
single color background correction
* `dir2metadata`: directory to which sample metadata of the dataset is stored
* `BMIQ_genome_loci`: Methylation beta values at each CpG
site (genome loci) aligned using hg19 genome
To lift the CpG sites from hg19 to hg38 genome assembly:

```sh
BMIQ_lift<-liftover(dir2chain, BMIQ_genome_loci)
```

* `dir2chain`: The directory to hg19 to hg38 chain files
* `BMIQ_genome_loci`: Methylation beta values at each CpG site
(genome loci) aligned using hg19 genome
* `BMIQ_lift`: Methylation beta values at each CpG site
(genome loci) aligned using hg38 genome

### Microarray data output

Finally, the beta values for each sample are written to the
right database directory (GEO/ENCODE/..).

## Sequencing data

### Install tools for sequnencing data processing

To install tools for running the RRBS and WGBS data processing scripts:

```sh
conda install -c bioconda sra-tools=2.10
conda install -c bioconda entrez-direct
conda install -c bioconda bismark
conda install -c bioconda trim-galore
conda install -c bioconda samtools openssl=1.0
```

### Overall processing of sequencing data

To process a list of sequencing data whose accessions are in /annotation/file database_sequencing_accession_(RRBS/WGBS).txt:

```sh
build-dataset-sequencing.sh [-i] [-d database] [-t datatype] [-c number of cores]
```

* `-i`: If true, ignore exisiting dataset. Default: False
* `-d`: Name of database that the datasets are from (GEO etc.)
* `-t`: Type of sequencing data (RRBS or WGBS)
* `-c`: Number of cores to use while processing

### Get metadata

To obtain the metadata of GEO sample and series:

```sh
get-metadata-geo-sequencing.R [accession_number]
```

* `accession_number`: Accession number of the GSE series

### Genome preparation

Before processing, to download and prepare the reference genome
for bisulfite sequencing alignment:

```sh
configure-reference-genome.sh
```

This step is only done once.

### Data downloading GEO

For RRBS data processing, Entrez-direct is first used to find the SRR accessions
that corresponds to a GEO sample and then download the raw sequence
file. To download the sequences:

```sh
download-files-GEO-sequencing.sh [-i] [-g series] [-s sample]
```

* `-i`: If true, ignore exisiting dataset. Default: False
* `-g`: Series accession of the dataset (GSE*)
* `-s`: Sample accession of the dataset (GSM*)

### Data processing RRBS data

To trim adaptor sequence, align to reference genome and extract methylation
percentage at each genome loci:

```sh
processing-geo-sequencing-RRBS.sh [-g series] [-s sample] [-c number of cores]
```

* `-g`: Series accession of the dataset (GSE*)
* `-s`: Sample accession of the dataset (GSM*)
* `-c`: Number of cores used during alignment and methylation percentage extraction
