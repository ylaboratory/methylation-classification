# Methylation-classification
This project uses whole-genome DNA methylation data to contruct a computational model for classifying human tissue/cell type or diseases.

## Installing R packages for preprocessing microarray data from GEO database
Set-up in RStudio container (in shell): 
```
$ podman exec -it [container-name or container-id] bash 
$ sudo apt-get update 
$ sudo apt-get install libbz2-dev 
$ sudo apt-get install -y liblzma-dev 
```
Set-up in RStudio (in R console):
```
>if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager") 
>BiocManager::install("minfi") 
>BiocManager::install("GEOquery") 
>BiocManager::install("wateRmelon") 
>BiocManager::install("IlluminaHumanMethylationEPICmanifest") 
>BiocManager::install("IlluminaHumanMethylation450kmanifest") 
```
## Generating genome liftover reference
### Download genomic information of 450k and 850k CpG site aligned to hg19
First, 450k or 850k (depending on the input data) manifest file that contains genomic information of microarray probe is needed in /local/usr directory for subsequent liftover. Run the main function of GEO-microarray.R, which takes in the series accession number (begins with GSE) as input. If the manifest file is not avilable, a prompt will be given and a text file containing the probe hg19 genomic information will be saved in /local/usr with file name 'hg19-[450/850]-coordinate.txt'

```
main.microarray.geo(accession_number)
```
### Liftover from hg19 to hg38
The hg19 coordinate text file needs to be uploaded to https://genome.ucsc.edu/cgi-bin/hgLiftOver for conversion to hg38, using options Minimum ratio of bases that must remap=0.95, original assembly: hg19, new assembly:hg38 \
The output files of converted loci (clicking on 'view conversion' tab and save in a directory) and deleted loci (clicking on 'Display faliure file' tab and save the information in a text file in the local/usr directory) should be downloaded. The deleted loci should be the same as 'hg19_450_deleted.txt' and 'hg19_850_deleted.txt' 
### Parse the converted loci
Then run liftover1 function in Parser-for-genome-after-liftover.R, which parse the downloaded file of converted loci and save it in the local/usr directory. 
```
liftover1(directory to converted loci file) 
```
The output is a text file containing the genomic location of CpG sites assayed on 450k or 850k platform lifted over to hg38 genome assembly. The output should be the same as 'hg38_450_converted_coordinate.txt' and 'hg38_850_converted_coordinate.txt' 

## Microarray data preprocessing into methylation beta values
Finally, to download and preprocess GEO data, run the main function of GEO-microarray.R, which takes in the series accession number (begins with GSE) as input. The raw data and processed methylation beta values are saved at /local/usr/GEO/microarray/accession_number 

```
main.microarray.geo(accession_number)
```

