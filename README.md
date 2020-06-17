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
### Download genomic information of 450k and 850k CpG sites aligned to hg19
To convert the hg19 genome loci to hg38 for each platform, first download the Rdata "Illumina450kCpGsites.rds" and "Illumina850kCpGsites.rds" to a user defined directory. These data contain the genomic location information of CpG sites assayed on Illumina 450k/850k beadchip platform. The data are extracted from the 'minfi', 'IlluminaHumanMethylation450kmanifest' and 'IlluminaHumanMethylationEPICmanifest' package installed from Bioconductor. 
### Convert the genomic information to text file
Then run liftover1 function in Genome-annotation-liftover.R, which takes in the directory where "Illumina450kCpGsites.rds" and "Illumina850kCpGsites.rds" are downloaded as input 
``` 
 liftover1(directory.to.CpGdata) 
```
The output is a text file containing the genomic location of CpG sites assayed on either 450k or 850k platform stored in directory 'local/usr'. The output should be the same as 'hg19_450_coordinate.txt' and 'hg19_850_coordinate.txt' 
### Liftover from hg19 to hg38
This text file needs to be uploaded to https://genome.ucsc.edu/cgi-bin/hgLiftOver with options Minimum ratio of bases that must remap=0.95, original assembly: hg19, new assembly:hg38 \
The output files of converted loci (clicking on 'view conversion' tab and save in a directory) and deleted loci (clicking on 'Display faliure file' tab and save the information in a text file in the local/usr directory) should be downloaded. The deleted loci should be the same as 'hg19_450_deleted.txt' and 'hg19_850_deleted.txt' 
### Parse the converted loci
Then run liftover2 function in Genome-annotation-liftover.R, which parse the downloaded file of converted loci and save it in the local/usr directory. 
```
liftover2(directory to converted loci file) 
```
The output is a text file containing the genomic location of CpG sites assayed on 450k or 850k platform lifted over to hg38 genome assembly. The output should be the same as 'hg38_450_converted_coordinate.txt' and 'hg38_850_converted_coordinate.txt' 

## Microarray data preprocessing into methylation beta values
Finally, to download and preprocess GEO data, run the main function of GEO-microarray.R, which takes in the series accession number (begins with GSE) as input. The raw data and processed methylation beta values are saved at /local/usr/GEO/microarray/accession_number
main.microarray.geo(accession_number)
