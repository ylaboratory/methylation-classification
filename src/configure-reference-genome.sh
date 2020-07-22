#!/bin/bash
src=$(dirname $0)
#download the reference genome 
if [ ! -f "$src/../annotation/hg38.fa" ]; then
	wget -O $src/../annotation/hg38.fa.gz ftp://ftp.ensembl.org/pub/release-100/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
	gunzip $src/../annotation/hg38.fa.gz
fi
#wget -O $src/../annotation/hg38.fa.gz ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
echo "preparing reference genome for bisulfite treatment in directory $src/../annotation"
bismark_genome_preparation --verbose $src/../annotation
