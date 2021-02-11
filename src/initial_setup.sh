#!/bin/bash
# this script sets up the directory structure and other necessary preliminaries
src_dir=$(dirname "$0")
data_type=(data raw processed)
database=(GEO TCGA ENCODE)
for i in "${data_type[@]}"; do
	        for j in "${database[@]}"; do
			                mkdir -p -v "$src_dir/../$i/$j"
					        done
					done
					mkdir -p -v "$src_dir"/../annotation

# download hg19 to hg38 chain files for liftover
annotation=$src_dir/../annotation/hg19ToHg38.over.chain.gz
if [[ ! -f "$annotation" ]]
then
	wget 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz' -O "$annotation"
	gunzip "$annotation"
fi
