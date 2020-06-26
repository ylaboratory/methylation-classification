#!/bin/bash
# This script sets up the directory structure
# relative path is based on src directory, where the location of this script is
src_dir=$(dirname $0)
data_type=(data raw)
database=(GEO TCGA ENCODE)
assay_type=(microarray sequencing)
for i in "${data_type[@]}"; do
	mkdir -p $src_dir/../$i
		for j in "${database[@]}"; do
			mkdir -p $src_dir/../$i/$j
		done
done
for j in "${database[@]}"; do
	for z in "${assay_type[@]}"; do
		 mkdir -p $src_dir/../processed/$j/$z
	 done
done
mkdir -p  $src_dir/../annotation
	

