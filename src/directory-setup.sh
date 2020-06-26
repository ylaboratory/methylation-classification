#!/bin/bash
# This script sets up the directory structure
# relative path is based on src directory, where the location of this script is
src_dir=$(dirname $0)
data_type=(data raw processed)
database=(GEO TCGA ENCODE)
for i in "${data_type[@]}"; do
	mkdir -p $src_dir/../$i
		for j in "${database[@]}"; do
			mkdir -p $src_dir/../$i/$j
		done
done
mkdir -p  $src_dir/../annotation
	

