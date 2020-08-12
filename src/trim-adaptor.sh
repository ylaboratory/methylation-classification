#!/bin/bash
src_dir=$(dirname $0)
read -p "enter the name and GSE number of raw sequencing file:" SRR GSE
filename='.fastq'
trim_galore --rrbs -j 6 -o $src_dir/../raw/GEO/$GSE  $src_dir/../raw/GEO/$GSE/$SRR$filename
