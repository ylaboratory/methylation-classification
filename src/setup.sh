#!/bin/bash
# This script sets up the directory structure
# relative path is based on src directory, where the location of this script is
src_dir=$(dirname $0)
data_type=(data raw processed)
database=(GEO TCGA ENCODE)
for i in "${data_type[@]}"; do
	for j in "${database[@]}"; do
		mkdir -p -v $src_dir/../$i/$j
	done
done
mkdir -p -v $src_dir/../annotation

# download hg19 to hg38 chain files for liftover
wget 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz' -O $src_dir/../annotation/hg19ToHg38.over.chain.gz
gunzip $src_dir/../annotation/hg19ToHg38.over.chain.gz

# read in the text file containing accession numbers of datasets to be processed
usage() {
	echo "script usage: $(basename $0) [-d database_name] "
}
ignore_exist="F"
filename='_microarray_accession.txt'
while getopts ":id:" Option; do
	case "$Option" in
		i)      ignore_exist="T";                                                                                                    	echo "ignoring existing dataset" ;;		
		d)	database=$OPTARG; echo "database $OPTARG is used";
			input="$src_dir/../annotation/$database$filename";
			echo $ignore_exist;
			while read -r line;do
				 Rscript $src_dir/build-microarray.R $database $line $ignore_exist
				 echo "processing done for $line"
			done < $input
					;;
		i)	ignore_exist=$OPTARG
			echo "ignoring existing dataset" ;;
		:)
			echo "Error: -$OPTARG requires an argument"; usage; exit -1;;
		?)
			echo "Error: unrecognized argument"; usage; exit -1;;
		*)	usage; exit -1;;
	esac
done

