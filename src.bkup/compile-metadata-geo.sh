#!/bin/bash
usage () {
		echo "script usage: $(basename $0) [ -f GSE file path]"
}
while getopts ":f:" Option; do
	case "$Option" in
		f) file=$OPTARG; echo "getting metadata for GEO series from list under $OPTARG" ;;
	esac
done
src_dir=$(dirname $0)
while IFS= read -r -u3 line
do
	echo $line
	GSE=$line
	GSE_prefix=${GSE%???}
	echo $GSE_prefix
	# get the matrix series file for the GSE from the ftp site
	wget -O $src_dir/../raw/GEO/"$GSE""_series_matrix.txt.gz" ftp://ftp.ncbi.nlm.nih.gov/geo/series/$GSE_prefix"nnn"/$GSE/matrix/"$GSE"*"_series_matrix.txt.gz"
	# extract the metadata required
	gunzip $src_dir/../raw/GEO/"$GSE""_series_matrix.txt.gz"
done 3< $src_dir/../annotation/$file
python extract-metadata-geo.py $src_dir/../annotation/$file

