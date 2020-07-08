#!/bin/bash
# relative path is based on src directory, where the location of this script is
src_dir=$(dirname $0)

# read in the text file containing accession numbers of datasets to be processed
usage() {
	echo "script usage: $(basename $0) [-i flag for ignoring existing dataset when downloading -d database_name] "
}
ignore_exist="F"
while getopts ":id:" Option; do
	case "$Option" in
		i)      ignore_exist="T"; echo "ignoring existing dataset" ;;		
		d)	database=$OPTARG; echo "database $OPTARG is used";
			Rscript $src_dir/build-microarray.R $database $ignore_exist;;
		:)
			echo "Error: -$OPTARG requires an argument"; usage; exit -1;;
		?)
			echo "Error: unrecognized argument"; usage; exit -1;;
		*)	usage; exit -1;;
	esac
done


