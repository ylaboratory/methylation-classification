#!/bin/bash
ignore_exist_flag=""
usage () {
	echo "script usage: $(basename $0) [-i flag for ignoring existing dataset when downloading -d database_name -c number of cores used -t sequencing data type]"
}
while getopts ":id:t:c:" Option; do
	case "$Option" in
		i)      ignore_exist_flag="-i";;
		t)      datatype=$OPTARG; echo "processing $OPTARG data" ;;
		d)	database=$OPTARG; echo "using $OPTARG database";;
		c)	core=$OPTARG; echo "using $OPTARG cores";;
		:)	echo "Error: -$OPTARG requires an argument"; usage; exit -1;;
		?)	echo "Error: unrecognized argument"; usage; exit -1;;
		*)	usage; exit -1;;
	esac
done
src_dir=$(dirname $0)
# Prepare reference genome if there is non existing
if [ -f $src_dir/../annotation/Bisulfite_Genome ]; then
	bash configure-reference-genome.sh
fi
accession_filename1="_sequencing_accession_"
accession_filename2=".txt"
while IFS= read -r -u3 line
do
	accession=( $line )
	GSE="${accession[0]}"
	GSM="${accession[1]}"
# Get metadata of the GEO series
	if [ ! -f $src_dir/../raw/GEO/$GSE ];then
		Rscript $src_dir/get-metadata-geo-sequencing.R $GSE
	fi
# Download SRR files corresponding to GSM accession
#	bash $src_dir/download-files-GEO-sequencing.sh -g $GSE -s $GSM $ignore_exist_flag
# Process the data	
#	if [ $database == "GEO" ] && [ $datatype == 'RRBS' ]; then
#	       echo "$GSM"	
#		bash $src_dir/download-files-GEO-sequencing.sh -g $GSE -s $GSM $ignore_exist_flag
 #       	bash $src_dir/processing-geo-sequencing-RRBS.sh -g $GSE -s $GSM -c $core
#	elif [ $database == "GEO" ] && [ $datatype == 'WGBS' ]; then
#		echo "$GSM"
#		bash $src_dir/download-files-GEO-sequencing.sh -g $GSE -s $GSM $ignore_exist_flag
#		bash $src_dir/processing-geo-sequencing-WGBS.sh -g $GSE -s $GSM -c $core
#	fi
	echo "$GSE"	
done 3< $src_dir/../annotation/$database$accession_filename1$datatype$accession_filename2
