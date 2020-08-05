#!/bin/bash
src=$(dirname $0)
ignore_exist="F"
usage() {
	        echo "script usage: $(basename $0) [-s sample_name -g series_name -i ignore_exist] "
	}
while getopts ":s:g:i" Option; do
case "$Option" in
	s)      GSM=$OPTARG; echo "GEO sample $OPTARG is used" ;;
	g)      GSE=$OPTARG; echo "GEO series $OPTARG is used ";;
	i)	ignore_exist="T"; echo "ignoring existing dataset";;
	:)      echo "Error: -$OPTARG requires an argument"; usage; exit -1;;
	?)      echo "Error: unrecognized argument"; usage; exit -1;;
	*)      usage; exit -1;;
esac
done

echo "getting SRA accession corresponding to $GSM"
filename=".fastq"
SRA=( $(esearch -db sra -query $GSM | efetch -format docsum | xtract -pattern DocumentSummary -element Run@acc) )
echo "downloading fastq file for ${SRA[@]}"
for i in ${SRA[@]}
do
	if [ ! -f "$src_dir/../raw/GEO/$GSE/$i$filename" ] || [ "$ignore_exist" == "T" ]; then
		exit_status=1
		while [ ! $exit_status -eq 0 ]
		do
			fasterq-dump -f -O $src/../raw/GEO/$GSE  $i
			exit_status=$?
		done
	fi
done
echo "finished downloading!"
