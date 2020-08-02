#!/bin/bash
# read in the accession numbers of datasets to be processed
usage() {
        echo "script usage: $(basename $0) [-s sample_name -g series_name -c number of cores] "
}
while getopts ":s:g:c:" Option; do
case "$Option" in 
	s)      GSM=$OPTARG; echo "Align sequence $OPTARG" ;;
	g)      GSE=$OPTARG; echo "GEO series $OPTARG is used ";;
	c)      core_count=$OPTARG; echo "$OPTARG cores are used ";;
	:)	echo "Error: -$OPTARG requires an argument"; usage; exit -1;;
	?)	echo "Error: unrecognized argument"; usage; exit -1;;
	*)      usage; exit -1;;
esac
done
# relative path is based on src directory, where the location of this script is
src_dir=$(dirname $0)
# get SRR accession from GSE and GSM
SRR=( $(esearch -db sra -query $GSM | efetch -format docsum | xtract -pattern DocumentSummary -element Run@acc) )
echo ${SRR[@]}
for i in "${SRR[@]}"
do
	filename_align='.fastq'
	filename_extract=".bam"
	count=$( ls -1q $src_dir/../raw/GEO/$GSE/$i??$filename_align | wc -l )
	echo $count
	if  [ $count -gt 1 ]; then
		echo "paired"
		bismark -o $src_dir/../processed/GEO/$GSE --multicore $core_count --temp_dir $src_dir/../processed/GEO/$GSE --genome $src_dir/../annotation -1 $src_dir/../raw/GEO/$GSE/$i"_1"$filename_align -2 $src_dir/../raw/GEO/$GSE/$i"_2"$filename_align
		deduplicate_bismark -p --output_dir $src_dir/../processed/GEO/$GSE -o $i $src_dir/../processed/GEO/$GSE/$i*$filename_extract
	bismark_methylation_extractor -p --no_overlap --comprehensive --multicore $core_count --zero_based --bedGraph --cutoff 20 -o $src_dir/../data/GEO/$GSE $src_dir/../processed/GEO/$GSE/$i*".deduplicated.bam"
	else
		echo "single"
		bismark -o $src_dir/../processed/GEO/$GSE --multicore $core_count --temp_dir $src_dir/../processed/GEO/$GSE --genome $src_dir/../annotation $src_dir/../raw/GEO/$GSE/$i$filename_align
		deduplicate_bismark -s --output_dir $src_dir/../processed/GEO/$GSE -o $i $src_dir/../processed/GEO/$GSE/$i*$filename_extract
		bismark_methylation_extractor -s --comprehensive --multicore $core_count --zero_based --bedGraph --cutoff 20 -o $src_dir/../data/GEO/$GSE $src_dir/../processed/GEO/$GSE/$i*".deduplicated.bam" 
	fi
done
