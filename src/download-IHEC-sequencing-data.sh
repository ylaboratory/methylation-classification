#!/bin/bash
# relative path is based on src directory, where the location of this script is
src_dir=$(dirname $0)
# read in the text file containing accession numbers of datasets to be processed
usage() {
		echo "script usage: $(basename $0) [-m manifest file name -d database_name -s sample_list file] "
}
while getopts ":m:d:s:" Option; do
	case "$Option" in
		m) manifest=$OPTARG; echo " using manifest file $manifest";;
		d) database=$OPTARG; echo "downloading data from $database";;
		s) sample=$OPTARG; echo "using sample list $sample";;
		:) echo "Error: -$OPTARG requires an argument"; usage; exit -1;;
		?) echo "Error: unrecognized argument"; usage; exit -1;;
		*) usage; exit -1;;
	esac
done
if [ ! -f $src_dir/bin/bigWigToBedGraph ];then
	wget -P $src_dir/bin http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigToBedGraph
	chmod 750 $src_dir/bin/bigWigToBedGraph
fi
mkdir -p $src_dir/../data/$database
#wget -P $src_dir/../data/$database -i $src_dir/../annotation/$manifest
Rscript rename-and-get-metadata-IHEC.R $manifest $sample $src_dir $database
for file in "$src_dir/../data/$database"/*;do
	if [[ $file == *_beta_values.bigWig ]];then
		base_name="bedGraph"
		$src_dir/bin/bigWigToBedGraph $file ${file/bigWig/$base_name}
	fi
done
