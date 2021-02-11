#!/bin/bash
ignore_exist_flag=false
accession_filename1="_sequencing_accession_"
accession_filename2=".txt"
datatype="WGBS"
dataset="geo"
src_dir=$(dirname "$0")
base_dir=$(cd "$src_dir"/../ || exit; pwd)
core=10
current_reference_url="ftp://ftp.ensembl.org/pub/release-100/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
all_ids_geo=""
all_ids_tcga=""
all_ids_encode=""

usage() {
  echo "
  **************************************************************************************
  $(basename $0)

  Master script for downloading and processing WGBS and RRBS
  methylation sequencing data.

  -h, --help              prints instructions
  -d <geo, encode>        build different datasets given one of the databases
                          [default: geo]
  -t <WGBS, RRBS>         sequencing data type [default: WGBS]

  Optional: 
  -c <integer>            number of cores to use [default: 10]
  -e                      flag when given will ignore existing datasets when downloading
  -i <file>               input file with list of accessions / ids to download
                          if this file is not provided this script will download 
                          all for a given database
  **************************************************************************************           
  "
}
while getopts ":d:t:c:i:e" option; do
	case "${option}" in
		e) ignore_exist_flag=true;;
		t) datatype=${OPTARG}; echo "processing ${OPTARG} data" ;;
		d) dataset=${OPTARG}; echo "using ${OPTARG} dataset";;
		c) core=${OPTARG}; echo "using ${OPTARG} cores";;
		i) input_file=${OPTARG};;
		:) echo "Error: -${OPTARG} requires an argument"; usage; exit 0;;
		*) usage; exit 0;;
	esac
done

# if dataset flag used check for valid entry
if [[ ! $dataset =~ ^(tcga|geo|encode)$ ]]; then
  echo "Error invalid argument."
  usage
  exit 0
fi

# if datatype flag used check for valid entry
if [[ ! $datatype =~ ^(WGBS|RRBS)$ ]]; then
  echo "Error invalid argument."
  usage
  exit 0
fi

# download GEO sequencing data 




# prepare reference genome if doesn't exist
if [ ! -d "$base_dir"/annotation/Bisulfite_Genome ]; then
	echo "[$(date +"%Y-%m-%d %T")] downloading reference bisulfite genome"
	if [ ! -f "$base_dir/annotation/hg38.fa" ]; then
		wget -O "$base_dir/annotation/hg38.fa.gz" "$current_reference_url"
		echo "[$(date +"%Y-%m-%d %T")] extracting reference"
		gunzip "$base_dir/annotation/hg38.fa.gz"
	fi
	echo "[$(date +"%Y-%m-%d %T")] preping ref genome in $base_dir/annotation"
	bismark_genome_preparation --verbose "$base_dir/annotation"
fi



# while IFS= read -r -u3 line
# do
# 	accession=( "$line" )
# 	GSE="${accession[0]}"
# 	GSM="${accession[1]}"
# # Get metadata of the GEO series
# 	if [ ! -f "$src_dir/../raw/GEO/$GSE" ];then
# 		Rscript "$src_dir/get-metadata-geo-sequencing.R $GSE $src_dir"
# 	fi
# # Download SRR files corresponding to GSM accession
# 	bash $src_dir/download-files-GEO-sequencing.sh -g $GSE -s $GSM $ignore_exist_flag
# # Process the data	
# 	if [ $database == "GEO" ] && [ $datatype == 'RRBS' ]; then
# 	       echo "$GSM"	
# 	       bash $src_dir/processing-geo-sequencing-RRBS.sh -g $GSE -s $GSM -c $core
# 	elif [ $database == "GEO" ] && [ $datatype == 'WGBS' ]; then
# 		echo "$GSM"
# 		bash $src_dir/processing-geo-sequencing-WGBS.sh -g $GSE -s $GSM -c $core
# 	fi
# 	echo "$GSE"	
# done 3< $src_dir/../annotation/$database$accession_filename1$datatype$accession_filename2
