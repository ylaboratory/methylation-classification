#!/bin/bash
# relative path is based on src directory, where the location of this script is
src_dir=$(dirname $0)
usage() {
	echo "script usage: $(basename $0) ( -t datatype ) "
	}
while getopts ":t:" Option; do
	case "$Option" in
		t) datatype=$OPTARG; echo "getting accession of $OPTARG data";;
	esac
done
# Get all WGBS accession from SRA
if [ ! -f $src_dir/../annotation/SRA-$datatype-accession.txt ]; then
	esearch -db sra -query "$datatype AND human[organism] AND bisulfite seq[Strategy] AND public[Access]" | efetch -format docsum | xtract -pattern DocumentSummary -element Experiment@acc >> $src_dir/../annotation/SRA-$datatype-accession.txt
fi
declare -a GSEarr
declare -a GSMarr
declare -a SRXarr 
while IFS= read -r -u3 line
do
	if [[ $line == SRX* ]]; then
		GSE=$( esearch -db gds -query $line | efetch -format docsum | xtract -pattern DocumentSummary -element GSE )
		GSM=$( esearch -db gds -query $line | efetch -format docsum | xtract -pattern DocumentSummary -element Accession )
		if [ ! -z $GSE ] && [ ${#GSM[@]} == 1 ]; then
			if [[ $GSE == *";"* ]]; then
				GSE=$( echo $GSE | cut -d ";" -f 1 )
			fi
			GSEarr=("GSE$GSE" "${GSEarr[@]}")
			GSMarr=("$GSM" "${GSMarr[@]}")
			SRXarr=("$line" "${SRXarr[@]}")
			echo $GSM
			echo $GSE
		fi
	fi
done 3< $src_dir/../annotation/SRA-$datatype-accession.txt
paste <(printf "%s\n" "${GSEarr[@]}") <(printf "%s\n" "${GSMarr[@]}") <(printf "%s\n" "${SRXarr[@]}") >> $src_dir/../annotation/GEO_sequencing_accession_"$datatype"_full.txt
