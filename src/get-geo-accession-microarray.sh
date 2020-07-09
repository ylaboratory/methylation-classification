#!/bin/bash
src_dir=$(dirname $0)
esearch -db gds -query "GPL13534[ACCN] OR GPL16304[ACCN] OR GPL21145[ACCN] OR GPL23976[ACCN] AND gse[E
> TYP] AND idat[suppFile] AND human[organism] NOT SuperSeries" | efetch -format docsum | xtract -pattern DocumentSummary -element GSE >> $src_dir/../annotation/GEO_microarray_accession_full.txt
sed -i -e 's/^/GSE/' $src_dir/../annotation/GEO_microarray_accession_full.txt
