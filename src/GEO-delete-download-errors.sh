#!/bin/bash
src_dir=$(dirname $0)
raw_dir=$src_dir/../raw/GEO

: ' by current geo download logic, the folders should be empty unless encountering an error while downloading the betavalues or the metadata - this script deletes those folders that are not empty to redownload during the next round of data download
'

for dir in $raw_dir/*/; do
    if [ -z "$(ls -A $dir)" ]; then
	echo "$dir Empty"
    else
	echo "$dir Not Empty"
	rm -r $dir
    fi
done    
