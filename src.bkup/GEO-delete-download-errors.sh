#!/bin/bash
src_dir=$(dirname $0)
raw_dir=$src_dir/../raw/GEO

for dir in $raw_dir/*/; do
    if [ -z "$(ls -A $dir)" ]; then
	echo "$dir Empty"
	rm -r $dir
    else
	echo "$dir Not Empty"
    fi
done    
