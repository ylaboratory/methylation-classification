#!/bin/bash
src_dir=$(dirname $0)
raw_dir=$src_dir/../raw/GEO

for dir in $raw_dir/*/; do
    if [ -z "$(ls -A $dir)" ]; then
	echo "$dir Empty"
    else
	echo "$dir Not Empty"
	rm -r $dir
    fi
done    
