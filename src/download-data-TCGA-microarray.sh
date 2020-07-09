# This script downloads the microarray data from TCGA database
#!/bin/bash
# Install the gdc_client if it is not yet available
src_dir=$(dirname $0)
gdc_client=$src_dir/gdc-client_v1.5.0_Ubuntu_x64.zip
if [ ! -f "$gdc_client" ]
then
	        wget https://gdc.cancer.gov/system/files/authenticated%20user/0/gdc-client_v1.5.0_Ubuntu_x64.zip -O $gdc_client
		unzip $gdc_client
fi
# Download all files from the manifest, here test data is used 
$src_dir/gdc-client download -d $src_dir/../raw/TCGA --no-annotations --http-chunk-size 9000 -m $src_dir/../annotation/TCGA_microarray_manifest_test.txt

