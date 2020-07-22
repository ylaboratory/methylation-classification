# This script downloads the microarray data from TCGA database
#!/bin/bash
usage() {
	        echo "script usage: $(basename $0) [-m filename of the manifest] "
	}
while getopts ":m:" Option; do
case "$Option" in
	m)      manifest=$OPTARG; echo "$OPTARG is used to download TCGA data" ;;
	:)	echo "Error: -$OPTARG requires an argument"; usage; exit -1;;
	?)	echo "Error: unrecognized argument"; usage; exit -1;;
	*)      usage; exit -1;;
esac
done
# Install the gdc_client if it is not yet available
src_dir=$(dirname $0)
mkdir -p $src_dir/bin
if [ ! -f "$src_dir/bin/gdc-client"* ]
then
	        wget https://gdc.cancer.gov/system/files/authenticated%20user/0/gdc-client_v1.5.0_Ubuntu_x64.zip -O $src_dir/bin/gdc-client.zip
fi
unzip -d $src_dir/bin $src_dir/bin/gdc-client.zip
# Download all files from the manifest, here test data is used 
$src_dir/bin/gdc-client download -d $src_dir/../raw/TCGA --no-annotations --http-chunk-size 9000 -m $src_dir/../annotation/$manifest

