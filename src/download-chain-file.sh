src=$(dirname $0)
wget 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz' -O $src/../annotation/hg19ToHg38.over.chain.gz
gunzip $src/../annotation/hg19ToHg38.over.chain.gz
