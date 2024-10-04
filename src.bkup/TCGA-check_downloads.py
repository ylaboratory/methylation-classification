import os
import pandas as pd

def read_id_list(file_path, delimiter='\t'):
    """Reads the list of IDs from a delimited text file."""
    id_list_df = pd.read_csv(file_path, delimiter="\t", header=0)
    id_list = id_list_df.iloc[:, 0].tolist()
    return id_list

def check_downloaded_ids(id_list, directory):
    """Checks which IDs from the list are present in the directory."""
    downloaded_ids = []
    for filename in os.listdir(directory):
        if filename in id_list:
            downloaded_ids.append(filename)
    return downloaded_ids

# Path to the file containing the list of IDs
id_list_file = '/grain/mk98/methyl/methylation-classification/annotation/TCGA_manifest_breast.txt'

# Directory containing downloaded files
download_directory = '/grain/mk98/methyl/methylation-classification/raw/TCGA/'

# Read the list of IDs
id_list = read_id_list(id_list_file)

# Check which IDs from the list are downloaded
downloaded_ids = check_downloaded_ids(id_list, download_directory)

output_file = '/grain/mk98/methyl/methylation-classification/annotation/TCGA-successful-downloads.txt'
with open(output_file, 'w') as file:
    for id in downloaded_ids:
        file.write(id + '\n')
