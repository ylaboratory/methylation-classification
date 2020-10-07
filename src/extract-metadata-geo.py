# This script parses the matrix series file from GEO
import argparse
import pandas as pd
import os
# read the list of accessions
parser = argparse.ArgumentParser(description='Parsing metadata from GEO ')
parser.add_argument('file_name', nargs='?', type=str, action='store', help='File name of the GEO accessions')
args = parser.parse_args()
file_path = args.file_name
accession_list = pd.read_csv(file_path, sep='\t', header=None, engine='python')


def separation(string):
    info = string.split('\t')
    for i in range(0, len(info)):
        info[i] = info[i].replace('"', '')
        info[i] = info[i].replace('\n', '')
    del info[0]
    return info


df_list = []
for i in range(0, accession_list.shape[0]):
    accession = accession_list.iloc[i, 0]
    metadata = open("./../raw/GEO/" + accession + "_series_matrix.txt", 'r')
    metadata_line = metadata.readlines()
    for j in metadata_line:
        info = j.split('\t')
        if info[0] == "!Sample_geo_accession":
            ids = separation(j)
        elif info[0] == "!Sample_source_name_ch1":
            source = separation(j)
        elif info[0] == "!Sample_characteristics_ch1":
            character = separation(j)
        elif info[0] == "!Sample_platform_id":
            platform = separation(j)
    d = {'SampleID': ids, 'Source': source, "Character": character, 'Platform': platform }
    df = pd.DataFrame(data=d)
    df_list.append(df)
summary_df = pd.concat(df_list)
summary_df.to_csv("./../data/GEO/Metadata_summary.txt", sep='\t', index=False, header=True)
