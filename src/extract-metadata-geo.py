# This script parses the matrix series file from GEO
import argparse
import os
# read the list of accessions
parser = argparse.ArgumentParser(description='Parsing metadata from GEO ')
parser.add_argument('file_name', nargs='?', type=str, action='store', help='File name of the GEO accessions')
args = parser.parse_args()
file_path = args.file_name
accession_list_file = open(file_path, 'r')
accession_list = accession_list_file.read().split('\n')


def separation(content):
    content = content.replace('"', '').replace('\n', '')
    info = content.split('\t')[1:]
    return info


ids_list = []
source_list = []
character_full_list = []
platform_list = []
for i in range(1, len(accession_list)-1):
    accession = accession_list[i]
    character_list = []
    metadata = open("../raw/GEO/" + accession + "_series_matrix.txt", 'r')
    metadata_line = metadata.readlines()
    for j in metadata_line:
        info = j.split('\t')
        if info[0] == "!Sample_geo_accession":
            ids = separation(j)
            ids_list = ids_list + ids
        elif info[0] == "!Sample_source_name_ch1":
            source = separation(j)
            source_list = source_list + source
        elif info[0] == "!Sample_characteristics_ch1":
            character = separation(j)
            character_list = character_list + character
        elif info[0] == "!Sample_platform_id":
            platform = separation(j)
            platform_list = platform_list + platform
# This block is used to reorganize the information of sample characters.
# Sometimes multiple rows in the series matrix file have the title of !Sample_characteristics_ch1
# We integrate all !Sample_characteristics_ch1 information for each sample into one block of text
    reorganize_character = [None]*len(ids)
    for k in range(0, len(ids)):
        character_concat = ""
        for z in range(0, len(character_list)):
            if z % len(ids) == k:
                character_concat = character_concat + " " + character_list[z]
        reorganize_character[k-1] = character_concat
    character_full_list = character_full_list + reorganize_character
f = open('../data/GEO/Metadata_summary.txt', 'w')
for i in range(0, len(ids_list)):
    f.write(ids_list[i] + '\t' + source_list[i] + '\t' + character_full_list[i] + '\t' + platform_list[i] + '\n')
f.close()
