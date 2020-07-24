#!/usr/bin/env python
# coding: utf-8
# This script finds the overlapping accession number of ENCODE and GEO database
# Input: accession numbers from ENCODE database
# Output: ENCODE accession numbers with no overlap with GEO database & mapping
# of overlapping ENCODE accession to GSM number
import requests
import json
import os
import argparse
# read the full encode accession list
parser = argparse.ArgumentParser(description='Get overlapping ENCODE and GEO accession ')
parser.add_argument('file_directory', nargs='?', type=str, action='store', default='./../annotation', help='directory containing the ENCODE accession list')
args = parser.parse_args()
file_dir = args.file_directory
if not os.path.exists(file_dir):
    raise ValueError("Path" + str(file_dir) + "does not exist")
text_file = open(file_dir+'/ENCODE_microarray_accession_full.txt', 'r')
encode_accession_full = text_file.read().split('\n')
encode_accession_full = [i for i in encode_accession_full if i is not ""]
# query individual ENCODE accession to get the external GEO accession if available
headers = {'accept': 'application/json'}
GSE_accession_list = []
GSM_accession_list = []
encode_accession_overlap = []
for j in range(0, len(encode_accession_full)):
    url = 'https://www.encodeproject.org/' + encode_accession_full[j] + '/?frame=object'
    response = requests.get(url, headers=headers)
    accession = response.json()
    external_accession = [i for i in accession["dbxrefs"] if "GEO" in i]
    if len(external_accession) > 0:
        for i in range(0, len(external_accession)):
            external_GEO = external_accession[i].split(':', 1)[1]
            if "GSE" in external_GEO:
                GSE_accession_list.append(external_GEO)
            elif "GSM" in external_GEO:
                GSM_accession_list.append(external_GEO)
        encode_accession_overlap.append(encode_accession_full[j])

# write the mapping of overlapped accession from ENCODE to GEO
f = open(file_dir+'/ENCODE_microarray_accession_GEO_overlap.txt', "w")
for i in range(0, len(encode_accession_overlap)):
    f.write(encode_accession_overlap[i] + '\t' + GSE_accession_list[i] + '\t' + GSM_accession_list [i] + '\n')
f.close()
# get the ENCODE accession after removing overlaps with GEO
encode_accession_GEO_removed = set(encode_accession_full) - set(encode_accession_overlap)
f = open(file_dir + '/ENCODE_microarray_accession_GEO_removed.txt', "w")
for i in encode_accession_GEO_removed:
    f.write(i + '\n')
f.close()
