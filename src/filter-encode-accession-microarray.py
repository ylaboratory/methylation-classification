#!/usr/bin/env python
# coding: utf-8
# This script finds the overlapping accession numbers of ENCODE and GEO database
# Input: accession numbers from ENCODE database
# Output: ENCODE accession numbers with no overlap with GEO database & mapping of overlapping ENCODE accession to GSM number
import requests, json
from collections import defaultdict
import argparse
# read the full encode accession list
parser = argparse.ArgumentParser(description='Get overlapping ENCODE and GEO accession ')
parser.add_argument('file_directory', type=str, action = 'store',
                    help='directory containing the ENCODE accession list')
args = parser.parse_args()
file_dir = args.file_directory
print(file_dir)
def read_encode(file_dir):
    text_file = open(file_dir+'/ENCODE_microarray_accession_full.txt', 'r')
    encode_accession_full = text_file.read().split('\n')
    return (list(filter(None, encode_accession_full)))
encode_accession_full = read_encode(file_dir)
# query individual ENCODE accession to get the external GEO accession if available
headers={'accept': 'application/json'}
overlap_dic = defaultdict(list)
for j in range(0,len(encode_accession_full)):
    url = 'https://www.encodeproject.org/' + encode_accession_full[j] + '/?frame=object'
    response = requests.get(url, headers=headers)
    accession = response.json()
    external_accession = [i for i in accession["dbxrefs"] if "GEO" in i]
    if len(external_accession) != 0:
        for i in range(0,len(external_accession)):
                external_GEO = external_accession[i].split(':', 1)[1]
                overlap_dic[encode_accession_full[j]].append(external_GEO)
# write the mapping of overlapped accession from ENCODE to GEO
def map_encode_geo(file_dir, overlap_dic):
    f = open(file_dir+'/ENCODE_microarray_accession_GEO_overlap.txt', "w")
    for key, value in overlap_dic.items():
        f.write( str(key) + '\t' + str(value) + '\n')
    f.close()
map_encode_geo(file_dir,overlap_dic)
# get the ENCODE accession after removing overlaps with GEO
encode_accession_GEO_removed = set(encode_accession_full) - set(overlap_dic.keys())
def write_encode(file_dir, encode_accession_GEO_removed):
    f = open(file_dir + '/ENCODE_microarray_accession_GEO_removed.txt', "w")
    for i in encode_accession_GEO_removed:
        f.write(i + '\n')
    f.close()
write_encode(file_dir,encode_accession_GEO_removed)
