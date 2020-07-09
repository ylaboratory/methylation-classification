#!/usr/bin/env python
# coding: utf-8

# In[16]:


import requests, json
import os
from collections import defaultdict


# In[53]:


text_file = open(os.getcwd()+'/../annotation/ENCODE_microarray_accession_full_raw.txt', 'r')
encode_accession = text_file.read().split('\n')
encode_accession = list(filter(None, encode_accession))


# In[55]:


headers={'accept': 'application/json'}
overlap_dic = defaultdict(list)
for j in range(0,len(encode_accession)):
    url='https://www.encodeproject.org/' + encode_accession[j] + '/?frame=object'
    response = requests.get(url, headers=headers)
    accession = response.json()
    overlap_element = [i for i in accession["dbxrefs"] if "GEO" in i]
    if len(overlap_element)!=0:
        for i in range(0,len(overlap_element)):
                sep = ':'
                overlap_GEO = overlap_element[i].split(sep, 1)[1]
                overlap_dic[encode_accession[j]].append(overlap_GEO)
            


# In[57]:


encode_accession_no_overlap = set(encode_accession) - set(overlap_dic.keys())


# In[61]:


f = open(os.getcwd()+'/../annotation/ENCODE_microarray_accession_full.txt', "w")
for i in encode_accession_no_overlap:
    f.write(i)
    f.write('\n')
f.close()

