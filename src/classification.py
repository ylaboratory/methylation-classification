from os.path import abspath, dirname
import pandas as pd
import numpy as np
import obonet
import networkx as nx

from sklearn.feature_selection import VarianceThreshold
from sklearn.model_selection import train_test_split
from sklearn.svm import SVC
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc, average_precision_score, f1_score
from sklearn.metrics import precision_recall_curve, classification_report

from collections import Counter

import warnings
warnings.filterwarnings('ignore')

import joblib
import sys
sys.modules['sklearn.externals.joblib'] = joblib

from pystruct.models import BinaryClf, MultiClassClf, MultiLabelClf
from pystruct.learners import OneSlackSSVM, LatentSSVM, NSlackSSVM, FrankWolfeSSVM,SubgradientSSVM
from sklearn.preprocessing import LabelEncoder, OneHotEncoder

import time

#import data (may take couple hours)
beta_ncit=pd.read_csv('./../data/GEO/all_betavalues.txt.gz', sep='\t')
meta_ncit=pd.read_csv('./../data/GEO/all_metadata_annotated.txt.gz',sep=',')

beta_ncit_all=beta_ncit.copy(deep=True)
meta_ncit_all=meta_ncit.copy(deep=True)
tissue_ontology_tree=obonet.read_obo('./../annotation/ncit.obo')

#remove X and Y chromosomes
index_remove=beta_ncit[beta_ncit['V1']=='chrX'].index
index_remove.append(beta_ncit[beta_ncit['V1']=='chrY'].index)
beta_ncit=beta_ncit.drop(index_remove)
beta_ncit.reset_index(drop=True)

#label features as: [chr# position]
feature=np.array(beta_ncit['V1'])+' '
feature=pd.DataFrame(np.array(feature)+np.array(list(map(str, beta_ncit['V2']))))
feature.index=beta_ncit.index

#Get samples with both methylation data and metadata annotation
beta_samples=list(beta_ncit.columns)[2:]
meta_samples=meta_ncit['sample_id']
both_beta_metadata=list(set(beta_samples).intersection(meta_samples))

#Remove overlapping samples
index_to_delete=['GSM3813550', 'GSM3813542', 'GSM3813555', 'GSM3813541', 'GSM3813540']
for i in index_to_delete:
    if i in both_beta_metadata: both_beta_metadata.remove(i)
        
#filter for those with beta values and metadata info
beta=beta_ncit[both_beta_metadata]
meta=meta_ncit[meta_ncit['sample_id'].isin(both_beta_metadata)]
meta.index=meta['sample_id']

#remove duplicated entries from multiple platforms other than microarray
meta=meta.drop(meta[meta['platform']=='GPL10558'].index)
meta=meta.drop(meta[meta['platform']=='GPL21135'].index)
meta=meta.drop_duplicates('sample_id', 'first')

#filter for normal (non-disease) tissue
normal_meta=meta[meta['disease_name']=='normal']
normal_beta=beta[normal_meta.index]

#typos in annotation
normal_meta.replace('Whoel Blood','Whole Blood', inplace=True)
normal_meta.replace('Breast Neoplasm','Breast', inplace=True)
normal_meta.replace('Embryonic Stemm Cell','Embryonic Stem Cell', inplace=True)
normal_meta.replace('T-Lymphocytes','T-Lymphocyte', inplace=True)
normal_meta.replace('Adipose','Adipose Tissue', inplace=True)
normal_meta.replace('Prostate','Prostate Gland', inplace=True)
normal_meta.replace('Interstinal Mucosa','Intestinal Mucosa',inplace=True)

#exclude mislabeled tissues/tissue with only one sample
#updated: Sept 2021
avoid=['Schizophrenia','Sigmoid Colon']
normal_sample_index=[i for i in range(len(normal_tissue_term)) if (normal_tissue_term.iloc[i] not in avoid)]
normal_meta=normal_meta.iloc[normal_sample_index]
normal_beta=beta[normal_meta.index]
normal_tissue_term=normal_meta['tissue_name']
normal_sample_value=normal_tissue_term.iloc[normal_sample_index]

#dictionary of names to id from NCIt
id_to_name = {id_:data.get('name') for id_, data in tissue_ontology_tree.nodes(data=True)}
name_to_id = {data.get('name'):id_ for id_, data in tissue_ontology_tree.nodes(data=True)}

#make sure there are no labels with one sample only 
only_one=[x for x in Counter(normal_meta['tissue_name']).keys() if Counter(normal_meta['tissue_name'])[x]==1]
if len(only_one)!=0:
    print(f"only one sample in {only_one}")
    exit()

#filter by variance and get mask
thresh=0.10
selector=VarianceThreshold(threshold=(thresh))
filtered_beta=pd.DataFrame(selector.fit_transform(normal_beta.T.values), index=normal_beta.columns)
filtered_meta=normal_meta
mask=selector.get_support(indices=True)

#filtered [chr# position] with mask
filtered_feature=feature.iloc[mask]

def get_indicator_matrix(node, beta, meta):
    #tissue in question
    id_node=name_to_id[node]
    tissue_name=meta['tissue_name']
    #tissue/organ/organism
    descendant_nodes=[gsm for gsm in beta.index if name_to_id[tissue_name.loc[gsm]] in nx.descendants(tissue_ontology_tree, id_node)]
    #cell/subtype
    ancestor_nodes=[gsm for gsm in beta.index if name_to_id[tissue_name.loc[gsm]] in nx.ancestors(tissue_ontology_tree, id_node)]

    #propagate label through cell/subtype levels
    positive_sample_index=[gsm for gsm in beta.index if (name_to_id[tissue_name[gsm]]==id_node or gsm in ancestor_nodes)]

    #control label as -1 for pystruct
    y_matrix=pd.DataFrame(-1, index=meta.index, columns=['bin'])
    y_matrix.loc[positive_sample_index]=1

    #get rid of nodes with conflicting tissue/organ/organism
    keep_nodes=[gsm for gsm in meta.index if gsm not in descendant_nodes]
    y_matrix_sample=y_matrix.loc[keep_nodes].values
    return(y_matrix_sample,keep_nodes,positive_sample_index)

#binary one vs. rest - no descendents, label propagate for ancestors
results=pd.DataFrame()

start_time=time.time()
for i in set(np.unique(filtered_meta['tissue_name'])):
    print(i)
    [y_matrix_single, keep_nodes, positive_sample]=get_indicator_matrix(i, filtered_beta, filtered_meta)
    x_matrix_single=filtered_beta.loc[keep_nodes].values
    
    #split train/test
    X_train,X_test,y_train,y_test=train_test_split(x_matrix_single, y_matrix_single, stratify=y_matrix_single,test_size=0.20, random_state=42) 

    clf=BinaryClf()
    svm=FrankWolfeSSVM(clf, verbose=1, tol=0.1)
    svm.fit(X_train, y_train)
    pred=svm.predict(X_test)
    acc=sum(1 for x,y in zip(pred,y_test) if x == y)/len(pred)
    results.loc['accuracy', i]=acc
print(results.T)
results.to_csv('./../results/binary.csv')    
