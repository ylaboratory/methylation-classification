# date updated: Apr 2024
# clean up and filter data from GEO saved in /data/GEO/compiled/ for disease samples
# save to /data/GEO/preprocessed/450K_Mvalues_disease

import pandas as pd
import numpy as np
import dill

from sklearn.preprocessing import StandardScaler
import utils

start_from_beginning = True

if start_from_beginning:
    print("starting from beginning...")
    
    print("...loading data...")
    #load metadata, sample_id of beta values, and beta values of those with both beta values and metadata
    meta_all=pd.read_csv('./../data/GEO/compiled/all_metadata_annotated_may2023.txt.gz', sep=',')
    beta_samples=pd.read_csv('./../data/GEO/compiled/450K_betavalues_probe.txt.gz', sep='\t', index_col=0, nrows=0).columns.tolist()
    probes=pd.read_csv('./../data/GEO/compiled/450K_betavalues_probe.txt.gz', sep='\t', usecols=['probe','probe.2'])
    beta_all=pd.read_csv('./../data/GEO/compiled/450K_betavalues_probe.txt.gz', sep='\t', 
                         usecols=list(set(beta_samples).intersection(set(meta_all['sample_id']))))
    beta_all.index=probes['probe']

    island=pd.read_csv('./../annotation/HM450.hg38.manifest.gencode.v36.tsv', sep='\t')
    island['position']=island['CpG_beg']+1
    island=island.drop(['CpG_beg','CpG_end', 'probe_strand'], axis=1)
    island.columns=['chr', 'probe', 'genesUniq','geneNames','transcriptTypes','transcriptIDs','disToTSS','CGI','CGI type','position']
    island.index=island['probe']
    island=island.drop(['probe'], axis=1)
    island=island.loc[probes['probe']]

    print(f"all meta shape: {meta_all.shape}")
    print(f"all beta shape: {beta_all.shape}")
    print()

    print("...removing X, Y, NA...")
    #remove X and Y chromosomes
    probe_remove=list(island[island['chr']=='chrX'].index)
    probe_remove.extend(list(island[island['chr']=='chrY'].index))
    probe_remove.extend(list(island[island['chr']=='chrM'].index))
    dropped_beta=beta_all.drop(index=probe_remove)
    dropped_island=island.drop(index=probe_remove)

    #remove probes with no info
    na_probe=list()
    for p in island.index:
        test_row=island.loc[p]
        if (pd.isna(test_row['chr']) and pd.isna(test_row['position'])):
            na_probe.extend([p])
    dropped_beta=dropped_beta.drop(index=na_probe)
    filtered_island=dropped_island.drop(na_probe, axis=0)

    #Get samples with both methylation data and metadata annotation
    beta_samples=list(dropped_beta.columns)
    meta_samples=meta_all['sample_id']
    both_beta_metadata=list(set(beta_samples).intersection(set(meta_samples)))

    #Remove overlapping samples
    index_to_delete=['GSM3813550', 'GSM3813542', 'GSM3813555', 'GSM3813541', 'GSM3813540']
    for i in index_to_delete:
        if i in both_beta_metadata: both_beta_metadata.remove(i)

    #filter for those with both beta values and metadata info
    beta=dropped_beta[both_beta_metadata]
    meta=meta_all[meta_all['sample_id'].isin(both_beta_metadata)]
    meta.index=meta['sample_id']
    meta=meta.drop(meta[meta['platform']=='GPL10558'].index)
    meta=meta.drop(meta[meta['platform']=='GPL21135'].index)
    meta=meta.drop_duplicates('sample_id', 'first')

    island=dropped_island.copy(deep=True)
    island
    
    print("...control, platform, disease filter...")

    #control filter
    control_meta=meta[meta['treatment']==False]
    control_beta=beta[control_meta.index].T

    #platform filter
    avoid_platform=['GPL21145', 
                    'GPL23976',
                    'GPL16304']
    drop_idx=list()
    drop_idx.extend(list(control_meta[control_meta['platform'].isin(avoid_platform)].index))
    control_meta=control_meta.drop(drop_idx)
    control_beta=control_beta.loc[control_meta.index]

    #disease filter

    normal_meta=control_meta[control_meta['normal']=='disease']
    normal_beta=control_beta.loc[normal_meta.index]

    #tissue names in uberon terms
    normal_meta['tissue_name']=normal_meta['tissue_name'].str.lower()
    normal_meta['tissue_name'].replace('umbilical blood','umbilical cord blood', inplace=True)
    normal_meta['tissue_name'].replace('t-lymphocyte','lymphocyte', inplace=True)
    normal_meta['tissue_name'].replace('spermatozoon','sperm', inplace=True)
    normal_meta['tissue_name'].replace('rectal','rectum', inplace=True)
    normal_meta['tissue_name'].replace('skin','skin epidermis', inplace=True)
    normal_meta['tissue_name'].replace('tracheal epithelium','trachea', inplace=True)
    normal_meta['tissue_name'].replace('islet of langerhans','islet of Langerhans', inplace=True)
    normal_meta['tissue_name'].replace('muscle','muscle tissue', inplace=True)
    normal_meta['tissue_name'].replace('hair','hair follicle', inplace=True)
    normal_meta['tissue_name'].replace('frontal lobe cortex','frontal lobe', inplace=True)
    normal_meta['tissue_name'].replace('left atrium','left cardiac atrium', inplace=True)
    normal_meta['tissue_name'].replace('nucleated red blood cell','blood', inplace=True)
    normal_meta['tissue_name'].replace('gastric mucosa', 'gastric gland', inplace=True)
    normal_meta['tissue_name'].replace('bronchial epithelium', 'bronchus', inplace=True)
    normal_meta['tissue_name'].replace('trachea smooth muscle tissue', 'smooth muscle of trachea', inplace=True)
    normal_meta['tissue_name'].replace('thymus gland', 'thymus', inplace=True)
    normal_meta['tissue_name'].replace('oropharyngeal tissue', 'oropharyngeal gland', inplace=True)
    normal_meta['tissue_name'].replace('ovarian surface epithelium', 'ovarian surface epithelial cell', inplace=True)
    normal_meta['tissue_name'].replace('prostatic tissue', 'prostate epithelium', inplace=True)
    normal_meta['tissue_name'].replace('airway','respiratory airway', inplace=True)
    normal_meta['tissue_name'].replace('cartilage','cartilage tissue', inplace=True)
    normal_meta['tissue_name'].replace('esophageal squamous epithelium','esophagus squamous epithelium', inplace=True)
    normal_meta['tissue_name'].replace('oral cavity epithelium','oral epithelium', inplace=True)
    
    print("...converting to Mvalues...")
    #convert to Mvalue
    normal_Mv=np.log2(np.divide((normal_beta),(1-normal_beta)))
    normal_Mv=pd.DataFrame(normal_Mv, columns=normal_beta.columns, index=normal_beta.index)
    island=island.loc[normal_Mv.columns]
    
    print(f"disease Mv shape: {normal_Mv.shape}")
    print(f"disease meta shape: {normal_meta.shape}")
    print()
    
    
    #save all as pickle file
    save_all_disease=True
    
    if save_all_disease:
        filename = f"./../data/GEO/preprocessed/450K_Mvalues_disease"
        print(f"...saving as {filename}...")
        with open(filename, "wb") as dill_file:
            dill.dump([normal_Mv, normal_meta, island], dill_file, protocol=4)