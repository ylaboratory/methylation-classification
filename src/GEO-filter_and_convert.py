#date updated: Sep 2024
#clean up and filter data from GEO saved in /data/GEO/compiled/ as all_Mv.txt.gz and all_metadata.txt.gz

date='sep2024'

import pandas as pd
import numpy as np
import pickle

from sklearn.preprocessing import StandardScaler
import utils

start_from_beginning = True

if start_from_beginning:
    print("starting from beginning...")
    
    print("...loading data...")
    #load metadata, sample_id of beta values, and beta values of those with both beta values and metadata
    meta_all=pd.read_csv(f'./../data/GEO/compiled/all_metadata_annotated_{date}.txt.gz', sep=',')
    beta_samples=pd.read_csv(f'./../data/GEO/compiled/450K_betavalues_probe_{date}.txt.gz', sep='\t', index_col=0, nrows=0).columns.tolist()
    probes=pd.read_csv(f'./../data/GEO/compiled/450K_betavalues_probe_{date}.txt.gz', sep='\t', usecols=['probe','probe.2'])
    beta_all=pd.read_csv(f'./../data/GEO/compiled/450K_betavalues_probe_{date}.txt.gz', sep='\t', 
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
    
    print("...control, platform, normal filter...")

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

    #normal filter

    normal_meta=control_meta[control_meta['normal']=='normal']
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
    
    print(f"normal Mv shape: {normal_Mv.shape}")
    print(f"normal meta shape: {normal_meta.shape}")
    print()
    
    
    #save all as pickle file
    save_all_normal=True
    
    if save_all_normal:
        filename = f"./../data/GEO/preprocessed/450K_Mvalues"
        print(f"...saving as {filename}...")
        with open(filename, "wb") as dill_file:
            pickle.dump([normal_Mv, normal_meta, island], dill_file, protocol=4)

else:
    filename = f"./../data/GEO/preprocessed/450K_Mvalues"
    print(f"loading from {filename}...")
    with open(filename, "rb") as dill_file:
        normal_Mv, normal_meta, island = pickle.load(dill_file)
        
    print(f"normal Mv shape: {normal_Mv.shape}")
    print(f"normal meta shape: {normal_meta.shape}")
    print()
        

atleast=2
atleast_sample_per_study = 4
atleast_sample_per_tissue = 10

print("...size filter...")

#multiple filter
mult_thresh=atleast

#sample size per study filter
large_enough_study = list()
for series in np.unique(normal_meta['series']):
    if normal_meta[normal_meta['series']==series].shape[0]>=atleast_sample_per_study:
        # [mult_idx.remove(idx) for idx in normal_meta[normal_meta['series']==series].index if idx in mult_idx]
        large_enough_study += [series]
mult_meta = normal_meta[normal_meta['series'].isin(large_enough_study)]

#filter for those in multiple datasets only
enough_study_in_tissue = list()
for tissue in np.unique(mult_meta['tissue_name']):
    num_gse=len(np.unique(mult_meta[mult_meta['tissue_name']==tissue]['series']))
    if num_gse>=mult_thresh:
        enough_study_in_tissue.append(tissue)
mult_meta = mult_meta[mult_meta['tissue_name'].isin(enough_study_in_tissue)]
        
#sample size per tissue filter
large_enough_tissue = list()
for tissue in np.unique(mult_meta['tissue_name']):
    if mult_meta[mult_meta['tissue_name']==tissue].shape[0]>=atleast_sample_per_tissue:
        large_enough_tissue += [tissue]
mult_meta = mult_meta[mult_meta['tissue_name'].isin(large_enough_tissue)]
        
# mult_meta=normal_meta.loc[mult_idx]
mult_Mv=normal_Mv.loc[mult_meta.index]

print(f"mult Mv shape: {mult_Mv.shape}")
print(f"mult meta shape: {mult_meta.shape}")
print()

#save all noscaling
save_all_noscaling = False

if save_all_noscaling:
    filename = f"./../data/GEO/preprocessed/450K_Mvalues_atleast{atleast}_noscaling"
    print(f"...saving as {filename}...")
    with open(filename, "wb") as dill_file:
        pickle.dump([mult_Mv, mult_meta, island], dill_file, protocol=4)
    print()
    
#save all scaled
save_all_scaled = True

if save_all_scaled:
    print(f"...samplewise scaling...")
    
    scaler = StandardScaler().fit(mult_Mv.transpose().values)
    mult_Mv_scaled = scaler.transform(mult_Mv.transpose().values)
    mult_Mv_scaled = pd.DataFrame(mult_Mv_scaled.transpose(), index=mult_Mv.index, columns=mult_Mv.columns)
    
    filename = f"./../data/GEO/preprocessed/450K_Mvalues_atleast{atleast}_samplewise"
    print(f"...saving as {filename}...")
    with open(filename, "wb") as dill_file:
        pickle.dump([mult_Mv_scaled, mult_meta, island], dill_file, protocol=4)
    print()

thresh = 2
print(f"...variance filter (threshold = {thresh})...")

highvar_features = utils.weighted_variance(mult_Mv, mult_meta, island, thresh)
highvar_Mv = mult_Mv[highvar_features]
highvar_meta = mult_meta
highvar_mapping = island.loc[highvar_features]

print(f"highvar Mv shape: {highvar_Mv.shape}")
print(f"highvar meta shape: {highvar_meta.shape}")
print()

#save highvar noscaling
save_thresh2_noscaling = False

if save_thresh2_noscaling:
    
    filename = f"./../data/GEO/preprocessed/450K_Mvalues_atleast{atleast}_thresh{thresh}_noscaling"
    print(f"...saving as {filename}...")
    with open(filename, "wb") as dill_file:
        pickle.dump([highvar_Mv, highvar_meta, highvar_mapping], dill_file, protocol=4)
    print()
    
#save scale first then var
save_scale_then_thresh = False

if save_scale_then_thresh:
    thresh = 0.1
    print(f"...scale first then variance filter (threshold = {thresh})...")

    scaled_highvar_features = utils.weighted_variance(mult_Mv_scaled, mult_meta, island, thresh)
    scaled_highvar_Mv = mult_Mv_scaled[scaled_highvar_features]
    scaled_highvar_meta = mult_meta
    scaled_highvar_mapping = island.loc[scaled_highvar_features]

    print(f"highvar Mv shape: {scaled_highvar_Mv.shape}")
    print(f"highvar meta shape: {scaled_highvar_meta.shape}")

    filename = f"./../data/GEO/preprocessed/450K_Mvalues_atleast{atleast}_samplewise_thresh{thresh}"
    print(f"...saving as {filename}...")
    with open(filename, "wb") as dill_file:
        pickle.dump([scaled_highvar_Mv, scaled_highvar_meta, scaled_highvar_mapping], dill_file, protocol=4)
    print()
    
#save highvar samplewise scaling
save_highvar_samplewise = False

if save_highvar_samplewise:
    scaler = StandardScaler().fit(highvar_Mv.transpose().values)
    highvar_Mv_scaled = scaler.transform(highvar_Mv.transpose().values)
    highvar_Mv_scaled = pd.DataFrame(highvar_Mv_scaled.transpose(), index=highvar_Mv.index, columns=highvar_Mv.columns)

    filename = f"./../data/GEO/preprocessed/450K_Mvalues_atleast{atleast}_thresh{thresh}_samplewise"
    print(f"...saving as {filename}...")
    with open(filename, "wb") as dill_file:
        pickle.dump([highvar_Mv_scaled, highvar_meta], dill_file, protocol=4)
    print()


save_highvar_samplewise_folds = False

if save_highvar_samplewise_folds:
    for fold in range(3):
        print(f"fold: {fold}")
        scaler = StandardScaler().fit(highvar_Mv.transpose().values)
        highvar_Mv_scaled = scaler.transform(highvar_Mv.transpose().values)
        highvar_Mv_scaled = pd.DataFrame(highvar_Mv_scaled.transpose(), index=highvar_Mv.index, columns=highvar_Mv.columns)

        holdout_Mv, holdout_meta, rest_Mv, rest_meta = utils.make_holdout(highvar_Mv_scaled, highvar_meta, seed=fold, disp=False)

        filename = f"./../data/GEO/preprocessed/450K_Mvalues_atleast{atleast}_thresh{thresh}_samplewise_fold{fold}"
        print(f"...saving as {filename}...")
        with open(filename, "wb") as dill_file:
            pickle.dump([rest_Mv, holdout_Mv], dill_file, protocol=4)

        # rest_Mv.to_csv(f"./../data/GEO/preprocessed/450K_Mvalues_atleast{atleast}_thresh{thresh}_samplewise_fold{fold}_train.csv.gz", compression='gzip')
        # holdout_Mv.to_csv(f"./../data/GEO/preprocessed/450K_Mvalues_atleast{atleast}_thresh{thresh}_samplewise_fold{fold}_holdout.csv.gz", compression='gzip')
        print()
