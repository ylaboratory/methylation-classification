import warnings
warnings.filterwarnings('ignore')
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import dill
dill.load_session(f'base.db')
import sys
sys.path.append('./../src/')
import pickle

from methylize import diff_meth_pos
from statsmodels.stats.multitest import multipletests

times = list()

import time
start_time = time.time()
fold_time = start_time

Mv_location = f"./../data/GEO/preprocessed/450K_Mvalues_atleast2_samplewise"
print(f"loading Mv, meta, mapping from {Mv_location}")
Mv, meta, mapping = dill.load(open(Mv_location, 'rb'))

filename = f"_default_{int(Mv.shape[1])}" 
le=LabelEncoder().fit(meta['tissue_name'].unique())
meta_le=le.transform(meta['tissue_name'].values)

print(meta_le.shape)
print(le.classes_)

with open('/grain/mk98/methyl/methylation-classification/data/GEO/sgkf/folds', 'rb') as f:
    fold_Mvs = dill.load(f)

for i in range(3):
    fold = i
    print(f"fold {i}:")
    rest_Mv, rest_meta, holdout_Mv, holdout_meta = fold_Mvs[i]
    
    print(f"train and validation shapes: {rest_Mv.shape, holdout_Mv.shape}")
    for tissue in meta.tissue_name.unique():
        tissue_rest_meta = rest_meta[rest_meta['tissue_name']==tissue]
        tissue_holdout_meta = holdout_meta[holdout_meta['tissue_name']==tissue]
        # if len(tissue_rest_meta['series'].unique())==0 or len(tissue_holdout_meta['series'].unique())==0:
            # print(f"tissue: {tissue}")
            # print(f"training: {tissue_rest_meta['series'].unique()}")
            # print(f"validation: {tissue_holdout_meta['series'].unique()}")
    print()
    
    res_per_tissue = dict()
    
    for tissue in sorted(holdout_meta['tissue_name'].unique()):
        print(f"tissue: {tissue}")
        
        rest_meta_tissue = rest_meta.copy()
        rest_meta_tissue['tissue'] = rest_meta_tissue['tissue_name']==tissue
        holdout_meta_tissue = holdout_meta.copy()
        holdout_meta_tissue['tissue'] = holdout_meta_tissue['tissue_name']==tissue
        
        res = diff_meth_pos(meth_data=rest_Mv, 
                            pheno_data=rest_meta_tissue, 
                            column='tissue',
                            regression_method="logistic", 
                            covariates='series',
                            export=False, 
                            verbose=False,
                           )
        
        interesting_probes1 = res[res['PValue'] <= 0.05].index
        print(interesting_probes1.shape)
        
        adjusted = multipletests(res.PValue, alpha=0.05)
        pvalue_cutoff_y = adjusted[3]
        interesting_probes2 = res[res['PValue'] <= pvalue_cutoff_y] #bonferoni correction for cutoff
        print(interesting_probes2.shape)
        
        res['PValue_cutoff'] = pvalue_cutoff_y
        
        # probes_per_tissue[tissue] = interesting_probes2 #if only want significant
        res_per_tissue[tissue] = res
        
    # with open(f'diffmeth_probes_fold{fold}.pkl', 'wb') as f:
    # with open(f'./../data/GEO/diffmeth/diffmeth_atleast2_fold{fold}', 'wb') as f:
    #     pickle.dump(res_per_tissue, f)
        
    # End time
    previous_time = fold_time
    fold_time = time.time()

    # Calculate execution time in seconds
    execution_time = fold_time - previous_time
    times.append(execution_time)
    print(f"fold {i} time: {execution_time:.4f}")
    
    print()
    
all_time = time.time()

# Calculate execution time in seconds
execution_time = all_time - start_time
times.append(execution_time)
print(f"all time: {execution_time:.4f}")

with open(f'diffmeth_times.txt', 'a') as f:
    f.write('\n'.join([str(x) for x in times]))
    f.write('\n')