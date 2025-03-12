import warnings
warnings.filterwarnings('ignore')
from mplearn.feature_selection._adaptive_stable_minipatch_selection import AdaSTAMPS
from mplearn.feature_selection.base_selector import DecisionTreeSelector, SupportVectorMachineSelector
    
from sklearn.preprocessing import LabelEncoder
import dill
import pandas as pd
import sys
sys.path.append('./../src/')
import numpy as np

from methylize import diff_meth_pos
from statsmodels.stats.multitest import multipletests

times = list()

import time
start_time = time.time()
fold_time = start_time
    
atleast = 2

Mv_location = f"./../data/GEO/preprocessed/450K_Mvalues_atleast{atleast}_samplewise"
print(f"loading Mv, meta, mapping from {Mv_location}")
Mv, meta, mapping = dill.load(open(Mv_location, 'rb'))

filename = f"_default_{int(Mv.shape[1])}" 

with open('/grain/mk98/methyl/methylation-classification/data/GEO/sgkf/folds', 'rb') as f:
    fold_Mvs = dill.load(f)
    
with open('/grain/mk98/methyl/methylation-classification/data/GEO/sgkf/folds_sampled_indices', 'rb') as f:
    sampled_dict = dill.load(f)

for num_samples in sampled_dict.keys():
    print(f"num_samples: {num_samples}")
    num_samples_times = [num_samples]
    for i in range(3):
        print(f"fold {i}:")
        rest_Mv, rest_meta, holdout_Mv, holdout_meta = fold_Mvs[i]
        rest_Mv = rest_Mv.loc[sampled_dict[num_samples][i]]
        rest_meta = rest_meta.loc[sampled_dict[num_samples][i]]
        
        print(f"train and validation shapes: {rest_Mv.shape, holdout_Mv.shape}")
        # for tissue in meta.tissue_name.unique():
        #     tissue_rest_meta = rest_meta[rest_meta['tissue_name']==tissue]
        #     tissue_holdout_meta = holdout_meta[holdout_meta['tissue_name']==tissue]
        #     if len(tissue_rest_meta['series'].unique())==0 or len(tissue_holdout_meta['series'].unique())==0:
        #         print(f"tissue: {tissue}")
        #         print(f"training: {tissue_rest_meta['series'].unique()}")
        #         print(f"validation: {tissue_holdout_meta['series'].unique()}")
        # print()
        
        res_per_tissue = dict()
    
        for tissue in sorted(rest_meta['tissue_name'].unique()):
            print(f"tissue: {tissue}")
            
            rest_meta_tissue = rest_meta.copy()
            rest_meta_tissue['tissue'] = rest_meta_tissue['tissue_name']==tissue
            
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
        
        # End time
        previous_time = fold_time
        fold_time = time.time()

        # Calculate execution time in seconds
        execution_time = fold_time - previous_time
        num_samples_times.append(execution_time)
        print(f"fold {i} time: {execution_time:.4f}")
        
        print()
    
    all_time = time.time()

    # Calculate execution time in seconds
    execution_time = all_time - start_time
    num_samples_times.append(execution_time)
    print(f"all time: {execution_time:.4f}")

    with open(f'diffmeth_times_sampled.txt', 'a') as f:
        f.write('\n'.join([str(x) for x in num_samples_times]))
        f.write('\n')
        f.write('\n')
