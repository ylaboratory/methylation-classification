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

times = list()

import time
start_time = time.time()
fold_time = start_time
    
atleast = 2

Mv_location = f"./../data/GEO/preprocessed/450K_Mvalues_atleast{atleast}_samplewise"
print(f"loading Mv, meta, mapping from {Mv_location}")
Mv, meta, mapping = dill.load(open(Mv_location, 'rb'))

filename = f"_default_{int(Mv.shape[1])}" 
le=LabelEncoder().fit(meta['tissue_name'].unique())
meta_le=le.transform(meta['tissue_name'].values)

print(meta_le.shape)
print(le.classes_)

mode = "ee"
E = 10
active_threshold = 0.5
gamma_min = 0.5
gamma_max = 0.5
gamma_len = 200
m_ratio = np.sqrt(Mv.shape[1])/Mv.shape[1]
n_ratio = np.sqrt(Mv.shape[0])/Mv.shape[0]

tau_u = 10
tau_l = 1
num_last_iterations = 10 
max_k = None
selection_frequency_threshold = 0.45

print("adaptive Exploitation & Exploration scheme...")
clf = DecisionTreeSelector(random_state=9)
selector = AdaSTAMPS(base_selector=clf,
                     minipatch_m_ratio=m_ratio,
                     minipatch_n_ratio=n_ratio,
                     random_state=9,
                     verbose=1,
                    sampling_options={'mode': mode,
                                     'E': E, #default
                                     'active_set_thr': active_threshold, #default
                                     'gamma_min': gamma_min, #default
                                     'gamma_max': gamma_max, #default
                                     'gamma_len': gamma_len}, #default
                    stopping_criteria_options={'tau_u':tau_u,
                                               'tau_l':tau_l,
                                               'num_last_iterations':num_last_iterations, #default = 100
                                               'max_k':max_k})

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
        
        rest_meta = le.transform(rest_meta['tissue_name'])   
        # holdout_meta = le.transform(holdout_meta['tissue_name'])    
        
        fitted_selector = selector.fit(rest_Mv.values, rest_meta)
        Mv_new = fitted_selector.transform(rest_Mv.values, pi_thr=selection_frequency_threshold)
        print(f"selected probes: {Mv_new.shape[1]}")
        
        # print(f"saving to ./../data/GEO/minipatch/minipatch{filename}_fold{i}_selector")
        # with open(f"./../data/GEO/minipatch/minipatch{filename}_fold{i}_selector", "wb") as dill_file:
        #     dill.dump(fitted_selector, dill_file, protocol=4)
        
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

    with open(f'minipatch_times_sampled.txt', 'a') as f:
        f.write('\n'.join([str(x) for x in num_samples_times]))
        f.write('\n')
        f.write('\n')
