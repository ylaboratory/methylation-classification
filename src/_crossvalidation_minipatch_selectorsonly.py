import warnings
warnings.filterwarnings('ignore')
from mplearn.feature_selection._adaptive_stable_minipatch_selection import AdaSTAMPS
from mplearn.feature_selection.base_selector import DecisionTreeSelector
    
from sklearn.preprocessing import LabelEncoder
import dill
import pandas as pd
import sys
sys.path.append('./../src/')
import utils
import numpy as np
import pickle

atleast = 2

Mv_location = f"./../data/GEO/preprocessed/450K_Mvalues_atleast{atleast}_samplewise"
print(f"loading Mv, meta, mapping from {Mv_location}")
Mv, meta, mapping = dill.load(open(Mv_location, 'rb'))

filename = f"_crossvalidation" 
le=LabelEncoder().fit(meta['tissue_name'].unique())
meta_le=le.transform(meta['tissue_name'].values)

print(meta_le.shape)
print(le.classes_)

m_ratio = np.sqrt(Mv.shape[1])/Mv.shape[1]
n_ratio = np.sqrt(Mv.shape[0])/Mv.shape[0]

print("adaptive Exploitation & Exploration scheme...")
clf = DecisionTreeSelector(random_state=9)
selector = AdaSTAMPS(base_selector=clf,
                     minipatch_m_ratio=m_ratio,
                     minipatch_n_ratio=n_ratio,
                     random_state=9,
                     verbose=1,
                     sampling_options=None,
                     stopping_criteria_options=None)

def load_fold_data():
    """Load fold data."""
    with open('./../data/GEO/preprocessed/450K_Mvalues_atleast2_samplewise_fold_Mvs', 'rb') as f:
        return pickle.load(f)
    
fold_Mvs = load_fold_data()
fold_selectors = dict()

for i, (rest_Mv, rest_meta, holdout_Mv, holdout_meta) in fold_Mvs.items():
    print(f"fold {i}:")
    print(f"pre-feature selection train and validation shapes: {rest_Mv.shape, holdout_Mv.shape}")
    # print(f"series overlap, if any: {set(rest_meta['series'].values).intersection(set(holdout_meta['series'].values))}")
    # print(f"if any tissue exists only in training or only in validation:")
    
    rest_meta = le.transform(rest_meta['tissue_name'])   
    holdout_meta = le.transform(holdout_meta['tissue_name'])    
    
    fitted_selector = selector.fit(rest_Mv.values, rest_meta)
    
print(f"saving to ./../data/GEO/minipatch/minipatch{filename}_selectors")
with open(f"./../data/GEO/minipatch/minipatch{filename}_selectors", "wb") as dill_file:
    dill.dump(fold_selectors, dill_file, protocol=4)
