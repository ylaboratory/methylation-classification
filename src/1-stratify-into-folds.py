"""
1-stratify-into-folds.py
This script sets up the data for analysis.
    - input: 
        - parquet files of preprocessed M-values and metadata
    - output: 
        - pickled files of [training M-values, training metadata]
        - pickled files of [full ontology, training ontology]
        
** To replicate published results, use split indices in ./../annotation/cv_splits.dill
"""
import dill
import os
from sklearn.model_selection import StratifiedGroupKFold
import random
import numpy as np
import pyarrow.parquet as pq
import pandas as pd
import networkx as nx

random.seed(9)

### Relative paths ###
DOWNLOAD_PATH = './../download' # Path to downloaded parquet files from huggingface
ANNOTATION_PATH = './../annotation' # Path to load/save the annotations/ontologies for easy load
PREPROCESSED_PATH = './../data/GEO/preprocessed' # Path to save the processed data for easy load
if not os.path.exists(PREPROCESSED_PATH):
    os.makedirs(PREPROCESSED_PATH, exist_ok=True)
    
### Ontologies ###
# Load ontologies using edgelist and save as dill
full_ontology = nx.read_edgelist(f'{ANNOTATION_PATH}/full_ontology.edgelist', create_using=nx.DiGraph)
training_ontology = nx.read_edgelist(f'{ANNOTATION_PATH}/training_ontology.edgelist', create_using=nx.DiGraph)

with open(f'{ANNOTATION_PATH}/ontologies.dill', 'wb') as f:
    dill.dump([full_ontology, training_ontology], f)

### Training data into cross-validation splits ###
# Load the data using pyarrow and save as dill for faster reading
files = {
    "Mv": f"{DOWNLOAD_PATH}/train.parquet",
    "meta": f"{DOWNLOAD_PATH}/metadata.parquet",
}
Mv, meta = [pq.read_table(path).to_pandas() for path in files.values()]
Mv = Mv.set_index('Sample') if 'Sample' in Mv.columns else Mv
meta = meta.set_index('Sample') if 'Sample' in meta.columns else meta
meta = meta.loc[Mv.index]

# Random shuffle while being grouped by Dataset
groups = meta['Dataset'].unique()
random.shuffle(groups)
meta = meta.reset_index().rename(columns={'index':'Sample'}).set_index("Dataset").loc[groups].reset_index().set_index("Sample")
Mv = Mv.loc[meta.index]
print(Mv.shape, meta.shape)


# Stratify into folds grouped by GSE
sgkf = StratifiedGroupKFold(n_splits=3, shuffle=False, random_state=None)
sgkf.get_n_splits(Mv, meta['training.ID'], meta['Dataset'])

fold_Mvs = dict()
fold_selectors = dict()

for i, (train_index, test_index) in enumerate(sgkf.split(X=Mv, y=meta['training.ID'], groups=meta['Dataset'])):
    print(f"fold {i}:")
    print(f"\ttrain: len={len(train_index)}, groups={len(meta['Dataset'][train_index].unique())}")
    print(f"\tvalidation:  len={len(test_index)}, groups={len(meta['Dataset'][test_index].unique())}")
    
    rest_Mv = Mv.iloc[train_index]
    rest_meta = meta.iloc[train_index]
    holdout_Mv = Mv.iloc[test_index]
    holdout_meta = meta.iloc[test_index]
    
    print(f"\ttrain and validation shapes: {rest_Mv.shape, holdout_Mv.shape}")
    
    fold_Mvs[i] = [rest_Mv, rest_meta, holdout_Mv, holdout_meta]

# OR replicate using cv_split.dill
print(f"\nusing cv_splits...")
cv_split = dill.load(open(f"{ANNOTATION_PATH}/cv_split.dill", "rb"))
Mv = Mv.loc[cv_split['all']]
meta = meta.loc[Mv.index]
for fold in range(3):
    train_index = cv_split[fold]['train']
    test_index = cv_split[fold]['validation']
    
    print(f"fold {fold}:")
    print(f"\ttrain: len={len(train_index)}, groups={len(meta.loc[train_index]['Dataset'].unique())}")
    print(f"\tvalidation:  len={len(test_index)}, groups={len(meta.loc[test_index]['Dataset'].unique())}")
    
    rest_Mv = Mv.loc[train_index]
    rest_meta = meta.loc[train_index]
    holdout_Mv = Mv.loc[test_index]
    holdout_meta = meta.loc[test_index]
    
    print(f"\ttrain and validation shapes: {rest_Mv.shape, holdout_Mv.shape}")
    
    fold_Mvs[fold] = [rest_Mv, rest_meta, holdout_Mv, holdout_meta]

# Double check for leakage 
print('\n...print if any overlap GSE between training and validation (should be none)...')
for fold, fold_data in fold_Mvs.items():
    print(f"fold{fold}: {set(fold_data[1]['Dataset']).intersection(fold_data[3]['Dataset'])}")

# Save data as dill for faster loading
with open(f'{PREPROCESSED_PATH}/training.dill', 'wb') as f:
    dill.dump([Mv, meta], f)
with open(f'{PREPROCESSED_PATH}/training_folds.dill', 'wb') as f:
    dill.dump(fold_Mvs, f)