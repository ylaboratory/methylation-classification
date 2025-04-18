import dill
import os
from sklearn.model_selection import StratifiedGroupKFold
import random
import numpy as np
import pyarrow.parquet as pq
import pandas as pd
import networkx as nx

# Path to downloaded parquet files from huggingface
DOWNLOAD_PATH = './../download'
# Path to load/save the annotations/ontologies for easy load
ANNOTATION_PATH = './../annotation'
# Path to save the processed data for easy load
PREPROCESSED_PATH = './../data/GEO/preprocessed'
if not os.path.exists(PREPROCESSED_PATH):
    os.makedirs(PREPROCESSED_PATH, exist_ok=True)

# Load the data using pyarrow and save as dill for faster reading
files = {
    "Mv": f"{DOWNLOAD_PATH}/training_mvalues.parquet",
    "meta": f"{DOWNLOAD_PATH}/training_meta.parquet",
}
Mv, meta = [pq.read_table(path).to_pandas() for path in files.values()]

random.seed(9)
groups = meta['Dataset'].unique()
random.shuffle(groups)

meta = meta.reset_index().rename(columns={'index':'Sample'}).set_index("Dataset").loc[groups].reset_index().set_index("Sample")
Mv = Mv.T.loc[meta.index]

print(Mv.shape, meta.shape)

with open(f'{PREPROCESSED_PATH}/training.dill', 'wb') as f:
    dill.dump([Mv, meta], f)


# Load ontologies using edgelist and save as dill
full_ontology = nx.read_edgelist(f'{DOWNLOAD_PATH}/full_ontology.edgelist', create_using=nx.DiGraph)
training_ontology = nx.read_edgelist(f'{DOWNLOAD_PATH}/training_ontology.edgelist', create_using=nx.DiGraph)

with open(f'{ANNOTATION_PATH}/ontologies.dill', 'wb') as f:
    dill.dump([full_ontology, training_ontology], f)


# Stratify into folds grouped by GSE
sgkf = StratifiedGroupKFold(n_splits=3, shuffle=False, random_state=None)
sgkf.get_n_splits(Mv, meta['training.ID'], meta['Dataset'])

fold_Mvs = dict()
fold_selectors = dict()

for i, (train_index, test_index) in enumerate(sgkf.split(Mv, meta['training.ID'], meta['Dataset'])):
    print(f"fold {i}:")
    print(f"\ttrain: len={len(train_index)}, groups={len(meta['Dataset'][train_index].unique())}")
    print(f"\tvalidation:  len={len(test_index)}, groups={len(meta['Dataset'][test_index].unique())}")
    
    rest_Mv = Mv.iloc[train_index]
    rest_meta = meta.iloc[train_index]
    holdout_Mv = Mv.iloc[test_index]
    holdout_meta = meta.iloc[test_index]
    
    print(f"\ttrain and validation shapes: {rest_Mv.shape, holdout_Mv.shape}")
    
    fold_Mvs[i] = [rest_Mv, rest_meta, holdout_Mv, holdout_meta]
    
print('\n...print if any overlap GSE between training and validation (should be none)...')
for fold, fold_data in fold_Mvs.items():
    print(f"fold{fold}: {set(fold_data[1]['Dataset']).intersection(fold_data[3]['Dataset'])}")

with open(f'{PREPROCESSED_PATH}/training_folds.dill', 'wb') as f:
    dill.dump(fold_Mvs, f)