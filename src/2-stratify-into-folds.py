import dill
from sklearn.model_selection import StratifiedGroupKFold
import random
import numpy as np
import pyarrow.parquet as pq
import pandas as pd

# Load the data using pyarrow for faster reading
file_path = './../preprocessed/training_Mvalues.parquet'
table = pq.read_table(file_path)

# Convert to pandas DataFrame if needed
Mv = table.to_pandas().set_index('probe')
Mv.columns = Mv.columns.str.split('_').str[0]
meta = pd.read_csv('./../annotation/training_meta.csv', header=0, index_col='Sample')

random.seed(9)
groups = meta['Dataset'].unique()
random.shuffle(groups)

meta = meta.reset_index().rename(columns={'index':'Sample'}).set_index("Dataset").loc[groups].reset_index().set_index("Sample")
Mv = Mv.T.loc[meta.index]

print(Mv.shape, meta.shape)

with open(f'./../data/GEO/preprocessed/training.dill', 'wb') as f:
    dill.dump([Mv, meta], f)

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

with open(f'./../data/GEO/preprocessed/training_folds.dill', 'wb') as f:
    dill.dump(fold_Mvs, f)