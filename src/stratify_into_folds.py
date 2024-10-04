import dill
import random
from sklearn.model_selection import StratifiedGroupKFold    

atleast=2
Mv, meta, mapping = dill.load(open(f'./../data/GEO/preprocessed/450K_Mvalues_atleast{atleast}_samplewise', 'rb'))

random.seed(9)
groups = meta['series'].unique()
random.shuffle(groups)
meta = meta.set_index("series").loc[groups].reset_index().set_index("sample_id")
# print(meta.head())
Mv = Mv.loc[meta.index]

# sgkf = StratifiedGroupKFold(n_splits=3, shuffle=True, random_state=9) # if shuffle=True, then some folds with no training for some labels
sgkf = StratifiedGroupKFold(n_splits=3, shuffle=False, random_state=None)
sgkf.get_n_splits(Mv, meta['tissue_name'], meta['series'])
print(sgkf)

fold_Mvs = dict()
fold_selectors = dict()

for i, (train_index, test_index) in enumerate(sgkf.split(Mv, meta['tissue_name'], meta['series'])):
    print(f"fold {i}:")
    # print(f"  Train: index={len(train_index)}, group={meta['series'][train_index].unique()}")
    # print(f"  Test:  index={len(test_index)}, group={meta['series'][test_index].unique()}")
    
    rest_Mv = Mv.iloc[train_index]
    rest_meta = meta.iloc[train_index]
    holdout_Mv = Mv.iloc[test_index]
    holdout_meta = meta.iloc[test_index]
    
    print(f"train and validation shapes: {rest_Mv.shape, holdout_Mv.shape}")
    
    fold_Mvs[i] = rest_Mv, rest_meta, holdout_Mv, holdout_meta
    
print('...print if any overlap GSE between training and holdout...')
for fold, fold_data in fold_Mvs.items():
    print(f"fold{fold}: {set(fold_data[1]['series']).intersection(fold_data[3]['series'])}")
    
with open(f'./../data/GEO/preprocessed/450K_Mvalues_atleast{atleast}_samplewise_fold_Mvs', 'wb') as f:
    dill.dump(fold_Mvs, f)