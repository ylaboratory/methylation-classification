import warnings
warnings.filterwarnings('ignore')
import joblib
import sys
sys.modules['sklearn.externals.joblib'] = joblib
sys.path.append('./../src/')
import dill
from sklearn.svm import SVC
from sklearn.multioutput import MultiOutputClassifier
from sklearn.metrics import accuracy_score, precision_score, jaccard_score, f1_score
from sklearn.preprocessing import MultiLabelBinarizer, LabelEncoder

import pandas as pd
import networkx as nx
import utils
import pickle
import numpy as np
import random

random.seed(9)
np.random.seed(9)

date = "sep2024"
atleast = 2
subtree = nx.read_multiline_adjlist(path=f"uberon_{date}_atleast{atleast}_adjlist",create_using=nx.MultiDiGraph)
mlb=MultiLabelBinarizer().fit([[utils.id_to_name[node] for node in subtree.nodes]])

Mv_location = f"./../data/GEO/preprocessed/450K_Mvalues_atleast{atleast}_samplewise"
print(f"loading Mv, meta, mapping from {Mv_location}")
Mv, meta, mapping = dill.load(open(Mv_location, 'rb'))

le=LabelEncoder().fit(meta['tissue_name'].unique())
meta_le=le.transform(meta['tissue_name'].values)

print(meta_le.shape)
print(le.classes_)

filename = "_crossvalidation"
minipatch_location = f"./../data/GEO/minipatch/minipatch{filename}_selectors"
print(f"loading minipatch selectors from {minipatch_location}")
fold_selectors = dill.load(open(minipatch_location, 'rb'))

clf_location = f"./../data/GEO/minipatch/multilabel{filename}_clfs"
print(f"loading minipatch selectors from {clf_location}")
fold_clfs = dill.load(open(clf_location, 'rb'))

pred_res = pd.DataFrame(columns = ['predict_proba', 'pred', 'true'])

def load_fold_data():
    """Load fold data."""
    with open('./../data/GEO/preprocessed/450K_Mvalues_atleast2_samplewise_fold_Mvs', 'rb') as f:
        return pickle.load(f)

def custom_macro_precision(y_true, y_pred):
    sample_score = precision_score(y_true, y_pred, average=None)
    sample_weight = y_true.sum(axis=0)
    valid_indices = sample_weight != 0
    return np.mean(sample_score[valid_indices])

def custom_tissue_f1_score(y_true, y_pred):
    sample_score = f1_score(y_true, y_pred, average=None)
    sample_weight = y_true.sum(axis=0)
    valid_indices = sample_weight != 0
    return np.mean(sample_score[valid_indices])
    
fold_Mvs = load_fold_data()
selection_frequency_threshold = 0.57
print(f"selection frequency threshold: {selection_frequency_threshold}")

for fold, (rest_Mv, rest_meta, holdout_Mv, holdout_meta) in fold_Mvs.items():
    print(f"fold {fold}:")
    
    print(f"pre-selection train and validation shapes: {rest_Mv.shape, holdout_Mv.shape}")
    print(f"series overlap, if any: {set(rest_meta['series'].values).intersection(set(holdout_meta['series'].values))}")
        
    fitted_selector = fold_selectors[fold]
    selection_freq = pd.DataFrame(fitted_selector.Pi_hat_last_k_, index=Mv.columns)
    minipatch_probes = list(selection_freq[selection_freq[0]>=selection_frequency_threshold].index)
    rest_Mv = rest_Mv[minipatch_probes]
    holdout_Mv = holdout_Mv[minipatch_probes]
    
    print(f"post-selection train and validation shapes: {rest_Mv.shape, holdout_Mv.shape}")
    
    rest_multi=utils.propagate_parent(subtree, rest_meta, outdict=False)
    rest_mlb=mlb.transform(rest_multi['tissue_name'].values)
    holdout_multi=utils.propagate_parent(subtree, holdout_meta, outdict=False)
    holdout_mlb=mlb.transform(holdout_multi['tissue_name'].values)

    clf = fold_clfs[fold]

    true=holdout_mlb
    pred=clf.predict(holdout_Mv.values)
    pred_prob = np.transpose([y_pred[:, 1] for y_pred in clf.predict_proba(holdout_Mv.values)])
    pred_set = mlb.inverse_transform(pred)
    true_set = mlb.inverse_transform(true)
    for i, gse in enumerate(holdout_meta.index):
        pred_res.loc[f'f{i}.{gse}'] = [pred_prob[i], pred_set[i], true_set[i]]

    accs = round(accuracy_score(true, pred),4)
    jaccs = round(jaccard_score(true, pred, average="samples"),4)    
    precs = round(precision_score(true, pred, average="samples"),4)
    prect = round(custom_macro_precision(true, pred), 4)
    f1 = round(custom_tissue_f1_score(true, pred),4)
    
    print(f'acc-samp: {accs}')
    print(f'jacc-samp: {jaccs}')
    print(f'prec-samp: {precs}')
    print(f'prec-tiss: {prect}')
    print(f'f1-tiss: {f1}')
    print()
    
    with open(f"./../data/GEO/minipatch/multilabel{filename}_results.txt", 'a') as f:
        f.write('\t'.join([str(accs), 
                          str(jaccs),
                          str(precs),
                          str(prect),
                          str(f1)]) + '\n')
        
    # fold_clfs[fold] = clf
        
# with open(f"./../data/GEO/minipatch/multilabel{filename}_clfs", 'wb') as f:
#     pickle.dump(fold_clfs, f)
    
with open(f"./../data/GEO/minipatch/multilabel{filename}_pred", 'wb') as f:
    pickle.dump(pred_res, f)