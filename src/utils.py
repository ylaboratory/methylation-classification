import warnings
warnings.filterwarnings('ignore')
import sklearn as sk
import matplotlib.pyplot as plt
import random
random.seed(9)
import networkx as nx
import numpy as np
np.random.seed(9)
import pandas as pd
from collections import Counter
from sklearn.metrics import f1_score
from sklearn.preprocessing import MultiLabelBinarizer
import dill
from sklearn.metrics import accuracy_score, precision_score, jaccard_score, f1_score

import obonet
# uberon=obonet.read_obo('/grain/mk98/methyl/methylation-classification/annotation/uberon_ext.obo')
# id_to_name = {id_:data.get('name') for id_, data in uberon.nodes(data=True)}
# name_to_id = {data.get('name'):id_ for id_, data in uberon.nodes(data=True)}
id_to_name = dill.load(open('./../annotation/id_to_name.dill', 'rb'))
name_to_id = dill.load(open('./../annotation/name_to_id.dill', 'rb'))

with open('./../annotation/ontologies.dill', 'rb') as f:
    _, subtree = dill.load(f)
mlb = MultiLabelBinarizer().fit([[node for node in subtree.nodes if node not in ['root', '']]])

def beta_to_mvalue(beta_df):
    if 'probe' in beta_df.columns:
        beta_values = beta_df.drop(columns=['probe'])
    else:
        beta_values = beta_df
    beta_values.replace(1, 1-1e-3, inplace=True)
    beta_values.replace(0, 1e-3, inplace=True)
    m_values = np.log2(beta_values / (1 - beta_values))
    if 'probe' in beta_df.columns: m_values.insert(0, 'probe', beta_df['probe'])
    return m_values

def generate_colors(targets):
    colors = list()
    for n in range(len(targets)):
        random.seed(n)
        color = "#%06x" % random.randint(0, 0xFFFFFF)
        colors+=[color]
    return colors

def custom_macro_precision(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    """Calculate custom macro precision score."""
    sample_score = precision_score(y_true, y_pred, average=None)
    sample_weight = y_true.sum(axis=0)
    valid_indices = sample_weight != 0
    return np.median(sample_score[valid_indices])

def custom_tissue_f1_score(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    """Calculate custom tissue F1 score."""
    sample_score = f1_score(y_true, y_pred, average=None)
    sample_weight = y_true.sum(axis=0)
    valid_indices = sample_weight != 0
    return np.median(sample_score[valid_indices])

def propagate_parent(subtree, meta, tissue_col='tissue_name', outdict=False):
    print(f"propagating with {tissue_col}")
    # print(meta.columns)
    multioutput=pd.DataFrame(index=meta.index, columns=[tissue_col])
    output_dict=dict()
    for idx in meta.index:
        sample=meta[tissue_col].loc[idx]
        # subtree_parents=[x for x in nx.ancestors(subtree, sample) if x not in ['root', '']]
        subtree_parents=[x for x in nx.ancestors(subtree, sample) if x not in ['root', '']]
        subtree_parents.append(sample)
        multioutput.at[idx, tissue_col]=set(subtree_parents)
        if sample not in output_dict.keys(): output_dict[sample]=set(subtree_parents)
    if outdict:
        return multioutput, output_dict
    else:
        return multioutput