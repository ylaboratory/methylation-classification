import warnings
warnings.filterwarnings('ignore')
import joblib
import sys
sys.modules['sklearn.externals.joblib'] = joblib
sys.path.append('./../src/')
import dill
# dill.load_session(f'base.db')
from sklearn.svm import SVC
from statsmodels.stats.multitest import multipletests

import pickle
import random
import numpy as np
import pandas as pd
import networkx as nx
from sklearn.preprocessing import MultiLabelBinarizer
from utils import id_to_name
random.seed(9)
np.random.seed(9)

date = "sep2024"
atleast = 2
subtree = nx.read_multiline_adjlist(path=f"uberon_{date}_atleast{atleast}_adjlist",create_using=nx.MultiDiGraph)
mlb=MultiLabelBinarizer().fit([[id_to_name[node] for node in subtree.nodes]])

import dill
atleast=2
    
Mv, meta, mapping = dill.load(open(f'./../data/GEO/preprocessed/450K_Mvalues_atleast{atleast}_samplewise', 'rb'))

import pandas as pd
import numpy as np

all_pred_res = pd.DataFrame(columns = meta['tissue_name'].unique())
        
def load_fold_data():
    """Load fold data."""
    with open('./../data/GEO/preprocessed/450K_Mvalues_atleast2_samplewise_fold_Mvs', 'rb') as f:
        return pickle.load(f)
    
fold_Mvs = load_fold_data()


for fold, (rest_Mv, rest_meta, holdout_Mv, holdout_meta) in fold_Mvs.items():
    print(f"fold: {fold}")

    with open(f"./../data/GEO/diffmeth/diffmeth_fold{fold}", 'rb') as f:
        res_per_tissue = pickle.load(f)
    
    probe_per_tissue = dict()
    for tissue in res_per_tissue.keys():
        res = res_per_tissue[tissue]
        adjusted = multipletests(res.PValue, alpha=0.05)
        pvalue_cutoff_y = adjusted[3]
        interesting_probes2 = list(res[res['PValue'] <= pvalue_cutoff_y].index) #bonferoni correction for cutoff
    
        probe_per_tissue[tissue] = interesting_probes2
    
    pred_res = pd.DataFrame(index= [f'f{fold}.{gse}' for gse in holdout_meta.index], 
                        columns=probe_per_tissue.keys()
                       )
    
    clfs = dict()
    for tissue, probe in probe_per_tissue.items():
        print(tissue)
        
        rest_meta['tissue_bool'] = rest_meta['tissue_name']==tissue
        holdout_meta['tissue_bool'] = holdout_meta['tissue_name']==tissue
        
        rest_probe = rest_Mv[probe]
        holdout_probe = holdout_Mv[probe]
        print(rest_probe.shape, holdout_probe.shape)

        clf = SVC(class_weight='balanced', kernel='linear', random_state=9, probability=True)
        clf.fit(rest_probe.values, rest_meta['tissue_bool'])

        true=holdout_meta['tissue_bool']
        pred=clf.predict(holdout_probe.values)
        pred_prob = clf.predict_proba(holdout_probe.values)
        
        pred_res[tissue] = [tuple(x) for x in pred_prob]
            
        clfs[tissue]=clf
  
    with open(f"./../data/GEO/diffmeth/clf_ovr_fold{fold}", 'wb') as f:
        pickle.dump(clfs, f)
        
    all_pred_res = pd.concat([all_pred_res, pred_res])
        
with open(f'./../data/GEO/diffmeth/diffmeth_ovr.pkl', 'wb') as f:
    pickle.dump(all_pred_res, f)
        
