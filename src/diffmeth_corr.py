import pickle
import dill
from statsmodels.stats.multitest import multipletests
atleast=2
    
Mv, meta, mapping = dill.load(open(f'./../data/GEO/preprocessed/450K_Mvalues_atleast{atleast}_samplewise', 'rb'))

import pandas as pd
import numpy as np
res = pd.DataFrame(columns = ['corr', 'pred', 'true', 'dict'])

import numpy as np
        
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
        tissue_res = res_per_tissue[tissue]
        adjusted = multipletests(tissue_res.PValue, alpha=0.05)
        pvalue_cutoff_y = adjusted[3]
        interesting_probes2 = list(tissue_res[tissue_res['PValue'] <= pvalue_cutoff_y].index) #bonferoni correction for cutoff
    
        probe_per_tissue[tissue] = interesting_probes2

    for i, sample in enumerate(holdout_Mv.index):
        if i%20==0: print(i)
        sample_Mv = holdout_Mv.loc[sample]
        sample_res = {tissue:np.corrcoef(sample_Mv[tissue_probes], rest_Mv[rest_meta['tissue_name']==tissue][tissue_probes])[1:,0].mean() for tissue, tissue_probes in probe_per_tissue.items()}
        sample_pred = max(sample_res, key=sample_res.get)

        res.loc[f"f{fold}.{sample}"] = [
            sample_res[sample_pred], 
            sample_pred, 
            holdout_meta.loc[sample]['tissue_name'],
            sample_res
            ]
        
with open(f'./../data/GEO/diffmeth/diffmeth_corr.pkl', 'wb') as f:
    pickle.dump(res, f)