import warnings
warnings.filterwarnings('ignore')
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import dill
# dill.load_session(f'base.db')
import sys
sys.path.append('./../src/')
# import utils
import pickle

from methylize import diff_meth_pos, volcano_plot, manhattan_plot
from statsmodels.stats.multitest import multipletests
from sklearn.preprocessing import MultiLabelBinarizer, LabelEncoder

atleast=2

Mv, meta, mapping = dill.load(open(f'./../data/GEO/preprocessed/450K_Mvalues_atleast{atleast}_samplewise', 'rb'))
le=LabelEncoder().fit(meta['tissue_name'].unique())
meta_le=le.transform(meta['tissue_name'].values)

print(meta_le.shape)
print(le.classes_)

def load_fold_data():
    """Load fold data."""
    with open('./../data/GEO/preprocessed/450K_Mvalues_atleast2_samplewise_fold_Mvs', 'rb') as f:
        return pickle.load(f)
    
fold_Mvs = load_fold_data()

for fold, (rest_Mv, rest_meta, holdout_Mv, holdout_meta) in fold_Mvs.items():
    print(f"fold: {fold}")
    
    # rest_Mv['series'] = rest_meta['series']
    res_per_tissue = dict()
    
    for tissue in sorted(holdout_meta['tissue_name'].unique()):
        print(f"tissue: {tissue}")
        
        rest_meta_tissue = rest_meta.copy()
        rest_meta_tissue['tissue'] = rest_meta_tissue['tissue_name']==tissue
        holdout_meta_tissue = holdout_meta.copy()
        holdout_meta_tissue['tissue'] = holdout_meta_tissue['tissue_name']==tissue
        
        res = diff_meth_pos(meth_data=rest_Mv, 
                            pheno_data=rest_meta_tissue, 
                            column='tissue',
                            regression_method="logistic", 
                            covariates='series',
                            export=False, 
                            verbose=False,
                           )
        # manhattan_plot(res, palette='default', save=False, array_type='450k', verbose=True)
        
        # clf = SVC(class_weight='balanced', kernel='linear', random_state=9)
        # clf.fit(rest_Mv.values, rest_meta_tissue)

        # true=holdout_meta_tissue.values
        # pred=clf.predict(holdout_Mv.values)

        # print(f'acc: {round(accuracy_score(true, pred),4)}')
        # print(f'prec: {round(precision_score(true, pred, average = "weighted"),4)}')

        # metrics.loc[tissue]['acc'] = round(accuracy_score(true, pred),4)
        # metrics.loc[tissue]['prec'] = round(precision_score(true, pred, average = "weighted"),4)
        # metrics.loc[tissue]['rec'] = round(recall_score(true, pred, average = "weighted"),4)
        
        interesting_probes1 = res[res['PValue'] <= 0.05].index
        print(interesting_probes1.shape)
        
        adjusted = multipletests(res.PValue, alpha=0.05)
        # pvalue_cutoff_y = -np.log10(adjusted[3])
        pvalue_cutoff_y = adjusted[3]
        # res['minuslog10value'] =  -np.log10(res['PValue'])
        # interesting_probes2 = res[res['minuslog10value'] >= pvalue_cutoff_y] #bonferoni correction for cutoff
        interesting_probes2 = res[res['PValue'] <= pvalue_cutoff_y] #bonferoni correction for cutoff
        print(interesting_probes2.shape)
        
        res['PValue_cutoff'] = pvalue_cutoff_y
        
        # probes_per_tissue[tissue] = interesting_probes2 #if only want significant
        res_per_tissue[tissue] = res
        
    # with open(f'diffmeth_probes_fold{fold}.pkl', 'wb') as f:
    with open(f'./../data/GEO/diffmeth/diffmeth_fold{fold}', 'wb') as f:
        pickle.dump(res_per_tissue, f)
        
