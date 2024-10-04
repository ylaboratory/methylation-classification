import warnings
warnings.filterwarnings('ignore')
import dill
# dill.load_session(f'base.db')
import sys
sys.path.append('./../src/')
import pickle

from methylize import diff_meth_pos
from statsmodels.stats.multitest import multipletests

from sklearn.preprocessing import MultiLabelBinarizer, LabelEncoder
import pandas as pd

atleast=2
thresh=2

Mv, meta, mapping = dill.load(open(f'./../data/GEO/preprocessed/450K_Mvalues_atleast{atleast}_samplewise', 'rb'))

le=LabelEncoder().fit(meta['tissue_name'].unique())
meta_le=le.transform(meta['tissue_name'].values)

from sklearn.preprocessing import StandardScaler
    
scaler = StandardScaler().fit(Mv.transpose().values)
Mv_scaled = scaler.transform(Mv.transpose().values)
Mv_scaled = pd.DataFrame(Mv_scaled.transpose(), index=Mv.index, columns=Mv.columns)

print(meta_le.shape)
print(le.classes_)

res_per_tissue = dict()

for tissue in sorted(meta['tissue_name'].unique()):
    print(f"tissue: {tissue}")
    
    meta_tissue = meta.copy()
    meta_tissue['tissue'] = meta_tissue['tissue_name']==tissue
    
    res = diff_meth_pos(meth_data=Mv_scaled, 
                        pheno_data=meta_tissue, 
                        column='tissue',
                        regression_method="logistic", 
                        covariates='series',
                        export=False, 
                        verbose=False,
                        )
    
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
with open(f'./../data/GEO/diffmeth/diffmeth_all', 'wb') as f:
    pickle.dump(res_per_tissue, f)
        
