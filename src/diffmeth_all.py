# Standard library
import sys
import random
from typing import Dict, Tuple, List

# Third-party
import warnings
warnings.filterwarnings('ignore')
import joblib
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from statsmodels.stats.multitest import multipletests
from methylize import diff_meth_pos

# Local
sys.path.append('./../src/')
sys.modules['sklearn.externals.joblib'] = joblib
import utils
import dill

# Configuration
RANDOM_SEED = 9
DATA_PATH = './../data/GEO'

# Set random seeds
random.seed(RANDOM_SEED)
np.random.seed(RANDOM_SEED)

def load_data() -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Load main methylation data."""
    Mv_location = f"{DATA_PATH}/preprocessed/training.dill"
    print(f"loading Mv, meta from {Mv_location}")
    return dill.load(open(Mv_location, 'rb'))

def scale_methylation_data(Mv: pd.DataFrame) -> pd.DataFrame:
    """Scale methylation values using StandardScaler."""
    print("scaling mvalues")
    scaler = StandardScaler().fit(Mv.transpose().values)
    Mv_scaled = scaler.transform(Mv.transpose().values)
    return pd.DataFrame(Mv_scaled.transpose(), index=Mv.index, columns=Mv.columns)

def analyze_tissue_methylation(Mv_scaled: pd.DataFrame, 
                             meta: pd.DataFrame, 
                             tissue: str) -> pd.DataFrame:
    """Analyze differential methylation for a specific tissue."""
    meta_tissue = meta.copy()
    meta_tissue['tissue'] = meta_tissue['training.ID'] == tissue
    
    return diff_meth_pos(meth_data=Mv_scaled,
                        pheno_data=meta_tissue,
                        column='tissue',
                        regression_method="logistic",
                        covariates='series',
                        export=False,
                        verbose=False)

def apply_multiple_testing(res: pd.DataFrame, alpha: float = 0.05) -> Tuple[pd.Index, float]:
    """Apply multiple testing correction and return significant probes."""
    adjusted = multipletests(res.PValue, alpha=alpha)
    pvalue_cutoff = adjusted[3]
    significant_probes = res[res['PValue'] <= pvalue_cutoff].index
    return significant_probes, pvalue_cutoff

def main():
    # Load and scale data
    Mv, meta = load_data()
    # Mv_scaled = scale_methylation_data(Mv)
    
    res_per_tissue = {}
    for tissue in sorted(meta['training.ID'].unique()):
        print(f"Tissue: {tissue}")
        
        # Analyze methylation
        res = analyze_tissue_methylation(Mv, meta, tissue)
        
        # Apply multiple testing correction
        significant_probes, pvalue_cutoff = apply_multiple_testing(res)
        print(f"Significant probes: {len(significant_probes)}")
        
        res['PValue_cutoff'] = pvalue_cutoff
        res_per_tissue[tissue] = res
    
    # Save results
    with open(f'{DATA_PATH}/diffmeth/diffmeth_all', 'wb') as f:
        dill.dump(res_per_tissue, f)

if __name__ == "__main__":
    main()
