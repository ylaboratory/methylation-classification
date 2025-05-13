# Standard library
import sys
import random
from typing import Dict, Tuple
import gc

# Third-party
import warnings
warnings.filterwarnings('ignore')
import joblib
import numpy as np
import pandas as pd
from methylize import diff_meth_pos
from statsmodels.stats.multitest import multipletests

# Local
sys.path.append('./../src/')
sys.modules['sklearn.externals.joblib'] = joblib
import dill

DATA_PATH = f'./../data/GEO'
# Path to easy load the preprocessed data
PREPROCESSED_PATH = f'{DATA_PATH}/preprocessed'
# Path to save the differential methylation results
DIFFMETH_PATH  = f"{DATA_PATH}/diffmeth"

# Configuration
RANDOM_SEED = 9

# Set random seeds
random.seed(RANDOM_SEED)
np.random.seed(RANDOM_SEED)

def load_fold_data() -> Dict:
    """Load cross-validation fold data."""
    try:
        with open(f'{PREPROCESSED_PATH}/training_folds.dill', 'rb') as f:
            return dill.load(f)
    except FileNotFoundError:
        raise FileNotFoundError("Training folds data file not found")
    except Exception as e:
        raise Exception(f"Error loading fold data: {str(e)}")

def analyze_tissue_methylation(rest_Mv: pd.DataFrame, 
                             rest_meta: pd.DataFrame, 
                             tissue: str) -> pd.DataFrame:
    """Analyze differential methylation for a specific tissue."""
    rest_meta_tissue = rest_meta.copy()
    rest_meta_tissue['tissue'] = rest_meta_tissue['training.ID'] == tissue
    
    return diff_meth_pos(meth_data=rest_Mv,
                        pheno_data=rest_meta_tissue,
                        column='tissue',
                        regression_method="logistic",
                        covariates='Dataset',
                        export=False,
                        verbose=False)

def apply_multiple_testing(res: pd.DataFrame, alpha: float = 0.05) -> Tuple[pd.Index, float]:
    """Apply multiple testing correction and return significant probes."""
    adjusted = multipletests(res.PValue, alpha=alpha)
    pvalue_cutoff = adjusted[3]
    significant_probes = res[res['PValue'] <= pvalue_cutoff].index
    return significant_probes, pvalue_cutoff

def main():
    fold_Mvs = load_fold_data()
    
    for fold, (rest_Mv, rest_meta, holdout_Mv, holdout_meta) in {k: fold_Mvs[k] for k in [0,1,2]}.items():
        print(f"Fold: {fold}")
        res_per_tissue = {}
        
        for tissue in sorted(rest_meta['training.ID'].unique()):
            print(f"Tissue: {tissue}")
            
            # Analyze methylation
            res = analyze_tissue_methylation(rest_Mv, rest_meta, tissue)
            
            # Apply multiple testing correction
            significant_probes, pvalue_cutoff = apply_multiple_testing(res)
            print(f"Significant probes: {len(significant_probes)}")
            
            res['PValue_cutoff'] = pvalue_cutoff
            res_per_tissue[tissue] = res
        
        # Save results
        with open(f'{DIFFMETH_PATH}/diffmeth_fold{fold}', 'wb') as f:
            dill.dump(res_per_tissue, f)
            
        gc.collect()

if __name__ == "__main__":
    main()
