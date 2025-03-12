# taskset -c 120-129 python -u diffmeth_time.py | tee diffmeth_time.txt

import sys
import random
import time
from typing import Tuple
import warnings

warnings.filterwarnings('ignore')
import joblib
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests
from methylize import diff_meth_pos
import dill

# Configuration
RANDOM_SEED = 9
DATA_PATH = './../data/GEO'
MAX_TIME = 24 * 3600  # 24 hours in seconds

random.seed(RANDOM_SEED)
np.random.seed(RANDOM_SEED)


def load_data() -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Load main methylation data."""
    Mv_location = f"{DATA_PATH}/preprocessed/training.dill"
    print(f"loading Mv, meta from {Mv_location}")
    return dill.load(open(Mv_location, 'rb'))


def analyze_tissue_methylation(Mv: pd.DataFrame, meta: pd.DataFrame, tissue: str) -> pd.DataFrame:
    """Analyze differential methylation for a specific tissue."""
    meta_tissue = meta.copy()
    meta_tissue['tissue'] = meta_tissue['training.ID'] == tissue

    return diff_meth_pos(
        meth_data=Mv,
        pheno_data=meta_tissue,
        column='tissue',
        regression_method="logistic",
        covariates='Dataset',
        export=False,
        verbose=False,
    )


def apply_multiple_testing(res: pd.DataFrame, alpha: float = 0.05) -> Tuple[pd.Index, float]:
    """Apply multiple testing correction and return significant probes."""
    adjusted = multipletests(res.PValue, alpha=alpha)
    pvalue_cutoff = adjusted[3]
    significant_probes = res[res['PValue'] <= pvalue_cutoff].index
    return significant_probes, pvalue_cutoff


def main():
    with open('./time/sampled_samples.dill', 'rb') as f:
        sampled_dict = dill.load(f)

    # Load data
    Mv, meta = load_data()
    diffmeth_times = pd.DataFrame(columns=sampled_dict.keys())
    try:
        for num_samples, idx in sampled_dict.items():
            print(f"num_samples: {num_samples}")
            start_time = time.time()
    
            sampled_Mv = Mv.loc[idx]
            sampled_meta = meta.loc[idx]
    
            res_per_tissue = {}
            timeout = False  # Flag to track timeout
    
            try:
                for tissue in sorted(sampled_meta['training.ID'].unique()):
                    if time.time() - start_time > MAX_TIME:
                        print(f"  Timeout reached for {num_samples} samples. Skipping...")
                        timeout = True
                        break
    
                    res = analyze_tissue_methylation(sampled_Mv, sampled_meta, tissue)
                    significant_probes, pvalue_cutoff = apply_multiple_testing(res)
                    res['PValue_cutoff'] = pvalue_cutoff
                    res_per_tissue[tissue] = res
    
                execution_time = time.time() - start_time
                execution_time_str = f"{MAX_TIME:.4f}" if timeout else f"{execution_time:.4f}"
    
                print(f"  time: {execution_time_str}")
                diffmeth_times.loc[0, num_samples] = execution_time_str
                diffmeth_times.to_csv('./../data/GEO/time/diffmeth_times.csv')
    
            except Exception as e:
                print(f"Error in processing {num_samples} samples: {e}")
        
    except Exception as e:
        print(f"Critical failure: {e}")
    
    finally:
        # **GUARANTEE SAVE BEFORE EXITING**
        diffmeth_times.to_csv('./../data/GEO/time/diffmeth_times.csv')

if __name__ == "__main__":
    main()
