import warnings
warnings.filterwarnings('ignore')
import joblib
import sys
sys.modules['sklearn.externals.joblib'] = joblib
sys.path.append('./../src/')
import dill
from statsmodels.stats.multitest import multipletests
import pandas as pd
import numpy as np
import random
from typing import Dict, Tuple, List

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

def load_fold_data() -> Dict:
    """Load cross-validation fold data."""
    try:
        with open(f'{DATA_PATH}/preprocessed/training_folds.dill', 'rb') as f:
            return dill.load(f)
    except FileNotFoundError:
        raise FileNotFoundError("Training folds data file not found")
    except Exception as e:
        raise Exception(f"Error loading fold data: {str(e)}")

def load_diffmeth_results(fold: str) -> Dict:
    """Load differential methylation results for a fold."""
    with open(f"{DATA_PATH}/diffmeth/diffmeth_fold{fold}", 'rb') as f:
        return dill.load(f)

def get_significant_probes(res_per_tissue: Dict) -> Dict[str, List[str]]:
    """Get significant probes for each tissue using multiple testing correction."""
    probe_per_tissue = {}
    for tissue, tissue_res in res_per_tissue.items():
        adjusted = multipletests(tissue_res.PValue, alpha=0.05)
        pvalue_cutoff = adjusted[3]
        probe_per_tissue[tissue] = list(tissue_res[tissue_res['PValue'] <= pvalue_cutoff].index)
    return probe_per_tissue

def calculate_correlation(sample_Mv: pd.Series, 
                        rest_Mv: pd.DataFrame, 
                        rest_meta: pd.DataFrame,
                        tissue: str, 
                        probes: List[str]) -> float:
    """Calculate correlation for a sample with tissue-specific probes."""
    tissue_data = rest_Mv[rest_meta['training ID'] == tissue][probes]
    return np.corrcoef(sample_Mv[probes], tissue_data)[1:, 0].mean()

def main():
    # Load data
    Mv, meta = load_data()
    fold_Mvs = load_fold_data()
    
    # Initialize results
    res = pd.DataFrame(columns=['corr', 'pred', 'true', 'dict'])
    
    for fold, (rest_Mv, rest_meta, holdout_Mv, holdout_meta) in fold_Mvs.items():
        print(f"Fold: {fold}")
        
        # Get significant probes for each tissue
        res_per_tissue = load_diffmeth_results(fold)
        probe_per_tissue = get_significant_probes(res_per_tissue)
        
        # Process each sample
        for i, sample in enumerate(holdout_Mv.index):
            if i % 20 == 0:
                print(i)
                
            sample_Mv = holdout_Mv.loc[sample]
            sample_res = {tissue: calculate_correlation(sample_Mv, rest_Mv, rest_meta, tissue, tissue_probes) 
                         for tissue, tissue_probes in probe_per_tissue.items()}
            
            sample_pred = max(sample_res, key=sample_res.get)
            res.loc[f"f{fold}.{sample}"] = [
                sample_res[sample_pred],
                sample_pred,
                holdout_meta.loc[sample]['training ID'],
                sample_res
            ]
    
    # Save results
    with open(f'{DATA_PATH}/diffmeth/diffmeth_corr.pkl', 'wb') as f:
        dill.dump(res, f)

if __name__ == "__main__":
    main()
