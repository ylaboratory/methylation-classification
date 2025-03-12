# Standard library
import sys
import random
from typing import Dict, Tuple, List

# Third-party
import warnings
import joblib
import numpy as np
import pandas as pd
import networkx as nx
from sklearn.svm import SVC
from statsmodels.stats.multitest import multipletests

# Local
sys.path.append('./../src/')
sys.modules['sklearn.externals.joblib'] = joblib
import utils
import dill

warnings.filterwarnings('ignore')

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

def get_significant_probes(res_per_tissue: Dict, tissue: str) -> List[str]:
    """Get significant probes for a tissue using multiple testing correction."""
    res = res_per_tissue[tissue]
    adjusted = multipletests(res.PValue, alpha=0.05)
    pvalue_cutoff = adjusted[3]
    return list(res[res['PValue'] <= pvalue_cutoff].index)

def train_tissue_classifier(rest_Mv: pd.DataFrame, 
                          rest_meta: pd.DataFrame,
                          probes: List[str],
                          tissue: str) -> Tuple[SVC, pd.DataFrame]:
    """Train classifier for a specific tissue."""
    rest_meta['tissue_bool'] = rest_meta['training.ID'] == tissue
    rest_probe = rest_Mv[probes]
    
    clf = SVC(class_weight='balanced', kernel='linear', random_state=RANDOM_SEED, probability=True)
    return clf.fit(rest_probe.values, rest_meta['tissue_bool']), rest_probe

def predict_tissue(clf: SVC, 
                  holdout_probe: pd.DataFrame, 
                  holdout_meta: pd.DataFrame,
                  tissue: str) -> np.ndarray:
    """Make predictions for a tissue."""
    holdout_meta['tissue_bool'] = holdout_meta['training.ID'] == tissue
    return clf.predict_proba(holdout_probe.values)

def main():
    # Load data
    Mv, meta = load_data()
    print(f"Mv, meta: {Mv.shape}, {meta.shape}")
    fold_Mvs = load_fold_data()
    
    # Initialize results dataframe
    all_pred_res = pd.DataFrame(columns=meta['training.ID'].unique())
    
    for fold, (rest_Mv, rest_meta, holdout_Mv, holdout_meta) in fold_Mvs.items():
        print(f"Fold: {fold}")
        
        # Load differential methylation results
        res_per_tissue = load_diffmeth_results(fold)
        
        # Process each tissue
        clfs = {}
        pred_res = pd.DataFrame(index=[f'f{fold}.{gse}' for gse in holdout_meta.index],
                              columns=res_per_tissue.keys())
        
        for tissue in res_per_tissue.keys():
            print(f"Tissue: {tissue}")
            
            # Get significant probes
            probes = get_significant_probes(res_per_tissue, tissue)
            print(f"Number of probes: {len(probes)}")
            
            # Train classifier
            clf, rest_probe = train_tissue_classifier(rest_Mv, rest_meta, probes, tissue)
            holdout_probe = holdout_Mv[probes]
            
            # Make predictions
            pred_prob = predict_tissue(clf, holdout_probe, holdout_meta, tissue)
            pred_res[tissue] = [tuple(x) for x in pred_prob]
            
            clfs[tissue] = clf
        
        # Save results
        with open(f"{DATA_PATH}/diffmeth/clf_ovr_fold{fold}", 'wb') as f:
            dill.dump(clfs, f)
            
        all_pred_res = pd.concat([all_pred_res, pred_res])
    
    # Save final predictions
    with open(f'{DATA_PATH}/diffmeth/diffmeth_ovr.pkl', 'wb') as f:
        dill.dump(all_pred_res, f)

if __name__ == "__main__":
    main()
