import warnings
warnings.filterwarnings('ignore')
import joblib
import sys
sys.modules['sklearn.externals.joblib'] = joblib
sys.path.append('./../src/')
import dill
from sklearn.metrics import accuracy_score, precision_score, jaccard_score, f1_score
import pandas as pd
import networkx as nx
import utils
import pickle
import numpy as np
import random
from typing import Dict, Tuple, List
import os

# Configuration
RANDOM_SEED = 9
DATA_PATH = './../data/GEO'
SELECTION_FREQ_RANGE = "[0.8,0.2]"
THRESHOLD = 0.65

# Set random seeds
random.seed(RANDOM_SEED)
np.random.seed(RANDOM_SEED)

def load_fold_data() -> Dict:
    """Load cross-validation fold data."""
    try:
        with open(f'{DATA_PATH}/preprocessed/training_folds.dill', 'rb') as f:
            return pickle.load(f)
    except FileNotFoundError:
        raise FileNotFoundError("Training folds data file not found")
    except Exception as e:
        raise Exception(f"Error loading fold data: {str(e)}")

def load_models() -> Tuple[Dict, Dict, Dict]:
    """Load trained models and thresholds."""
    minipatch_location = f"{DATA_PATH}/minipatch/crossvalidation_selectors_{SELECTION_FREQ_RANGE}"
    clf_location = f"{DATA_PATH}/minipatch/crossvalidation_clfs_{SELECTION_FREQ_RANGE}"

    threshold = THRESHOLD
    fold_thresholds = {0:threshold, 1:threshold, 2:threshold}
    
    fold_selectors = dill.load(open(minipatch_location, 'rb'))
    fold_clfs = dill.load(open(clf_location, 'rb'))
    
    return fold_selectors, fold_clfs, fold_thresholds

def evaluate_fold(fold: str, 
                 rest_Mv: pd.DataFrame, 
                 rest_meta: pd.DataFrame,
                 holdout_Mv: pd.DataFrame, 
                 holdout_meta: pd.DataFrame,
                 selector, 
                 clf, 
                 threshold: float) -> Tuple[np.ndarray, List, List]:
    """Evaluate a single fold."""
    selection_freq = pd.DataFrame(selector.Pi_hat_last_k_, index=rest_Mv.columns)
    minipatch_probes = list(selection_freq[selection_freq[0] >= threshold].index)
    
    if len(minipatch_probes) == 0:
        print(f"Warning: No probes selected for fold {fold}")
        return None, None, None
        
    holdout_Mv_selected = holdout_Mv[minipatch_probes]
    
    holdout_multi = utils.propagate_parent(utils.subtree, holdout_meta, tissue_col='training.ID',outdict=False)
    holdout_mlb = utils.mlb.transform(holdout_multi['training.ID'].values)
    
    pred_prob = np.array([y_pred[:, 1] for y_pred in clf.predict_proba(holdout_Mv_selected.values)]).T
    pred_proba_bin = (pred_prob >= 0.5).astype(int)
    pred_proba_pred = [np.where(row)[0] for row in pred_proba_bin]
    pred_proba_pred = [tuple(utils.mlb.classes_[idx]) for idx in pred_proba_pred]
    pred = utils.mlb.transform(pred_proba_pred)
    
    return pred_prob, pred_proba_pred, utils.mlb.inverse_transform(holdout_mlb)

def main():
    # Load data and models
    fold_selectors, fold_clfs, fold_thresholds = load_models()
    fold_Mvs = load_fold_data()
    
    pred_res = pd.DataFrame(columns=['predict_proba', 'pred', 'true'])
    
    for fold, (rest_Mv, rest_meta, holdout_Mv, holdout_meta) in fold_Mvs.items():
        print(f"Fold {fold}:")
        print(f"Pre-selection shapes: {rest_Mv.shape, holdout_Mv.shape}")
        
        threshold = fold_thresholds[fold]
        print(f"Frequency threshold: {threshold}")
        
        pred_prob, pred_set, true_set = evaluate_fold(
            fold, rest_Mv, rest_meta, holdout_Mv, holdout_meta,
            fold_selectors[fold], fold_clfs[fold][threshold], threshold
        )
        
        if pred_prob is not None:
            for i, gse in enumerate(holdout_meta.index):
                pred_res.loc[f'f{fold}.{gse}'] = [pred_prob[i], pred_set[i], true_set[i]]
            
            # Calculate metrics
            holdout_multi = utils.propagate_parent(utils.subtree, holdout_meta, tissue_col='training.ID', outdict=False)
            holdout_mlb = utils.mlb.transform(holdout_multi['training.ID'].values)
            pred = utils.mlb.transform(pred_set)
            
            metrics = {
                'acc-samp': accuracy_score(holdout_mlb, pred),
                'jacc-samp': jaccard_score(holdout_mlb, pred, average="samples"),
                'prec-samp': precision_score(holdout_mlb, pred, average="samples"),
                'prec-tiss': utils.custom_macro_precision(holdout_mlb, pred),
                'f1-tiss': utils.custom_tissue_f1_score(holdout_mlb, pred)
            }
            
            print('\n'.join(f'{k}: {round(v,4)}' for k,v in metrics.items()))
            
            with open(f"{DATA_PATH}/minipatch/crossvalidation_results.txt", 'a') as f:
                f.write('\t'.join(str(round(v,4)) for v in metrics.values()) + '\n')
    
    with open(f"{DATA_PATH}/minipatch/crossvalidation_pred", 'wb') as f:
        pickle.dump(pred_res, f)

if __name__ == "__main__":
    main()
