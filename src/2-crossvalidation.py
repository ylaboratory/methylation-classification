# Standard library
import sys, os
import time
import random
from typing import Dict, Tuple, List
import warnings
warnings.filterwarnings('ignore')

# Third-party
import joblib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.svm import SVC
from sklearn.multioutput import MultiOutputClassifier
from sklearn.metrics import accuracy_score, precision_score, jaccard_score
from sklearn.preprocessing import LabelEncoder
from mplearn.feature_selection._adaptive_stable_minipatch_selection import AdaSTAMPS
from mplearn.feature_selection.base_selector import DecisionTreeSelector

# Local
sys.path.append('./../src/')
sys.modules['sklearn.externals.joblib'] = joblib
import utils
import dill

DATA_PATH = f'./../data/GEO'
# Path to easy load the preprocessed data
PREPROCESSED_PATH = f'{DATA_PATH}/preprocessed'
# Path to save the minipatch results
MINIPATCH_PATH  = f"{DATA_PATH}/minipatch"
if not os.path.exists(MINIPATCH_PATH):
    os.makedirs(MINIPATCH_PATH, exist_ok=True)

# Configuration
RANDOM_SEED = 9
SELECTION_RANGE = {
    'begin': 0.8,
    'end': 0.2,
    'step': -0.1
}
KNEE_DROP_THRESHOLD = 0.05

# Set random seeds
random.seed(RANDOM_SEED)
np.random.seed(RANDOM_SEED)

def timestamped_print(msg: str):
    """Helper function to print messages with a timestamp."""
    print(f"{time.strftime('%Y-%m-%d %H:%M:%S')} - {msg}")

def load_data() -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Load main methylation data."""
    Mv_location = f'{PREPROCESSED_PATH}/training.dill'
    timestamped_print(f"loading Mv, meta from {Mv_location}")
    return dill.load(open(Mv_location, 'rb'))

def load_fold_data() -> Dict:
    """Load cross-validation fold data."""
    try:
        with open(f'{PREPROCESSED_PATH}/training_folds.dill', 'rb') as f:
            return dill.load(f)
    except FileNotFoundError:
        raise FileNotFoundError("Training folds data file not found")
    except Exception as e:
        raise Exception(f"Error loading fold data: {str(e)}")

def initialize_selector(rest_Mv: pd.DataFrame) -> AdaSTAMPS:
    """Initialize the AdaSTAMPS selector."""
    m_ratio = np.sqrt(rest_Mv.shape[1])/rest_Mv.shape[1]
    n_ratio = np.sqrt(rest_Mv.shape[0])/rest_Mv.shape[0]
    
    clf = DecisionTreeSelector(random_state=RANDOM_SEED)
    return AdaSTAMPS(base_selector=clf,
                     minipatch_m_ratio=m_ratio,
                     minipatch_n_ratio=n_ratio,
                     random_state=RANDOM_SEED,
                     verbose=1,
                     )

def evaluate_threshold(fitted_selector, threshold, rest_Mv, holdout_Mv, rest_mlb, holdout_mlb):
    """Evaluate performance for a given threshold with tissue-specific metrics."""
    timestamped_print('  Evaluating thresholds...')
    selection_freq = pd.DataFrame(fitted_selector.Pi_hat_last_k_, index=rest_Mv.columns)
    minipatch_probes = list(selection_freq[selection_freq[0] >= threshold].index)
    
    if len(minipatch_probes) == 0:
        return 0.0, 0, None
        
    rest_Mv_selected = rest_Mv[minipatch_probes]
    holdout_Mv_selected = holdout_Mv[minipatch_probes]
    
    clf = MultiOutputClassifier(SVC(class_weight='balanced', kernel='linear', random_state=RANDOM_SEED, probability=True))
    clf.fit(rest_Mv_selected.values, rest_mlb)
    # pred = clf.predict(holdout_Mv_selected.values)
    pred_prob = np.array([y_pred[:, 1] for y_pred in clf.predict_proba(holdout_Mv_selected.values)]).T
    pred_proba_bin = (pred_prob >= 0.5).astype(int)
    pred_proba_pred = [np.where(row)[0] for row in pred_proba_bin]
    pred_proba_pred = [tuple(utils.mlb.classes_[idx]) for idx in pred_proba_pred]
    pred = utils.mlb.transform(pred_proba_pred)
    
    # Calculate tissue-specific precision
    tissue_precision = precision_score(holdout_mlb, pred, average=None)
    tissue_names = utils.mlb.classes_
    
    # Replace zero precision with None where no samples exist
    sample_counts = holdout_mlb.sum(axis=0)
    tissue_metrics = {tissue: None if count == 0 else score 
                     for tissue, score, count in zip(tissue_names, tissue_precision, sample_counts)}
    
    # Add overall metrics
    tissue_metrics.update({
        'num_probes': len(minipatch_probes),
        'accuracy': accuracy_score(holdout_mlb, pred),
        'jaccard': jaccard_score(holdout_mlb, pred, average="samples"),
        'custom_f1': utils.custom_tissue_f1_score(holdout_mlb, pred)
    })

    timestamped_print(f"  acc-samp {round(accuracy_score(holdout_mlb, pred),4)}|jacc-samp {round(jaccard_score(holdout_mlb, pred, average='samples'),4)}|prec-tiss {round(np.median(tissue_precision),4)}|f1-tissue {round(utils.custom_tissue_f1_score(holdout_mlb, pred),4)}")
    
    return tissue_metrics, len(minipatch_probes), clf

def find_all_knee_points(fold_results, fold, thresholds, drop_threshold: float = KNEE_DROP_THRESHOLD):
    """Find all knee points in performance curve."""
    knees = []
    f1_scores = [fold_results[fold][t]['custom_f1'] for t in thresholds]
    
    for i in range(1, len(f1_scores)):
        if f1_scores[i-1] - f1_scores[i] >= drop_threshold:
            knee_threshold = thresholds[i-1]
            knee_f1 = f1_scores[i-1]
            knees.append((knee_threshold, knee_f1, "knee"))
        elif f1_scores[i-1] - f1_scores[i] <= -drop_threshold:
            knee_threshold = thresholds[i]
            knee_f1 = f1_scores[i]
            knees.append((knee_threshold, knee_f1, "knee"))
    
    return knees

def determine_optimal_knee_by_voting(all_knee_points):
    threshold_votes = {}
    
    # Count votes for each threshold
    for knee in all_knee_points:
        threshold, _, _ = knee
        if threshold not in threshold_votes:
            threshold_votes[threshold] = 0
        threshold_votes[threshold] += 1

    if len(threshold_votes) == 0:
        timestamped_print(f"WARNING: no knee points")
        return None
    
    # Find the maximum vote count
    max_votes = max(threshold_votes.values())
    optimal_thresholds = [threshold for threshold, votes in threshold_votes.items() if votes == max_votes]
    
    return optimal_thresholds


def plot_results(fold_results: Dict, knee_points: Dict, thresholds: List[float], selection_freq_range: str):
    
    fig, ax1 = plt.subplots(1, 1, figsize=(15, 10))
    
    # Plot overall F1 scores
    for fold in fold_results:
        f1_scores = [fold_results[fold][t]['custom_f1'] for t in thresholds]
        ax1.plot(thresholds, f1_scores, label=f'Fold {fold}')
        
        # Mark knee points for each fold
        for knee in knee_points[fold]:
            threshold, f1_score, _ = knee
            ax1.plot(threshold, f1_score, 'ro')  # Red dot for knee points
    
    ax1.set_xlabel('Threshold')
    ax1.set_ylabel('F1 Score')
    ax1.legend()
    
    plt.tight_layout()
    plt.savefig(f"{MINIPATCH_PATH}/crossvalidation_plot_{selection_freq_range}.pdf")
    plt.show()

# Main code
def main():
    Mv, meta = load_data()
    fold_Mvs = load_fold_data()

    # Generate thresholds
    selection_frequency_thresholds = [round(x, 2) for x in np.arange(
        SELECTION_RANGE['begin'], 
        SELECTION_RANGE['end'] + SELECTION_RANGE['step'], 
        SELECTION_RANGE['step']
    )]
    selection_freq_range = f"[{SELECTION_RANGE['begin']},{SELECTION_RANGE['end']}]"
    timestamped_print(f"Testing thresholds: {selection_frequency_thresholds}")

    fold_results = {}
    fold_selectors = {}
    fold_clfs = {}
    fold_best_thresholds = {}
    knee_points = {}
    probe_sets_results = {}

    for fold, (rest_Mv, rest_meta, holdout_Mv, holdout_meta) in fold_Mvs.items():
        timestamped_print(f"\nFold {fold}:")
        fold_results[fold] = {}
        fold_clfs[fold] = {}
        
        rest_multi = utils.propagate_parent(utils.subtree, rest_meta, tissue_col='training.ID', outdict=False)
        rest_mlb = utils.mlb.transform(rest_multi['training.ID'].values)
        holdout_multi = utils.propagate_parent(utils.subtree, holdout_meta, tissue_col='training.ID', outdict=False)
        holdout_mlb = utils.mlb.transform(holdout_multi['training.ID'].values)

        timestamped_print('Fitting minipatch learning...')
        selector = initialize_selector(rest_Mv)
        le = LabelEncoder().fit(meta['training.ID'].unique())
        fitted_selector = selector.fit(rest_Mv.values, le.transform(rest_meta['training.ID']))
        fold_selectors[fold] = fitted_selector
        
        for threshold in selection_frequency_thresholds:
            selection_freq = pd.DataFrame(fitted_selector.Pi_hat_last_k_, index=rest_Mv.columns)
            selected_probes = frozenset(selection_freq[selection_freq[0] >= threshold].index)
            
            if selected_probes in probe_sets_results:
                metrics, num_probes, clf = probe_sets_results[selected_probes]
            else:
                timestamped_print(f"  Threshold: {threshold} | Selected probes: {len(selected_probes)}")
                metrics, num_probes, clf = evaluate_threshold(
                    fitted_selector, threshold, rest_Mv, holdout_Mv, rest_mlb, holdout_mlb
                )
                probe_sets_results[selected_probes] = (metrics, num_probes, clf)
            
            fold_results[fold][threshold] = metrics
            fold_clfs[fold][threshold] = clf
            
    timestamped_print("\nOptimal Threshold Analysis:")
    all_knee_points = []

    for fold in fold_results:
        knees = find_all_knee_points(fold_results, fold, selection_frequency_thresholds)
        knee_points[fold] = knees
        all_knee_points.extend(knees)
        
        timestamped_print(f"\nFold {fold}: {knees}")
        for knee in knees:
            threshold, f1_score, method = knee
            timestamped_print(f"  Method: {method}")
            timestamped_print(f"  Threshold: {threshold:.3f}")
            timestamped_print(f"  F1 score: {f1_score:.3f}")
            timestamped_print(f"  Number of probes: {fold_results[fold][threshold]['num_probes']}")

    timestamped_print(f"\nSaving...")

    with open(f"{MINIPATCH_PATH}/crossvalidation_selectors_{selection_freq_range}", 'wb') as f:
        dill.dump(fold_selectors, f)
    with open(f"{MINIPATCH_PATH}/crossvalidation_clfs_{selection_freq_range}", 'wb') as f:
        dill.dump(fold_clfs, f)
    with open(f"{MINIPATCH_PATH}/crossvalidation_results_{selection_freq_range}", 'wb') as f:
        dill.dump(fold_results, f)
        
    # Plot and save results
    plot_results(fold_results, knee_points, selection_frequency_thresholds, selection_freq_range)

    # Determine the optimal knee point by voting
    optimal_threshold = determine_optimal_knee_by_voting(all_knee_points)
    timestamped_print(f"\nOptimal Threshold by Voting: {optimal_threshold}")

if __name__ == "__main__":
    main()
