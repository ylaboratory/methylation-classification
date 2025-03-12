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
import networkx as nx
from sklearn.svm import SVC
from sklearn.preprocessing import LabelEncoder
from sklearn.multioutput import MultiOutputClassifier
from mplearn.feature_selection._adaptive_stable_minipatch_selection import AdaSTAMPS
from mplearn.feature_selection.base_selector import DecisionTreeSelector

# Local
sys.path.append('./../src/')
sys.modules['sklearn.externals.joblib'] = joblib
import utils
import dill

# Configuration
RANDOM_SEED = 9
DATA_PATH = './../data/GEO'
SELECTION_THRESHOLD = 0.65

# Set random seeds
random.seed(RANDOM_SEED)
np.random.seed(RANDOM_SEED)

def load_data() -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Load main methylation data."""
    Mv_location = f"{DATA_PATH}/preprocessed/training.dill"
    print(f"loading Mv, meta from {Mv_location}")
    return dill.load(open(Mv_location, 'rb'))

def initialize_selector(Mv: pd.DataFrame) -> AdaSTAMPS:
    """Initialize the AdaSTAMPS selector."""
    m_ratio = np.sqrt(Mv.shape[1])/Mv.shape[1]
    n_ratio = np.sqrt(Mv.shape[0])/Mv.shape[0]
    
    clf = DecisionTreeSelector(random_state=RANDOM_SEED)
    return AdaSTAMPS(base_selector=clf,
                     minipatch_m_ratio=m_ratio,
                     minipatch_n_ratio=n_ratio,
                     random_state=RANDOM_SEED,
                     verbose=1)

def train_classifier(Mv_new: np.ndarray, meta_mlb: np.ndarray) -> MultiOutputClassifier:
    """Train the multilabel classifier."""
    clf = MultiOutputClassifier(SVC(class_weight='balanced', 
                                  kernel='linear', 
                                  random_state=RANDOM_SEED, 
                                  probability=True))
    return clf.fit(Mv_new, meta_mlb)

def main():
    # Load data
    Mv, meta = load_data()
    
    # Initialize label encoder
    le = LabelEncoder().fit(meta['training.ID'].unique())
    meta_le = le.transform(meta['training.ID'].values)
    
    print("Adaptive Exploitation & Exploration scheme...")
    selector = initialize_selector(Mv)
    fitted_selector = selector.fit(Mv.values, meta_le)
    Mv_new = fitted_selector.transform(Mv.values, pi_thr=SELECTION_THRESHOLD)
    
    # Save selector and frequencies
    selection_frequencies = pd.DataFrame(fitted_selector.Pi_hat_last_k_, index=Mv.columns)
    with open(f"{DATA_PATH}/minipatch/minipatch_whole_selector", "wb") as f:
        dill.dump(fitted_selector, f)
    with open(f"{DATA_PATH}/minipatch/minipatch_whole_frequencies", "wb") as f:
        dill.dump(selection_frequencies, f)
    
    print(f"Transformed data shape: {Mv_new.shape}")
    
    print("Multilabel classification...")
    meta_multi = utils.propagate_parent(utils.subtree, meta, tissue_col='training.ID', outdict=False)
    meta_mlb = utils.mlb.transform(meta_multi['training.ID'].values)
    
    clf = train_classifier(Mv_new, meta_mlb)
    
    # Save classifier
    location = f"{DATA_PATH}/minipatch/multilabel_whole_clf"
    print(f"Saving to {location}...")
    with open(location, "wb") as f:
        dill.dump(clf, f)

if __name__ == "__main__":
    main()
