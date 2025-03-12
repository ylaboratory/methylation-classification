# taskset -c 110-119 python -u minipatch_time.py | tee minipatch_time.txt  

# Standard library
import sys
import random
from typing import Dict, Tuple, List
import time

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
    with open('./time/sampled_samples.dill', 'rb') as f:
        sampled_dict = dill.load(f)
    
    # Load data
    Mv, meta = load_data()
    le = LabelEncoder().fit(meta['training.ID'].unique())

    minipatch_times = pd.DataFrame(columns=sampled_dict.keys())

    for num_samples, idx in sampled_dict.items():
        print(f"num_samples: {num_samples}")
        start_time = time.time()

        sampled_Mv = Mv.loc[idx]
        sampled_meta_le = le.transform(meta.loc[idx,'training.ID'].values)

        selector = initialize_selector(sampled_Mv)
        fitted_selector = selector.fit(sampled_Mv.values, sampled_meta_le)
        Mv_new = fitted_selector.transform(sampled_Mv.values, pi_thr=SELECTION_THRESHOLD)
        print(f"  transformed data shape: {Mv_new.shape}")
          
        # Calculate execution time in seconds
        execution_time = time.time() - start_time
        print(f"  time: {execution_time:.4f}")

        minipatch_times.loc[0, num_samples] = execution_time
        minipatch_times.to_csv('./../data/GEO/time/minipatch_times.csv')

    minipatch_times.to_csv('./../data/GEO/time/minipatch_times.csv')

if __name__ == "__main__":
    main()
