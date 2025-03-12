import warnings
warnings.filterwarnings('ignore')
import time
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from joblib import Parallel, delayed
import numpy as np
import pandas as pd
from utils import name_to_id, id_to_name
import dill
from tqdm import tqdm
from datetime import datetime

DATA_PATH = './../data/GEO'
FOLD = 0  # Define the fold

def log(msg):
    """Prints a message with a timestamp."""
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}")

# Load data
def load_fold_data():
    """Load cross-validation fold data."""
    return dill.load(open(f'{DATA_PATH}/preprocessed/training_folds.dill', 'rb'))

def load_diffmeth_results(fold):
    """Load differential methylation results for a fold."""
    return dill.load(open(f"{DATA_PATH}/diffmeth/diffmeth_fold{fold}", 'rb'))

def get_significant_probes(res_per_tissue):
    """Get significant probes for each tissue using multiple testing correction."""
    from statsmodels.stats.multitest import multipletests
    return {
        tissue: list(tissue_res[multipletests(tissue_res.PValue, alpha=0.05)[0]].index)
        for tissue, tissue_res in res_per_tissue.items()
    }

# PCA Computation
def compute_pca(tissue, holdout_id):
    """Computes PCA projection for a given training tissue."""
    start_time = time.time()
    train_id = name_to_id[tissue][0]
    probes = probe_per_tissue.get(train_id, [])

    if not probes:
        log(f"Skipping {tissue} (no significant probes)")
        return tissue, None, None, None

    pca = PCA(n_components=2).fit(training_mv[probes])
    train_pca = pca.transform(training_mv[probes])
    holdout_pca = pca.transform(holdout_mv.loc[holdout_meta['training.ID'] == holdout_id, probes])

    train_idx = training_meta['training.ID'] == train_id
    elapsed = time.time() - start_time
    log(f"PCA done for {tissue} ({len(probes)} probes) in {elapsed:.2f}s")
    return tissue, train_pca, holdout_pca, train_idx

# PCA Plotting
def plot_all_pca(holdout_tissue, grid_shape=(10, 6)):
    start_time = time.time()
    holdout_id = name_to_id[holdout_tissue][0]
    tissues = [id_to_name[x][0] for x in training_meta['training.ID'].unique()]

    log(f"Starting PCA computation for {len(tissues)} tissues...")

    # Compute PCA projections in parallel with progress bar
    results = Parallel(n_jobs=-1)(
        delayed(compute_pca)(tissue, holdout_id) for tissue in tqdm(tissues, desc="Computing PCA")
    )

    # Filter out failed PCA computations
    results = [r for r in results if r[1] is not None]

    # Plot results
    fig, axes = plt.subplots(*grid_shape, figsize=(grid_shape[1] * 3, grid_shape[0] * 3))
    axes = axes.flatten()

    for i, (tissue, train_pca, holdout_pca, train_idx) in enumerate(results):
        ax = axes[i]
        ax.scatter(*train_pca[~train_idx].T, color='gray', alpha=0.4, label='Other Training')
        ax.scatter(*train_pca[train_idx].T, color='blue', alpha=0.6, label=tissue)
        ax.scatter(*holdout_pca.T, color='red', alpha=0.8, label=holdout_tissue, edgecolors='k')

        ax.set_title(tissue, fontsize=10)
        ax.set_xticks([]), ax.set_yticks([])

    plt.tight_layout()
    plt.savefig(f'pca_all_holdout_{holdout_tissue}.png', dpi=300)
    plt.show()

    log(f"Total execution time: {time.time() - start_time:.2f}s")

# Load Data
log("Loading data...")
folds = load_fold_data()
training_mv, training_meta, holdout_mv, holdout_meta = folds[FOLD]
res_per_tissue = load_diffmeth_results(FOLD)
probe_per_tissue = get_significant_probes(res_per_tissue)
log("Data loaded.")

# Run PCA and plot
plot_all_pca(holdout_tissue='airway')
