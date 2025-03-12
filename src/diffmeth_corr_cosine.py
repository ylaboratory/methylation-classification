import warnings
import dill
import pandas as pd
import numpy as np
import random
import logging
from tqdm import tqdm
from statsmodels.stats.multitest import multipletests

# Configure logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

warnings.filterwarnings('ignore')

# Configuration
RANDOM_SEED = 9
DATA_PATH = './../data/GEO'
random.seed(RANDOM_SEED)
np.random.seed(RANDOM_SEED)

def load_data():
    """Load methylation data."""
    Mv_location = f"{DATA_PATH}/preprocessed/training.dill"
    logger.info(f"Loading Mv and meta from {Mv_location}")
    return dill.load(open(Mv_location, 'rb'))

def load_fold_data():
    """Load cross-validation fold data."""
    return dill.load(open(f'{DATA_PATH}/preprocessed/training_folds.dill', 'rb'))

def load_diffmeth_results(fold):
    """Load differential methylation results for a fold."""
    return dill.load(open(f"{DATA_PATH}/diffmeth/diffmeth_fold{fold}", 'rb'))

def get_significant_probes(res_per_tissue):
    """Get significant probes for each tissue using multiple testing correction."""
    return {
        tissue: list(tissue_res[multipletests(tissue_res.PValue, alpha=0.05)[0]].index)
        for tissue, tissue_res in res_per_tissue.items()
    }

def batch_correlation(holdout_Mv, rest_Mv, rest_meta, probe_per_tissue):
    """Compute correlations in large batches using matrix operations and track progress per tissue."""
    corr_results = {}

    # Progress bar for tissue processing
    for tissue, probes in tqdm(probe_per_tissue.items(), desc="Processing Tissues", unit="tissue"):
        if not probes or tissue not in rest_meta['training.ID'].values:
            continue  # Skip if no significant probes

        # Get relevant training samples for the tissue
        tissue_samples = rest_meta.index[rest_meta['training.ID'] == tissue]
        if tissue_samples.empty:
            continue

        tissue_data = rest_Mv.loc[tissue_samples, probes]  # (N_tissue_samples, P)
        if tissue_data.empty:
            continue
        
        # Compute correlations in bulk
        holdout_data = holdout_Mv[probes]  # (N_holdout_samples, P)
        logger.info(f"{tissue}: {holdout_data.shape} x {tissue_data.shape}...")

        # Compute correlation without standardization (np.corrcoef standardizes globally, letting context change things)
        holdout_norms = np.sqrt((holdout_data.to_numpy()**2).sum(axis=1, keepdims=True))  # (N_holdout, 1)
        tissue_norms = np.sqrt((tissue_data.to_numpy()**2).sum(axis=1, keepdims=True))   # (N_tissue, 1)
        
        corr_matrix = (holdout_data @ tissue_data.T) / (holdout_norms * tissue_norms.T)  # (N_holdout, N_tissue)
        
        logger.info(f"{tissue}: {holdout_norms.shape} x {tissue_norms.shape}: {corr_matrix.shape}...")
        
        mean_corrs = corr_matrix.mean(axis=1)

        # Store results
        for i, sample in enumerate(holdout_Mv.index):
            if sample not in corr_results:
                corr_results[sample] = {}
            corr_results[sample][tissue] = mean_corrs[i]

    return corr_results

def main():
    # Load data
    Mv, meta = load_data()
    logger.info(f"{Mv.shape}, {meta.shape}")
    fold_Mvs = load_fold_data()
    
    # Store results
    res_dict = {}

    for fold, (rest_Mv, rest_meta, holdout_Mv, holdout_meta) in fold_Mvs.items():
        logger.info(f"Processing Fold {fold}")

        # Load differential methylation results
        res_per_tissue = load_diffmeth_results(fold)
        probe_per_tissue = get_significant_probes(res_per_tissue)

        # Compute batch correlations
        sample_corrs = batch_correlation(holdout_Mv, rest_Mv, rest_meta, probe_per_tissue)

        # Process results
        for sample, sample_corr in sample_corrs.items():
            if sample_corr:
                pred_tissue = max(sample_corr, key=sample_corr.get)
                res_dict[f"f{fold}.{sample}"] = [sample_corr[pred_tissue], pred_tissue, holdout_meta.loc[sample]['training.ID'], sample_corr]

        # Save fold results
        fold_res_df = pd.DataFrame.from_dict(res_dict, orient='index', columns=['corr', 'pred', 'true', 'dict'])
        fold_res_df.to_pickle(f'{DATA_PATH}/diffmeth/diffmeth_corr_f{fold}.pkl')

    # Save final results
    final_res_df = pd.DataFrame.from_dict(res_dict, orient='index', columns=['corr', 'pred', 'true', 'dict'])
    final_res_df.to_pickle(f'{DATA_PATH}/diffmeth/diffmeth_corr.pkl')

if __name__ == "__main__":
    main()
