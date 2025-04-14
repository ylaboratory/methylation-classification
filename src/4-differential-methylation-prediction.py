import warnings
warnings.filterwarnings('ignore')
import dill
import pandas as pd
import numpy as np
import random
import logging
from tqdm import tqdm
from statsmodels.stats.multitest import multipletests
from utils import name_to_id, id_to_name
additional = f"allprobes"

# Configure logging
logging.basicConfig(
    filename=f"diffmeth_corr_{additional}.log",
    level=logging.INFO, 
    format="%(asctime)s - %(levelname)s - %(message)s",
    filemode="w"
)
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
    import sys    
    if not (rest_Mv.index == rest_meta.index).all():
        sys.exit("Error: Index mismatch between rest_Mv and rest_meta. Exiting.")
    
    logger.info(f"correlating: {holdout_Mv.shape} x {rest_Mv.shape}...")

    corr_matrix = np.corrcoef(holdout_Mv, rest_Mv, rowvar=True)
    n_holdout = holdout_Mv.shape[0] #axis change
    logger.info(f"{n_holdout}, {corr_matrix.shape}...")
    corr_matrix = corr_matrix[:n_holdout, n_holdout:]
    logger.info(f"{n_holdout}, {corr_matrix.shape}...")
    mean_corrs = {
        tissue: corr_matrix[:, rest_meta['training.ID'] == tissue].mean(axis=1)
        for tissue in probe_per_tissue.keys()
    }

    # Store results in a structured format: {sample: {tissue: mean_corr}}
    corr_results = {
        sample: {tissue: mean_corrs[tissue][i] for tissue in probe_per_tissue.keys()}
        for i, sample in enumerate(holdout_Mv.index)
    }

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
        all_probes = set().union(*probe_per_tissue.values())

        sample_corrs = batch_correlation(holdout_Mv[all_probes], rest_Mv[all_probes], rest_meta, probe_per_tissue)

        # Process results
        for sample, sample_corr in sample_corrs.items():
            if sample_corr:
                pred_tissue = max(sample_corr, key=sample_corr.get)
                res_dict[f"f{fold}.{sample}"] = [sample_corr[pred_tissue], pred_tissue, holdout_meta.loc[sample]['training.ID'], sample_corr]

        # Save fold results
        fold_res_df = pd.DataFrame.from_dict(res_dict, orient='index', columns=['corr', 'pred', 'true', 'dict'])
        fold_res_df.to_pickle(f'{DATA_PATH}/diffmeth/diffmeth_corr_f{fold}_{additional}.pkl')

    # Save final results
    final_res_df = pd.DataFrame.from_dict(res_dict, orient='index', columns=['corr', 'pred', 'true', 'dict'])
    final_res_df.to_pickle(f'{DATA_PATH}/diffmeth/diffmeth_corr_{additional}.pkl')

if __name__ == "__main__":
    main()
