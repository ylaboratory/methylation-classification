import warnings
import joblib
import sys
import dill
import pandas as pd
import numpy as np
import random
from typing import Dict, List, Tuple
from concurrent.futures import ProcessPoolExecutor
from statsmodels.stats.multitest import multipletests
import logging
from tqdm import tqdm
from multiprocessing import Manager

# Configure logging
# Configure logging to both console and a log file
logging.basicConfig(level=logging.INFO, 
                    format='%(asctime)s - %(levelname)s - %(message)s',
                    handlers=[
                        logging.StreamHandler(),  # For console output
                        logging.FileHandler("diffmeth_corr_optimized_log.txt")  # For logging to a file
                    ])
logger = logging.getLogger(__name__)

warnings.filterwarnings('ignore')
sys.modules['sklearn.externals.joblib'] = joblib
sys.path.append('./../src/')

# Configuration
RANDOM_SEED = 9
DATA_PATH = './../data/GEO'

# Set random seeds
random.seed(RANDOM_SEED)
np.random.seed(RANDOM_SEED)

def load_data() -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Load main methylation data."""
    Mv_location = f"{DATA_PATH}/preprocessed/training.dill"
    logger.info(f"Loading Mv and meta from {Mv_location}")
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

def calculate_correlation(sample_Mv: pd.Series, rest_Mv: pd.DataFrame, rest_meta: pd.DataFrame,
                          tissue: str, probes: List[str]) -> float:
    """Calculate correlation for a sample with tissue-specific probes."""
    tissue_data = rest_Mv[rest_meta['training.ID'] == tissue][probes]
    return np.corrcoef(sample_Mv[probes], tissue_data)[1:, 0].mean()

def process_samples(batch_samples, fold, rest_Mv, rest_meta, probe_per_tissue, holdout_Mv, holdout_meta):
    batch_res = {}
    for sample in batch_samples:
        sample_Mv = holdout_Mv.loc[sample]
        sample_res = {tissue: calculate_correlation(sample_Mv, rest_Mv, rest_meta, tissue, tissue_probes) 
                     for tissue, tissue_probes in probe_per_tissue.items()}
        
        sample_pred = max(sample_res, key=sample_res.get)
        # Store results for each sample in the batch
        batch_res[sample] = [
            sample_res[sample_pred],
            sample_pred,
            holdout_meta.loc[sample]['training.ID'],
            sample_res
        ]
    return batch_res

def main():
    # Load data
    Mv, meta = load_data()
    logger.info(f"Mv, meta: {Mv.shape}, {meta.shape}")
    fold_Mvs = load_fold_data()
    
    # Initialize results
    res = pd.DataFrame(columns=['corr', 'pred', 'true', 'dict'])

    # Use a manager to safely share the 'res' DataFrame across processes
    with Manager() as manager:
        shared_res = manager.dict()  # Shared dictionary for res

        with ProcessPoolExecutor() as executor:
            for fold, (rest_Mv, rest_meta, holdout_Mv, holdout_meta) in fold_Mvs.items():
                logger.info(f"Processing Fold {fold}")
                
                # Get significant probes for each tissue
                res_per_tissue = load_diffmeth_results(fold)
                probe_per_tissue = get_significant_probes(res_per_tissue)
                
                # Show progress bar for sample processing
                with tqdm(total=len(holdout_Mv), desc=f"Fold {fold} Samples", unit="sample") as pbar:
                    # Process samples in batches
                    batch_size = 4
                    tasks = []
                    for i in range(0, len(holdout_Mv), batch_size):
                        batch_samples = holdout_Mv.index[i:i + batch_size]
                        task = executor.submit(process_samples, batch_samples, fold, rest_Mv, rest_meta, probe_per_tissue, holdout_Mv, holdout_meta)
                        tasks.append(task)

                    # Wait for all tasks to complete and aggregate the results
                    for task in tasks:
                        batch_res = task.result()
                        # Update shared dictionary with the batch results
                        for sample, sample_result in batch_res.items():
                            shared_res[f"f{fold}.{sample}"] = sample_result
                        
                        pbar.update(len(batch_res))

                # Convert the shared dictionary back to a DataFrame
                fold_res_df = pd.DataFrame.from_dict(shared_res, orient='index', columns=['corr', 'pred', 'true', 'dict'])

                # Save fold results
                logger.info(f"Saving fold {fold} results...")
                with open(f'{DATA_PATH}/diffmeth/diffmeth_corr_f{fold}.pkl', 'wb') as f:
                    dill.dump(fold_res_df, f)

        # Optionally, save the entire results at the end
        logger.info("Saving final results...")
        final_res_df = pd.DataFrame.from_dict(shared_res, orient='index', columns=['corr', 'pred', 'true', 'dict'])
        with open(f'{DATA_PATH}/diffmeth/diffmeth_corr.pkl', 'wb') as f:
            dill.dump(final_res_df, f)

if __name__ == "__main__":
    main()
