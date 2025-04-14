import random
import numpy as np
import pandas as pd
import dill
import networkx as nx
from sklearn.metrics import precision_score, f1_score
from sklearn.preprocessing import MultiLabelBinarizer
import obonet

# Set random seeds for reproducibility
random.seed(9)
np.random.seed(9)

# Load mappings and ontologies
id_to_name = dill.load(open('./../annotation/id_to_name.dill', 'rb'))
name_to_id = dill.load(open('./../annotation/name_to_id.dill', 'rb'))

with open('./../annotation/ontologies.dill', 'rb') as f:
    _, subtree = dill.load(f)

# Fit MultiLabelBinarizer for ontology nodes
mlb = MultiLabelBinarizer().fit([[node for node in subtree.nodes if node not in ['root', '']]])

def custom_macro_precision(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    """
    Calculate custom macro precision score.
    Args:
        y_true (np.ndarray): Ground truth binary labels.
        y_pred (np.ndarray): Predicted binary labels.
    Returns:
        float: Median precision score across valid classes.
    """
    sample_score = precision_score(y_true, y_pred, average=None)
    sample_weight = y_true.sum(axis=0)
    valid_indices = sample_weight != 0
    return np.median(sample_score[valid_indices])

def custom_tissue_f1_score(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    """
    Calculate custom tissue F1 score.
    Args:
        y_true (np.ndarray): Ground truth binary labels.
        y_pred (np.ndarray): Predicted binary labels.
    Returns:
        float: Median F1 score across valid classes.
    """
    sample_score = f1_score(y_true, y_pred, average=None)
    sample_weight = y_true.sum(axis=0)
    valid_indices = sample_weight != 0
    return np.median(sample_score[valid_indices])

def propagate_parent(subtree, meta, tissue_col='tissue_name', outdict=False):
    """
    Propagate parent nodes in the ontology for each sample.
    Args:
        subtree (networkx.DiGraph): Ontology subtree.
        meta (pd.DataFrame): Metadata containing tissue information.
        tissue_col (str): Column name for tissue labels.
        outdict (bool): Whether to return an output dictionary.
    Returns:
        pd.DataFrame or tuple: DataFrame with propagated parents, optionally with a dictionary.
    """
    print(f"Propagating with {tissue_col}")
    multioutput = pd.DataFrame(index=meta.index, columns=[tissue_col])
    output_dict = {}

    for idx in meta.index:
        sample = meta.at[idx, tissue_col]
        subtree_parents = [x for x in nx.ancestors(subtree, sample) if x not in ['root', '']]
        subtree_parents.append(sample)
        multioutput.at[idx, tissue_col] = set(subtree_parents)
        if sample not in output_dict:
            output_dict[sample] = set(subtree_parents)

    return (multioutput, output_dict) if outdict else multioutput