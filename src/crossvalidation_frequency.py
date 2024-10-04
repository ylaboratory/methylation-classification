import warnings
warnings.filterwarnings('ignore')
import joblib
import sys
sys.modules['sklearn.externals.joblib'] = joblib
sys.path.append('./../src/')
import dill
from sklearn.svm import SVC
from sklearn.multioutput import MultiOutputClassifier
from sklearn.metrics import f1_score
from sklearn.preprocessing import MultiLabelBinarizer, LabelEncoder
import matplotlib.pyplot as plt
from mplearn.feature_selection._adaptive_stable_minipatch_selection import AdaSTAMPS
from mplearn.feature_selection.base_selector import DecisionTreeSelector

import pandas as pd
import networkx as nx
import utils
import pickle
import numpy as np
import random
import io

random.seed(9)
np.random.seed(9)

date = "sep2024"
atleast = 2
subtree = nx.read_multiline_adjlist(path=f"uberon_{date}_atleast{atleast}_adjlist",create_using=nx.MultiDiGraph)
mlb = MultiLabelBinarizer().fit([[utils.id_to_name[node] for node in subtree.nodes]])

Mv_location = f"./../data/GEO/preprocessed/450K_Mvalues_atleast{atleast}_samplewise"
print(f"loading Mv, meta, mapping from {Mv_location}")
Mv, meta, mapping = dill.load(open(Mv_location, 'rb'))

le = LabelEncoder().fit(meta['tissue_name'].unique())
meta_le = le.transform(meta['tissue_name'].values)

print(meta_le.shape)
# print(le.classes_)

selection_frequency_thresholds = [round(x,2) for x in np.arange(0.5, 0.61, 0.01)]
selection_freq_range = "[0.5,0.6]_dropknee"
print(f"testing {selection_frequency_thresholds}")

def load_fold_data():
    """Load fold data."""
    with open('./../data/GEO/preprocessed/450K_Mvalues_atleast2_samplewise_fold_Mvs', 'rb') as f:
        return pickle.load(f)

def custom_tissue_f1_score(y_true, y_pred):
    sample_score = f1_score(y_true, y_pred, average=None)
    sample_weight = y_true.sum(axis=0)
    valid_indices = sample_weight != 0
    return np.mean(sample_score[valid_indices])

def find_elbow(x_actual, y_actual):
    # Normalize the data
    x = (x_actual - np.min(x_actual)) / (np.max(x_actual) - np.min(x_actual))
    y = (y_actual - np.min(y_actual)) / (np.max(y_actual) - np.min(y_actual))
    
    # Calculate the angle for each point
    npoints = len(x)
    angles = np.zeros(npoints - 2)
    for i in range(1, npoints - 1):
        x1, y1 = x[i] - x[i-1], y[i] - y[i-1]
        x2, y2 = x[i+1] - x[i], y[i+1] - y[i]
        angles[i-1] = np.arctan2(x1*y2 - y1*x2, x1*x2 + y1*y2)
    
    # Find the point of maximum curvature
    elbow_index = np.argmax(angles) + 1
    return x_actual[elbow_index], y_actual[elbow_index]

def find_knee(x_actual, y_actual):
    # Normalize the data
    x = (x_actual - np.min(x_actual)) / (np.max(x_actual) - np.min(x_actual))
    y = (y_actual - np.min(y_actual)) / (np.max(y_actual) - np.min(y_actual))
    
    first = np.array([x[0], y[0]])
    last = np.array([x[-1], y[-1]])
    line_vec = last - first
    point_vec = np.array([x, y]).T - first
    distances = np.abs(np.cross(line_vec, point_vec)) / np.linalg.norm(line_vec)
    
    # The knee is the point with the maximum distance
    knee_index = np.argmax(distances)
    return x_actual[knee_index], y_actual[knee_index]

def find_knee_performance_drop(x, y, drop_threshold=0.01):
    for i in range(1, len(y)):
        if y[i-1] - y[i] > drop_threshold:
            return x[i-1], y[i-1]
    return None, None

fold_Mvs = load_fold_data()
fold_best_thresholds = {}
fold_selectors = {}
fold_clfs = {}
fold_results = {}

for fold, (rest_Mv, rest_meta, holdout_Mv, holdout_meta) in fold_Mvs.items():
    print(f"Fold {fold}:")
    fold_results[fold] = {}
    
    rest_multi = utils.propagate_parent(subtree, rest_meta, outdict=False)
    rest_mlb = mlb.transform(rest_multi['tissue_name'].values)
    holdout_multi = utils.propagate_parent(subtree, holdout_meta, outdict=False)
    holdout_mlb = mlb.transform(holdout_multi['tissue_name'].values)
    
    m_ratio = np.sqrt(rest_Mv.shape[1])/rest_Mv.shape[1]
    n_ratio = np.sqrt(rest_Mv.shape[0])/rest_Mv.shape[0]
    
    clf = DecisionTreeSelector(random_state=9)
    selector = AdaSTAMPS(base_selector=clf,
                         minipatch_m_ratio=m_ratio,
                         minipatch_n_ratio=n_ratio,
                         random_state=9,
                         verbose=0)
    
    fitted_selector = selector.fit(rest_Mv.values, le.transform(rest_meta['tissue_name']))
    fold_selectors[fold] = fitted_selector
    
    best_threshold = 0
    best_f1_score = 0
    
    for threshold in selection_frequency_thresholds:
        selection_freq = pd.DataFrame(fitted_selector.Pi_hat_last_k_, index=rest_Mv.columns)
        minipatch_probes = list(selection_freq[selection_freq[0] >= threshold].index)
        
        rest_Mv_selected = rest_Mv[minipatch_probes]
        holdout_Mv_selected = holdout_Mv[minipatch_probes]
        
        clf = MultiOutputClassifier(SVC(class_weight='balanced', kernel='linear', random_state=9, probability=True))
        clf.fit(rest_Mv_selected.values, rest_mlb)
        
        pred = clf.predict(holdout_Mv_selected.values)
        f1 = custom_tissue_f1_score(holdout_mlb, pred)
        
        # Inside the main loop, after calculating f1 score for each threshold:
        fold_results[fold][threshold] = {
            'num_probes': len(minipatch_probes),
            'f1_score': f1
        }

# In your main code:
elbow_points = {}
knee_points = {}
plt.figure(figsize=(10, 6))
for fold in fold_results:
    thresholds = list(fold_results[fold].keys())
    f1_scores = [fold_results[fold][t]['f1_score'] for t in thresholds]
    num_probes = [fold_results[fold][t]['num_probes'] for t in thresholds]
    
    # Find the knee point
    knee_x, knee_y = find_knee_performance_drop(np.array(thresholds), np.array(f1_scores))
    if knee_x is not None:
        plt.plot(knee_x, knee_y, 'b*', markersize=15)
        print(f"Knee for fold {fold}: threshold = {knee_x}, F1 score = {knee_y}")
    else:
        print(f"No knee found for fold {fold}")

    knee_points[fold] = knee_x
    
    # Update the best threshold for this fold
    best_threshold = knee_x
    fold_best_thresholds[fold] = best_threshold
    
    # Grab fold rest_Mv
    rest_Mv = fold_Mvs[fold][0]
    rest_meta = fold_Mvs[fold][1]
    rest_multi = utils.propagate_parent(subtree, rest_meta, outdict=False)
    rest_mlb = mlb.transform(rest_multi['tissue_name'].values)
    
    # Save the best classifier for this fold
    fitted_selector = fold_selectors[fold]
    selection_freq = pd.DataFrame(fitted_selector.Pi_hat_last_k_, index=rest_Mv.columns)
    minipatch_probes = list(selection_freq[selection_freq[0] >= best_threshold].index)
    rest_Mv_selected = rest_Mv[minipatch_probes]
    clf = MultiOutputClassifier(SVC(class_weight='balanced', kernel='linear', random_state=9, probability=True))
    clf.fit(rest_Mv_selected.values, rest_mlb)
    fold_clfs[fold] = clf
    
    print(f"Best threshold for fold {fold}: {best_threshold} with {len(minipatch_probes)} probes")
    
    # Plot the results
    plt.plot(thresholds, f1_scores, marker='o', label=f'Fold {fold}')
    # plt.plot(elbow_x, elbow_y, 'r*', markersize=15, label=f'Elbow' if fold == '2' else '')
    plt.plot(knee_x, knee_y, 'b*', markersize=15, label=f'Knee' if fold == '2' else '')
    
plt.xticks(thresholds)
plt.xlabel('Threshold')
plt.ylabel('F1 Score')
plt.title('F1 Score vs Threshold for Each Fold')
plt.legend()
plt.grid(True)

# Add number of probes for each fold below x-axis
for i, threshold in enumerate(thresholds):
    probe_counts = []
    for fold in fold_results:
        num_probes = fold_results[fold][threshold]['num_probes']
        probe_counts.append(f"{num_probes}")
    
    plt.text(threshold, plt.ylim()[0] - 0.02, '\n'.join(probe_counts), 
             ha='center', va='top', fontsize=8, rotation=90)
    
plt.savefig(f'./../data/GEO/minipatch/minipatch_crossvalidation_frequency_{selection_freq_range}.png')
plt.show()

# Vote on the best overall threshold based on the knee points
best_overall_threshold = max(set(fold_best_thresholds.values()), key=list(fold_best_thresholds.values()).count)
print(f"Best overall threshold (based on knee): {best_overall_threshold}")

with open(f"./../data/GEO/minipatch/minipatch_crossvalidation_selectors_{selection_freq_range}", 'wb') as f:
    pickle.dump(fold_selectors, f)

with open(f"./../data/GEO/minipatch/multilabel_crossvalidation_clfs_{selection_freq_range}", 'wb') as f:
    pickle.dump(fold_clfs, f)
        
with io.open(f'./../data/GEO/minipatch/minipatch_crossvalidation_frequency_{selection_freq_range}.txt', 'w', encoding='utf-8') as f:
    f.write(f"Best overall threshold (performance drop > 0.01): {best_overall_threshold}\n\n")
    for fold, thresholds in fold_results.items():
        f.write(f"Fold {fold}:\n")
        f.write(f"Best threshold (performance drop > 0.01): {fold_best_thresholds[fold]}\n\n")
        f.write("All thresholds:\n")
        for threshold, results in thresholds.items():
            f.write(f"  Threshold: {threshold}\n")
            f.write(f"  Number of probes: {results['num_probes']}\n")
            f.write(f"  F1 score: {results['f1_score']}\n\n")
        f.write("\n")
        
        
        
        
        
        
        
        
        

# import re

# # Read the data from the file
# with open('./../data/GEO/minipatch/frequency_results.txt', 'r') as f:
#     data = f.read()

# # Parse the data
# folds = re.findall(r'Fold (\d+):(.*?)(?=Fold|\Z)', data, re.DOTALL)

# def find_elbow(x_actual, y_actual):
#     # Normalize the data
#     x = (x_actual - np.min(x_actual)) / (np.max(x_actual) - np.min(x_actual))
#     y = (y_actual - np.min(y_actual)) / (np.max(y_actual) - np.min(y_actual))
    
#     # Calculate the angle for each point
#     npoints = len(x)
#     angles = np.zeros(npoints - 2)
#     for i in range(1, npoints - 1):
#         x1, y1 = x[i] - x[i-1], y[i] - y[i-1]
#         x2, y2 = x[i+1] - x[i], y[i+1] - y[i]
#         angles[i-1] = np.arctan2(x1*y2 - y1*x2, x1*x2 + y1*y2)
    
#     # Find the point of maximum curvature
#     elbow_index = np.argmax(angles) + 1
#     return x_actual[elbow_index], y_actual[elbow_index]

# def find_knee(x_actual, y_actual):
#     # Normalize the data
#     x = (x_actual - np.min(x_actual)) / (np.max(x_actual) - np.min(x_actual))
#     y = (y_actual - np.min(y_actual)) / (np.max(y_actual) - np.min(y_actual))
    
#     first = np.array([x[0], y[0]])
#     last = np.array([x[-1], y[-1]])
#     line_vec = last - first
#     point_vec = np.array([x, y]).T - first
#     distances = np.abs(np.cross(line_vec, point_vec)) / np.linalg.norm(line_vec)
    
#     # The knee is the point with the maximum distance
#     knee_index = np.argmax(distances)
#     return x_actual[knee_index], y_actual[knee_index]

# # In your main code:
# elbow_points = {}
# knee_points = {}

# plt.figure(figsize=(12, 8))

# for fold, fold_data in folds:
#     thresholds = []
#     f1_scores = []
    
#     # Extract threshold and F1 score for each threshold in the fold
#     for match in re.finditer(r'Threshold: ([\d.]+).*?F1 score: ([\d.]+)', fold_data, re.DOTALL):
#         threshold, f1_score = match.groups()
#         thresholds.append(float(threshold))
#         f1_scores.append(float(f1_score))
        
#     # Find the elbow point
#     elbow_x, elbow_y = find_elbow(np.array(thresholds), np.array(f1_scores))
#     elbow_points[fold] = elbow_x
    
#     # Find the knee point
#     knee_x, knee_y = find_knee(np.array(thresholds), np.array(f1_scores))
#     knee_points[fold] = knee_x
    
#     # Plot the data for this fold
#     plt.plot(thresholds, f1_scores, marker='o', label=f'Fold {fold}')
    
#     # Mark the elbow point
#     elbow_x = elbow_points[fold]
#     elbow_y = f1_scores[thresholds.index(elbow_x)]
#     plt.plot(elbow_x, elbow_y, 'r*', markersize=15, label=f'Elbow' if fold == '2' else '')
    
#     plt.plot(knee_x, knee_y, 'b*', markersize=15, label=f'Knee' if fold == '2' else '')
    
#     # Mark the best threshold
#     best_threshold = thresholds[np.argmax(f1_scores)]
#     best_f1 = max(f1_scores)
#     plt.plot(best_threshold, best_f1, 'g^', markersize=15, label=f'Best performance' if fold == '2' else '')



# plt.xlabel('Threshold')
# plt.xticks(thresholds)
# plt.ylabel('F1 Score')
# plt.title('F1 Score vs Threshold for Each Fold (with Elbow Points)')
# plt.legend()
# plt.grid(True)

# # Add number of probes for each fold below x-axis
# for i, threshold in enumerate(thresholds):
#     probe_counts = []
#     for fold, fold_data in folds:
#         match = re.search(f'Threshold: {threshold}.*?Number of probes: (\d+)', fold_data, re.DOTALL)
#         if match:
#             probe_counts.append(f"Fold {fold}: {match.group(1)}")
#     plt.text(threshold, plt.ylim()[0] - 0.01, '\n'.join(probe_counts), 
#              ha='center', va='top', fontsize=8, rotation=90)

# plt.show()
