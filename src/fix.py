import numpy as np
import matplotlib.pyplot as plt
import re

def find_knees_moving_average(x, y, window=3, threshold=0.01):
    ma = np.convolve(y, np.ones(window), 'valid') / window
    ma = np.pad(ma, (window-1, 0), mode='edge')
    diff = y - ma
    knee_indices = np.where((diff[1:-1] > diff[:-2]) & (diff[1:-1] > diff[2:]))[0] + 1
    significant_knees = [i for i in knee_indices if diff[i] > threshold * np.max(diff)]
    return [(x[i], y[i]) for i in significant_knees]

def find_knee_performance_drop(x, y, drop_threshold=0.01):
    for i in range(1, len(y)):
        if y[i-1] - y[i] > drop_threshold:
            return x[i-1], y[i-1]
    return None, None

selection_freq_range = "[0.1,0.8]"
filename = f'./../data/GEO/minipatch/minipatch_crossvalidation_frequency_{selection_freq_range}.txt'

with open(filename, 'r') as f:
    data = f.read()

folds = re.findall(r'Fold (\d+):(.*?)(?=Fold|\Z)', data, re.DOTALL)

plt.figure(figsize=(12, 8))

for fold, fold_data in folds:
    thresholds = []
    f1_scores = []
    probe_counts = []
    
    for match in re.finditer(r'Threshold: ([\d.]+).*?Number of probes: (\d+).*?F1 score: ([\d.]+)', fold_data, re.DOTALL):
        threshold, num_probes, f1_score = match.groups()
        thresholds.append(float(threshold))
        f1_scores.append(float(f1_score))
        probe_counts.append(int(num_probes))
    
    plt.plot(thresholds, f1_scores, marker='o', label=f'Fold {fold}')
    
    knee_x, knee_y = find_knee_performance_drop(np.array(thresholds), np.array(f1_scores))
    if knee_x is not None:
        plt.plot(knee_x, knee_y, 'b*', markersize=15)
        print(f"Knee for fold {fold}: threshold = {knee_x}, F1 score = {knee_y}")
    else:
        print(f"No knee found for fold {fold}")
    
    best_threshold = knee_x if knee_x is not None else thresholds[np.argmax(f1_scores)]

    print(f"Best threshold for fold {fold}: {best_threshold}")

    for i, threshold in enumerate(thresholds):
        plt.text(threshold, plt.ylim()[0] - 0.01, f"{probe_counts[i]}", 
                 ha='center', va='top', fontsize=8, rotation=90)

plt.xlabel('Threshold')
plt.xticks(thresholds)
plt.ylabel('F1 Score')
plt.title(f'F1 Score vs Threshold for Each Fold (Frequency Range: {selection_freq_range})')
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.savefig(f'./../data/GEO/minipatch/minipatch_crossvalidation_frequency_{selection_freq_range}_knees.png')
plt.show()
