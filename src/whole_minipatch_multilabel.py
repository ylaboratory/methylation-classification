import warnings
warnings.filterwarnings('ignore')
from mplearn.feature_selection._adaptive_stable_minipatch_selection import AdaSTAMPS
from mplearn.feature_selection.base_selector import DecisionTreeSelector
    
from sklearn.preprocessing import LabelEncoder, MultiLabelBinarizer
from sklearn.multioutput import MultiOutputClassifier
from sklearn.svm import SVC
import networkx as nx
import dill
import pandas as pd
import sys
sys.path.append('./../src/')
import utils
import numpy as np
    
atleast = 2
Mv_location = f"./../data/GEO/preprocessed/450K_Mvalues_atleast{atleast}_samplewise"
print(f"loading Mv, meta, mapping from {Mv_location}")
Mv, meta, mapping = dill.load(open(Mv_location, 'rb'))

filename = f"_whole" 
le=LabelEncoder().fit(meta['tissue_name'].unique())
meta_le=le.transform(meta['tissue_name'].values)

print(meta_le.shape)
print(le.classes_)

m_ratio = np.sqrt(Mv.shape[1])/Mv.shape[1]
n_ratio = np.sqrt(Mv.shape[0])/Mv.shape[0]

selection_frequency_threshold = 0.57

print("adaptive Exploitation & Exploration scheme...")
clf = DecisionTreeSelector(random_state=9)
selector = AdaSTAMPS(base_selector=clf,
                     minipatch_m_ratio=m_ratio,
                     minipatch_n_ratio=n_ratio,
                     random_state=9,
                     verbose=1,
                     sampling_options=None,
                     stopping_criteria_options=None)

fitted_selector = selector.fit(Mv.values, le.transform(meta['tissue_name'].values))
Mv_new = fitted_selector.transform(Mv.values, pi_thr=selection_frequency_threshold)
with open(f"./../data/GEO/minipatch/minipatch{filename}_selector", "wb") as dill_file:
    dill.dump(fitted_selector, dill_file, protocol=4)
    
selection_frequencies = pd.DataFrame(fitted_selector.Pi_hat_last_k_, index=Mv.columns)
with open(f"./../data/GEO/minipatch/minipatch{filename}_frequencies", "wb") as dill_file:
    dill.dump(selection_frequencies, dill_file, protocol=4)
    
print(Mv_new.shape)

date = "sep2024"
atleast = 2
subtree = nx.read_multiline_adjlist(path=f"uberon_{date}_atleast{atleast}_adjlist",create_using=nx.MultiDiGraph)
mlb=MultiLabelBinarizer().fit([[utils.id_to_name[node] for node in subtree.nodes]])

print("multilabel classification...")
meta_multi=utils.propagate_parent(subtree, meta, outdict=False)
meta_mlb=mlb.transform(meta_multi['tissue_name'].values)

clf = MultiOutputClassifier(SVC(class_weight='balanced', kernel='linear', random_state=9, probability=True))
clf.fit(Mv_new, meta_mlb)

location = f"./../data/GEO/minipatch/multilabel{filename}_clf"
print(f"saving to {location}...")
with open(location, "wb") as dill_file:
    dill.dump(clf, dill_file, protocol=4)
