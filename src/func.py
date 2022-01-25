from os.path import abspath, dirname
import importlib as imp
import pandas as pd
import numpy as np
import obonet
import networkx as nx
import matplotlib.pyplot as plt
import time
from sklearn.feature_selection import VarianceThreshold
from sklearn.model_selection import StratifiedKFold
from sklearn.ensemble import AdaBoostClassifier, RandomForestClassifier, StackingClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.svm import LinearSVC

import warnings
warnings.filterwarnings('ignore')
from networkx.drawing.nx_agraph import graphviz_layout
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay
from collections import Counter

# from combat.pycombat import pycombat
# from sklearn.decomposition import PCA
import seaborn as sns

from sklearn.metrics import jaccard_score, precision_score, recall_score
from sklearn.multioutput import MultiOutputClassifier
from sklearn.preprocessing import LabelBinarizer, MultiLabelBinarizer

import func

tissue_ontology_tree=obonet.read_obo('./../annotation/ncit.obo')
#dictionary of names to id from ncit
id_to_name = {id_:data.get('name') for id_, data in tissue_ontology_tree.nodes(data=True)}
name_to_id = {data.get('name'):id_ for id_, data in tissue_ontology_tree.nodes(data=True)}

def make_holdout(beta, meta, disp=False):
    holdout=[]
    for tissue in np.unique(meta['tissue_name']):
        gses=np.unique(meta[meta['tissue_name']==tissue]['series'])
        gse=gses[-1]
        if gse not in holdout:
            holdout.append(gse)

    holdout_meta=meta.loc[meta['series'].isin(holdout)]
    holdout_beta=beta.loc[meta['series'].isin(holdout)]

    rest_meta=meta.loc[set(beta.index)-set(holdout_beta.index)]
    rest_beta=beta.loc[set(beta.index)-set(holdout_beta.index)]
    
    return holdout_beta, holdout_meta, rest_beta, rest_meta

def propagate_parent(tree, subtree_label, actual_label):
    multioutput=pd.DataFrame(index=subtree_label.index, columns=['labels'])
    for sample in subtree_label.index:
        actual=actual_label.loc[sample]
        current=subtree_label.loc[sample]
        actual_parents=set([id_to_name[x] for x in nx.descendants(tissue_ontology_tree, name_to_id[actual])])
        if isinstance(current, pd.Series):
            subtree_parents=set()
            for label in current.values:
                if label!='Other':
                    label_parents=set([id_to_name[x] for x in nx.descendants(tissue_ontology_tree, name_to_id[label])])
                    subtree_parents=subtree_parents.union(label_parents)
        else:
            subtree_parents=set([id_to_name[x] for x in nx.descendants(tissue_ontology_tree, name_to_id[current])])
        shared=set(subtree_parents).union(actual_parents)
        shared.add(current)
        shared.add(actual)
        multioutput.at[sample, 'labels']=shared
    return multioutput

def create_subtree_root(ontology_tree, tissues, root, display=False):
    nodes=[name_to_id[x] for x in tissues]
    all_higher_level_nodes=set()
    all_higher_level_nodes.add(name_to_id[root])
    for node in nodes:
        if display: print(f"node: {id_to_name[node]}")
        higher_level_nodes=set(nx.descendants(ontology_tree, node))
        
        if "Anatomic Structure, System, or Substance" not in [id_to_name[x] for x in higher_level_nodes]:
            print(f"node: {id_to_name[node]}")
            print(f"higher nodes: {[id_to_name[x] for x in higher_level_nodes]}")
        
        higher_level_nodes_root=set()
        if name_to_id[root] in [x for x in ontology_tree.successors(node)]:
            higher_level_nodes_root.add(name_to_id[root])
            if display: print(f"1: {[id_to_name[x] for x in higher_level_nodes_root]}")
        else:
            higher_level_nodes_root.update([x for x in higher_level_nodes if name_to_id[root] in nx.descendants(ontology_tree, x)])
            if display: print(f"2: {[id_to_name[x] for x in higher_level_nodes_root]}")

        if len(higher_level_nodes_root)!=0:
            all_higher_level_nodes.update(higher_level_nodes_root)
            all_higher_level_nodes.update([node])
            if display: print(f"3: {[id_to_name[x] for x in all_higher_level_nodes]}")
            
    return(ontology_tree.subgraph(all_higher_level_nodes))

def next_level(subtree, root, beta, meta):
    next_meta=pd.DataFrame()
    next_beta=pd.DataFrame()
    print(f"sublabels: {[id_to_name[label] for label in subtree.predecessors(name_to_id[root])]}")

    for x in list(subtree.predecessors(name_to_id[root])):
        sublabel=[id_to_name[label] for label in nx.ancestors(subtree, x)]
        sublabel.append(id_to_name[x])
        sublabel_idx=meta.index[meta['tissue_name'].isin(sublabel)]
        sub_meta=meta.loc[sublabel_idx]
        sub_beta=beta.loc[sublabel_idx]

        sub_meta['tissue_name']=id_to_name[x]

        next_meta=pd.concat([next_meta, sub_meta])
        next_beta=pd.concat([next_beta, sub_beta])
    return next_beta, next_meta

def next_level_other(subtree, root, beta, meta, disp=False):
    next_meta=pd.DataFrame()
    next_beta=pd.DataFrame()
    sublabels_name=[id_to_name[label] for label in subtree.predecessors(name_to_id[root])]
    sublabels_id=subtree.predecessors(name_to_id[root])
    if disp: print(f"sublabels: {sublabels_name}")

    for x in list(sublabels_id):
        sublabel=[id_to_name[label] for label in nx.ancestors(subtree, x)]
        sublabel.append(id_to_name[x])
        sublabel_idx=meta.index[meta['tissue_name'].isin(sublabel)]
        sub_meta=meta.loc[sublabel_idx]
        sub_beta=beta.loc[sublabel_idx]

        sub_meta['tissue_name']=id_to_name[x]

        next_meta=pd.concat([next_meta, sub_meta])
        next_beta=pd.concat([next_beta, sub_beta])
    other_beta=beta.loc[~beta.index.isin([next_beta.index])]
    other_meta=meta.loc[~meta.index.isin([next_meta.index])]
    other_meta['tissue_name']='Other'
    next_beta=pd.concat([next_beta, other_beta])
    next_meta=pd.concat([next_meta, other_meta])
    return next_beta, next_meta

def svm_noparent(holdout_beta, holdout_meta, mlb, rest_beta, rest_meta, meta, tree, name_to_id, print_result=False):
    start_time=time.time()
    num_folds=5
     #check if enough data compared to num_folds
    for label in np.unique(rest_meta['tissue_name']):
        if len(rest_meta[rest_meta['tissue_name']==label])<5:
            num_folds=len(rest_meta[rest_meta['tissue_name']==label])
            break
        elif len(rest_meta[rest_meta['tissue_name']==label])<2:
            print(f"not enough samples in {label} with {len(rest_meta[rest_meta['tissue_name']==label])}")
            return
    skf = StratifiedKFold(n_splits=num_folds)
    test_ave=0
    holdout_ave=0
    test_prec_ave=0
    holdout_prec_ave=0
    test_rec_ave=0
    holdout_rec_ave=0
   
    # clf=MultiOutputClassifier(LinearSVC(multi_class='crammer_singer'))
    clf=LinearSVC(multi_class='crammer_singer')
    holdout_df = pd.DataFrame(mlb.transform(holdout_meta['tissue_name']),columns=mlb.classes_, index=holdout_meta.index)
    
    for train_idx, test_idx in skf.split(rest_beta, rest_meta['tissue_name']):
        train_beta=rest_beta.iloc[train_idx]
        train_meta=rest_meta.iloc[train_idx]
        test_beta=rest_beta.iloc[test_idx]
        test_meta=rest_meta.iloc[test_idx]
        
        train_df = pd.DataFrame(mlb.transform(train_meta['tissue_name']),columns=mlb.classes_, index=train_meta.index)
        test_df = pd.DataFrame(mlb.transform(test_meta['tissue_name']),columns=mlb.classes_, index=test_meta.index)
        if print_result:
            print(f"train beta shape: {train_beta.shape}")
            print(f"test beta shape: {test_beta.shape}")
            print(f"holdout beta shape: {holdout_beta.shape}")
            print_result=False

        clf.fit(train_beta, train_meta['tissue_name'])

        #look at misclassified labels
        test_pred=clf.predict(test_beta)
        holdout_pred=clf.predict(holdout_beta)

        test_correct=[(x,y) for x,y in zip(test_pred,test_meta['tissue_name'].values) if (x == y)]
        test_score=len(test_correct)/len(test_pred)

        holdout_correct=[(x,y) for x,y in zip(holdout_pred,holdout_meta['tissue_name'].values) if (x == y)]
        holdout_score=len(holdout_correct)/len(holdout_pred)

        test_ave+=test_score
        holdout_ave+=holdout_score
        
        test_prec=precision_score(test_pred, test_meta['tissue_name'], average='weighted')
        test_rec=recall_score(test_pred, test_meta['tissue_name'], average='weighted')
        
        holdout_prec=precision_score(holdout_pred, holdout_meta['tissue_name'], average='weighted')
        holdout_rec=recall_score(holdout_pred, holdout_meta['tissue_name'], average='weighted')
        
        test_prec_ave+=test_prec
        test_rec_ave+=test_rec
        holdout_prec_ave+=holdout_prec
        holdout_rec_ave+=holdout_rec
        
        if print_result:
            print(f"test acc: {test_score}")
            print(f"test prec: {test_prec}")
            print(f"test rec: {test_rec}")
            print()
            print(f"holdout acc: {holdout_score}")
            print(f"holdout prec: {holdout_prec}")
            print(f"holdout rec: {holdout_rec}")
            print()
            
    print(f"test acc average: {test_ave/num_folds}")
    print(f"test prec average: {test_prec_ave/num_folds}")
    print(f"test rec average: {test_rec_ave/num_folds}")
    print()
    print(f"holdout acc average: {holdout_ave/num_folds}")
    print(f"holdout prec average: {holdout_prec_ave/num_folds}")
    print(f"holdout rec average: {holdout_rec_ave/num_folds}")

    print()
    print("--- %s seconds ---" % (time.time() - start_time))

    return test_beta, test_pred, test_df, holdout_beta, holdout_pred, holdout_df, clf

def find_next_parent_with_children(tree, parent, beta, meta):
    global clf_dict
    for label_below in [id_to_name[x] for x in tree.predecessors(parent)]:
        children=list(subtree.predecessors(name_to_id[label_below]))
        if ((len(children)>1) & (label_below not in clf_dict.keys())):
            recursive_svm(subtree, label_below, beta, meta)
        elif len(children)==1:
            find_next_parent_with_children(tree, name_to_id[label_below], beta, meta)   

def find_next_parent_with_children_other(tree, parent, beta, meta):
    global clf_dict
    for label_below in [id_to_name[x] for x in tree.predecessors(parent)]:
        children=list(subtree.predecessors(name_to_id[label_below]))
        if ((len(children)>1) & (label_below not in clf_dict.keys())):
            recursive_svm_other(subtree, label_below, beta, meta)
        elif len(children)==1:
            find_next_parent_with_children_other(tree, name_to_id[label_below], beta, meta)  
            
def recursive_svm(subtree, current_parent, beta, meta):
    global clf_dict
    print(f"parent node: {current_parent}")
    parent_id=name_to_id[current_parent]
    
    next_beta, next_meta = next_level(subtree, current_parent, beta, meta, False)
    holdout_beta, holdout_meta, rest_beta, rest_meta = make_holdout(next_beta, next_meta)
    lb=LabelBinarizer().fit(meta['tissue_name'])
    try:
        test_beta, test_pred, test_df, holdout_beta, holdout_pred, holdout_df, clf=svm_noparent(holdout_beta, holdout_meta,
                                                                                           lb,
                                                                                           rest_beta, rest_meta,
                                                                                           subtree_meta, tissue_ontology_tree, name_to_id, True)
        clf_dict[current_parent]=clf
    except ValueError:
        print("VALUE ERROR")
    print()
    
    if len(list(subtree.predecessors(parent_id)))!=0:
        find_next_parent_with_children(subtree, parent_id, beta, meta)
        
def recursive_svm_other(subtree, current_parent, beta, meta):
    global clf_dict
    print(f"parent node: {current_parent}")
    parent_id=name_to_id[current_parent]
    
    next_beta, next_meta = next_level_other(subtree, current_parent, beta, meta, False)
    holdout_beta, holdout_meta, rest_beta, rest_meta = make_holdout(next_beta, next_meta)
    lb=LabelBinarizer().fit(meta['tissue_name'])
    try:
        test_beta, test_pred, test_df, holdout_beta, holdout_pred, holdout_df, clf=svm_noparent(holdout_beta, holdout_meta,
                                                                                           lb,
                                                                                           rest_beta, rest_meta,
                                                                                           subtree_meta, tissue_ontology_tree, name_to_id, True)
        clf_dict[current_parent]=clf
        print(Counter(lb.inverse_transform(test_df.values)))
        print(Counter(holdout_meta['tissue_name']))
    except ValueError:
        print("VALUE ERROR")
    print()
    
    if len(list(subtree.predecessors(parent_id)))!=0:
        find_next_parent_with_children_other(subtree, parent_id, beta, meta)
        
