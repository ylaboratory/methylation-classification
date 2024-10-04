import sklearn as sk
import pandas as pd
import matplotlib.pyplot as plt
import random
random.seed(9)
import networkx as nx
import numpy as np
np.random.seed(9)
from collections import Counter
import obonet

uberon=obonet.read_obo('/grain/mk98/methyl/methylation-classification/annotation/uberon_ext.obo')
id_to_name = {id_:data.get('name') for id_, data in uberon.nodes(data=True)}
name_to_id = {data.get('name'):id_ for id_, data in uberon.nodes(data=True)}

def generate_colors(targets):
    colors = list()
    for n in range(len(targets)):
        random.seed(n)
        color = "#%06x" % random.randint(0, 0xFFFFFF)
        colors+=[color]
    return colors

def plot_pca(Mv, meta, return_pcs=False):
    standardized_Mv = sk.preprocessing.StandardScaler().fit_transform(Mv)
    pca = sk.decomposition.PCA(n_components=2)
    principalComponents = pca.fit_transform(standardized_Mv)
    principalDf = pd.DataFrame(data = principalComponents
                               , columns = ['principal component 1', 'principal component 2'])
    principalDf.index=meta.index
    fig = plt.figure(figsize = (10,10))
    ax = fig.add_subplot(1,1,1) 
    ax.set_xlabel('Principal Component 1', fontsize = 15)
    ax.set_ylabel('Principal Component 2', fontsize = 15)
    ax.set_title(
        f'PCA: {standardized_Mv.shape[1]}, {[round(x, 4) for x in pca.explained_variance_ratio_]}'
        , fontsize = 20)
    targets = meta['tissue_name'].unique()
    colors = generate_colors(targets)
    for target, color in zip(targets,colors):
        indicesToKeep = meta[meta['tissue_name'] == target].index
        ax.scatter(principalDf.loc[indicesToKeep, 'principal component 1']
                   , principalDf.loc[indicesToKeep, 'principal component 2']
                   , c = color
                   , s = 10)
    ax.legend(targets)
    ax.grid()
    
    if return_pcs:
        return principalDf
    
def weighted_variance(Mv, meta, island, threshold):
    from sklearn.model_selection import train_test_split
    # holdout_Mv, holdout_meta, rest_Mv, rest_meta = make_holdout(Mv, meta, seed=fold, disp=False)
    train_Mv, test_Mv, train_meta, test_meta = train_test_split(Mv, meta, random_state=9)
    tissue_count=dict(Counter(train_meta['tissue_name']))
    train_Mv_mean=train_Mv.mean()
    train_Mv_diff_squared=(train_Mv-train_Mv_mean)**2
    tissue_weights=[1/tissue_count[train_meta['tissue_name'].loc[x]] for x in train_Mv_diff_squared.index]
    train_Mv_weighted_diff_sum_squared=train_Mv_diff_squared.mul(tissue_weights, axis=0)
    train_Mv_weighted_var=train_Mv_weighted_diff_sum_squared.sum()/np.sum(tissue_weights)
    top_train_Mv_weighted_var=train_Mv_weighted_var[train_Mv_weighted_var>threshold]

    top_weighted_var_island=list(island.loc[top_train_Mv_weighted_var.index].index)
    print(f'high variance features: {len(top_weighted_var_island)}')
    
    return top_weighted_var_island
    
def weighted_variance_folds(Mv, meta, island, threshold):
    all_filtered_island=list()
    for fold in range(3):
        holdout_Mv, holdout_meta, rest_Mv, rest_meta = make_holdout(Mv, meta, seed=fold, disp=False)
        tissue_count=dict(Counter(rest_meta['tissue_name']))
        rest_Mv_mean=rest_Mv.mean()
        rest_Mv_diff_squared=(rest_Mv-rest_Mv_mean)**2
        tissue_weights=[1/tissue_count[rest_meta['tissue_name'].loc[x]] for x in rest_Mv_diff_squared.index]
        rest_Mv_weighted_diff_sum_squared=rest_Mv_diff_squared.mul(tissue_weights, axis=0)
        rest_Mv_weighted_var=rest_Mv_weighted_diff_sum_squared.sum()/np.sum(tissue_weights)
        top_rest_Mv_weighted_var=rest_Mv_weighted_var[rest_Mv_weighted_var>threshold]

        top_weighted_var_island=island.loc[top_rest_Mv_weighted_var.index]
        all_filtered_island.append(list(top_weighted_var_island.index))
    
    shared_island_weighted=sorted(set(all_filtered_island[0]).intersection(set(all_filtered_island[1])).intersection(set(all_filtered_island[2])))
    print(f'shared high variance features: {len(shared_island_weighted)}')
    
    return shared_island_weighted
    
def high_weight_features(multioutput, mlb, probes, weight_threshold, tissues=None, scaled=True, disp=True):
    if tissues==None: tissues=mlb.classes_
    if scaled: print("weights will be scaled")
    for tissue in tissues:
        if tissue not in mlb.classes_: print("tissue not in classes"); continue
        tissue_to_probes=dict()
        tissue_to_values=dict()
        tissue_to_values_signed=dict()
        for i, tissue in enumerate(tissues):
            #if weight threshold by value
            if weight_threshold<1: 
                weights=abs(multioutput.estimators_[i].coef_)
                high_weights=list(weights>weight_threshold)[0]
                high_weight_probes = [probe for (probe, high_weight) in zip(probes, high_weights) if high_weight]
                high_weight_values = [weight for (weight, high_weight) in zip(weights, high_weights) if high_weight]
                if scaled:
                    if disp: print(f'tissue: {tissue}, mean: {np.mean(weights)}, std: {np.std(weights)}')
                    high_weight_values= [weight-np.mean(weights)/np.std(weights) for weight in high_weight_values]
                if disp: print(f'tissue: {tissue}, num. important features: {len(high_weight_probes)}')
            #if top # weights wanted
            else:
                weights=abs(multioutput.estimators_[i].coef_)
                weights_signed=multioutput.estimators_[i].coef_
                high_weight_probes=[probes[x] for x in np.argsort(list(weights)[0])[::-1][:weight_threshold]]
                high_weight_values=[list(weights)[0][x] for x in np.argsort(list(weights)[0])[::-1][:weight_threshold]]
                high_weight_values_signed=[list(weights_signed)[0][x] for x in np.argsort(list(weights)[0])[::-1][:weight_threshold]]
                if scaled:
                    if disp: print(f'tissue: {tissue}, mean: {np.mean(weights)}, std: {np.std(weights)}')
                    # high_weight_values= [(weight-np.mean(weights))/np.std(weights) for weight in high_weight_values] #divide by std
                    high_weight_values= [(weight-np.mean(weights)) for weight in high_weight_values] #no dividing by std
                    high_weight_values_signed= [(weight-np.mean(weights_signed)) for weight in high_weight_values_signed] #not abs
                if disp: print(f'tissue: {tissue}, average weight: {np.average(high_weight_values):0.02e}')
            tissue_to_probes[tissue]=high_weight_probes
            tissue_to_values[tissue]=high_weight_values
            tissue_to_values_signed[tissue]=high_weight_values_signed
        return tissue_to_probes, tissue_to_values, tissue_to_values_signed

def create_subtree_root(tree, tissues, root, display=False):
    nodes=[name_to_id[x] for x in tissues if x in name_to_id.keys()]
    nodes.extend([x for x in tissues if x in id_to_name.keys()])
    root_node=name_to_id[root]
    root_children=nx.ancestors(tree,root_node)
    root_children.add(root_node)
    all_parents_under_root=[root_node]
    for num, node in enumerate(nodes):
        if node in tree.nodes:
            all_parents=nx.descendants(tree,node)
            parents_under_root=[x for x in all_parents if x in root_children]
            all_parents_under_root.append(node)
            [all_parents_under_root.append(x) for x in parents_under_root]
            if display: print([id_to_name[x] for x in all_parents_under_root])
        else: 
            print(f"{id_to_name[node]} not in tree")
            print()
    return tree.subgraph(set(all_parents_under_root))

def make_holdout(beta, meta, seed=0, disp=False):
    holdout=[]
    # random.seed(seed)
    for tissue in meta['tissue_name'].unique():
        gses=np.unique(meta[meta['tissue_name']==tissue]['series'])
        # if len(set(gses).intersection(set(holdout)))>0: 
        #     if disp: 
        #         print(f"tissue: {tissue}")
        #         print(f"num GSE: {len(gses)}")
        #         print(f"holdout GSE: already in holdout")
        #     continue
            
        if seed<len(gses): 
            gse=gses[seed]
            gse_meta = meta[meta['series'] == gse]
            holdout+=list(gse_meta[gse_meta['tissue_name']==tissue]['sample_id'])
        else: 
            # gse=gses[seed%len(gses)]
            # print(tissue)
            # print(gses)
            continue
            
        # if gse not in holdout:
        #     holdout.append(gse)
        # if disp:
        #     print(f"tissue: {tissue}")
        #     print(f"num GSE: {len(gses)}")
        #     print(f"holdout GSE: {gse}")
    if disp: 
        print(holdout)

    holdout_meta=meta.loc[meta['sample_id'].isin(holdout)]
    holdout_beta=beta.loc[holdout_meta.index]

    rest_meta=meta.loc[set(beta.index)-set(holdout_beta.index)]
    rest_beta=beta.loc[rest_meta.index]
    
    return holdout_beta, holdout_meta, rest_beta, rest_meta

def propagate_parent(subtree, meta, outdict=False):
    multioutput=pd.DataFrame(index=meta.index, columns=['tissue_name'])
    output_dict=dict()
    for idx in meta.index:
        sample=meta['tissue_name'].loc[idx]
        subtree_parents=[id_to_name[x] for x in nx.descendants(subtree, name_to_id[sample])]
        subtree_parents.append(sample)
        multioutput.at[idx, 'tissue_name']=set(subtree_parents)
        if sample not in output_dict.keys(): output_dict[sample]=set(subtree_parents)
    if outdict:
        return multioutput, output_dict
    else:
        return multioutput

def collapse_tree(tree, root, tissues):
    tree_copy=copy.deepcopy(tree)
    tree_copy=nx.MultiDiGraph(tree_copy)
    for node in sorted(tree.nodes):
        if id_to_name[node] not in tissues:
            parents=list(tree_copy.successors(node))
            children=list(tree_copy.predecessors(node))
            if len(parents)==1 and len(children)==1:
                tree_copy.remove_node(node)
                tree_copy.add_edge(children[0], parents[0])
    return tree_copy
