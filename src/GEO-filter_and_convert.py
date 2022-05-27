#date updated: Apr 2022
#clean up and filter data from GEO saved in /data/GEO/ as all_Mv.txt.gz and all_metadata.txt.gz

import pandas as pd
import numpy as np

#load metadata, sample_id of beta values, and beta values of those with both beta values and metadata
meta_all=pd.read_csv('./../data/GEO/all_metadata_annotated_may2022.txt.gz', sep=',')
beta_samples=pd.read_csv('./../data/GEO/all_betavalues.txt.gz', sep='\t', index_col=0, nrows=0).columns.tolist()
beta_all=pd.read_csv('./../data/GEO/all_betavalues.txt.gz', sep='\t', 
                     usecols=list(set(beta_samples).intersection(set(meta_all['sample_id']))))

print(f"all meta shape: {meta_all.shape}")
print(f"all beta shape: {beta_all.shape}")

#load site indices of beta values i.e. "chr1" "18295" stored in beta value file under V1 and V2
beta_vs=pd.read_csv('./../data/GEO/all_betavalues.txt.gz', sep='\t', usecols=['V1','V2'])

#remove X and Y chromosomes
index_remove=beta_vs[beta_vs['V1']=='chrX'].index
index_remove=index_remove.append(beta_vs[beta_vs['V1']=='chrY'].index)
beta_all=beta_all.drop(index=index_remove)
beta_all=beta_all.reset_index(drop=True)
beta_vs=beta_vs.drop(index=index_remove)
beta_vs=beta_vs.reset_index(drop=True)

#save features as "chr# ####"
feature=np.array(beta_vs['V1'])+' '
feature=pd.DataFrame(np.array(feature)+np.array(list(map(str, beta_vs['V2']))))
feature.index=beta_vs.index

#Get samples with both methylation data and metadata annotation
beta_samples=list(beta_all.columns)
meta_samples=meta_all['sample_id']
both_beta_metadata=list(set(beta_samples).intersection(set(meta_samples)))

#Remove overlapping samples
index_to_delete=['GSM3813550', 'GSM3813542', 'GSM3813555', 'GSM3813541', 'GSM3813540']
for i in index_to_delete:
    if i in both_beta_metadata: both_beta_metadata.remove(i)
        
#filter for those with both beta values and metadata info
beta=beta_all[both_beta_metadata]
meta=meta_all[meta_all['sample_id'].isin(both_beta_metadata)]
meta.index=meta['sample_id']
meta=meta.drop(meta[meta['platform']=='GPL10558'].index)
meta=meta.drop(meta[meta['platform']=='GPL21135'].index)
meta=meta.drop_duplicates('sample_id', 'first')

#filter for normal, non-cancerous tissue with no treatment
normal_meta=meta[meta['disease_name']=='normal']
normal_meta=normal_meta[normal_meta['treatment']==False]
normal_beta=beta[normal_meta.index].T

#metadata typos
#typos
normal_meta.replace('Whoel Blood','Whole Blood', inplace=True)
normal_meta.replace('Breast Neoplasm','Breast', inplace=True)
normal_meta.replace('Embryonic Stemm Cell','Embryonic Stem Cell', inplace=True)
normal_meta.replace('T-Lymphocytes','T-Lymphocyte', inplace=True)
normal_meta.replace('Adipose','Adipose Tissue', inplace=True)
normal_meta.replace('Prostate','Prostate Gland', inplace=True)
normal_meta.replace('Interstinal Mucosa','Intestinal Mucosa',inplace=True)
normal_meta.replace('Whole Blood','Blood', inplace=True)
normal_meta.replace('Peripheral Blood','Blood', inplace=True)
normal_meta.replace('Umbilical Cord Blood','Umbilical Blood', inplace=True)
normal_tissue_term=normal_meta['tissue_name']

#exclude mislabeled tissues
avoid=['Schizophrenia','Sigmoid Colon',
       'Bronchoalveolar Lavage Fluid',
       'Mesenchymal Stem Cell','Mesenchymal',
       'Stem Cell','Embryonic Stem Cell', 'HL60']
avoid_series=['GSE79185']

drop_idx=list(mult_meta[mult_meta['tissue_name'].isin(avoid)].index)
drop_idx.extend(list(mult_meta[mult_meta['series'].isin(avoid_series)].index))

normal_meta=normal_meta.drop(drop_idx)
normal_beta=normal_beta.loc[normal_meta.index]

mult_thresh=0

#filter for those in multiple datasets only
mult_gse=list()
for tissue in np.unique(normal_meta['tissue_name']):
    num_gse=len(np.unique(normal_meta[normal_meta['tissue_name']==tissue]['series']))
    if num_gse>mult_thresh:
        mult_gse.append(tissue)
        
mult_idx=normal_meta.index[normal_meta['tissue_name'].isin(mult_gse)]
mult_meta=normal_meta.loc[mult_idx]
mult_beta=normal_beta.loc[mult_idx]

#convert to Mvalue
Mv=np.log2(np.divide((mult_beta),(1-mult_beta)))
Mv=pd.DataFrame(Mv, columns=mult_beta.columns, index=mult_beta.index)

# #perform PCA with n_components = num.samples
# from sklearn.decomposition import PCA
# pca=PCA(svd_solver='full', n_components=0.999)
# Mv_pca=pd.DataFrame(pca.fit_transform(Mv.values), index=Mv.index)
# import pickle
# with open(f'./../data/GEO_new/pca_pkl', 'wb') as files:
#     pickle.dump(pca, files, protocol=4)

#save
Mv.to_csv(f'./../data/GEO_may/all_Mv.txt.gz',sep='\t', compression='gzip')
mult_meta.to_csv(f'./../data/GEO_may/all_metadata.txt.gz', sep='\t', compression='gzip')
feature.to_csv(f'./../data/GEO_may/all_feature.txt.gz',sep='\t', compression='gzip')