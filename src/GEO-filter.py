#Feb 2022
#clean up and filter data from GEO saved in /data/GEO/ as _all_betavalues.txt.gz and _all_metadata.txt.gz

import pandas as pd

#load metadata, sample_id of beta values, and beta values of those with both beta values and metadata
meta_all=pd.read_csv('./../data/GEO/all_metadata_annotated_new.txt.gz', sep=',')
beta_new_samples=pd.read_csv('./../data/GEO/all_betavalues.txt.gz', sep='\t', index_col=0, nrows=0).columns.tolist()
beta_all=pd.read_csv('./../data/GEO/all_betavalues.txt.gz', sep='\t', usecols=list(set(beta_new_samples).intersection(set(meta_new['sample_id']))))

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
normal_beta=beta[normal_meta.index]

#metadata typos
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
normal_meta.replace('T-Lymphocytes', 'T-Lymphocyte', inplace=True)
normal_tissue_term=normal_meta['tissue_name']

#exclude mislabeled tissues
avoid=['Schizophrenia','Sigmoid Colon',
       'Bronchoalveolar Lavage Fluid',
       'Mesenchymal Stem Cell','Mesenchymal',
       'Stem Cell','Embryonic Stem Cell']
normal_sample_index=[i for i in range(len(normal_tissue_term)) \
                     if (normal_tissue_term.iloc[i] not in avoid)]

normal_meta=normal_meta.iloc[normal_sample_index]
normal_beta=beta[normal_meta.index]

#check beta value distribution
import matplotlib.pyplot as plt

plt.hist(np.array(normal_beta.values).flatten(), range=(0,1))
plt.title('beta value distribution')
plt.show()
plt.hist(np.array(normal_beta.values).var(1))
plt.title('beta value variance distribution')
plt.show()

#filter by variance threshold from beta value variance distibution
thresh=0.12

#filter by variance and get mask
selector=VarianceThreshold(threshold=(thresh))
selector.fit(normal_beta.T.values)
mask=selector.get_support(indices=True)

#filtered [chr# position] with mask
filtered_feature=feature.iloc[mask]
filtered_beta=pd.DataFrame(selector.transform(normal_beta.T.values), 
                           index=normal_beta.columns,
                           columns=filtered_feature.values)
filtered_meta=normal_meta 

#filter for those in multiple datasets only
mult_gse=list()
studies_by_tissue=dict()
for tissue in np.unique(filtered_meta['tissue_name']):
    num_gse=len(np.unique(filtered_meta[filtered_meta['tissue_name']==tissue]['series']))
    if num_gse>1:
        mult_gse.append(tissue)
        studies_by_tissue[tissue]=num_gse
        
mult_idx=filtered_meta.index[filtered_meta['tissue_name'].isin(mult_gse)]
mult_meta=filtered_meta.loc[mult_idx]
mult_beta=filtered_beta.loc[mult_idx]

#save
mult_beta.to_csv(f'./../data/GEO/all_betavalue_filtered{thresh}_mult.txt.gz',sep='\t', compression='gzip')
mult_meta.to_csv(f'./../data/GEO/all_metadata_filtered{thresh}_mult.txt.gz', sep='\t', compression='gzip')
filtered_feature.to_csv(f'./../data/GEO/all_feature_filtered{thresh}_mult.txt.gz',sep='\t', compression='gzip')