# Methylation-classification

Methylation classification is a pan-tissue dataset and analysis paper from publicaly available DNA methylation microarrays (450K) from the Gene Expression Omnibus. 

## Getting started

```
git clone https://github.com/ylaboratory/methylation-classification.git
cd methylation-classification
conda env create -f env.yml --name methyl-classify
conda activate methyl-classify
mkdir download
cd download
huggingface-cli download ylab/methyl-classification
```

Additional libraries required may include:  

-[minipatch-learning](https://github.com/DataSlingers/minipatch-learning)  
-[mlcm](https://github.com/mrh110/mlcm)  

## Data and usage

Analysis code, annotation, data files, and skeleton preprocessing code (if running starting from raw files) are located in this repo and [huggingface](https://huggingface.co/datasets/ylab/methyl-classification/tree/main).  

Huggingface annotation/meta files contain the following:  
- Sample ID  
- training.ID: UBERON ID used for training  
- training.Name: tissue name of ID used for training  
- Dataset: GSE study ID  
- Original.ID: annotated tissue's most descriptive UBERON ID  
- Original.Name: initially annotated tissue    
  
Analysis code are numbered according to order:  
- src/1-stratify-into-folds.py: stratify training dataset and save as formats used in following code  
- src/2-crossvalidation.py: perform crossvalidation for minipatch learning threshold selection  
- src/2-crossvalidation-prediction.py: save prediction results from the chosen threshold  
- src/3-differential-methylation.py: run baseline differential methylation  
- src/3-differential-methylation-prediction.py: save correlation-based prediction results from differential methylation  
- src/4-final-model.py: fit a final multioutput SVM using entire training data  
- src/5-figures-training.ipynb: training and crossvalidation performance  
- src/5-figures-label-transfer.ipynb: label transfer performance  
- src/0-preprocess.R: preprocessing if starting from raw idat files  
  
- utils.py: analysis functions  
- process-microarray.R: preprocesssing functions  
  
Additional files, genomic annotations, and ID to text/color for visualizations are in annotation/. Each dataset's basename files contain more detailed filenames for preprocessing, if user chooses to start from raw idat files self-downloaded from the Gene Expression Omnibus.  
- File: detailed basenames  
- FileSeries: GSE the data file might be under (occasionally different from Dataset)  

## Cite

Please read and cite our paper for more information and usage:  
coming soon

## Contact

Mirae Sunny Kim - [mk98@rice.edu](mk98@rice.edu)  
Project Link: [https://github.com/ylaboratory/methylation-classification](https://github.com/ylaboratory/methylation-classification)




