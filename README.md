# Ontology-aware DNA methylation classification of healthy human tissues and cell types

This repository contains code for assembling, analyzing, and performing multi-label tissue and cell type classification
using Illumina 450k DNA metnylation data obtained from the Gene Expression Omnibus (GEO). The full processed dataset used for
training and evaluation can be found on [HuggingFace](https://huggingface.co/datasets/ylab/methyl-classification).

## Citation

> Ontology-aware DNA methylation classification with a curated atlas of human tissues and cell types.
Kim M, Dannenfelser R, Cui Y, Allen G, and Yao V. (2025) bioRxiv. https://doi.org/10.1101/2025.04.18.649618

## Getting started

Conda is recommended for easy installation of relevant dependencies.
Once the repo is cloned create a new environment with the dependencies
specified in `methylation-classification/env.yml`.

```sh
conda env create -f env.yml --name methyl-classify
```

Activate the environment:

```sh
conda activate methyl-classify
```

The 450k DNA methylation atlas is stored on
[HuggingFace](https://huggingface.co/datasets/ylab/methyl-classification/tree/main).
One way to retrieve these files is to download the files using
the HuggingFace command line interface. Storing these files within the `downloads` directory
of `methylation-classification` will enable seamless with the existing path structure for
downstream analysis and classification scripts.

```
cd methylation-classification
mkdir download
cd download
huggingface-cli download ylab/methyl-classification
```

To run the feature selection part of our pipeline we require
another python library, [minipatch](https://github.com/DataSlingers/minipatch-learning).
We refer to user to their [GitHub](https://github.com/DataSlingers/minipatch-learning)
for installation instructions.

## Code structure and usage

All relevant source code is stored in `src`, supporting files containing the
ontology structure and other relevant metadata are in `annotation`. Scripts are numbered
roughly in by the order in which they should be run as part of a larger pipeline.
We provide a further breakdown of scripts in the following sections.

### Atlas generation

For the sake of brevity we only include our preprocessing script for normalizing
downloaded idat files from GEO (`src/0-preprocess.R`). Run this script by pointing
the `RAW_DATA_DIR` file path to point to a directory of downloaded idats. We provide
additional mappings in `annotation` if users wish to start directly from these files. 

We recommend starting directly with the already preprocessed atlas data
in the next steps.

### Classification preprocessing

 `src/1-stratify-into-folds.py` starts with the 450k atlas data already downloaded
 from our HuggingFace repo. This initial script divides the training data into folds
 with respect to dataset and tissue / cell labels.

### Cross-validation and final model training

These scripts require successful installation of Minipatch as described above
and that the stratification script has been run above. In this step we first
run `src/2-crossvalidation.py` to perform crossvalidation with Minipatch at different
thresholds. Once a threshold is selected for the optimal feature set we then run 
`src/2-crossvalidation-prediction.py` to calculate prediction results using the chosen threshold. 

Once the parameters are chosen we can then run `src/4-final-model.py` to save a final 
multi-output SVM using the entire training dataset.
 
### Differential methylation

In the paper we compare our ontology-aware classification method with results from
differential methylation applied to the entire atlas. These scripts calculate
baseline differential methylation (`src/3-differential-methylation.py`) and evaluate
using differential methylated probes (`src/3-differential-methylation-prediction.py`). Be
aware that these scripts are slow and require large amounts of memory.

### Additional analyses and evaluations

We generate additional figures evaluating the performance and composition of the
altas in the following jupyter notebooks: `src/5-figures-training.ipynb` and
`src/5-figures-label-transfer.ipynb`. The former evaluates primarily on the 
training set tissues while the latter assesses the generalizability of our
classifier on unseen tissue and cell type labels.

  




