# correlation test
library(data.table)
library(ggplot2)
library(arrow)

airway_meta = fread('/grain/mk98/methyl/methylation-classification/extra/holdout_airway_meta.csv')
airway_samples = airway_meta$Sample
rm(airway_meta)

training_subset_meta = fread('/grain/mk98/methyl/methylation-classification/extra/training_subset_meta.csv')
training_subset_meta = training_subset_meta[,c('Sample', 'Display.Name', 'training.ID')]
names(training_subset_meta) = c('sample', 'name', 'id')

probes_per_tissue = fread('/grain/mk98/methyl/methylation-classification/extra/tissue_to_probes_fold0.csv', fill=T, sep=",", header=F)
rm(probes_per_tissue)
airway_probes = probes_per_tissue[V1 == 'UBERON:0001005']
airway_probes = t(airway_probes)
airway_probes = as.data.table(airway_probes)
airway_probes = airway_probes[-1]
airway_probes = airway_probes[V1 != ""]
airway_probes = airway_probes$V1

# subset the m.val data to another one for training for correlation
setwd('/grain/rad4/quick')
annots = fread('methyl-classification/data/meta_training_feb26.txt')
annots[, V1 := NULL] # remove the index column
names(annots) = c("sample", "dataset", 'toss1', 'tissue.name', 'tissue.id', 'name', 'display.name', 'mergeID', 'trainID', 'toss2', 'toss3')
annots[, toss1 := NULL]
annots[, toss2 := NULL]
annots[, toss3 := NULL]
# get color scheme
id_to_name = fread('methyl-classification/id_display_mappings.txt')
id_to_name = rbind(id_to_name, data.table(ID = 'root', color = 'grey', display.name = ''))
annots = merge(annots[,c('sample', 'dataset', 'trainID')], id_to_name, by.x = 'trainID', by.y = 'ID')

# load in all m-vals (slow-ish)
parq = '/srv/risotto.cs.rice.edu/scratch/mk98/methyl/data/GEO/preprocessed/training_Mvalues.parquet'
m.val = as.data.table(read_parquet(parq))
setnames(m.val, sub("_.*", "", names(m.val)))

subset.m.val <- lapply(unique(annots$trainID), function(tissue) {
  df = m.val[, c(annots[trainID == tissue]$sample, "probe"), with = FALSE]
  return(list(tissue = tissue, data = df))
})

names(subset.m.val) <- unique(annots$trainID)
rm(m.val)

# subset to airway probes
training.subset.m <- lapply(unique(training_subset_meta$id), function(tis) {
  dt = subset.m.val[[tis]]$data[probe %in% airway_probes, ]
  return(list(tissue = tis, data = dt))
})
names(training.subset.m) <- unique(training_subset_meta$id)

# get airway heldout samples
sams_for_cor = training.subset.m[['UBERON:0001005']]$data
sams_for_cor = sams_for_cor[,c(airway_samples, 'probe'), with=F]

# filter and smush together training data into one dt
training_dt <- lapply(unique(training_subset_meta$id), function(tis) {
  tmp = training.subset.m[[tis]]$data
  dt = tmp[, colnames(tmp) %in% c(training_subset_meta$sample, 'probe'), with = FALSE]
  return(dt)
})
training_dt <- Reduce(function(x, y) merge(x, y, by = "probe", all = TRUE), training_dt)

# calculate pairwise corr for all samples
all.cor <- merge(sams_for_cor, training_dt, by='probe')
all.cor[, probe := NULL]
cor_matrix <- cor(as.matrix(all.cor), use = "pairwise.complete.obs")
nms <- row.names(cor_matrix)
cor_matrix = as.data.table(cor_matrix)
cor_matrix[, rown := nms]
cm <- melt(cor_matrix)
names(cm) <- c('sam1', 'sam2', 'cor')
cm = cm[!(sam1 == sam2)]
cm = cm[sam1 %in% airway_samples]
cm = merge(cm, training_subset_meta, by.x='sam2', by.y='sample')

top_labs <- cm[, .SD[which.max(cor)], by = sam1]
top_labs <- top_labs[,.N,by='name']

