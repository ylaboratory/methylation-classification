#date updated: May 2022

library(data.table)
library(R.utils)
meta_dir='data/GEO/'
database_type='GEO'

annot<-function(total_metadata, meta_dir, sample, normal_vec, tissue_vec, 
                treatment_vec=c(rep('False',length(normal_vec))),
                disease_vec=c(rep('Unknown',length(normal_vec)))){
  sample_matrix<-read.csv(paste0(meta_dir, sample,'_sample_metadata.txt'), sep='\t')
  if (dim(sample_matrix)[1]!=length(disease_vec) || 
      dim(sample_matrix)[1]!=length(tissue_vec) || 
      dim(sample_matrix)[1]!=length(treatment_vec)){
    print(sample)
    stop("vector size not consistent")
  }
  
  if (!(sample %in% total_metadata$series)){
    total_metadata=rbind(total_metadata, 
                         data.table(sample_id=sample_matrix$Samples, 
                                    series=sample_matrix$Series, 
                                    platform=sample_matrix$Platform,
                                    normal=normal_vec, 
                                    tissue_name=tissue_vec, 
                                    treatment=treatment_vec,
                                    disease=disease_vec))
  }
  return(total_metadata)
}

append_to_remove<-function(list, item){
  if (!(item %in% list)){list<-append(list, item)}
  return(list)
}

total_metadata<-data.table(sample_id=character(), 
                           series=character(), 
                           platform=character(), 
                           normal=character(), 
                           tissue_name=character(), 
                           treatment=character(),
                           disease=character())
to_remove<-list()

# total_metadata<-annot(total_metadata, meta_dir, 'GSE102031', c(rep('normal',16)), c(rep('Embryonic Stem Cell', 16)), c(rep('False', 16)))
# total_metadata<-annot(total_metadata, meta_dir, 'GSE100780', c(rep('normal', 4)), c(rep('HL60', 4)), c(rep('False',2), rep('True',2))) 
# total_metadata<-annot(total_metadata, meta_dir, 'GSE103006', c(rep('normal', 43)), c(rep('Hematopoietic Stem Cell', 43)), c(rep('False',43)))
to_remove<-append_to_remove(to_remove, 'GSE103027')
to_remove<-append_to_remove(to_remove, 'GSE103271')
to_remove<-append_to_remove(to_remove, 'GSE103279')
to_remove<-append_to_remove(to_remove, 'GSE103280')
to_remove<-append_to_remove(to_remove, 'GSE103287')
to_remove<-append_to_remove(to_remove, 'GSE103328')
to_remove<-append_to_remove(to_remove, 'GSE104359')
total_metadata<-annot(total_metadata, meta_dir, 'GSE106360', 
                      c(rep('disease', 47)), 
                      c(rep('Breast', 47)), 
                      c(rep('False', 47)), 
                      c(rep('Breast Carcinoma', 47)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE106437', 
                      c(rep('disease', 3)), 
                      c(rep('Colon',3)), 
                      c('True', rep('False', 2)),
                      c(rep('Colorectal Adenocarcinoma', 3)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE106438', 
                      c(rep('disease', 3)), 
                      c(rep('Colon', 3)), 
                      c('True', rep('False', 2)),
                      c(rep('Colon Carcinoma', 3)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE106556', 
                      c(rep(list('disease','normal'), 10)), 
                      c(rep('Rectal', 18), rep('Sigmoid Colon',2)), 
                      c(rep('False', 20)),
                      c(rep(c('Rectal Neoplasm', 'normal'), 9), 'Colon Neoplasm','normal'))
total_metadata<-annot(total_metadata, meta_dir, 'GSE106727', 
                      c(rep('disease', 18)), 
                      c(rep('Brain',18)), 
                      c(rep('False', 18)),
                      c(rep('Brain Neoplasm',18)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE108123', 
                      c(rep('disease', 87)), 
                      c(rep('Lung', 87)), 
                      c(rep('False', 87)),
                      c(rep('Lung Squamous Cell Carcinoma', 87)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE108143', 
                      c(rep('normal', 2)), 
                      c(rep('Embryonic Stem Cell',2)), 
                      c('False','True'))
total_metadata<-annot(total_metadata, meta_dir, 'GSE108202', 
                      c(rep('normal', 3), rep('disease',3), rep('normal',2)), 
                      c(rep('Fallopian Tube',8)), 
                      c(rep('False',8)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE108785', 
                      c(rep('disease',3), rep('normal',3)), 
                      c(rep('Blood', 6)), 
                      c(rep('False',6)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE108982', 
                      c(rep('disease',16)), 
                      c(rep('Breast',16)), 
                      c(rep('False',16)),
                      c(rep('Breast Neoplasm',16)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE109330', 
                      c(rep('disease',27)), 
                      c(rep('Brain',27)), 
                      c(rep('False',27)),
                      c(rep('Glioblastoma',27)))
to_remove<-append_to_remove(to_remove, 'GSE109364')
total_metadata<-annot(total_metadata, meta_dir, 'GSE109541', 
                      c(rep('disease',4)), 
                      c(rep('Stomach',4)), 
                      c(rep('False',4)),
                      c(rep('Gastric Adenocarcinoma',4)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE109904', 
                      c(rep('disease',6)), 
                      c(rep(c('Bone','Blood'),3)), 
                      c(rep('False',6)),
                      c(rep(c('Sarcoma','Leukemia'),3)))
to_remove<-append_to_remove(to_remove, 'GSE110184')
to_remove<-append_to_remove(to_remove, 'GSE110697')
total_metadata<-annot(total_metadata, meta_dir, 'GSE100197', 
                      c(rep('disease',51),rep('normal',51)),
                      c(rep('Placenta',102)), 
                      c(rep('False',102)),
                      c(rep('Eclampsia',22),rep('normal',11),rep('Eclampsia',18),rep('normal',51)))
to_remove<-append_to_remove(to_remove, 'GSE100249')
total_metadata<-annot(total_metadata, meta_dir, 'GSE100386', 
                      c(rep('normal',46)), 
                      c(rep('Peripheral Blood Mononuclear Cell',46)),
                      c(rep('False',5), rep('True',5),rep(c(rep('False',6), rep('True',6)), 3)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE100503', 
                      c(rep('disease',13)), 
                      c(rep('Breast',13)), 
                      c(rep('False',13)),
                      c(rep('Breast Ductal Carcinoma In Situ',13)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE100561', 
                      c(rep(c(rep('normal',3), rep('disease',3)),2)), 
                      c(rep('Monocyte',12)),
                      c(rep('False',12)))
# total_metadata<-annot(total_metadata, meta_dir, 'GSE100653', 
#                       c(rep('normal',18)), 
#                       c(rep('Breast',18)), 
#                       c(rep('False',18))) #wrong platform
total_metadata<-annot(total_metadata, meta_dir, 'GSE100940', 
                      c(rep('normal',24)), 
                      c(rep('Blood',24)), 
                      c(rep('False',24)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE101443', 
                      c(rep(c('disease','normal'),4)), 
                      c(rep('Breast',8)), 
                      c(rep('False',8)),
                      c(rep(c('Breast Neoplasm','normal'),4)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE101641', 
                      c(rep('normal',16),rep('disease',32)), 
                      c(rep('Nasal Cavity Epithelium',48)),
                      c(rep('False',48)),
                      c(rep('Normal',16),rep('Cystic Fibrosis',32)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE101658', 
                      c(rep('disease',15)), 
                      c(rep('Hippocampus',15)),
                      c(rep('False',15)),
                      c(rep('Multiple Sclerosis',15)))
to_remove<-append_to_remove(to_remove, 'GSE101673')   
to_remove<-append_to_remove(to_remove, 'GSE101733')
total_metadata<-annot(total_metadata, meta_dir, 'GSE101840', 
                      c(rep('normal',6)), 
                      c(rep('Umbilical Blood',6)),
                      c(rep('False',6)))                      
total_metadata<-annot(total_metadata, meta_dir, 'GSE101961', 
                      c(rep('normal',121)), 
                      c(rep('Breast',121)), 
                      c(rep('False',121)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE102119', 
                      c(rep('disease',146)), 
                      c(rep('Ovarian Carcinoma',146)), 
                      c(rep('False',146)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE102177', 
                      c(rep('normal',18), rep('disease',18)), 
                      c(rep('Blood',36)),
                      c(rep('False',36)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE102504', 
                      c(rep('disease',25)), 
                      c(rep('Peripheral Blood Mononuclear Cell',25)), 
                      c(rep('False',25)),
                      c(rep('Chronic Fatigue Syndrome',25)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE102970', 
                      c(rep('normal',48)), 
                      c(rep('Spermatozoon',48)), 
                      c(rep('False',48)))
to_remove<-append_to_remove(to_remove, 'GSE102994')
total_metadata<-annot(total_metadata, meta_dir, 'GSE103006', 
                      c(rep('normal',24)),
                      c(rep('Umbilical Blood',24)), 
                      c(rep('True',24)))
to_remove<-append_to_remove(to_remove, 'GSE103010')
total_metadata<-annot(total_metadata, meta_dir, 'GSE103413', c(rep('normal', 67)), 
                      c(rep('Brain',9),rep('Liver',10),rep('Breast',7),rep('Leukocyte',36),rep('Melanocyte',3),rep('Chorionic Villus',2)),
                      c(rep('False',67)))
to_remove<-append_to_remove(to_remove, 'GSE103502')
total_metadata<-annot(total_metadata, meta_dir, 'GSE103659', 
                      c(rep('disease',181)),
                      c(rep('Brain',181)),
                      c(rep('False',181)),
                      c(rep('Brain Neoplasm',181)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE103768', 
                      c(rep('normal',57)), 
                      c(rep('Adipose Tissue',57)), 
                      c(rep('False',57))) #obesity study over time
total_metadata<-annot(total_metadata, meta_dir, 'GSE103911', 
                      c(rep('disease', 27),rep('normal',27), rep('disease', 11)), 
                      c(rep('T-Lymphocyte',65)), 
                      c(rep('False',54), rep('True', 11)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE104087', 
                      c(rep('disease',40)),
                      c(rep('Nasal Cavity Epithelium',40)),
                      c(rep('False',40)),
                      c(rep('Asthma',40)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE104287', 
                      c(rep('normal', 48)), 
                      c(rep('Peripheral Blood Mononuclear Cell', 32), rep('Natural Killer Cell',16)),
                      c(rep('True',20), 'False', rep('True', 3), rep('False', 4), 'True', rep('False',3), rep('True',2), 'False', 'True', rep('False', 3),
                        rep('True', 2), rep('False',2), rep('True',2), 'False', 'True', 'False'))
total_metadata<-annot(total_metadata, meta_dir, 'GSE104376', 
                      c(rep('normal',45)), 
                      c(rep('Umbilical Blood',45)), 
                      c(rep('False',45))) #pregnancy anxiety low or high
to_remove<-append_to_remove(to_remove, 'GSE104471')
total_metadata<-annot(total_metadata, meta_dir, 'GSE104728', 
                      c(rep('disease',44)), 
                      c(rep('Brain',44)), 
                      c(rep('False',44)),
                      c(rep('Medulloblastoma',44)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE104770', 
                      c(rep('disease',14)), 
                      c(rep('B-Lymphocyte',14)), 
                      c(rep('False',14)),
                      c(rep('Plasma Cell Leukemia',14)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE104812', 
                      c(rep('normal',48)), 
                      c(rep('Blood',48)), 
                      c(rep('False',48)))
to_remove<-append_to_remove(to_remove, 'GSE105066')
total_metadata<-annot(total_metadata, meta_dir, 'GSE105123', 
                      c(rep('normal',108)), 
                      c(rep('Peripheral Blood Mononuclear Cell',108)), 
                      c(rep('False',19), rep('True', 89)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE105260', 
                      c(rep('disease',26),rep('normal',9),rep('disease',9)),
                      c(rep('Renal',44)),
                      c(rep('False',44)),
                      c(rep('Renal Cell Carcinoma',26),rep('normal',9),rep('Renal Cell Carcinoma',9)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE106089', 
                      c(rep('normal',84)),
                      c(rep('Placenta',84)), 
                      c(rep('False',84))) #low gestational age newborn study
total_metadata<-annot(total_metadata, meta_dir, 'GSE106099', 
                      c(rep('normal',12), rep('disease', 3), rep('normal', 6), rep('disease', 9)),
                      c(rep('endothelial cell',30)),
                      c(rep('False',30))) #placenta endothelial cells
total_metadata<-annot(total_metadata, meta_dir, 'GSE106360', 
                      c(rep('disease',47)),
                      c(rep('Breast',47)),
                      c(rep('False',47)),
                      c(rep('Breast Neoplasm',47)))
# total_metadata<-annot(total_metadata,meta_dir, 'GSE106556', c(rep('normal',6),rep('disease',6),rep('normal',5),rep('disease',6)),
                     # c(rep('Hematopoietic Tissue',6),rep('Myeloid Leukemia',6),rep('Hematopoietic Tissue',5),rep('Myeloid Leukemia',6)),
                     # c(rep('False',23)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE107038', 
                      c(rep('normal',40)),
                      c(rep('Liver',40)),
                      c(rep('False',40)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE107211', 
                      c(rep('normal',15)),
                      c(rep('Blood',15)),
                      c(rep('False',15)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE107226', 
                      c(rep('disease',8), rep('normal',4)),
                      c(rep('Fibroblast',12)),
                      c(rep('False',12)),
                      c(rep('Pulmonary Fibrosis',8),rep('normal',4)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE107351', 
                      c(rep('normal',41),rep('disease',75)),
                      c(rep('Blood',116)),
                      c(rep('False',116)),
                      c(rep('normal',41),rep('Lynch Syndrome',61),rep('MLH1 Gene Mutation',14)))
to_remove<-append_to_remove(to_remove, 'GSE107352')
total_metadata<-annot(total_metadata, meta_dir, 'GSE107737', 
                      c(rep('normal',12),rep('disease',12)), 
                      c(rep('Blood',24)),
                      c(rep('False',24)),
                      c(rep('normal',12),rep('Hypopituitarism',12)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE108058', 
                      c(rep('normal',30)),
                      c(rep('Spermatozoon',30)),
                      c(rep('True',10), rep('False', 20)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE108423', 
                      c(rep('disease',9), rep('normal', 12)),
                      c(rep('Blood',21)),
                      c(rep('False',21)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE108562', 
                      c(rep('normal',6)),
                      c(rep('Mesenchymal Stem Cell',6)),
                      c(rep('False',3), rep('True', 3)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE108567', 
                      c(rep('normal',59)),
                      c(rep('Chorionic Villus',59)),
                      c(rep('False',6), rep('True', 3), rep('False',2), rep('True', 1), rep('False', 1), 'True', 'False', 'True', 'False',
                        rep('True', 8), 'False', rep('True', 6), rep('False', 8), 'True', 'False', 'True', rep('False',3), rep('True',3),
                        rep('False',4), rep('True',5), 'False')) 
                      #two different processing methods, one with more samples currently used
to_remove<-append_to_remove(to_remove, 'GSE109042')
total_metadata<-annot(total_metadata, meta_dir, 'GSE109096', c(rep('disease',11)),c(rep('Cardiac Ventricle',11)), c(rep('False',11)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE109430', 
                      c(rep('normal',12),rep('disease',24)), 
                      c(rep('Leukocyte',36)), 
                      c(rep('False',36)),
                      c(rep('normal',12),rep('Kawasaki Disease',24)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE109446', 
                      c(rep(c('disease','normal'),29)),
                      c(rep('Nasal Cavity Epithelium', 58)), 
                      c(rep('False',58)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE109905', 
                      c(rep('disease',38), rep('normal',31)),
                      c(rep('Blood',69)),
                      c(rep('False',69)),
                      c(rep('Atrial Septal Defect', 38), rep('normal',31)))
to_remove<-append_to_remove(to_remove, 'GSE109914')
to_remove<-append_to_remove(to_remove,'GSE110607')
# total_metadata<-annot(total_metadata, meta_dir, 'GSE110607', c(rep('normal',21)),c(rep('Adipose Tissue',15),rep('Blood',6)), c(rep('False',21)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE110776', 
                      c(rep('disease',24)),
                      c(rep('Bronchoalveolar Lavage Fluid', 24)), 
                      c(rep('False',24)))
to_remove<-append_to_remove(to_remove,'GSE110778')
to_remove<-append_to_remove(to_remove, 'GSE111396')
total_metadata<-annot(total_metadata, meta_dir, 'GSE111428', 
                      c(rep('disease',6)), 
                      c(rep('Ependymoma',6)), 
                      c(rep('False',6)),
                      c(rep('Ependymoma',6))) 
total_metadata<-annot(total_metadata, meta_dir, 'GSE111632', 
                      c(rep('normal',6), rep('disease',6)), 
                      c(rep('Adipose Tissue',12)),
                      c(rep('False',12)))
to_remove<-append_to_remove(to_remove, 'GSE111933')
to_remove<-append_to_remove(to_remove,'GSE112012')
total_metadata<-annot(total_metadata, meta_dir, 'GSE112047', 
                      c(rep('disease',47)),
                      c(rep('Prostatic Tissue',47)),
                      c(rep('False',47)),
                      c(rep('Prostate Neoplasm',47)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE112306', 
                      c(rep('normal',6)),
                      c(rep('Epithelial Cell',3),rep('Mesenchymal',3)),
                      c(rep('True',6)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE112314', 
                      c(rep('normal',22)),
                      c(rep('Saliva',22)),
                      c(rep('False',22))) #range of childhood stress
total_metadata<-annot(total_metadata, meta_dir, 'GSE112696', 
                      c(rep('normal',6),rep('disease',6)), 
                      c(rep('T-Lymphocyte',12)),
                      c(rep('False',12)),
                      c(rep('normal',6),rep('Myasthenia Gravis',6)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE112987', 
                      c(rep('normal',64),rep('disease',39)), 
                      c(rep('Blood',103)),
                      c(rep('False',103)),
                      c(rep('normal',64),rep('Fetal Alcohol Spectrum Disorder',39)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE113012', 
                      c(rep('disease',7),rep('normal',28)),
                      c(rep('Blood',35)),
                      c(rep('False',35)),
                      c(rep('Fetal Alcohol Spectrum Disorder',7),rep('normal',28)))
to_remove<-append_to_remove(to_remove, 'GSE113061')
to_remove<-append_to_remove(to_remove, 'GSE113775')
to_remove<-append_to_remove(to_remove, 'GSE113061')
to_remove<-append_to_remove(to_remove, 'GSE113775')
total_metadata<-annot(total_metadata, meta_dir, 'GSE113967', 
                      c(rep('disease',134)), 
                      c(rep('Blood',134)),
                      c(rep('False',134)))
to_remove<-append_to_remove(to_remove,'GSE114676')
to_remove<-append_to_remove(to_remove, 'GSE114683')
total_metadata<-annot(total_metadata, meta_dir, 'GSE114935', 
                      c(rep('normal',47)), 
                      c(rep('Blood',47)),
                      c(rep('True',47))) #pregnancy
to_remove<-append_to_remove(to_remove, 'GSE115399')
to_remove<-append_to_remove(to_remove, 'GSE115783')
total_metadata<-annot(total_metadata, meta_dir, 'GSE115797', 
                      c(rep(c('normal','disease'),24)), 
                      c(rep('Skin',48)),
                      c(rep('False',48)),
                      c(rep(c('normal','Psoriasis'),24)))
to_remove<-append_to_remove(to_remove,'GSE115852')
total_metadata<-annot(total_metadata, meta_dir, 'GSE115920', 
                      c(rep('normal',6)), 
                      c(rep('Spermatozoon',6)), 
                      c(rep('False',6))) #different age groups
total_metadata<-annot(total_metadata, meta_dir, 'GSE116300', 
                      c(rep('normal',9), rep('disease',35)), 
                      c(rep('Blood',44)),
                      c(rep('False',44)))
to_remove<-append_to_remove(to_remove, 'GSE116754')
to_remove<-append_to_remove(to_remove, 'GSE116924')
total_metadata<-annot(total_metadata, meta_dir, 'GSE117050', 
                      c(rep('normal',19), rep('disease',19)),
                      c(rep('T-Lymphocyte',38)),
                      c(rep('False',38)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE117439', 
                      c(rep('disease',52)),
                      c(rep('Breast',52)), 
                      c(rep('False',52)),
                      c(rep('Breast Carcinoma',52)))
to_remove<-append_to_remove(to_remove, 'GSE117448')
total_metadata<-annot(total_metadata, meta_dir, 'GSE117852', 
                      c(rep('disease',32)),
                      c(rep('Pancreas',32)),
                      c(rep('False',32)),
                      c(rep('Pancreatic Neoplasm',32)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE118132', 
                      c(rep('disease',18)),
                      c(rep('Lung',18)),
                      c(rep('False',18)),
                      c(rep('Lung Carcinoma',18)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE118250', 
                      c(rep('normal',12)),
                      c(rep('Saliva',12)),
                      c(rep('False',6), rep('True',6)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE118260', 
                      c(rep('normal',20)),
                      c(rep('Intestinal Mucosa',10),rep('Saliva',10)),
                      c(rep('False',20)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE118469', 
                      c(rep('normal',6),rep('disease',15)), 
                      c(rep('Peripheral Blood Mononuclear Cell',21)),
                      c(rep('False',21)),
                      c(rep('normal',6), rep('Tuberculosis', 15)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE118570', 
                      c(rep('disease',43)), #drug addiction study
                      c(rep('T-Lymphocyte',43)),
                      c(rep('False',43)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE119078', 
                      c(rep('normal', 4), rep('diseaese', 25), rep('normal', 17), rep('disease', 2), rep('normal', 2), rep('disease',4), rep('normal', 5)), 
                      c(rep('Saliva',59)), 
                      c(rep('False',59)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE119778', 
                      c(rep('disease',34),rep('normal',34)), 
                      c(rep('Blood',68)),
                      c(rep('False',68)),
                      c(rep('Williams Syndrome',34),rep('normal',34)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE119846', 
                      c(rep('normal', 60)), 
                      c(rep('Hematopoietic Stem Cell',60)), 
                      c(rep('True',60)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE120062', 
                      c(rep('normal', 9), rep('disease',29)), 
                      c(rep('Placenta',38)), 
                      c(rep('False',38)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE120307', 
                      c(rep('normal',2), rep('disease',3), rep('normal', 6), rep('disease', 3), rep('normal',2), 'disease', rep('normal',5), 
                        rep('disease', 4), 'normal', 'disease', rep('normal',2), 'disease', rep('normal',2), 'disease'), 
                      c(rep('Blood',34)), 
                      c(rep('False',34)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE124367', c(rep('normal',12)), c(rep('Breast',12)), c(rep('False',12)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE124565', 
                      c(rep('disease',10),rep('normal',12)), 
                      c(rep('Neutrophil',22)),
                      c(rep('False',22)),
                      c(rep('Antiphospholipid Syndrome',10),rep('normal',12)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE124567', 
                      c(rep('normal', 54)), 
                      c(rep('Retina',51),rep('Embryonic Stem Cell',3)), 
                      c(rep('False',29), rep('True',25)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE125605', 
                      c(rep('normal',13), 'disease', 'normal', 'disease', 'normal', rep('disease',2), 'normal', rep('disease',13), 'normal',
                        rep('disease', 1), 'normal', rep('disease', 3), 'normal', 'disease', 'normal'), 
                      c(rep('Placenta',42)), 
                      c(rep('False',42)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE126017', 
                      c(rep('normal',18), rep('disease',36)),
                      c(rep('Spermatozoon',54)), 
                      c(rep('False',54)),
                      c(rep('normal',18), rep('Psoriasis',36)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE128126', 
                      c(rep('normal',15)), 
                      c(rep('Embryonic Stem Cell',15)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE127824', 
                      c(rep('normal',24)), 
                      c(rep('Umbilical Blood',24)), 
                      c(rep('False',24)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE128801', 
                      c(rep('disease',12)), 
                      c(rep('Blood',12)), 
                      c(rep('False',12)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE129266', 
                      c(rep('normal',12)), 
                      c(rep('Mesenchymal Stem Cell',12)), 
                      c(rep('False',12))) #timeseries
total_metadata<-annot(total_metadata, meta_dir, 'GSE130354', 
                      c(rep('disease',12)), 
                      c(rep('Melanoma',12)), 
                      c(rep('False',12)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE131350', 
                      c(rep('disease',48), rep('normal',12)), 
                      c(rep('Adrenal Gland',60)),
                      c(rep('False',60)),
                      c(rep('Adenocarcinoma',48), rep('normal',12)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE132181', 
                      c(rep('normal', 392)), 
                      c(rep(c('Umbilical Blood','Peripheral Blood Mononuclear Cell'),196)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE132399', 
                      c(rep('normal',26), rep('disease',28)), 
                      c(rep('Liver', 54)),
                      c(rep('True',5), rep('False',49)),
                      c(rep('normal', 26), rep('Hepatoblastoma',28)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE137716', 
                      c(rep('normal',16),rep('disease',17)), 
                      c(rep('Tracheal Epithelium',33)),
                      c(rep('False',33)),
                      c(rep('normal',16),rep('Asthma',17)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE138279', 
                      c(rep('normal',65)), 
                      c(rep('Saliva',65)), 
                      c(rep('False',65))) #childhood trauma
total_metadata<-annot(total_metadata, meta_dir, 'GSE139307', 
                      c(rep('normal',37)), 
                      c(rep('Spermatozoon',37)) , 
                      c(rep('False',37))) #dioxin exposure
total_metadata<-annot(total_metadata, meta_dir, 'GSE141338', 
                      c(rep('normal',6),rep('disease',42)), 
                      c(rep('Breast',48)),
                      c(rep('False',48)),
                      c(rep('normal',6),rep('Breast Neoplasm',42)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE142554', 
                      c(rep('disease',6), rep('disease',6)), 
                      c(rep('Blood',12)), 
                      c(rep('False',12)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE143209', 
                      c(rep('normal',64)), 
                      c(rep('Islet of Langerhans',64)), 
                      c(rep('False',64)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE144894', 
                      c(rep('disease', 120)), 
                      c(rep('Peripheral Blood Mononuclear Cell',120)),
                      c(rep('False',120)),
                      c(rep('Chronic Lymphocytic Leukemia',120)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE144977', 
                      c(rep('normal',89)),
                      c(rep('Placenta',89)),
                      c(rep('False',89)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE147740', 
                      c(rep('normal',1129)), 
                      c(rep('Blood',1129))) #wide variety
total_metadata<-annot(total_metadata, meta_dir, 'GSE150901', 
                      c(rep('normal',48)), 
                      c(rep('Buccal Mucosa', 48)),
                      c('False', rep(c(rep('True',6), rep('False',2)), 5), rep('True', 6), 'False')) #assistive reproduction technologies
total_metadata<-annot(total_metadata, meta_dir, 'GSE151278', 
                      c(rep('disease',70)), 
                      c(rep('Blood',70)),
                      c(rep('True',70)),
                      c(rep('Psoriasis',70)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE151042', 
                      c(rep('normal',492)), 
                      c(rep('Umbilical Blood',492)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE151407', 
                      c(rep('normal',78)), 
                      c(rep('Vastus Lateralis',78)),
                      c(rep(c('False', 'True', 'True'),25), rep('True',3)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE151485', 
                      c(rep('normal', 100)), 
                      c(rep('Saliva',100)),
                      c('False', rep('True',2), rep('False',2), rep('True',2), rep('False',1), rep('True',2), rep('False',2), rep('True',2), 'False',
                        'True', 'False', rep(c('True','True','False'),6), 'False', rep('True',2), 'False', rep('True',2), 
                        rep(c(rep('False',3), rep('True',2)),2), rep(c('False', 'True', 'True'),6), 'False', rep(c('False', 'True', 'True'),2),
                        'False', 'True', rep(c('False', 'True', 'True'),3), 'False', 'True', 'False', rep('True',4), rep('False',2), 
                        'True', rep('False',2), 'True'))
total_metadata<-annot(total_metadata, meta_dir, 'GSE151600', 
                      c(rep('normal', 16)), 
                      c(rep('Skin',16)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE152380', 
                      c(rep('normal',90)), 
                      c(rep('Umbilical Blood',90)),
                      c(rep('True',18), rep('False', 17), rep('True',55)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE47915', 
                      c(rep('disease',8)), 
                      c(rep('Prostatic Tissue',8)), 
                      c(rep('False',8)),
                      c(rep('Prostate Neoplasm',8)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE48472', 
                      c(rep('normal',56)), 
                      c(rep('Blood',6),rep('Liver',5),rep('Muscle',6),rep('visceral fat',6),rep('Pancreas',4),
                        rep('subcutaneous adipose tissue',6),rep('Spleen',3),rep('Blood',5),rep('Buccal Mucosa',5),
                        rep('Hair',5),rep('Saliva',5)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE49618', 
                      c(rep('normal',21)), 
                      c(rep('Bone Marrow',21))) #different stem/immune celltypes
total_metadata<-annot(total_metadata, meta_dir, 'GSE52025', 
                      c(rep('normal',62)), 
                      c(rep('Fibroblast',62)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE59038', 
                      c(rep('disease',24)), 
                      c(rep('Colon',24)),
                      c(rep('True',24)),
                      c(rep('Colon Adenocarcinoma', 24)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE59524', 
                      c(rep('normal',24)), 
                      c(rep('visceral fat',24))) #visceral
total_metadata<-annot(total_metadata, meta_dir, 'GSE60655', 
                      c(rep('normal',36)), 
                      c(rep('Vastus Lateralis', 36)),
                      c(rep(c('False','True'), 18)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE61107', 
                      c(rep('disease',23), rep('normal',24), 'disease'), 
                      c(rep('Frontal Lobe Cortex', 48))) #post mortem
total_metadata<-annot(total_metadata, meta_dir, 'GSE61278', 
                      c(rep('normal', 66), rep('disease',44)), 
                      c(rep('Liver', 110))) #mix of fetal and adult
total_metadata<-annot(total_metadata, meta_dir, 'GSE61446', 
                      c(rep('normal',67)), 
                      c(rep('Liver',67))) #obesity
total_metadata<-annot(total_metadata, meta_dir, 'GSE61450', 
                      c(rep('normal',71)), 
                      c(rep('subcutaneous adipose tissue',71))) #obesity
total_metadata<-annot(total_metadata, meta_dir, 'GSE61452', 
                      c(rep('normal',60)), 
                      c(rep('Muscle',60))) #obesity
total_metadata<-annot(total_metadata, meta_dir, 'GSE61453', 
                      c(rep('normal',71)), 
                      c(rep('visceral fat',71))) #obesity
total_metadata<-annot(total_metadata, meta_dir, 'GSE62727', 
                      c(rep('disease', 7), rep('normal',4)), 
                      c(rep('Left Atrium',11)),
                      c(rep('False', 11)),
                      c(rep('Atrial Fibrillation', 7),rep('normal',4)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE63179', 
                      c(rep('normal',8)), 
                      c(rep('Cerebellum',8)),
                      c('False','True','False','True','False','False','True','True'))                      
total_metadata<-annot(total_metadata, meta_dir, 'GSE64096', 
                      c(rep('normal',40)), 
                      c(rep('Spermatozoon',40))) #two different density
total_metadata<-annot(total_metadata, meta_dir, 'GSE66351', 
                      c(rep('disease', 190)), 
                      c(rep('Frontal Lobe Cortex', 190)),
                      c(rep('False',190)),
                      c(rep('Alzheimer\'s Disease', 190)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE72245', 
                      c(rep('disease',118)), 
                      c(rep('Breast',118)),
                      c(rep('True',118)),
                      c(rep('Breast Carcinoma',118)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE72251', 
                      c(rep('disease',119)), 
                      c(rep('Breast',119)),
                      c(rep('True',119)),
                      c(rep('Breast Carcinoma',119)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE72254', 
                      c(rep('disease',58)), 
                      c(rep('Breast',58)),
                      c(rep('True',58)),
                      c(rep('Breast Carcinoma',58)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE74167', 
                      c(rep('normal',21), rep('disease',21)), 
                      c(rep('Mesenchymal Stem Cell',42)),
                      c(rep('False',42)),
                      c(rep('Mucous Connective Tissue',42)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE74214', 
                      c(rep('normal',18)), 
                      c(rep('Breast',18)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE75133', 
                      c(rep('normal',15)), 
                      c(rep('Mesenchymal Stem Cell',15)), 
                      c(rep('True',15)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE75405', 
                      c(rep('normal',24)), 
                      c(rep('T-Lymphocyte',24)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE75443', 
                      c(rep('disease',12)), 
                      c(rep('Glioblastoma',12)),
                      c(rep('False',12)),
                      c(rep('Glioblastoma',12)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE75704', 
                      c(rep('normal',72), rep('disease',94)), 
                      c(rep('Frontal Lobe Cortex',166)),
                      c(rep('False',166)),
                      c(rep('normal',72), rep('Progressive Supranuclear Palsy', 94))) #post mortem
total_metadata<-annot(total_metadata, meta_dir, 'GSE76372', 
                      c(rep('normal',9)), 
                      c(rep('Lymphocyte',2),'Induced Pluripotent Stem Cell', rep('Lymphocyte',3), rep('Induced Pluripotent Stem Cell',2), 'Lymphocyte'),
                      c('False', rep('True',3), 'False', rep('True',3), 'False'))
total_metadata<-annot(total_metadata, meta_dir, 'GSE76503', 
                      c(rep('normal',48)), 
                      c(rep('Blood',48)), 
                      c(rep('False',12), rep('True',12),rep('False',12), rep('True',12)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE77269', 
                      c(rep('disease',60)), 
                      c(rep('Lung',60)),
                      c(rep('False',60)),
                      c(rep('Hepatocellular Carcinoma',60)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE78732', 
                      c(rep('normal',10),rep('disease',19)),
                      c(rep('Liver',29)),
                      c(rep('False',29)),
                      c(rep('normal',10),rep('Hepatoblastoma',19)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE79009', 
                      c(rep('disease',125)),
                      c(rep('Schwannoma',125)),
                      c(rep('False',125)),
                      c(rep('Schwannoma',125)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE79329', 
                      c(rep('normal',34)), 
                      c(rep('Leukocyte',34))) #pollutant study
total_metadata<-annot(total_metadata, meta_dir, 'GSE80685', 
                      c(rep('disease',172)), 
                      c(rep('Prostatic Tissue',172)),
                      c(rep('False',172)),
                      c(rep('Prostate Neoplasm',172)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE82084', 
                      c(rep('disease',16), rep('normal',20)), 
                      c(rep('Umbilical Blood', 16), rep(c('Granulocyte','Monocyte','Nucleated Red Blood Cell','T-Lymphocyte'),5)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE82234', 
                      c(rep('normal',6)), 
                      c(rep('Umbilical Blood',6))) #different passages
total_metadata<-annot(total_metadata, meta_dir, 'GSE83917', 
                      c(rep('disease',160)), 
                      c(rep('Prostatic Tissue',160)),
                      c(rep('False',160)),
                      c(rep('Prostate Neoplasm',160)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE83933', 
                      c(rep('disease',39)), 
                      c(rep('Meninges',39)),
                      c(rep('True',20), rep('False',19)),
                      c(rep('Meningioma',39)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE83944', 
                      c(rep('normal',48)), 
                      c(rep('Neutrophil',48))) #same sample, different time of day
total_metadata<-annot(total_metadata, meta_dir, 'GSE84274', 
                      c('normal',rep('disease',6), 'normal','disease','normal',rep('disease',2),'normal',rep('disease',2),rep('normal',2),
                        rep('disease',7)), 
                      c(rep('Ascending Aorta', 24))) 
total_metadata<-annot(total_metadata, meta_dir, 'GSE85042', 
                      c(rep('normal',71)), 
                      c(rep('Umbilical Blood',71))) #gestational stimuli study
total_metadata<-annot(total_metadata, meta_dir, 'GSE85210', 
                      c(rep('normal',253)), 
                      c(rep('Blood',253))) #smoking study
total_metadata<-annot(total_metadata, meta_dir, 'GSE85506', 
                      c(rep('disease',47)), 
                      c(rep('Blood',47)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE85845', 
                      c(rep('disease',8), rep('normal',8)),
                      c(rep('Lung',16)),
                      c(rep('False',16)), 
                      c(rep('Lung Adenocarcinoma',8), rep('normal',8)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE86078', 
                      c(rep('disease',149)), 
                      c(rep('Colonic Mucosa', 149)),
                      c(rep('False',149)),
                      c(rep('Colorectal Carcinoma', 149)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE86258', 
                      c(rep('disease',7),rep('normal',7)), 
                      c(rep('Fibroblast',14)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE87582', 
                      c(rep('disease',21)), 
                      c(rep('Blood',21)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE87640', 
                      c(rep('disease', 240)), 
                      c(rep('Leukocyte',240))) #complicated, some healthy mixed
total_metadata<-annot(total_metadata, meta_dir, 'GSE88883', 
                      c(rep('normal',100)), 
                      c(rep('Breast',100)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE89218', 
                      c(rep('normal',163)), 
                      c(rep('Blood',163))) #PTSD study
total_metadata<-annot(total_metadata, meta_dir, 'GSE89251', 
                      c(rep('disease',136)), 
                      c(rep('T-Lymphocyte',136)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE89472', 
                      c(rep(c('disease','normal'),5)), 
                      c(rep('Blood',10)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE89803', 
                      c(rep('disease',138), rep('normal',4)), 
                      c(rep('Bile Duct', 142)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE89852', 
                      c(rep('normal',37), rep('disease',37)), 
                      c(rep('Liver', 74)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE89925', 
                      c(rep('normal',3),rep('disease',3)), 
                      c(rep('Fibroblast',6)),
                      c(rep('False',6)),
                      c(rep('normal',3), rep('Teratoma',3)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE90060', 
                      c(rep('normal', 34)), 
                      c(rep('Endometrium',34))) #differen part of implanting
total_metadata<-annot(total_metadata, meta_dir, 'GSE92577', 
                      c(rep('disease',12)), 
                      c(rep('Brain',12)),
                      c(rep('False',12)),
                      c(rep('Brain Neoplasm',12)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE92909', 
                      c(rep('disease',6)), 
                      c(rep('Breast', 6)),
                      c(rep('False',6)),
                      c(rep('Breast Carcinoma',6)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE93208', 
                      c(rep('normal', 19)), 
                      c(rep('Chorionic Villus', 19)),
                      c(rep(c('True','False'),6), rep('False',3), rep('True',4)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE94326', 
                      c(rep('normal', 8)), 
                      c(rep('Brain',4),rep('Prostatic Tissue',4)), 
                      c(rep('False',2),rep('True',6)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE95486', 
                      c(rep('normal',21), rep('disease',3)), 
                      c(rep('Blood',24)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE97784', 
                      c(rep(list('disease','normal'),12)), 
                      c(rep('Buccal Mucosa',24)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE98056', 
                      c(rep('normal',69)),
                      c(rep('Leukocyte',69)),
                      c(rep('False',35),rep('True',34)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE98203', 
                      c(rep('normal',88)), 
                      c(rep('Neuron',88)),
                      c(rep('True', 2), rep('False',3), 'True', rep('False',3), rep('True',3), 'False',rep('True',6), 'False', rep('True',5), 'False',
                        rep('True',3), 'False', rep('True',3), 'False','True','False', rep('True',2),'False',rep('True',5),'False',
                        rep('True',3), 'False','True','False',rep('True',6), 'False',rep('True',4), rep('False',2),
                        'True','False','True',rep('False',3),'True','False',rep('True',2),'False',rep('True',5), rep('False',4), 
                        rep('True',4))) #postmortem. heroine, suicide
total_metadata<-annot(total_metadata, meta_dir, 'GSE99553', 
                      c('disease', rep('normal',5), rep('disease',3), 'normal',rep('diseae',5), 'normal','normal','disease', 'normal','disease','disease',
                      rep('normal',3), 'disease','normal','disease',rep('normal',2), rep('disease',2),'normal',rep('disease',2),rep('normal',5), 'disease',
                      rep('normal',2),'disease',rep('normal',4), 'disease', rep('normal',2), rep('disease',5),'normal',rep('disease',2), rep('normal',2), 
                      'disease', rep('normal',2), rep('disease',4), 'normal','disease','normal',rep('disease',2), 'normal',rep('disease',2),rep('normal',2),
                      rep('disease',2), rep('normal',2), rep('disease',2), 'normal'),
                      c(rep('Gastric Mucosa',84)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE99624', 
                      c(rep('disease',32),rep('normal',16)), 
                      c(rep('Blood',48)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE99755', 
                      c(rep('normal',67)), 
                      c(rep('Blood',67))) #police PTSD

##Nov 2021
total_metadata<-annot(total_metadata, meta_dir, 'GSE103502', 
                      c(rep(c('normal', 'disease'), 4)), 
                      c(rep('Lung', 8)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE104471', 
                       c(rep('disease',3), rep('normal', 2), rep('disease',2), 
                         'normal', 'disease', rep('normal',4), rep('disease',5), 
                         'normal', 'disease', rep('normal',4), rep('disease', 3),
                         rep('normal',2), rep('disease',2), 'normal', 'disease', rep('normal',4),
                         rep('disease',5), 'normal', 'disease', rep('normal',4), rep('disease',2), 
                         'disease', rep('normal',2), rep('disease',2),'normal','disease',rep('normal',4),rep('disease',5),
                         'normal','disease',rep('normal',4)),
                       c(rep('Peripheral Blood Mononuclear Cell', 24), rep('Nasal Cavity Epithelium', 24), rep('Bronchial Epithelium', 24)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE110724', 
                      c(rep('normal',21)),
                      c(rep('subcutaneous adipose tissue', 6), rep('visceral fat', 9), rep('Leukocyte',6)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE111396', 
                      c(rep('disease', 4), 'normal', rep('disease', 6), rep('normal', 3), rep('disease', 2), 'normal', 'disease', 
                        rep('normal',2), rep('disease',3), rep('normal', 9), rep('disease',2), 'normal', 'disease', 'normal', 
                        rep('disease', 4), 'normal', rep('disease', 2),'normal', rep('disease',4),rep('normal',2), 
                        'disease','normal', 'disease','normal', rep('disease',2), rep('normal',2), rep('disease',2)),
                      c(rep('Fibroblast', 10), rep('Airway', 15), rep('Fibroblast', 36)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE111942',
                      c(rep('disease', 25), rep('normal', 18)),
                      c(rep('Peripheral Blood Mononuclear Cell', 43)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE112012',
                      c('disease', 'normal'),
                      c(rep('Esophagus', 2)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE113061',
                      c(rep('normal',3), rep('disease', 5)),
                      c(rep('Trachea Smooth Muscle Tissue', 8)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE113775', c(rep(c('disease','normal'),1)), c(rep('Esophagus', 2)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE114753', 
                      c(rep('normal',156)), 
                      c(rep('Spermatozoon', 156))) #smoking study
total_metadata<-annot(total_metadata, meta_dir, 'GSE115783', 
                      c(rep('disease', 15), rep('normal', 2), rep('disease', 3), 
                        'normal', rep('disease',3), rep('normal',3), rep('disease',13)),
                      c(rep('Pituitary Gland', 40)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE116375',
                      c(rep('normal', 16)),
                      c(rep('Umbilical Vein', 4), rep('Bone Marrow', 12)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE116754',
                      c(rep('normal', 38)),
                      c(rep('Embryonic Stem Cell', 9),
                        rep('Aorta Smooth Muscle Tissue', 1),
                        rep('Brain', 10),
                        rep('Pancreas', 4),
                        rep('Islet of Langerhans', 5),
                        rep('Umbilical Cord', 3),
                        rep('Small Intestine', 3), 
                        'Thymus Gland', 'Diaphragm', 'Duodenum'))
total_metadata<-annot(total_metadata, meta_dir,'GSE118696',
                      c(rep('normal', 24)),
                      c(rep(c('Monocyte',rep('Macrophage',7)), 2),
                        'Monocyte',rep('Macrophage',6), 'Monocyte'),
                      c(rep('False',2), rep('True',6), rep('False',2), rep('True',6), rep('False',3), rep('True',4), 'False'))
total_metadata<-annot(total_metadata, meta_dir, 'GSE120250',
                      c(rep('normal',88)),
                      c(rep('Placenta',88)),
                      c(rep('True',44), rep('False',44)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE121192',
                      c(rep('normal',4),
                        rep('disease', 10),
                        rep('normal', 6),
                        rep('disease', 10),
                        rep('normal', 6),
                        rep('disease', 10)),
                      c(rep('T-Lymphocyte', 30),
                        rep('Monocyte', 16)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE123678',
                      c(rep('normal', 8), rep('disease', 70)),
                      c(rep('Corpus Callosum', 5), rep('Frontal Lobe Cortex', 3), rep('Glioma', 70)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE124366',
                      c(rep('normal',215)),
                      c(rep('Buccal Mucosa', 20),
                        rep('Peripheral Blood Mononuclear Cell', 13),
                        rep('Buccal Mucosa', 5),
                        rep('Peripheral Blood Mononuclear Cell', 6),
                        rep('Buccal Mucosa', 48),
                        rep('Peripheral Blood Mononuclear Cell', 86),
                        rep('Buccal Mucosa', 37)
                        ))
total_metadata<-annot(total_metadata, meta_dir, 'GSE124368',
                      c(rep('disease',4)),
                      c(rep('Breast',4)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE127857',
                      c(rep('disease', 24), rep('normal', 34), rep('disease', 53)),
                      c(rep('Tumor Tissue', 11), rep('Intestine', 13), rep('Tissue', 11), rep('Intestine', 13), rep('Mucosa', 10), rep('Unknown',53)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE130029',
                      c(rep('disease', 20), rep('normal', 11)),
                      c(rep('T-Lymphocyte', 31)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE135097',
                      c(rep('normal', 6), rep('disease', 6)),
                      c(rep('Fibroblast', 12)),
                      c(rep(c(rep('False', 2), rep('True', 4)), 2)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE136704',
                      c(rep(c('normal','disease'), 11), rep('disease', 2), 
                        rep(c('normal', 'disease'), 5), rep('disease', 3), 'normal', rep('disease', 6)),
                      c(rep('Oropharyngeal Tissue', 44)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE137134',
                      c(rep('normal',6), rep('disease', 6)),
                      c(rep('Fibroblast', 12)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE146552',
                      c(rep('disease', 20), rep('normal', 18)),
                      c(rep('Ovary', 20), rep('Fallopian Tube', 8), rep('Ovary', 2), rep('Ovarian Surface Epithelium', 8)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE148663',
                      c(rep('disease', 22), rep('normal', 10)),
                      c(rep('Leukocyte', 32)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE149395',
                      c(rep('normal', 11)),
                      c(rep('Pancreas',11)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE152204',
                      c(rep('disease', 41), rep('normal', 48)),
                      c(rep('Blood', 89)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE154971',
                      c(rep('disease', 16), rep('normal', 8)),
                      c(rep('Lymphocyte', 24)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE155353',
                      c(rep('normal',34), rep('disease',82)),
                      c(rep('Pancreas', 116)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE156669',
                      c(rep('normal', 5), rep('disease', 7)),
                      c(rep('Buccal Mucosa', 12)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE156792',
                      c(rep('disease', 7), rep('normal', 6), rep('disease',6), rep('normal', 6), rep('disease', 18), 
                        rep(c(rep('normal', 6), rep('disease',6)),4), rep('normal', 4), 'disease', 'normal', rep('disease',5),
                        rep('normal', 7)),
                      c(rep('T-Lymphocyte', 109)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE156994',
                      c(rep('normal', 105), rep('disease', 114)),
                      c(rep('Blood', 219)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE157272',
                      c(rep('normal',10), rep('disease',35), rep('normal',3)),
                      c(rep('Prostatic Tissue',48)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE157651',
                      c(rep('normal',23), rep('disease',17)),
                      c(rep('Fibroblast',15), rep('Airway',17), rep('Fibroblast',8)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE162484',
                      c(rep(c('disease', 'normal'),4), 'disease', 'disease', 'normal', rep('disease', 4)),
                      c(rep('Cartilage', 15)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE163372',
                      c('normal', rep('disease',3), 'normal', rep('disease',2), 'normal', 'disease', 'normal', 
                        rep('disease', 2), 'normal', rep('disease',3),'normal',rep('disease',2), 'normal', 'disease'),
                      c(rep('Kidney',21)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE164988',
                      c(rep(c('disease','normal'),12)),
                      c(rep('Stomach', 24)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE166611',
                      c(rep('normal',32)),
                      c(rep('Blood', 32)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE166652',
                      c(rep('normal',28), rep('disease', 28)),
                      c(rep('Myoblast',56)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE166787',
                      c(rep('normal',28), rep('disease', 28)),
                      c(rep('Myoblast',56)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE174422',
                      c(rep('normal',256)),
                      c(rep('Blood', 256)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE178212',
                      c(rep(c('normal','disease'),11), 'disease', 'normal', rep('disease',3), 'normal', 
                        rep('disease',2), 'normal', 'disease', 'normal', rep('disease',2), 'normal', rep('disease',4)),
                      c(rep('Esophageal Squamous Epithelium', 40)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE178216',
                      c('normal', rep('disease',2), rep(c('normal','disease'),3), rep('disease',4), 'normal', 
                        rep('disease',2), 'normal','disease','normal',rep('disease',3)),
                      c(rep('Oral Cavity Epithelium', 22)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE74071',
                      c('disease', 'normal', 'disease', rep('normal', 2), rep('disease', 9), rep(c('normal', 'disease'), 2), 
                        rep(c('disease', 'normal'), 3), rep('disease',4)),
                      c(rep('Pancreas', 28)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE79185',
                      c(rep('disease', 61)),
                      c(rep('Breast',5), rep('Brain', 6), rep('Colon', 7), rep('Leukemia',6), rep('Lung',9), rep('Melanoma', 10), 
                        rep('Ovary', 7), rep('Prostatic Tissue', 3), rep('Kidney', 8)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE81846',
                      c(rep('disease',10), rep('normal',6)),
                      c(rep('Leukocyte',16)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE83691',
                      c('normal', rep('disease',5), rep('normal',2), rep('disease',9), 'normal', rep('disease',8)),
                      c(rep('Liver',26)))
# total_metadata<-annot(total_metadata, meta_dir, 'GSE85566',
#                       c(rep('disease',11), rep('normal',8), rep(c('disease','normal'),2), rep('disease',6), 'normal','disease',
#                         rep('normal',3),rep(c('disease','normal'),2), rep('disease',2), 'normal', 'disease', rep('normal',3), 
#                         rep(c('disease','normal'),2), rep('disease',2), rep(c('normal','disease'),2), rep('normal',2),
#                         'disease', 'normal', rep('disease',3), 'normal', rep('disease',5), rep(c('normal', 'disease'),2),
#                         rep('disease',2), 'normal', rep('disease',6), 'normal', rep('disease',8), rep('normal',2),
#                         rep('disease',2), 'normal','disease', rep('normal',4), 'disease','normal', rep('disease',2), 
#                         rep(c('normal','disease'),2), rep('disease',2), rep('normal',2), rep('disease',3),
#                         rep('normal',2), rep('disease',5), 'normal', rep('disease',3)),
#                       c(rep('Tracheal Epithelium',115)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE116338',
                      c(rep('normal',10), rep('disease',38)),
                      c(rep('Prostatic Tissue',48)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE130030',
                      c(rep('disease',14), rep('normal',14)),
                      c(rep('T-Lymphocyte', 28)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE157341',
                      c(rep('disease', 239), rep('normal', 35)),
                      c(rep('Liver',274)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE166069',
                      c(rep('normal',106)),
                      c(rep('Melanocyte',106)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE178218',
                      c(rep(c('disease', 'normal'),2), rep('disease',3), rep(c('disease', 'normal'),3), rep('disease', 2), rep('normal', 2), 
                        rep('disease',3), 'normal', 'disease', rep('normal',2), rep('disease',4), 'normal', rep('disease',2)),
                      c(rep('Larynx',31)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE67170',
                      c(rep('normal',20), rep('disease',69)),
                      c(rep('T-Lymphocyte',10), rep('Peripheral Blood Mononuclear Cell', 10), rep('Peripheral Blood Mononuclear Cell',69)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE81211',
                      c(rep('normal',3), rep('disease',9)),
                      c(rep('Colon',12)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE85938',
                      c(rep('normal', 31), rep('disease',2)),
                      c(rep('Astrocyte', 33)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE87095',
                      c(rep('disease',3), rep('normal', 10), rep('disease',3), rep('normal',10), rep('disease',9), rep('normal', 11), rep('disease', 5),
                        rep('normal',10), rep('disease', 7), rep('normal',7), rep('disease',6), rep('normal',6), rep('disease',6), rep('normal',6),
                        rep('disease',6), rep('normal',6), 'disease', 'normal', 'disease', 'normal', 'disease', rep('normal',5),'disease'),
                      c(rep('B-Lymphocyte',122)))

total_metadata<-annot(total_metadata, meta_dir, 'GSE93266',
                      c(rep(c(rep('disease',6), rep('normal',6)),4), rep('disease',6), rep('normal',2), rep('disease',19)),
                      c(rep('Peripheral Blood Mononuclear Cell', 75)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE94462',
                      c(rep('disease', 10), rep('normal',4), rep('disease',2)),
                      c(rep('Cornea', 16)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE95761',
                      c(rep('normal', 6)),
                      c(rep('Embryonic Stem Cell', 6)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE99511',
                      c(rep('normal', 28), rep('disease', 40)),
                      c(rep('Cervix Squamous Epithelium', 68)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE101908',
                      c(rep('normal',2), rep(c('disease','normal'),3), rep('disease',2), rep('normal',3), rep('disease',2), rep('normal',2),
                        rep('disease',2), rep('normal',5), 'disease', rep('normal', 3), rep('disease', 5), 'normal', rep('disease',3), 'normal',
                        rep('disease',2), 'normal', 'disease'),
                      c(rep('Prostatic Tissue', 42)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE103541',
                      c(rep('normal',145)),
                      c(rep(c('B-Lymphocyte', rep('T-Lymphocyte',2), 'Granulocyte', 'Monocyte'), 29)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE112905',
                      c(rep('normal',4), rep(c('disease','normal'), 4), rep('diease',3), 'normal', rep('disease',2), 'normal', 'disease'),
                      c(rep('Peripheral Blood Mononuclear Cell', 20)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE114763',
                      c(rep('normal',40)),
                      c(rep('Vastus Lateralis', 40)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE115382',
                      c(rep('normal',4), rep('disease', 4)),
                      c(rep('Astrocyte', 8)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE119617',
                      c(rep('normal',5), rep('disease',21)),
                      c(rep('Bone Marrow', 26)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE120854',
                      c(rep('normal',10), rep('disease',24)),
                      c(rep('Myometrium',34)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE122148',
                      c(rep('normal',5), rep('disease',3)),
                      c(rep('Endometrial Tissue', 8)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE122408',
                      c(rep('normal', 180)),
                      c(rep('Peripheral Blood Mononuclear Cell',180)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE123003',
                      c(rep('disease',8), rep('normal',8)),
                      c(rep('T-Lymphocyte',16)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE123914',
                      c(rep('normal', 69)),
                      c(rep('Blood', 69)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE123995',
                      c(rep('normal',56)),
                      c(rep('Hepatocyte',56)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE125367',
                      c(rep('normal', 27), rep('disease', 17)),
                      c(rep('Blood', 44)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE130748',
                      c(rep('normal',2), rep('disease',2), rep('normal',2), rep('disease',2), rep('normal',2), rep('disease',2), rep('normal',4),
                        rep('disease',4), rep('normal',4), rep('disease',2), rep('normal',6), rep('normal',2), rep('disease',3)),
                      c(rep('Leukocyte',37)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE132547',
                      c(rep('normal',24)),
                      c(rep('Bronchoalveolar Lavage Fluid',24))) #weird lung holdout set
total_metadata<-annot(total_metadata, meta_dir, 'GSE132866',
                      c(rep('disease', 5), rep('normal', 3)),
                      c(rep('Leukocyte', 8)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE133062',
                      c(rep('normal', 70)),
                      c(rep('Alveolar Cell', 70)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE133774',
                      c(rep('normal',6), rep('disease', 4)),
                      c(rep('Blood', 10)))

# total_metadata<-annot(total_metadata, meta_dir, '',
#                       c(),
#                       c())


##combine and save##
df <- apply(total_metadata,2,as.character)
write.csv(df, file=paste0('data/', database_type, '/all_metadata_annotated_sep2022.txt'),
          row.names=F,
          sep="\t",
          quote=F)
          
gzip(filename=paste0('data/', database_type, '/all_metadata_annotated_sep2022.txt'), 
     destname=paste0('data/', database_type, '/all_metadata_annotated_sep2022.txt.gz'), overwrite=TRUE, remove=TRUE)

