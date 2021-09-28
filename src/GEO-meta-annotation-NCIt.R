library(data.table)
meta_dir='data/GEO/'
database_type='GEO'

annot<-function(total_metadata, meta_dir, sample, disease_vec, tissue_vec, 
                treatment_vec=c(rep('False',length(disease_vec)))){
  sample_matrix<-read.csv(paste0(meta_dir, sample,'_sample_metadata.txt.gz'), sep='\t')
  if (!(sample %in% total_metadata$series)){
    total_metadata=rbind(total_metadata, 
                         data.table(sample_id=sample_matrix$Samples, 
                                    series=sample_matrix$Series, 
                                    platform=sample_matrix$Platform,
                                    disease_name=disease_vec, 
                                    tissue_name=tissue_vec, 
                                    treatment=treatment_vec))
  }
  return(total_metadata)
}

append_to_remove<-function(list, item){
  if (!(item %in% list)){list<-append(list, item)}
  return(list)
}

total_metadata<-data.table(sample_id=character(), series=character(), platform=character(), disease_name=character(), tissue_name=character(), treatment=character())
to_remove<-list()

total_metadata<-annot(total_metadata, meta_dir, 'GSE102031', c(rep('normal',16)), c(rep('Embryonic Stem Cell', 16)), c(rep('False', 16)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE100780', c(rep('normal', 4)), c(rep('HL60', 4)), c(rep('False',2), rep('True',2))) 
total_metadata<-annot(total_metadata, meta_dir, 'GSE103006', c(rep('normal', 43)), c(rep('Hematopoietic Stem Cell', 43)), c(rep('False',43)))
to_remove<-append_to_remove(to_remove, 'GSE103027')
to_remove<-append_to_remove(to_remove, 'GSE103271')
to_remove<-append_to_remove(to_remove, 'GSE103279')
to_remove<-append_to_remove(to_remove, 'GSE103280')
to_remove<-append_to_remove(to_remove, 'GSE103287')
to_remove<-append_to_remove(to_remove, 'GSE103328')
to_remove<-append_to_remove(to_remove, 'GSE104359')
total_metadata<-annot(total_metadata, meta_dir, 'GSE106360', c(rep('disease', 47)), c(rep('Breast Carcinoma', 47)), c(rep('False', 47)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE106437', c(rep('disease', 3)), c(rep('Colorectal Adenocarcinoma', 3)), c('True', rep('False', 2)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE106438', c(rep('disease', 3)), c(rep('Colon Carcinoma', 3)), c('True', rep('False', 2)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE106556', 
                      c(rep(list('disease','normal'), 10)), 
                      c(rep(list('Rectal Neoplasm', 'Rectal'), 9), 'Colon Neoplasm','Sigmoid Colon'), 
                      c(rep('False', 20)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE106727', c(rep('disease', 18)), c(rep('Brain Neoplasm',18)), c(rep('False', 18)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE108123', c(rep('disease', 87)), c(rep('Lung Squamous Cell Carcinoma', 87)), c(rep('False', 87)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE108143', c(rep('normal', 2)), c(rep('Embryonic Stem Cell',2)), c('False','True'))
total_metadata<-annot(total_metadata, meta_dir, 'GSE108202', c(rep('normal', 3), rep('disease',3), rep('normal',2)), c(rep('Fallopian Tube',8)), 
                      c(rep('False',3)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE108785', c(rep('disease',3), rep('normal',3)), c(rep('Whole Blood', 6)), c(rep('False',6)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE108982', c(rep('disease',16)), c(rep('Breast Neoplasm',16)), c(rep('False',16)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE109330', c(rep('disease',27)), c(rep('Glioblastoma',27)), c(rep('False',27)))
to_remove<-append_to_remove(to_remove, 'GSE109364')
total_metadata<-annot(total_metadata, meta_dir, 'GSE109541', c(rep('disease',4)), c(rep('Gastric Adenocarcinoma',4)), c(rep('False',4)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE109904', c(rep('disease',6)), c(rep(list('Sarcoma','Leukemia',3))), c(rep('False',6)))
to_remove<-append_to_remove(to_remove, 'GSE110184')
to_remove<-append_to_remove(to_remove, 'GSE110697')
total_metadata<-annot(total_metadata, meta_dir, 'GSE100197', 
                      c(rep('dissease',51),rep('normal',51)),
                      c(rep('Eclampsia',22),rep('Placenta',11),rep('Eclampsia',18),rep('Placenta',51)), 
                      c(rep('False',102)))
to_remove<-append_to_remove(to_remove, 'GSE100249')
total_metadata<-annot(total_metadata, meta_dir, 'GSE100386', c(rep('normal',46)), c(rep('Peripheral Blood Mononuclear Cell',46)),c(rep('False',46)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE100503', c(rep('disease',13)), c(rep('Breast Ductal Carcinoma In Situ',13)), c(rep('False',13)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE100561', c(rep('normal',12)), c(rep('Peripheral Blood Mononuclear Cell',12)),c(rep('False',12)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE100653', c(rep('normal',18)), c(rep('Breast',18)), c(rep('False',18)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE100940', c(rep('normal',24)), c(rep('Blood',24)), c(rep('False',24)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE101443', c(rep(list('disease','normal'),4)), c(rep(list('Breast Neoplasm','Breast'),4)), 
                      c(rep('False',8)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE101641', c(rep('normal',16),rep('disease',32)), 
                      c(rep('Nasal Cavity Epithelium',16),rep('Cystic Fibrosis',32)),
                      c(rep('False',48)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE101658', c(rep('disease',15)), c(rep('Multiple Sclerosis',15)),c(rep('False',15)))
to_remove<-append_to_remove(to_remove, 'GSE101673')   
to_remove<-append_to_remove(to_remove, 'GSE101733')
total_metadata<-annot(total_metadata, meta_dir, 'GSE101840', c(rep('normal',5)), c(rep('Whole Blood',5)),c(rep('False',5)))                      
total_metadata<-annot(total_metadata, meta_dir, 'GSE101961', c(rep('normal',121)), c(rep('Breast',121)), c(rep('False',121)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE102119', c(rep('disease',146)), c(rep('Ovarian Carcinoma',146)), c(rep('False',146)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE102177', c(rep('normal',36)), c(rep('Peripheral Blood',36)),c(rep('False',36)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE102504', c(rep('disease',25)), c(rep('Chronic Fatigue Syndrome',25)), c(rep('False',25)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE102970', c(rep('normal',48)), c(rep('Spermatozoon',48)), c(rep('False',48)))
to_remove<-append_to_remove(to_remove, 'GSE102994')
total_metadata<-annot(total_metadata, meta_dir, 'GSE103006', c(rep('normal',24)),c(rep('Umbilical Cord Blood',24)), c(rep('False',24)))
to_remove<-append_to_remove(to_remove, 'GSE103010')
total_metadata<-annot(total_metadata, meta_dir, 'GSE103413', c(rep('normal', 67)), 
                      c(rep('Brain',9),rep('Liver',7),rep('Breast',7),rep('Leukocyte',36),rep('Melanocyte',3),rep('Villus',2)),
                      c(rep('False',67)))
to_remove<-append_to_remove(to_remove, 'GSE103502')
total_metadata<-annot(total_metadata, meta_dir, 'GSE103659', c(rep('disease',181)),c(rep('Brain Neoplasm',181)),c(rep('False',181)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE103768', c(rep('normal',57)), c(rep('Adipose Tissue',57)), c(rep('False',57)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE103911', c(rep('normal',65)), c(rep('T-Lymphocyte',65)), c(rep('False',65)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE104087', c(rep('disease',40)),c(rep('Asthma',40)),c(rep('False',40)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE104287', c(rep('normal', 48)), 
                      c(rep('Peripheral Blood Mononuclear Cell', 32), rep('Natural Killer Cell',16)),
                      c(rep('False',48)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE104376', c(rep('normal',45)), c(rep('Umbilical Cord Blood',45)), c(rep('False',45)))
to_remove<-append_to_remove(to_remove, 'GSE104471')
total_metadata<-annot(total_metadata, meta_dir, 'GSE104728', c(rep('disease',44)), c(rep('Medulloblastoma',44)), c(rep('False',44)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE104770', c(rep('disease',14)), c(rep('Plasma Cell Leukemia',14)), c(rep('False',14)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE104812', c(rep('normal',48)), c(rep('Whole Blood',48)), c(rep('False',48)))
to_remove<-append_to_remove(to_remove, 'GSE105066')
total_metadata<-annot(total_metadata, meta_dir, 'GSE105123', c(rep('normal',108)), c(rep('Peripheral Blood Mononuclear Cell',108)), c(rep('False',108)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE105260', c(rep('disease',26),rep('normal',9),rep('disease',9)),
                      c(rep('Renal Cell Carcinoma',26),rep('Renal',9),rep('Renal Cell Carcinoma',9)),
                      c(rep('False',44)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE106089', c(rep('normal',84)),c(rep('Placenta',84)), c(rep('False',84)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE106099', c(rep('normal',30)),c(rep('Endothelial Cell',30)),c(rep('False',30)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE106360', c(rep('disease',47)),c(rep('Breast Neoplasm',47)),c(rep('False',47)))
total_metadata<-annot(total_metadata,meta_dir, 'GSE106556', c(rep('normal',6),rep('disease',6),rep('normal',5),rep('disease',6)),
                     c(rep('Hematopoietic Tissue',6),rep('Myeloid Leukemia',6),rep('Hematopoietic Tissue',5),rep('Myeloid Leukemia',6)),
                     c(rep('False',23)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE107038', c(rep('normal',40)),c(rep('Liver',40)),c(rep('False',40)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE107211', c(rep('normal',15)),c(rep('Whole Blood',15)),c(rep('False',15)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE107226', c(rep('disease',8), rep('normal',4)),c(rep('Pulmonary Fibrosis',8),rep('Fibroblast',4)),
                      c(rep('False',12)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE107351', 
                      c(rep('noraml',41),rep('disease',75)),c(rep('Blood',41),rep('Lynch Syndrome',61),rep('MLH1 Gene Mutation',14)),
                      c(rep('False',116)))
to_remove<-append_to_remove(to_remove, 'GSE107352')
total_metadata<-annot(total_metadata, meta_dir, 'GSE107737', c(rep('normal',12),rep('disease',12)), c(rep('Whole Blood',12),rep('Hypopituitarism',12)),
                      c(rep('False',24)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE108058', c(rep('normal',30)),c(rep('Spermatozoon',30)),c(rep('False',30)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE108423', c(rep('normal',21)),c(rep('Peripheral Blood',20)),c(rep('False',20)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE108562', c(rep('normal',6)),c(rep('Umbilical Cord Blood',6)),c(rep('False',6)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE108567', c(rep('normal',59)),c(rep('Chorionic Villus',59)),c(rep('False',59)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE108567', c(rep('disease',89)),
                      c(rep('Breast Carcinoma',30),rep('Brain Neoplasm',4),rep('Lung Carcinoma',18),rep('Melanoma',37)), c(rep('False',89)))
to_remove<-append_to_remove(to_remove, 'GSE109042')
total_metadata<-annot(total_metadata, meta_dir, 'GSE109096', c(rep('normal',11)),c(rep('Cardiac Ventricle',11)), c(rep('False',11)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE109430', c(rep('normal',12),rep('disease',24)), c(rep('Leukocyte',12),rep('Kawasaki Disease',24)), 
                      c(rep('False',36)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE109446', c(rep('normal',58)),c(rep('Olfactory Epithelial Cell', 58)), c(rep('False',58)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE109905', c(rep('disease',38), rep('normal',31)),
                      c(rep('Atrial Septal Defect', 38), rep('Whole Blood',31)),
                      c(rep('False',69)))
to_remove<-append_to_remove(to_remove, 'GSE109914')
to_remove<-append_to_remove(to_remove,'GSE110607')
total_metadata<-annot(total_metadata, meta_dir, 'GSE110607', c(rep('normal',21)),c(rep('Adipose Tissue',15),rep('Blood',6)), c(rep('False',21)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE110776', c(rep('normal',24)),c(rep('Bronchoalveolar Lavage Fluid', 24)), c(rep('False',21)))
to_remove<-append_to_remove(to_remove,'GSE110778')
to_remove<-append_to_remove(to_remove, 'GSE111396')
total_metadata<-annot(total_metadata, meta_dir, 'GSE111428', c(rep('disease',6)), c(rep('Ependymoma',6)), c(rep('False',6))) 
total_metadata<-annot(total_metadata, meta_dir, 'GSE111632', c(rep('normal',12)), c(rep('Adipose Tissue',12)),c(rep('False',12)))
to_remove<-apend_to_remove(to_remove, 'GSE111933')
total_metadata<-annot(total_metadata, meta_dir, 'GSE111933', c(rep('disease', 25), rep('normal',18)),
                      c(rep('Rheumatoid Arthritis',25),rep('Peripheral Blood Mononuclear Cell',18)), c(rep('False',43)))
to_remove<-append_to_remove(to_remove,'GSE112012')
total_metadata<-annot(total_metadata, meta_dir, 'GSE112047', c(rep('disease',47)),c(rep('Prostate Neoplasm',47)),c(rep('False',47)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE112306', c(rep('normal',6)),c(rep('Epithelial Cell',3),rep('Mesenchymal',3)),c(rep('False',6)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE112314', c(rep('normal',22)),c(rep('Saliva',22)),c(rep('False',22)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE112696', c(rep('normal',6),rep('disease',6)), c(rep('Peripheral Blood',6),rep('Myasthenia Gravis',6)),
                      c(rep('False',12)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE112987', c(rep('normal',64),rep('disease',39)), 
                      c(rep('Blood',64),rep('Fetal Alcohol Spectrum Disorder',39)),
                      c(rep('False',103)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE113012', c(rep('disease',7),rep('noraml',28)),
                      c(rep('Fetal Alcohol Spectrum Disorder',7),rep('Blood',28)),
                      c(rep('False',35)))
to_remove<-append_to_remove(to_remove, 'GSE113061')
to_remove<-append_to_remove(to_remove, 'GSE113775')
to_remove<-append_to_remove(to_remove, 'GSE113061')
to_remove<-append_to_remove(to_remove, 'GSE113775')
total_metadata<-annot(total_metadata, meta_dir, 'GSE113967', c(rep('normal',134)), c(rep('Whole Blood',134)),c(rep('False',134)))
to_remove<-append_to_remove(to_remove,'GSE114676')
to_remove<-append_to_remove(to_remove, 'GSE114683')
total_metadata<-annot(total_metadata, meta_dir, 'GSE114935', c(rep('normal',47)), c(rep('Whole Blood',47)),c(rep('False',47)))
to_remove<-append_to_remove(to_remove, 'GSE115399')
to_remove<-append_to_remove(to_remove, 'GSE115783')
total_metadata<-annot(total_metadata, meta_dir, 'GSE115797', c(rep(list('normal','disease'),24)), c(rep(list('Skin','Psoriasis'),24)),
                      c(rep('False',48)))
to_remove<-append_to_remove(to_remove,'GSE115852')
total_metadata<-annot(total_metadata, meta_dir, 'GSE115920', c(rep('normal',6)), c(rep('Spermatozoon',6)), c(rep('False',6)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE116057', c(rep('disease',12)),c(rep('Acute Lymphoblastic Leukemia', 12)), c(rep('False',12)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE116300', c(rep('normal',44)), c(rep('Blood',44)),c(rep('False',44)))
to_remove<-append_to_remove(to_remove, 'GSE116754')
to_remove<-append_to_remove(to_remove, 'GSE116924')
total_metadata<-annot(total_metadata, meta_dir, 'GSE117050', c(rep('normal',38)),c(rep('T-Lymphocyte',38)),c(rep('False',38)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE117439', c(rep('disease',52)),c(rep('Breast Carcinoma',52)), c(rep('False',52)))
to_remove<-append_to_remove(to_remove, 'GSE117448')
total_metadata<-annot(total_metadata, meta_dir, 'GSE117852', c(rep('disease',32)),c(rep('Pancreatic Neoplasm',32)),c(rep('False',32)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE118132', c(rep('disease',18)),c(rep('Lung Carcinoma',18)),c(rep('False',18)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE118250', c(rep('normal',12)),c(rep('Saliva',12)),c(rep('False',12)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE118260', c(rep('normal',20)),c(rep('Intestinal Mucosa',10),rep('Saliva',10)),c(rep('False',20)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE118469', c(rep('normal',6),rep('disease',15)), 
                      c(rep('Peripheral Blood Mononuclear Cell',6), rep('Tuberculosis', 15)),
                      c(rep('False',21)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE118570', c(rep('normal',43)),c(rep('T-Lymphocytes',43)),c(rep('False',43)))
to_remove<-append_to_remove(to_remove, 'GSE118696')
total_metadata<-annot(total_metadata, meta_dir, 'GSE119078', c(rep('normal', 59)), c(rep('Saliva',59)), c(rep('False',59)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE119778', c(rep('disease',34),rep('normal',34)), c(rep('Williams Syndrome',34),rep('Whole Blood',34)),
                      c(rep('False',68)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE119846', c(rep('normal', 60)), c(rep('Hematopoietic Stem Cell',60)), c(rep('False',60)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE120062', c(rep('normal', 38)), c(rep('Placenta',38)), c(rep('False',38)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE120307', c(rep('normal',34)), c(rep('Peripheral Blood',34)), c(rep('False',34)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE124367', c(rep('normal',12)), c(rep('Breast',12)), c(rep('False',12)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE124565', c(rep('disease',10),rep('normal',12)), 
                      c(rep('Antiphospholipid Syndrome',10),rep('Neutrophil',12)),
                      c(rep('False',22)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE124567', c(rep('normal', 54)), c(rep('Retina',51),rep('Embryonic Stem Cell',3)), c(rep('False',54)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE125465', c(rep('normal',12)), c(rep('Stem Cell',12)),c(rep('False',12)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE125605', c(rep('normal',42)), c(rep('Placenta',42)), c(rep('False',42)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE126017', c(rep('normal',18), rep('disease',36)),
                      c(rep('Spermatozoon',18), rep('Psoriasis',36)), c(rep('False',54)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE128126', c(rep('normal',15)), c(rep('Embryonic Stem Cell',15)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE127824', c(rep('normal',24)), c(rep('Umbilical Cord Blood',24)), c(rep('False',24)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE128801', c(rep('normal',12)), c(rep('Blood',12)), c(rep('False',12)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE129266', c(rep('normal',12)), c(rep('Mesenchymal Stem Cell',12)), c(rep('False',12)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE130354', c(rep('disease',12)), c(rep('Melanoma',12)), c(rep('False',12)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE131350', c(rep('disease',48), rep('normal',12)), c(rep('Adenocarcinoma',48), rep('Adrenal Gland',12)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE132181', c(rep('normal', 392)), rep(c('Umbilical Cord Blood','Peripheral Blood Mononuclear Cell'),196))
total_metadata<-annot(total_metadata, meta_dir, 'GSE132399', c(rep('normal',26), rep('disease',28)), c(rep('Liver', 26), rep('Hepatoblastoma',28)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE137716', c(rep('normal',16),rep('disease',17)), c(rep('Airway',16),rep('Asthma',17)),
                      c(rep('False',33)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE138279', c(rep('normal',65)), c(rep('Saliva',65)), c(rep('False',65)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE139307', c(rep('normal',37)), c(rep('Spermatozoon',37)) , c(rep('False',37)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE141338', c(rep('normal',6),rep('disease',42)), c(rep('Breast',6),rep('Breast Neoplasm',42)),
                      c(rep('False',48)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE142554', c(rep('normal',12)), c(rep('Blood',12)), c(rep('False',12)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE143209', c(rep('normal',64)), c(rep('Islet of Langerhans',64)), c(rep('False',64)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE144894', c(rep('disease', 120)), c(rep('Chronic Lymphocytic Leukemia',120)),
                      c(rep('False',120)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE144977', c(rep('normal',89)),c(rep('Placenta',89)),c(rep('False',89)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE147740', c(rep('normal',1129)), c(rep('Blood',1129)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE150901', c(rep('normal',48)), c(rep('Buccal Mucosa', 48)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE151278', c(rep('disease',70)), c(rep('Psoriasis',70)),c(rep('False',70)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE151042', c(rep('normal',492)), c(rep('Umbilical Cord Blood',492)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE151407', c(rep('normal',78)), c(rep('Vastus Lateralis',78)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE151485', c(rep('normal', 100)), c(rep('Saliva',100)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE151600', c(rep('normal', 16)), c(rep('Skin',16)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE152380', c(rep('normal',90)), c(rep('Umbilical Cord Blood',90)),c(rep('False',90)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE47915', c(rep('disease',8)), c(rep('Prostate Neoplasm',8)), c(rep('False',8)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE48472', c(rep('normal',56)), 
                      c(rep('Blood',6),rep('Liver',5),rep('Muscle',6),rep('Omentum',6),rep('Pancreas',4),
                        rep('Adipose',6),rep('Spleen',3),rep('Blood',5),rep('Buccal Mucosa',5),
                        rep('Hair',5),rep('Saliva',5)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE49618', c(rep('normal',21)), c(rep('Bone Marrow',21)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE52025', c(rep('normal',62)), c(rep('Fibroblast',62)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE59038', c(rep('disease',24)), c(rep('Colon Adenocarcinoma',24)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE59524', c(rep('normal',24)), c(rep('Adipose Tissue',24)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE60655', c(rep('normal',36)), c(rep('Vastus Lateralis', 36)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE61107', c(rep('disease',23), rep('normal',24), 'disease'), c(rep('Schizophrenia', 23), rep('Frontal Lobe Cortex', 24), 'Schizophrenia'))
total_metadata<-annot(total_metadata, meta_dir, 'GSE61278', c(rep('normal', 66), rep('disease',44)), c(rep('Liver', 66), rep('Adult Liver Carcinoma', 44)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE61446', c(rep('normal',67)), c(rep('Liver',67))) #obesity
total_metadata<-annot(total_metadata, meta_dir, 'GSE61450', c(rep('normal',71)), c(rep('Adipose Tissue',71))) #obesity
total_metadata<-annot(total_metadata, meta_dir, 'GSE61452', c(rep('normal',60)), c(rep('Muscle',60))) #obesity
total_metadata<-annot(total_metadata, meta_dir, 'GSE61453', c(rep('normal',71)), c(rep('Adipose Tissue',71))) #obesity
total_metadata<-annot(total_metadata, meta_dir, 'GSE62727', c(rep('disease', 7), rep('normal',4)), c(rep('Atrial Fibrillation', 7),rep('Left Atrium',4)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE63179', c(rep('normal',8)), c(rep('Cerebellum',8)))                      
total_metadata<-annot(total_metadata, meta_dir, 'GSE64096', c(rep('normal',40)), c(rep('Spermatozoon',40)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE66351', c(rep('disease', 190)), c(rep('Alzheimer\'s Disease', 190)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE72245', c(rep('disease',118)), c(rep('Breast Carcinoma',118)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE72251', c(rep('disease',119)), c(rep('Breast Carcinoma',119)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE72254', c(rep('disease',58)), c(rep('Breast Carcinoma',58)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE74167', c(rep('normal',42)), c(rep('Mesenchymal Stem Cell',42)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE74214', c(rep('normal',18)), c(rep('Breast',18)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE75133', c(rep('normal',15)), c(rep('Mesenchymal Stem Cell',15)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE75405', c(rep('normal',24)), c(rep('Peripheral Blood',24)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE75443', c(rep('disease',12)), c(rep('Glioblastoma',12)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE75704', c(rep('normal',72), rep('disease',94)), c(rep('Frontal Lobe Cortex',72), rep('Progressive Supranuclear Palsy')))
total_metadata<-annot(total_metadata, meta_dir, 'GSE76372', c(rep('noraml',9)), c(rep('Lymphocyte',2),'Induced Pluripotent Stem Cell', rep('Lymphocyte',3), rep('Induced Pluripotent Stem Cell',2), 'Lymphocyte'))
total_metadata<-annot(total_metadata, meta_dir, 'GSE76503', c(rep('normal',48)), c(rep('Whole Blood',48)), c(rep('False',36), rep('True',12)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE77269', c(rep('disease',60)), c(rep('Hepatocellular Carcinoma',60)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE78732', c(rep('normal',10),rep('disease',19)),c(rep('Liver',10),rep('Hepatoblastoma',19)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE79009', c(rep('disease',125)),c(rep('Schwannoma',125)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE79329', c(rep('normal',34)), c(rep('Leukocyte',34)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE80685', c(rep('disease',172)), c(rep('Prostate Neoplasm',172)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE82084', c(rep('normal',36)), c(rep('Umbilical Cord Blood',36)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE82234', c(rep('normal',6)), c(rep('Umbilical Cord Blood',6)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE83917', c(rep('disease',160)), c(rep('Prostate Neoplasm',160)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE83933', c(rep('disease',39)), c(rep('Meningioma',39)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE83944', c(rep('normal',48)), c(rep('Neutrophil',48)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE84274', c(rep('normal',18), rep('disease',18)), c(rep('Ascending Aorta', 18))) 
total_metadata<-annot(total_metadata, meta_dir, 'GSE85042', c(rep('normal',71)), c(rep('Umbilical Cord Blood',71)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE85210', c(rep('normal',253)), c(rep('Whole Blood',253)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE85506', c(rep('normal',47)), c(rep('Peripheral Blood',47)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE85845', c(rep('disease',8), rep('normal',8)),c(rep('Lung Adenocarcinoma',8), rep('Lung',8)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE86078', c(rep('disease',149)), c(rep('Colorectal Carcinoma', 149)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE86258', c(rep('disease',7),rep('normal',7)), c(rep('Fibroblast',14)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE87582', c(rep('normal',21)), c(rep('Peripheral Blood',21)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE87640', c(rep('normal', 240)), c(rep('Leukocyte',240)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE88883', c(rep('normal',100)), c(rep('Breast',100)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE89218', c(rep('normal',163)), c(rep('Whole Blood',163)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE89251', c(rep('normal',136)), c(rep('CD4-Positive T-Lymphocyte',136)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE89472', c(rep(list('disease','normal'),5)), c(rep('Blood',10)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE89803', c(rep('disease',138), rep('normal',4)), c(rep('Bile Duct', 142)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE89852', c(rep('normal',37), rep('disease',27)), c(rep('Liver', 74)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE89925', c(rep('normal',3),rep('disease',3)), c(rep('Fibroblast',3),rep('Teratoma',3)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE90060', c(rep('normal', 34)), c(rep('Endometrium',34)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE92577', c(rep('disease',12)), c(rep('Brain Neoplasm',12)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE92909', c(rep('disease',6)), c(rep('Breast Carcinoma', 6)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE93208', c(rep('normal', 19)), c(rep('Placenta', 19)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE94326', c(rep('normal', 8)), c(rep('Brain',4),rep('Prostate',4)), c(rep('False',4),rep('True',4)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE95486', c(rep('normal',21), rep('disease',3)), c(rep('Blood',24)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE97784', c(rep(list('disease','normal'),12)), c(rep('Oral Mucosa',24)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE98056', c(rep('normal',69)),c(rep('Leukocyte',69)),c(rep('False',35),rep('True',34)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE98203', c(rep('normal',88)), c(rep('Neuron',88)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE99553', c(rep('normal',84)), c(rep('Gastric Mucosa',84)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE99624', c(rep('disease',32),rep('normal',16)), c(rep('Blood',48)))
total_metadata<-annot(total_metadata, meta_dir, 'GSE99755', c(rep('normal',67)), c(rep('Whole Blood',67)))

df <- apply(total_metadata,2,as.character)
write.csv(df, file=paste0('data/', database_type, '/all_metadata_annotated.txt'),
          row.names=F,
          sep="\t",
          quote=F)
          
gzip(filename=paste0('data/', database_type, '/all_metadata_annotated.txt'), 
     destname=paste0('data/', database_type, '/all_metadata_annotated.txt.gz'), overwrite=TRUE, remove=TRUE)
