import torch
import scanpy as sc
import scvi
import re
import pandas as pd
import numpy as np
import anndata as ad
import matplotlib.pyplot as plt
from collections import Counter
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import ShuffleSplit, GridSearchCV, train_test_split
from sklearn.metrics import confusion_matrix, classification_report, accuracy_score, roc_auc_score, accuracy_score, recall_score, f1_score, precision_score
from sklearn.preprocessing import LabelEncoder
from scipy.sparse import issparse
from skorch import NeuralNetClassifier
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from torch.nn.parallel import DistributedDataParallel as DDP
from torch.utils.data import DataLoader, DistributedSampler

def table(vector):
    count_dict = Counter(vector)   
    elements = list(count_dict.keys())
    counts = list(count_dict.values())  
    table_data = {'element': elements, 'count': counts}
    df = pd.DataFrame(table_data)    
    return df

merge_sn = sc.read_h5ad("/gpfs/hpc/home/chenchao/ranm/SCvsSN/SN_Datasets/counts_merge_sn.h5ad")
meta = pd.read_csv("/gpfs/hpc/home/chenchao/ranm/SCvsSN/SN_Datasets/merge_sn/counts_merge_sn_obs.csv", index_col=0)
merge_sn = merge_sn[merge_sn.obs.index.isin(meta.index)]
RG_del = pd.read_csv('/gpfs/hpc/home/chenchao/ranm/SCvsSN/figures/monocle/RG_del.csv')
RG_del = RG_del['x'].tolist()
merge_sn = merge_sn[~merge_sn.obs.index.isin(RG_del)]
Diff_cell = ['RG', 'Ast', 'Glioblast']
diff_gene = pd.read_excel("/gpfs/hpc/home/chenchao/ranm/SCvsSN/SN_Datasets/RandomForest/Gene_importance_3Celltype_233gene_mean.xlsx", 
                         header=0, 
                         index_col=0)
NPC_marker = diff_gene.index
merge_sn = merge_sn[:, merge_sn.var_names.isin(NPC_marker)]

results = train_random_forest_model(merge_sn, 'Zhu', Diff_cell, NPC_marker)

####count to log2(CPM)
def convert_to_cpm(adata):
    counts_matrix = adata.X
    if issparse(counts_matrix):
        counts_matrix = counts_matrix.toarray()
    total_counts_per_cell = counts_matrix.sum(axis=1)
    total_counts_per_cell = np.maximum(total_counts_per_cell, 1)
    cpm_matrix = (counts_matrix / total_counts_per_cell[:, np.newaxis]) * 1e6
    log2_cpm_matrix = np.log2(cpm_matrix+1)
    return adata.__class__(log2_cpm_matrix, obs=adata.obs, var=adata.var, dtype='float32')
	

Zhu_adata = merge_sn[merge_sn.obs['Dataset'] == 'Zhu']
Zhu_adata = Zhu_adata[Zhu_adata.obs['Celltype'].isin(Diff_cell)]
Zhu_adata = convert_to_cpm(Zhu_adata)
Zhu_df = Zhu_adata.to_df()
Zhu_celltype = np.array(Zhu_adata.obs['Celltype'])
NPC_marker_index = np.where(Zhu_adata.var_names.isin(NPC_marker))[0]
Zhu_NPC_marker = Zhu_df.iloc[:, NPC_marker_index]
Zhu_celltype = np.array(Zhu_adata.obs['Celltype'])
seed_value = 12345   
test_size = 0.3         
Zhu_X_train, Zhu_X_test, Zhu_y_train, Zhu_y_test = train_test_split(Zhu_NPC_marker, Zhu_celltype, test_size=test_size, random_state=seed_value)
label_encoder = LabelEncoder()
Zhu_y_train = label_encoder.fit_transform(Zhu_y_train)
Zhu_y_test = label_encoder.fit_transform(Zhu_y_test)

param_test1 = {"n_estimators":range(1000,10001,500)}
gsearch1 = GridSearchCV(estimator=RandomForestClassifier(n_jobs=-1),param_grid=param_test1,
                        scoring='accuracy',cv=10,error_score='raise')
gsearch1.fit(Zhu_X_train, Zhu_y_train)
print(gsearch1.best_score_)
print(gsearch1.best_estimator_)
print("best accuracy:%f" % gsearch1.best_score_)
print(gsearch1.cv_results_)
Zhu_numbers = re.findall(r'\d+|\-\d+', str(gsearch1.best_estimator_))
Zhu_best_estimator = [int(num) if num.lstrip('-').isdigit() else float(num) for num in Zhu_numbers][0]

maxdepth = {'max_depth':range(5,21,1)}
gsearch2 = GridSearchCV(estimator = RandomForestClassifier(n_estimators = Zhu_best_estimator,n_jobs=-1),param_grid = maxdepth,
                        scoring = 'accuracy',cv = 10,error_score='raise')
gsearch2.fit(Zhu_X_train, Zhu_y_train)
print(gsearch2.best_score_)
print(gsearch2.best_estimator_)
Zhu_numbers = re.findall(r'\d+|\-\d+', str(gsearch2.best_estimator_))
Zhu_best_depth = [int(num) if num.lstrip('-').isdigit() else float(num) for num in Zhu_numbers][0]

minsamples = {'min_samples_split':range(2,101,2)}
gsearch3 = GridSearchCV(estimator = RandomForestClassifier(n_estimators=Zhu_best_estimator, max_depth=Zhu_best_depth,n_jobs=-1),param_grid = minsamples,
                        scoring = 'accuracy',cv = 10, error_score='raise')
gsearch3.fit(Zhu_X_train,Zhu_y_train)
print(gsearch3.best_score_)
print(gsearch3.best_estimator_)
Zhu_numbers = re.findall(r'\d+|\-\d+', str(gsearch3.best_estimator_))
Zhu_best_sample = [int(num) if num.lstrip('-').isdigit() else float(num) for num in Zhu_numbers][1]
Zhu_model = RandomForestClassifier(n_estimators=Zhu_best_estimator, max_depth=Zhu_best_depth, min_samples_split=Zhu_best_sample, random_state=seed_value)
Zhu_model.fit(Zhu_X_train, Zhu_y_train)
Zhu_y_train_pred = Zhu_model.predict(Zhu_X_train)
Zhu_y_test_pred = Zhu_model.predict(Zhu_X_test)
Zhu_score_train = Zhu_model.score(Zhu_X_test, Zhu_y_test)
Zhu_score_test = accuracy_score(Zhu_y_test_pred, Zhu_y_test)
y_pred_proba = Zhu_model.predict_proba(Zhu_X_test)[:, 1] 
Zhu_importances = Zhu_model.feature_importances_
Gene_labels = Zhu_NPC_marker.columns[0:]
indices = np.argsort(Zhu_importances)[::-1] 
Zhu_gene_importance  = [(f + 1, Gene_labels[indices[f]], Zhu_importances[indices[f]])
                          for f in range(len(Gene_labels))]
Zhu_gene_importance = pd.DataFrame(Zhu_gene_importance, columns=['Rank', 'Gene', 'Importance'])

cameron_adata = merge_sn[merge_sn.obs['Dataset'] == 'cameron']
cameron_adata = cameron_adata[cameron_adata.obs['Celltype'].isin(Diff_cell)]
cameron_adata = convert_to_cpm(cameron_adata)
cameron_df = cameron_adata.to_df()
cameron_celltype = np.array(cameron_adata.obs['Celltype'])
NPC_marker_index = np.where(cameron_adata.var_names.isin(NPC_marker))[0]
cameron_NPC_marker = cameron_df.iloc[:, NPC_marker_index]
cameron_celltype = np.array(cameron_adata.obs['Celltype'])
seed_value = 12345     
test_size = 0.3         
cameron_X_train, cameron_X_test, cameron_y_train, cameron_y_test = train_test_split(cameron_NPC_marker, cameron_celltype, test_size=test_size, random_state=seed_value)
label_encoder = LabelEncoder()
cameron_y_train = label_encoder.fit_transform(cameron_y_train)
cameron_y_test = label_encoder.fit_transform(cameron_y_test)

param_test1 = {"n_estimators":range(1000,10001,500)}
gsearch1 = GridSearchCV(estimator=RandomForestClassifier(n_jobs=-1),param_grid=param_test1,
                        scoring='accuracy',cv=10,error_score='raise')
gsearch1.fit(cameron_X_train, cameron_y_train)
print(gsearch1.best_score_)
print(gsearch1.best_estimator_)
print("best accuracy:%f" % gsearch1.best_score_)
print(gsearch1.cv_results_)
cameron_numbers = re.findall(r'\d+|\-\d+', str(gsearch1.best_estimator_))
cameron_best_estimator = [int(num) if num.lstrip('-').isdigit() else float(num) for num in cameron_numbers][0]

maxdepth = {'max_depth':range(5,21,1)}
gsearch2 = GridSearchCV(estimator = RandomForestClassifier(n_estimators = cameron_best_estimator,n_jobs=-1),param_grid = maxdepth,
                        scoring = 'accuracy',cv = 10,error_score='raise')
gsearch2.fit(cameron_X_train, cameron_y_train)
print(gsearch2.best_score_)
print(gsearch2.best_estimator_)
cameron_numbers = re.findall(r'\d+|\-\d+', str(gsearch2.best_estimator_))
cameron_best_depth = [int(num) if num.lstrip('-').isdigit() else float(num) for num in cameron_numbers][0]

minsamples = {'min_samples_split':range(2,101,2)}
gsearch3 = GridSearchCV(estimator = RandomForestClassifier(n_estimators=cameron_best_estimator, max_depth=cameron_best_depth,n_jobs=-1),param_grid = minsamples,
                        scoring = 'accuracy',cv = 10, error_score='raise')
gsearch3.fit(cameron_X_train,cameron_y_train)
print(gsearch3.best_score_)
print(gsearch3.best_estimator_)
cameron_numbers = re.findall(r'\d+|\-\d+', str(gsearch3.cv_results_['params'][gsearch3.best_index_]))
cameron_best_sample = int([int(num) if num.lstrip('-').isdigit() else float(num) for num in cameron_numbers][0])


cameron_model = RandomForestClassifier(n_estimators=cameron_best_estimator, max_depth=cameron_best_depth, min_samples_split=cameron_best_sample, random_state=seed_value)
cameron_model.fit(cameron_X_train, cameron_y_train)
cameron_y_train_pred = cameron_model.predict(cameron_X_train)
cameron_y_test_pred = cameron_model.predict(cameron_X_test)
cameron_score_train = cameron_model.score(cameron_X_test, cameron_y_test)
cameron_score_test = accuracy_score(cameron_y_test_pred, cameron_y_test)
y_pred_proba = cameron_model.predict_proba(cameron_X_test)[:, 1] 
cameron_importances = cameron_model.feature_importances_
Gene_labels = cameron_NPC_marker.columns[0:]
indices = np.argsort(cameron_importances)[::-1]
cameron_gene_importance  = [(f + 1, Gene_labels[indices[f]], cameron_importances[indices[f]])
                          for f in range(len(Gene_labels))]
cameron_gene_importance = pd.DataFrame(cameron_gene_importance, columns=['Rank', 'Gene', 'Importance'])

Velmeshev_adata = merge_sn[merge_sn.obs['Dataset'] == 'Velmeshev']
Velmeshev_adata = Velmeshev_adata[Velmeshev_adata.obs['Celltype'].isin(Diff_cell)]
Velmeshev_adata = convert_to_cpm(Velmeshev_adata)
Velmeshev_df = Velmeshev_adata.to_df()
Velmeshev_celltype = np.array(Velmeshev_adata.obs['Celltype'])
NPC_marker_index = np.where(Velmeshev_adata.var_names.isin(NPC_marker))[0]
Velmeshev_NPC_marker = Velmeshev_df.iloc[:, NPC_marker_index]
Velmeshev_celltype = np.array(Velmeshev_adata.obs['Celltype'])
seed_value = 12345    
test_size = 0.3       
Velmeshev_X_train, Velmeshev_X_test, Velmeshev_y_train, Velmeshev_y_test = train_test_split(Velmeshev_NPC_marker, Velmeshev_celltype, test_size=test_size, random_state=seed_value)
label_encoder = LabelEncoder()
Velmeshev_y_train = label_encoder.fit_transform(Velmeshev_y_train)
Velmeshev_y_test = label_encoder.fit_transform(Velmeshev_y_test)

param_test1 = {"n_estimators":range(1000,10001,500)}
gsearch1 = GridSearchCV(estimator=RandomForestClassifier(n_jobs=-1),param_grid=param_test1,
                        scoring='accuracy',cv=10,error_score='raise')
gsearch1.fit(Velmeshev_X_train, Velmeshev_y_train)
print(gsearch1.best_score_)
print(gsearch1.best_estimator_)
print("best accuracy:%f" % gsearch1.best_score_)
print(gsearch1.cv_results_)
Velmeshev_numbers = re.findall(r'\d+|\-\d+', str(gsearch1.best_estimator_))
Velmeshev_best_estimator = [int(num) if num.lstrip('-').isdigit() else float(num) for num in Velmeshev_numbers][0]

maxdepth = {'max_depth':range(5,21,1)}
gsearch2 = GridSearchCV(estimator = RandomForestClassifier(n_estimators = Velmeshev_best_estimator,n_jobs=-1),param_grid = maxdepth,
                        scoring = 'accuracy',cv = 10,error_score='raise')
gsearch2.fit(Velmeshev_X_train, Velmeshev_y_train)
print(gsearch2.best_score_)
print(gsearch2.best_estimator_)
Velmeshev_numbers = re.findall(r'\d+|\-\d+', str(gsearch2.best_estimator_))
Velmeshev_best_depth = [int(num) if num.lstrip('-').isdigit() else float(num) for num in Velmeshev_numbers][0]

minsamples = {'min_samples_split':range(2,101,2)}
gsearch3 = GridSearchCV(estimator = RandomForestClassifier(n_estimators=Velmeshev_best_estimator, max_depth=Velmeshev_best_depth,n_jobs=-1),param_grid = minsamples,
                        scoring = 'accuracy',cv = 10, error_score='raise')
gsearch3.fit(Velmeshev_X_train,Velmeshev_y_train)
print(gsearch3.best_score_)
print(gsearch3.best_estimator_)
Velmeshev_numbers = re.findall(r'\d+|\-\d+', str(gsearch3.cv_results_['params'][gsearch3.best_index_]))
Velmeshev_best_sample = int([int(num) if num.lstrip('-').isdigit() else float(num) for num in Velmeshev_numbers][0])

Velmeshev_model = RandomForestClassifier(n_estimators=Velmeshev_best_estimator, max_depth=Velmeshev_best_depth, min_samples_split=Velmeshev_best_sample, random_state=seed_value)
Velmeshev_model.fit(Velmeshev_X_train, Velmeshev_y_train)
Velmeshev_y_train_pred = Velmeshev_model.predict(Velmeshev_X_train)
Velmeshev_y_test_pred = Velmeshev_model.predict(Velmeshev_X_test)
Velmeshev_score_train = Velmeshev_model.score(Velmeshev_X_test, Velmeshev_y_test)
Velmeshev_score_test = accuracy_score(Velmeshev_y_test_pred, Velmeshev_y_test)
y_pred_proba = Velmeshev_model.predict_proba(Velmeshev_X_test)[:, 1]  
#Velmeshev_Auc = roc_auc_score(Velmeshev_y_test, y_pred_proba)
#print(f"Velmeshev ROC AUC: {Velmeshev_Auc:.3f}")
Velmeshev_importances = Velmeshev_model.feature_importances_
Gene_labels = Velmeshev_NPC_marker.columns[0:]
indices = np.argsort(Velmeshev_importances)[::-1] 
Velmeshev_gene_importance  = [(f + 1, Gene_labels[indices[f]], Velmeshev_importances[indices[f]])
                          for f in range(len(Gene_labels))]
Velmeshev_gene_importance = pd.DataFrame(Velmeshev_gene_importance, columns=['Rank', 'Gene', 'Importance'])

Herring_adata = merge_sn[merge_sn.obs['Dataset'] == 'Herring']
Herring_adata = Herring_adata[Herring_adata.obs['Celltype'].isin(Diff_cell)]
Herring_adata = convert_to_cpm(Herring_adata)
Herring_df = Herring_adata.to_df()
Herring_celltype = np.array(Herring_adata.obs['Celltype'])
NPC_marker_index = np.where(Herring_adata.var_names.isin(NPC_marker))[0]
Herring_NPC_marker = Herring_df.iloc[:, NPC_marker_index]
Herring_celltype = np.array(Herring_adata.obs['Celltype'])
seed_value = 12345     
test_size = 0.3        
Herring_X_train, Herring_X_test, Herring_y_train, Herring_y_test = train_test_split(Herring_NPC_marker, Herring_celltype, test_size=test_size, random_state=seed_value)
label_encoder = LabelEncoder()
Herring_y_train = label_encoder.fit_transform(Herring_y_train)
Herring_y_test = label_encoder.fit_transform(Herring_y_test)

param_test1 = {"n_estimators":range(1000,10001,500)}
gsearch1 = GridSearchCV(estimator=RandomForestClassifier(n_jobs=-1),param_grid=param_test1,
                        scoring='accuracy',cv=10,error_score='raise')
gsearch1.fit(Herring_X_train, Herring_y_train)
print(gsearch1.best_score_)
print(gsearch1.best_estimator_)
print("best accuracy:%f" % gsearch1.best_score_)
print(gsearch1.cv_results_)
Herring_numbers = re.findall(r'\d+|\-\d+', str(gsearch1.best_estimator_))
Herring_best_estimator = [int(num) if num.lstrip('-').isdigit() else float(num) for num in Herring_numbers][0]
maxdepth = {'max_depth':range(5,21,1)}
gsearch2 = GridSearchCV(estimator = RandomForestClassifier(n_estimators = Herring_best_estimator,n_jobs=-1),param_grid = maxdepth,
                        scoring = 'accuracy',cv = 10,error_score='raise')
gsearch2.fit(Herring_X_train, Herring_y_train)
print(gsearch2.best_score_)
print(gsearch2.best_estimator_)
Herring_numbers = re.findall(r'\d+|\-\d+', str(gsearch2.best_estimator_))
Herring_best_depth = [int(num) if num.lstrip('-').isdigit() else float(num) for num in Herring_numbers][0]
minsamples = {'min_samples_split':range(2,101,2)}
gsearch3 = GridSearchCV(estimator = RandomForestClassifier(n_estimators=Herring_best_estimator, max_depth=Herring_best_depth,n_jobs=-1),param_grid = minsamples,
                        scoring = 'accuracy',cv = 10, error_score='raise')
gsearch3.fit(Herring_X_train,Herring_y_train)
print(gsearch3.best_score_)
print(gsearch3.best_estimator_)
Herring_numbers = re.findall(r'\d+|\-\d+', str(gsearch3.cv_results_['params'][gsearch3.best_index_]))
Herring_best_sample = int([int(num) if num.lstrip('-').isdigit() else float(num) for num in Herring_numbers][0])

Herring_model = RandomForestClassifier(n_estimators=Herring_best_estimator, max_depth=Herring_best_depth, min_samples_split=Herring_best_sample, random_state=seed_value)
Herring_model.fit(Herring_X_train, Herring_y_train)
Herring_y_train_pred = Herring_model.predict(Herring_X_train)
Herring_y_test_pred = Herring_model.predict(Herring_X_test)
Herring_score_train = Herring_model.score(Herring_X_test, Herring_y_test)
Herring_score_test = accuracy_score(Herring_y_test_pred, Herring_y_test)
y_pred_proba = Herring_model.predict_proba(Herring_X_test)[:, 1]  
#Herring_Auc = roc_auc_score(Herring_y_test, y_pred_proba)
#print(f"Herring ROC AUC: {Herring_Auc:.3f}")
Herring_importances = Herring_model.feature_importances_
Gene_labels = Herring_NPC_marker.columns[0:]
indices = np.argsort(Herring_importances)[::-1] # 下标排序
Herring_gene_importance  = [(f + 1, Gene_labels[indices[f]], Herring_importances[indices[f]])
                          for f in range(len(Gene_labels))]
Herring_gene_importance = pd.DataFrame(Herring_gene_importance, columns=['Rank', 'Gene', 'Importance'])

Trevino_adata = merge_sn[merge_sn.obs['Dataset'] == 'Trevino']
Trevino_adata = Trevino_adata[Trevino_adata.obs['Celltype'].isin(Diff_cell)]
Trevino_adata = convert_to_cpm(Trevino_adata)
Trevino_df = Trevino_adata.to_df()
Trevino_celltype = np.array(Trevino_adata.obs['Celltype'])
NPC_marker_index = np.where(Trevino_adata.var_names.isin(NPC_marker))[0]
Trevino_NPC_marker = Trevino_df.iloc[:, NPC_marker_index]
Trevino_celltype = np.array(Trevino_adata.obs['Celltype'])
seed_value = 12345      
test_size = 0.3         
Trevino_X_train, Trevino_X_test, Trevino_y_train, Trevino_y_test = train_test_split(Trevino_NPC_marker, Trevino_celltype, test_size=test_size, random_state=seed_value)
label_encoder = LabelEncoder()
Trevino_y_train = label_encoder.fit_transform(Trevino_y_train)
Trevino_y_test = label_encoder.fit_transform(Trevino_y_test)


param_test1 = {"n_estimators":range(1000,10001,500)}
gsearch1 = GridSearchCV(estimator=RandomForestClassifier(n_jobs=-1),param_grid=param_test1,
                        scoring='accuracy',cv=10,error_score='raise')
gsearch1.fit(Trevino_X_train, Trevino_y_train)
print(gsearch1.best_score_)
print(gsearch1.best_estimator_)
print("best accuracy:%f" % gsearch1.best_score_)
print(gsearch1.cv_results_)
Trevino_numbers = re.findall(r'\d+|\-\d+', str(gsearch1.best_estimator_))
Trevino_best_estimator = [int(num) if num.lstrip('-').isdigit() else float(num) for num in Trevino_numbers][0]

maxdepth = {'max_depth':range(5,21,1)}
gsearch2 = GridSearchCV(estimator = RandomForestClassifier(n_estimators = Trevino_best_estimator,n_jobs=-1),param_grid = maxdepth,
                        scoring = 'accuracy',cv = 10,error_score='raise')
gsearch2.fit(Trevino_X_train, Trevino_y_train)
print(gsearch2.best_score_)
print(gsearch2.best_estimator_)
Trevino_numbers = re.findall(r'\d+|\-\d+', str(gsearch2.best_estimator_))
Trevino_best_depth = [int(num) if num.lstrip('-').isdigit() else float(num) for num in Trevino_numbers][0]

minsamples = {'min_samples_split':range(2,101,2)}
gsearch3 = GridSearchCV(estimator = RandomForestClassifier(n_estimators=Trevino_best_estimator, max_depth=Trevino_best_depth,n_jobs=-1),param_grid = minsamples,
                        scoring = 'accuracy',cv = 10, error_score='raise')
gsearch3.fit(Trevino_X_train,Trevino_y_train)
print(gsearch3.best_score_)
print(gsearch3.best_estimator_)
Trevino_numbers = re.findall(r'\d+|\-\d+', str(gsearch3.cv_results_['params'][gsearch3.best_index_]))
#Trevino_best_sample = 2
Trevino_best_sample = int([int(num) if num.lstrip('-').isdigit() else float(num) for num in Trevino_numbers][0])


Trevino_model = RandomForestClassifier(n_estimators=Trevino_best_estimator, max_depth=Trevino_best_depth, min_samples_split=Trevino_best_sample, random_state=seed_value)
Trevino_model.fit(Trevino_X_train, Trevino_y_train)
Trevino_y_train_pred = Trevino_model.predict(Trevino_X_train)
Trevino_y_test_pred = Trevino_model.predict(Trevino_X_test)
Trevino_score_train = Trevino_model.score(Trevino_X_test, Trevino_y_test)
Trevino_score_test = accuracy_score(Trevino_y_test_pred, Trevino_y_test)
y_pred_proba = Trevino_model.predict_proba(Trevino_X_test)[:, 1]  
#Trevino_Auc = roc_auc_score(Trevino_y_test, y_pred_proba)
#print(f"Trevino ROC AUC: {Trevino_Auc:.3f}")
Trevino_importances = Trevino_model.feature_importances_
Gene_labels = Trevino_NPC_marker.columns[0:]
indices = np.argsort(Trevino_importances)[::-1]
Trevino_gene_importance  = [(f + 1, Gene_labels[indices[f]], Trevino_importances[indices[f]])
                          for f in range(len(Gene_labels))]
Trevino_gene_importance = pd.DataFrame(Trevino_gene_importance, columns=['Rank', 'Gene', 'Importance'])


Ramos_adata = merge_sn[merge_sn.obs['Dataset'] == 'Ramos']
Ramos_adata = Ramos_adata[Ramos_adata.obs['Celltype'].isin(Diff_cell)]
Ramos_adata = convert_to_cpm(Ramos_adata)
Ramos_df = Ramos_adata.to_df()
Ramos_celltype = np.array(Ramos_adata.obs['Celltype'])
NPC_marker_index = np.where(Ramos_adata.var_names.isin(NPC_marker))[0]
Ramos_NPC_marker = Ramos_df.iloc[:, NPC_marker_index]
Ramos_celltype = np.array(Ramos_adata.obs['Celltype'])
seed_value = 12345     
test_size = 0.3        
Ramos_X_train, Ramos_X_test, Ramos_y_train, Ramos_y_test = train_test_split(Ramos_NPC_marker, Ramos_celltype, test_size=test_size, random_state=seed_value)
label_encoder = LabelEncoder()
Ramos_y_train = label_encoder.fit_transform(Ramos_y_train)
Ramos_y_test = label_encoder.fit_transform(Ramos_y_test)

param_test1 = {"n_estimators":range(1000,10001,500)}
gsearch1 = GridSearchCV(estimator=RandomForestClassifier(n_jobs=-1),param_grid=param_test1,
                        scoring='accuracy',cv=10,error_score='raise')
gsearch1.fit(Ramos_X_train, Ramos_y_train)
print(gsearch1.best_score_)
print(gsearch1.best_estimator_)
print("best accuracy:%f" % gsearch1.best_score_)
print(gsearch1.cv_results_)
Ramos_numbers = re.findall(r'\d+|\-\d+', str(gsearch1.best_estimator_))
Ramos_best_estimator = [int(num) if num.lstrip('-').isdigit() else float(num) for num in Ramos_numbers][0]

maxdepth = {'max_depth':range(5,21,1)}
gsearch2 = GridSearchCV(estimator = RandomForestClassifier(n_estimators = Ramos_best_estimator,n_jobs=-1),param_grid = maxdepth,
                        scoring = 'accuracy',cv = 10,error_score='raise')
gsearch2.fit(Ramos_X_train, Ramos_y_train)
print(gsearch2.best_score_)
print(gsearch2.best_estimator_)
Ramos_numbers = re.findall(r'\d+|\-\d+', str(gsearch2.best_estimator_))
Ramos_best_depth = [int(num) if num.lstrip('-').isdigit() else float(num) for num in Ramos_numbers][0]

minsamples = {'min_samples_split':range(2,101,2)}
gsearch3 = GridSearchCV(estimator = RandomForestClassifier(n_estimators=Ramos_best_estimator, max_depth=Ramos_best_depth,n_jobs=-1),param_grid = minsamples,
                        scoring = 'accuracy',cv = 10, error_score='raise')
gsearch3.fit(Ramos_X_train,Ramos_y_train)
print(gsearch3.best_score_)
print(gsearch3.best_estimator_)
Ramos_numbers = re.findall(r'\d+|\-\d+', str(gsearch3.cv_results_['params'][gsearch3.best_index_]))
Ramos_best_sample = int([int(num) if num.lstrip('-').isdigit() else float(num) for num in Ramos_numbers][0])


Ramos_model = RandomForestClassifier(n_estimators=Ramos_best_estimator, max_depth=Ramos_best_depth, min_samples_split=Ramos_best_sample, random_state=seed_value,n_jobs=-1)
Ramos_model.fit(Ramos_X_train, Ramos_y_train)
Ramos_y_train_pred = Ramos_model.predict(Ramos_X_train)
Ramos_y_test_pred = Ramos_model.predict(Ramos_X_test)
Ramos_score_train = Ramos_model.score(Ramos_X_test, Ramos_y_test)
Ramos_score_test = accuracy_score(Ramos_y_test_pred, Ramos_y_test)
y_pred_proba = Ramos_model.predict_proba(Ramos_X_test)[:, 1]  
#Ramos_Auc = roc_auc_score(Ramos_y_test, y_pred_proba)
#print(f"Ramos ROC AUC: {Ramos_Auc:.3f}")
#confusion_matrix(Ramos_y_test, Ramos_y_test_pred)
#classification_report(Ramos_y_test, Ramos_y_test_pred)
Ramos_importances = Ramos_model.feature_importances_
Gene_labels = Ramos_NPC_marker.columns[0:]
indices = np.argsort(Ramos_importances)[::-1] # 下标排序
Ramos_gene_importance = [(f + 1, Gene_labels[indices[f]], Ramos_importances[indices[f]])
                          for f in range(len(Gene_labels))]
Ramos_gene_importance = pd.DataFrame(Ramos_gene_importance, columns=['Rank', 'Gene', 'Importance'])


Glia_NPC_marker = pd.concat([Zhu_NPC_marker, cameron_NPC_marker, Velmeshev_NPC_marker,
                             Herring_NPC_marker, Trevino_NPC_marker, Ramos_NPC_marker], ignore_index=True)
Glia_celltype = np.concatenate((Zhu_celltype, cameron_celltype, Velmeshev_celltype,
                                Herring_celltype, Trevino_celltype, Ramos_celltype))
seed_value = 12345      
test_size = 0.3         
Glia_X_train, Glia_X_test, Glia_y_train, Glia_y_test = train_test_split(Glia_NPC_marker, Glia_celltype, test_size=test_size, random_state=seed_value)
label_encoder = LabelEncoder()
Glia_y_train = label_encoder.fit_transform(Glia_y_train)
Glia_y_test = label_encoder.fit_transform(Glia_y_test)

param_test1 = {"n_estimators":range(1000,10001,500)}
gsearch1 = GridSearchCV(estimator=RandomForestClassifier(n_jobs=-1),param_grid=param_test1,
                        scoring='accuracy',cv=10,error_score='raise')
gsearch1.fit(Glia_X_train, Glia_y_train)
print(gsearch1.best_score_)
print(gsearch1.best_estimator_)
print("best accuracy:%f" % gsearch1.best_score_)
print(gsearch1.cv_results_)
Glia_numbers = re.findall(r'\d+|\-\d+', str(gsearch1.best_estimator_))
Glia_best_estimator = [int(num) if num.lstrip('-').isdigit() else float(num) for num in Glia_numbers][0]

maxdepth = {'max_depth':range(5,21,1)}
gsearch2 = GridSearchCV(estimator = RandomForestClassifier(n_estimators = Glia_best_estimator,n_jobs=-1),param_grid = maxdepth,
                        scoring = 'accuracy',cv = 10,error_score='raise')
gsearch2.fit(Glia_X_train, Glia_y_train)
print(gsearch2.best_score_)
print(gsearch2.best_estimator_)
Glia_numbers = re.findall(r'\d+|\-\d+', str(gsearch2.best_estimator_))
Glia_best_depth = [int(num) if num.lstrip('-').isdigit() else float(num) for num in Glia_numbers][0]

minsamples = {'min_samples_split':range(2,101,2)}
gsearch3 = GridSearchCV(estimator = RandomForestClassifier(n_estimators=Glia_best_estimator, max_depth=Glia_best_depth,n_jobs=-1),param_grid = minsamples,
                        scoring = 'accuracy',cv = 10, error_score='raise')
gsearch3.fit(Glia_X_train,Glia_y_train)
print(gsearch3.best_score_)
print(gsearch3.best_estimator_)
Glia_numbers = re.findall(r'\d+|\-\d+', str(gsearch3.cv_results_['params'][gsearch3.best_index_]))
Glia_best_sample = int([int(num) if num.lstrip('-').isdigit() else float(num) for num in Glia_numbers][0])

Glia_model = RandomForestClassifier(n_estimators=Glia_best_estimator, max_depth=Glia_best_depth, min_samples_split=Glia_best_sample, random_state=seed_value)
Glia_model.fit(Glia_X_train, Glia_y_train)
Glia_y_train_pred = Glia_model.predict(Glia_X_train)
Glia_y_test_pred = Glia_model.predict(Glia_X_test)
Glia_score_train = Glia_model.score(Glia_X_test, Glia_y_test)
Glia_score_test = accuracy_score(Glia_y_test_pred, Glia_y_test)
y_pred_proba = Glia_model.predict_proba(Glia_X_test)[:, 1]  
#Glia_Auc = roc_auc_score(Glia_y_test, y_pred_proba)
#print(f"Glia ROC AUC: {Glia_Auc:.3f}")
Glia_importances = Glia_model.feature_importances_
Gene_labels = Glia_NPC_marker.columns[0:]
indices = np.argsort(Glia_importances)[::-1]
Glia_gene_importance  = [(f + 1, Gene_labels[indices[f]], Glia_importances[indices[f]])
                          for f in range(len(Gene_labels))]
Glia_gene_importance = pd.DataFrame(Glia_gene_importance, columns=['Rank', 'Gene', 'Importance'])

Glia_gene_importance.to_csv('/gpfs/hpc/home/chenchao/ranm/SCvsSN/SN_Datasets/RandomForest/Merge_Gene_importance_linear200.txt', sep='\t', index=False)
dill.dump_session('/gpfs/hpc/home/chenchao/ranm/SCvsSN/SN_Datasets/RandomForest/RandomForest_6Dataset_3Celltype_linear200.pkl')

dataframes = [
    pd.DataFrame(Zhu_gene_importance), 
    pd.DataFrame(cameron_gene_importance),
    pd.DataFrame(Velmeshev_gene_importance),
    pd.DataFrame(Herring_gene_importance),
    pd.DataFrame(Trevino_gene_importance),
    pd.DataFrame(Ramos_gene_importance)
]

with pd.ExcelWriter('/gpfs/hpc/home/chenchao/ranm/SCvsSN/SN_Datasets/RandomForest/Gene_importance_3Celltype_233gene.xlsx') as writer:
    for i, df in enumerate(dataframes):
        df.to_excel(writer, sheet_name=f'Sheet{i+1}', index=False)
