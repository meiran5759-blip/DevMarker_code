import torch
import astrorfnet
import torch
import dill
import os
import scanpy as sc
import scvi
import re
import pandas as pd
import numpy as np
import anndata as ad
import matplotlib.pyplot as plt
from collections import Counter
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

Diff_NPC_marker = pd.concat([Zhu_NPC_marker, cameron_NPC_marker, Velmeshev_NPC_marker,
                             Herring_NPC_marker, Trevino_NPC_marker, Ramos_NPC_marker], ignore_index=True)
Diff_celltype = np.concatenate((Zhu_celltype, cameron_celltype, Velmeshev_celltype,
                                Herring_celltype, Trevino_celltype, Ramos_celltype))
Mean_importance= pd.read_excel('/gpfs/hpc/home/chenchao/ranm/SCvsSN/SN_Datasets/RandomForest/Gene_importance_3Celltype_233gene_mean.xlsx', header=0, index_col=None)
RF_marker = Mean_importance.Gene[0:50].values
Diff_NPC_marker = Diff_NPC_marker.iloc[:,Diff_NPC_marker.columns.isin(RF_marker)]
Diff_NPC_marker_Torch = Diff_NPC_marker.to_numpy()
Diff_NPC_marker_Torch = torch.from_numpy(Diff_NPC_marker_Torch).float()  
label_encoder = LabelEncoder()
Diff_celltype_Torch = label_encoder.fit_transform(Diff_celltype)
Diff_celltype_Torch = torch.tensor(Diff_celltype_Torch, dtype=torch.int64)    
seed_value = 12345    
test_size = 0.3         
Diff_X_train, Diff_X_test, Diff_y_train, Diff_y_test = train_test_split(Diff_NPC_marker_Torch, Diff_celltype_Torch, test_size=test_size, random_state=seed_value)

class NeuralNet(nn.Module):
    def __init__(self, input_size, hidden_size1, hidden_size2, hidden_size3,output_size):
        super(NeuralNet, self).__init__()
        self.fc1 = nn.Linear(input_size, hidden_size1)
        self.fc2 = nn.Linear(hidden_size1, hidden_size2)
        self.fc3 = nn.Linear(hidden_size2, hidden_size3)
        self.fc4 = nn.Linear(hidden_size3, output_size)
    def forward(self, x):
        x = F.relu(self.fc1(x))  
        x = F.relu(self.fc2(x)) 
        x = F.relu(self.fc3(x))  
        x = self.fc4(x)        
        x = F.softmax(x, dim=1)  
        return x

net = NeuralNetClassifier(
    NeuralNet,
    iterator_train__shuffle=True,
    device='cpu',
    optimizer=None
)
params = {
    'lr': [0.1],
    'max_epochs': range(1000,1501,100),
    'module__input_size': [10],
    'module__hidden_size1': [400,350,300],
    'module__hidden_size2': [200,150,100],
    'module__hidden_size3': [100,50,25],
    'module__output_size': [3],
    'optimizer': [torch.optim.SGD, torch.optim.Adam]
}
gs_50 = GridSearchCV(net, params, refit=False, cv=10, 
                  scoring='accuracy',n_jobs=-1,error_score='raise')
gs_50.fit(Diff_X_train, Diff_y_train)
print(gs_50.best_score_)
print(gs_50.best_params_)
pd.DataFrame(gs_50.cv_results_).to_csv('/gpfs/hpc/home/chenchao/ranm/SCvsSN/SN_Datasets/RandomForest/gridsearch_50.csv', index=False)
net = NeuralNetClassifier(
    NeuralNet,
    iterator_train__shuffle=True,
    device='cpu',
    optimizer=torch.optim.SGD,
    lr = 0.1,
    max_epochs = 1800,
    module__input_size = 50,
    module__hidden_size1 = 400,
    module__hidden_size2 =200,
    module__hidden_size3 = 50,
    module__output_size = 3,
    batch_size = 500
)
net.fit(Diff_X_train, Diff_y_train)
