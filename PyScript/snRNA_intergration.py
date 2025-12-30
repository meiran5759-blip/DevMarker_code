import dill
import os
import scanpy as sc
import scvi
import pandas as pd
import numpy as np
import anndata
import matplotlib.pyplot as plt
import seaborn as sns
import scrublet as scr
import scanpy.external as sce
import loompy
from collections import Counter
import cellhint

####Merge six datasets and retain the original celltype annotations for each dataset
velmeshev_counts = sc.read_10x_mtx('/gpfs/hpc/home/chenchao/ranm/2nd_3rd_exp/velmeshev/', var_names='gene_symbols', cache=True)
meta = pd.read_csv(r'/gpfs/hpc/home/chenchao/ranm/2nd_3rd_exp/metadata.tsv',delimiter='\t')
indices_fetal = meta[meta['age'].str.contains('GW|ga')].index
meta = meta.iloc[indices_fetal]
replacement_rules = pd.DataFrame({'Original': table(meta.age).element,
                                  'Replacement': [33,37,18,17,21,20,32,36,35,23,26,25,16,26,18,23,28,14,15,22,19,39,31,30,20,22,32]})
meta['PCW'] = meta['age'].replace(dict(zip(replacement_rules['Original'], replacement_rules['Replacement'])))


velmeshev_counts = velmeshev_counts[indices_fetal, :]
velmeshev_counts.obs['Celltype'] = meta.lineage.to_numpy()
velmeshev_counts.obs['Dataset'] = meta.dataset.to_numpy()
velmeshev_counts.obs['PCW'] = meta.PCW.to_numpy()
velmeshev_counts.obs['chemistry'] = meta.chemistry.to_numpy()
velmeshev_counts.obs['individual'] = meta['individual'].to_numpy()
del_celltype_index = [5,6,8]
del_celltype = velmeshev_counts.obs.Celltype.unique()[del_celltype_index]
velmeshev_counts = velmeshev_counts[~velmeshev_counts.obs['Celltype'].isin(del_celltype)]
velmeshev_counts = velmeshev_counts[:, ~velmeshev_counts.var_names.duplicated()]
velmeshev_counts.var.rename(columns={'gene_ids': 'Gene'}, inplace=True)

####2022_Ramos
directory_path = "/gpfs/hpc/home/chenchao/ranm/SCvsSN/SN_Datasets/2022_Ramos/metadata"
csv_files = [os.path.join(directory_path, file) for file in os.listdir(directory_path) if file.endswith('.csv')]
Ramos_metadata = pd.concat([pd.read_csv(file) for file in csv_files], axis=0)
Ramos_metadata.columns = ['CellID'] + list(Ramos_metadata.columns[1:]) ####修改第一列的列名
Ramos_metadata = Ramos_metadata.iloc[:, [0, 1, 11]]

Ramos_metadata['CellID_first_part'] = Ramos_metadata['CellID'].str.split('_').str[0]
Ramos_metadata['CellID'] = 'Ramos_' + Ramos_metadata['sample'] + '_' + Ramos_metadata['CellID_first_part']
Ramos_metadata.drop(columns=['CellID_first_part'], inplace=True)
Ramos_metadata = Ramos_metadata.drop_duplicates(subset=['CellID'])
Ramos_cellID = velmeshev_counts.obs[velmeshev_counts.obs['Dataset'] == 'Ramos'].index
positions = []
for ID in Ramos_cellID:
    index = int(np.where(Ramos_metadata['CellID'] == ID)[0])
    positions.append(index)

##positions = np.where(Ramos_metadata['CellID'].isin(Ramos_cellID[]))[0]
Ramos_metadata = Ramos_metadata.iloc[positions,:]
Ramos_index = np.where(velmeshev_counts.obs['Dataset'] == 'Ramos')[0]
replacement_rules = pd.DataFrame({'Original': table(Ramos_metadata['celltypes']).element,
                                  'Replacement': ['EN','EN','IN','IPC','EN','Ast','Unknow','EN','gIPC','OPC','CycPro','EN','Mic','VSMC','IPC','gIPC','gIPC','EN','RG','OPC','RG','RG','EPD','EN','IN']})
Ramos_metadata['Celltype'] = Ramos_metadata['celltypes'].replace(dict(zip(replacement_rules['Original'], replacement_rules['Replacement'])))

####2022_Herring
Herring_cellID = velmeshev_counts.obs[velmeshev_counts.obs['Dataset'] == 'Herring'].index
Herring_metadata = pd.read_csv('/gpfs/hpc/home/chenchao/ranm/SCvsSN/SN_Datasets/2022_Herring/metadata/RNA-all_BCs-meta-data.csv')
Herring_metadata['cell_part1'] = Herring_metadata['cell'].str.split('-').str[0]
Herring_metadata['cell'] = 'Herring_' + Herring_metadata['donor'] + '_' + Herring_metadata['cell_part1'] + '-1'
Herring_metadata.drop(columns=['cell_part1'], inplace=True)
Herring_metadata = Herring_metadata[Herring_metadata.donor.isin(['RL2103','RL2107','RL2121'])]
positions = []
for ID in Herring_cellID:
    index = int(np.where(Herring_metadata['cell'] == ID)[0])
    positions.append(index)

Herring_metadata = Herring_metadata.iloc[positions,:]
Herring_index = np.where(velmeshev_counts.obs['Dataset'] == 'Herring')[0]
replacement_rules = pd.DataFrame({'Original': table(Herring_metadata['major_clust']).element,
                                  'Replacement': ['EN','IN','EN','Ast','EN','EN','IN','EN','IN','OPC','Mic','IN','IN','Olig','Poor-Quality','IN','IN','IN']})
Herring_metadata['Celltype'] = Herring_metadata['major_clust'].replace(dict(zip(replacement_rules['Original'], replacement_rules['Replacement'])))

####2021_Trevino
Trevino_cellID = velmeshev_counts.obs[velmeshev_counts.obs['Dataset'] == 'Trevino'].index
Trevino_metadata = pd.read_csv('/gpfs/hpc/home/chenchao/ranm/SCvsSN/SN_Datasets/2021_Trevino/GSE162170_multiome_cell_metadata.txt', sep='\t')
Trevino_metadata['cell'] = 'Trevino_' + Trevino_metadata['cell']
positions = []
for ID in Trevino_cellID:
    index = int(np.where(Trevino_metadata['cell'] == ID)[0])
    positions.append(index)

Trevino_metadata = Trevino_metadata.iloc[positions,:]
Trevino_index = np.where(velmeshev_counts.obs['Dataset'] == 'Trevino')[0]
replacement_rules = pd.DataFrame({'Original': table(Trevino_metadata['celltype']).element,
                                  'Replacement': ['EN','IN','EN','IN','EN','EN','IN','RG','EN','EN','CycPro','gIPC']})
Trevino_metadata['Celltype'] = Trevino_metadata['celltype'].replace(dict(zip(replacement_rules['Original'], replacement_rules['Replacement'])))

####2023_Cameron
metadata_path = '/gpfs/hpc/home/chenchao/ranm/SCvsSN/SN_Datasets/2023_Cameron/cameron_5region_metadata'
txt_files = [f for f in os.listdir(metadata_path) if f.endswith('.txt')]
data_frames = []
for txt_file in txt_files:
    file_path = os.path.join(metadata_path, txt_file)
    df = pd.read_csv(file_path, sep='\t')  # 假设文件以制表符分隔
    data_frames.append(df)

metadata_cameron = pd.concat(data_frames, ignore_index=True)
print(f"Combined DataFrame shape: {metadata_cameron.shape}")
metadata_cameron['chemistry'] = 'V3'
metadata_cameron['individual'] = metadata_cameron['sample'].str.split('_').str[0]
metadata_cameron['Celltype'] = metadata_cameron['cellIDs'].str.split('-').str[1]
#metadata_cameron = metadata_cameron[~metadata_cameron['Celltype'].isin(['N', 'Endo'])]
print(metadata_cameron.shape)

count_path = '/gpfs/hpc/home/chenchao/ranm/SCvsSN/SN_Datasets/2023_Cameron/cameron_5region_count'
txt_files = [f for f in os.listdir(count_path) if f.endswith('.txt')]
data_frames = []
for txt_file in txt_files:
    file_path = os.path.join(count_path, txt_file)
    df = pd.read_csv(file_path, sep=' ', index_col=0)  # 假设文件以制表符分隔，第一列是基因名
    data_frames.append(df)
    print(f"{txt_file}: {df.shape[0]} rows, {df.shape[1]} columns")

merge_count = data_frames[0]
for df in data_frames[1:]:
    merge_count = merge_count.join(df, how='inner')

merge_count_transposed = merge_count.T
metadata_cameron.set_index('cells', inplace=True)
metadata_cameron = metadata_cameron.sort_index()
merge_count_transposed = merge_count_transposed.sort_index()
cameron_adata = sc.AnnData(merge_count_transposed, obs=metadata_cameron)
dataset_values = ['Cameron'] * len(cameron_adata.obs)
cameron_adata.obs['Dataset'] = dataset_values
cameron_adata.obs = cameron_adata.obs.drop(columns=['cellIDs','sample'])
PCW_values = ['14'] * len(cameron_adata.obs)
cameron_adata.obs['PCW'] = PCW_values
cameron_adata = cameron_adata[~cameron_adata.obs['Celltype'].isin(['CR', 'N'])]
cameron_adata = cameron_adata[~cameron_adata.obs_names.str.contains('Cer')]
replacement_rules = pd.DataFrame({'Original': table(cameron_adata.obs['Celltype']).element,
                                  'Replacement': ['RG','EN','IN','Mic','CycPro','Endo','OPC','IPC']})
cameron_adata.obs['Celltype'] = cameron_adata.obs['Celltype'].replace(dict(zip(replacement_rules['Original'], replacement_rules['Replacement'])))

####2023_Zhu
folder = '/gpfs/hpc/home/chenchao/ranm/SCvsSN/SN_Datasets/2023_Zhu/ZHU'
Zhu_2023 = sc.read_10x_mtx(folder, var_names='gene_symbols')
Zhu_metadata = pd.read_csv('/gpfs/hpc/home/chenchao/ranm/SCvsSN/SN_Datasets/2023_Zhu/scp_metadata.csv')
Zhu_metadata = Zhu_metadata.drop(columns=['age_group', 'nCount_RNA','nFeature_RNA', 'nCount_ATAC', 'nFeature_ATAC', 'TSS_percentile',
                                  'nucleosome_signal', 'percent_mt', 'biosample_id', 'disease','disease__ontology_label', 'library_preparation_protocol',
                                  'library_preparation_protocol__ontology_label', 'organ','organ__ontology_label', 'sex', 'species', 'species__ontology_label'])
Zhu_metadata = Zhu_metadata.drop(index=Zhu_metadata.index[0])
Zhu_metadata.index = Zhu_metadata['NAME']
Zhu_metadata = Zhu_metadata.drop(columns='NAME')
Zhu_metadata['age'] = Zhu_metadata['age'].astype(str)
Zhu_metadata = Zhu_metadata[Zhu_metadata['age'].str.contains('gw')]
replacement_rules = pd.DataFrame({'Original': table(Zhu_metadata['age']).element,
                                  'Replacement': [21,22,17,16]})
Zhu_metadata['PCW'] = Zhu_metadata['age'].replace(dict(zip(replacement_rules['Original'], replacement_rules['Replacement'])))
Zhu_2023 = Zhu_2023[:len(Zhu_metadata.index), :]
Zhu_2023.obs['chemistry'] = '10x multiome'
Zhu_2023.obs['Dataset'] = 'Zhu'
Zhu_2023.obs['PCW'] = Zhu_metadata.PCW.to_numpy()
Zhu_2023.obs['individual'] = Zhu_metadata.donor_id.to_numpy()
Zhu_2023.obs['Celltype'] = Zhu_metadata.celltype.to_numpy()
Zhu_2023.var.rename(columns={'gene_ids': 'Gene'}, inplace=True)
Zhu_2023 = Zhu_2023[~Zhu_2023.obs['Celltype'].isin(['VSMC', 'Pericytes', 'EN'])]
replacement_rules = pd.DataFrame({'Original': table(Zhu_2023.obs['Celltype']).element,
                                  'Replacement': ['EN','IN','Endo','RG','OPC','IN','Mic','IN','EN','Ast','IPC']})
Zhu_2023.obs['Celltype'] = Zhu_2023.obs['Celltype'].replace(dict(zip(replacement_rules['Original'], replacement_rules['Replacement'])))

common_genes = set(cameron_adata.var_names) & set(velmeshev_counts.var_names) & set(Zhu_2023.var_names)
Zhu_2023= Zhu_2023[:, list(common_genes)]
velmeshev_counts = velmeshev_counts[:, list(common_genes)]
cameron_adata = cameron_adata[:, list(common_genes)]

merge_sn_count = anndata.concat([Zhu_2023, cameron_adata, velmeshev_counts], join='outer', index_unique='_')
merge_sn_count.obs['chemistry'] = merge_sn_count.obs['chemistry'].replace({"10x multiome": "multiome"})
merge_sn_count.obs.index = [name[:-2] for name in merge_sn_count.obs.index]
merge_sn_count.obs['Celltype'] = merge_sn_count.obs['Celltype'].astype(str)
merge_sn_count.obs['PCW'] = merge_sn_count.obs['PCW'].astype(str)
merge_sn_count.obs['individual'] = merge_sn_count.obs['individual'].astype(str)
Ramos_index = np.where(merge_sn_count.obs['Dataset'] == 'Ramos')[0]
Trevino_index = np.where(merge_sn_count.obs['Dataset'] == 'Trevino')[0]
Herring_index = np.where(merge_sn_count.obs['Dataset'] == 'Herring')[0]
merge_sn_count.obs['Celltype'][Herring_index] = Herring_metadata['Celltype']
merge_sn_count.obs['Celltype'][Trevino_index] = Trevino_metadata['Celltype']
merge_sn_count.obs['Celltype'][Ramos_index] = Ramos_metadata['Celltype']
celltypes_to_remove = ['VSMC', 'Unknow', 'EPD', 'Olig', 'Poor-Quality','Endo']
merge_sn_count = merge_sn_count[~merge_sn_count.obs['Celltype'].isin(celltypes_to_remove)].copy()
sc.pp.filter_cells(merge_sn_count,min_genes=500)
sc.pp.filter_genes(merge_sn_count,min_cells=3)
merge_sn_count.var['mt'] = merge_sn_count.var_names.str.startswith('MT-') 
merge_sn_count.var['hb'] = merge_sn_count.var_names.str.contains('^HB[^P]')
merge_sn_count.var['ribo'] = merge_sn_count.var_names.str.startswith('RPS','RPL')
sc.pp.calculate_qc_metrics(merge_sn_count, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
merge_sn_count = merge_sn_count[merge_sn_count.obs['pct_counts_mt'] < 10, :]
merge_sn_count = merge_sn_count[merge_sn_count.obs['n_genes_by_counts'] < 6000, :]
merge_sn_count = merge_sn_count[merge_sn_count.obs['n_genes_by_counts'] > 500,:]

replacement_rules = pd.DataFrame({'Original': table(merge_sn_count.obs['Celltype']).element,
                                  'Replacement': ['EN','IN','RG','OPC','Mic','Ast','IPC','CycPro','EN','Ast','Mic','GlioProg','gIPC']})
merge_sn_count.obs['Celltype'] = merge_sn_count.obs['Celltype'].replace(dict(zip(replacement_rules['Original'], replacement_rules['Replacement'])))
merge_sn_count.write_h5ad("/gpfs/hpc/home/chenchao/ranm/SCvsSN/SN_Datasets/merge_sn/merge_sn_count.h5ad")

