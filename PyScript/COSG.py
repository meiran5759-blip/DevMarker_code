import scanpy as sc
import scvi
import pandas as pd
import numpy as np
import anndata
import matplotlib.pyplot as plt
import cosg as cosg


merge_sn = sc.read_h5ad("/gpfs/hpc/home/chenchao/ranm/SCvsSN/SN_Datasets/merge_sn_count.h5ad")
meta = pd.read_csv("/gpfs/hpc/home/chenchao/ranm/SCvsSN/SN_Datasets/merge_sn/merge_sn_count_obs.csv", index_col=0)
merge_sn = merge_sn[merge_sn.obs.index.isin(meta.index)]
merge_sn.obs['Celltype'] = merge_sn.obs['Celltype'].astype('category')
merge_sn.obs = meta.loc[merge_sn.obs_names]
Ramos_cell = merge_sn[merge_sn.obs['Dataset'] == 'Ramos']
Zhu_cell = merge_sn[merge_sn.obs['Dataset'] == 'Zhu']
cameron_cell = merge_sn[merge_sn.obs['Dataset'] == 'cameron']
Velmeshev_cell = merge_sn[merge_sn.obs['Dataset'] == 'Velmeshev']
Herring_cell = merge_sn[merge_sn.obs['Dataset'] == 'Herring']
Trevino_cell = merge_sn[merge_sn.obs['Dataset'] == 'Trevino']
cosg.cosg(Zhu_cell,
          key_added='cosg',
          mu=1,
          n_genes_user=100,
          groupby='Celltype')
cosg.cosg(cameron_cell,
          key_added='cosg',
          mu=1,
          n_genes_user=100,
          groupby='Celltype')
cosg.cosg(Velmeshev_cell,
          key_added='cosg',
          mu=1,
          n_genes_user=100,
          groupby='Celltype')
cosg.cosg(Herring_cell,
          key_added='cosg',
          mu=1,
          n_genes_user=100,
          groupby='Celltype')           
cosg.cosg(Trevino_cell,
          key_added='cosg',
          mu=1,
          n_genes_user=100,
          groupby='Celltype') 
cosg.cosg(Ramos_cell,
          key_added='cosg',
          mu=1,
          n_genes_user=100,
          groupby='Celltype') 
Ramos_cosgtop100 = Ramos_cell.uns['cosg']['names']
Ramos_cosgtop100 = pd.DataFrame(Ramos_cosgtop100)
Ramos_cosgtop100.to_csv("/gpfs/hpc/home/chenchao/ranm/SCvsSN/SN_Datasets/marker_identify/COSG_top100/Ramos_cosg_top100.txt", sep='\t')
Zhu_cosgtop100 = Zhu_cell.uns['cosg']['names']
Zhu_cosgtop100 = pd.DataFrame(Zhu_cosgtop100)
Zhu_cosgtop100.to_csv("/gpfs/hpc/home/chenchao/ranm/SCvsSN/SN_Datasets/marker_identify/COSG_top100/Zhu_cosg_top100.txt", sep='\t')
cameron_cosgtop100 = cameron_cell.uns['cosg']['names']
cameron_cosgtop100 = pd.DataFrame(cameron_cosgtop100)
cameron_cosgtop100.to_csv("/gpfs/hpc/home/chenchao/ranm/SCvsSN/SN_Datasets/marker_identify/COSG_top100/cameron_cosg_top100.txt", sep='\t')
Velmeshev_cosgtop100 = Velmeshev_cell.uns['cosg']['names']
Velmeshev_cosgtop100 = pd.DataFrame(Velmeshev_cosgtop100)
Velmeshev_cosgtop100.to_csv("/gpfs/hpc/home/chenchao/ranm/SCvsSN/SN_Datasets/marker_identify/COSG_top100/Velmeshev_cosg_top100.txt", sep='\t')
Herring_cosgtop100 = Herring_cell.uns['cosg']['names']
Herring_cosgtop100 = pd.DataFrame(Herring_cosgtop100)
Herring_cosgtop100.to_csv("/gpfs/hpc/home/chenchao/ranm/SCvsSN/SN_Datasets/marker_identify/COSG_top100/Herring_cosg_top100.txt", sep='\t')
Trevino_cosgtop100 = Trevino_cell.uns['cosg']['names']
Trevino_cosgtop100 = pd.DataFrame(Trevino_cosgtop100)
Trevino_cosgtop100.to_csv("/gpfs/hpc/home/chenchao/ranm/SCvsSN/SN_Datasets/marker_identify/COSG_top100/Trevino_cosg_top100.txt", sep='\t')
