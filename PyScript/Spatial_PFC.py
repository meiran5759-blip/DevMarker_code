import squidpy as sq
import scanpy as sc
import scvi
import pandas as pd
import numpy as np
import anndata
import matplotlib.pyplot as plt
import seaborn as sns
import scrublet as scr
import scanpy.external as sce
import plotly


PFC_GW24 = sq.read.vizgen(path="/gpfs/hpc/home/chenchao/ranm/SCvsSN/Spatial/Merfish",
                          counts_file="cell_by_gene.csv",
                          meta_file="cell_metadata.csv",
                          transformation_file="micron_to_mosaic_pixel_transform.csv",
                          library_id="spatial")
PFC_GW24 = PFC_GW24[:, [gene for gene in PFC_GW24.var_names if "Blank" not in gene]].copy()
all_metadata = pd.read_excel('/gpfs/hpc/home/chenchao/ranm/SCvsSN/Spatial/Merfish/allcell_metadata.xlsx', index_col=0)
PFC_meta = all_metadata[all_metadata.Sample_ID.isin(['ARKFrozen-62-PFC'])]
#PFC_meta = PFC_meta.set_index("Cell_ID")
PFC_GW24 = PFC_GW24[PFC_GW24.obs_names.isin(PFC_meta.index)]
PFC_meta = PFC_meta.rename(columns={"Subclass": "Celltype"})
columns_to_add = ['Sample_ID', 'Volume (µm3)', 'Celltype','Niche','nCount_Vizgen', 'nFeature_Vizgen']  # 替换为实际的列名
PFC_GW24.obs = PFC_GW24.obs.join(PFC_meta[columns_to_add])
replacement_rules = pd.DataFrame({'Original': table(PFC_GW24.obs.Celltype).element,
                                  'Replacement': ['IN','EN','Ast','OPC','IPC-EN','RG','IPC-Glia','Vasc','Mic','Oligo']})
PFC_GW24.obs['Celltype'] = PFC_GW24.obs.Celltype.replace(dict(zip(replacement_rules['Original'], replacement_rules['Replacement'])))
PFC_GW24 = PFC_GW24[~PFC_GW24.obs.Celltype.isin(['Vasc','Oligo'])]
PFC_GW24.obs["Celltype"] = pd.Categorical(PFC_GW24.obs["Celltype"], categories=["RG", "IPC-EN", "EN", "IN", "IPC-Glia", "OPC", "Ast", "Mic"], ordered=True)
PFC_GW24.obs["Niche"] = PFC_GW24.obs["Niche"].replace({
    "Dorsal lateral ganglionic eminence": "DLGE",
    "Intermediate zone": "IZ",
    "Ventricular zone and subventricular zone": "VZ_SVZ",
    "Cortical layer 6 and subplate": "L6_SP",
    "Cortical layer 4": "L4",
    "Cortical upper layer in development": "UL",
    "Cortical layer 5": "L5",
    "Cortical layer 1": "L1",
    "White matter and Meninge": "WM_M",
    "Cortical layer 2 and 3": "L2_L3"
})
PFC_GW24.obs["Niche"] = pd.Categorical(PFC_GW24.obs["Niche"], categories=["L1", "L2_L3", "L4", "L5", "L6_SP", "UL", "VZ_SVZ", "IZ", "DLGE", "WM_M"], ordered=True)
celltype_color = {'IN':'#FFAA1D',
  'EN':'#00AD43',
  'IPC-Glia':'#61bada',
  'Ast':'#da6f6d',
  'OPC':'#20B2AA',
  'IPC-EN':'#96C8A2',
  'RG':'#FF6700',
  'Mic':'#C154C1'}
Niche_color = {
  "L1": "#00FFFF", 
  "L2_L3": "#d25774",  
  "L4": "#2ca02c", 
  "L5": "#df5734",  
  "L6_SP": "#D99A6C", 
  "UL": "#d8a0c0", 
  "VZ_SVZ": "#c8c7e1",  
  "IZ": "#9d3b62",  
  "DLGE": "#A57164", 
  "WM_M": "#c4daec"  
}
celltype_color = [celltype_color[ct] for ct in np.sort(PFC_GW24.obs["Celltype"].unique())]
celltype_color = mcolors.ListedColormap(celltype_color)
Niche_color = [Niche_color[ct] for ct in np.sort(PFC_GW24.obs["Niche"].unique())]
Niche_color = mcolors.ListedColormap(Niche_color)
sc.pl.embedding(
    PFC_GW24,
    basis='spatial', 
    color='Celltype',  
    palette=celltype_color, 
    size=5  
)
plt.savefig('/gpfs/hpc/home/chenchao/ranm/SCvsSN/Spatial/Figure/Scatter_celltype.pdf', bbox_inches='tight')
sc.pl.embedding(
    PFC_GW24,
    basis='spatial', 
    color='Niche',  
    palette=Niche_color,  
    size=5  
)
plt.savefig('/gpfs/hpc/home/chenchao/ranm/SCvsSN/Spatial/Figure/Scatter_Niche.pdf', bbox_inches='tight')
