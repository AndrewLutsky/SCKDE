import scanpy as sc
from scipy.stats import gaussian_kde
from anndata import AnnData
import pandas as pd
import seaborn as sb
import numpy as np

def _sckde_uni(adata: AnnData, key:str) -> list[f64]:
    """ Function to generate gaussian kde from cells by gene expression matrix,
    and a list of relevant genes to query the anndata """
    
    return None

def _sckde_multi(adata:AnnData, keys:list[str]):#-> gaussian_kde:
    """ Function to generate gaussian kde from cells by gene expression matrix,
    and a list of relevant genes to query the anndata """
    
    # Check for keys in the adata
    for key in keys:
        if key not in adata.var.index:
            raise KeysNotFound(key)
   
    adata = adata[:, keys].copy()

    # Check for umap key in obsm
    if 'X_umap' not in adata.obsm:
        raise KeysNotFound('X_umap') 
    # Grab the first two UMAP axes
    um = adata.obsm['X_umap'][:, 0:2]

    kdes = []
    for key in keys:
        norm_const = np.sum(adata[:, key].X)
        if norm_const == 0.:
            norm_const = 1
        weights = np.array((adata[:, key].X / norm_const).toarray()).T.reshape(-1)
        kdes.append(gaussian_kde(um.T, weights =weights))
   
    # Grab the UMAP coordinates 
    um = adata.obsm['X_umap'][:, 0:2].T

    # Instantiate a density dictionary mapping key to a density
    # distribution. 
    density_map = {}
    
    # Calculate the density for each key
    for kde, key in zip(kdes, keys):
        density_map[key] = kde(um)
        print(um)
        print(np.sum(kde(um)))

    # Now let us find the joint distribution over all the keys.
    prod = np.prod(list(density_map.values()), axis = 0)
    return(prod)

    if len(keys) == 0:
        "Given keys is empty, utilizing umap-coordinate density!"
        return gaussian_kde(um.T)





# Error Class
class KeysNotFound(Exception):
    def __init__(self, key:str):
        super().__init__(f"{key} not found in the anndata!")

