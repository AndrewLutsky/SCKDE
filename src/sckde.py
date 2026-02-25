"""
Author: Andrew Lutsky
Main module to compute the sckde weighted
gaussian kernel density estimate.
"""

from scipy.stats import gaussian_kde
from anndata import AnnData
import numpy as np


def sckde(
    adata: AnnData, keys, basis="X_umap", dimensions_embedding=2, scale=True
):
    """
    Function to generate gaussian kde from cells by gene expression matrix,
    and a list of relevant genes to query the anndata.


    """
    if isinstance(keys, str):
        keys = [keys]
    elif isinstance(keys, list[str]):
        pass
    else:
        raise ValueError("Keys must be a string or a list of strings!")


    # Check for keys in the adata
    for key in keys:
        if key not in adata.var.index:
            raise KeysNotFound(key)

    
    adata = adata[:, keys].copy()

    # Check for umap key in obsm
    if basis not in adata.obsm:
        raise KeysNotFound(basis)
    # Grab the first two UMAP axes
    um = adata.obsm[basis][:, 0:dimensions_embedding]

    kdes = []
    for key in keys:
        norm_const = np.sum(adata[:, key].X)
        if norm_const == 0.0:
            norm_const = 1
        X = adata[:, key].X / norm_const
        weights = np.asarray(X).reshape(-1)
        weights = weights.astype(float)
        print(weights)
        kdes.append(gaussian_kde(um.T, weights=weights))

    # Grab the UMAP coordinates
    um = adata.obsm[basis][:, 0:dimensions_embedding].T

    if len(keys) == 0:
        print("Given keys is empty, utilizing umap-coordinate density!")
        return gaussian_kde(um.T)
    # Instantiate a density dictionary mapping key to a density
    # distribution.
    density_map = {}

    # Calculate the density for each key
    for kde, key in zip(kdes, keys):
        density_map[key] = kde(um)

    # Now let us find the joint distribution over all the keys.
    prod = np.prod(list(density_map.values()), axis=0)

    # Now let us scale the values to sum to 1.
    if scale:
        prod /= np.sum(prod)

    return prod


def sckde_trajectory(adata: AnnData, keys: list[str], trajectory_basis="trajectory"):
    """
    This is a wrapper function around sckde, specifically
    for computing the weighted gaussian kde along the
    pseudotime path.


    """
    return sckde(adata, keys, basis=trajectory_basis, dimensions_embedding=1)


# Error Class
class KeysNotFound(Exception):
    """
    An Exception thrown if specific keys are not available
    in the given AnnData layer.
    """

    def __init__(self, key: str):
        super().__init__(f"{key} not found in the anndata!")
