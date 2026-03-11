"""
Author: Andrew Lutsky
Main module to compute the sckde weighted
gaussian kernel density estimate.
"""

from scipy.stats import gaussian_kde
from KDEpy import FFTKDE
from anndata import AnnData
import numpy as np
from typing import Union
from sckde.utils import _todo
from sckde.exceptions import TooManyDimensions
import logging

def sckde(adata: AnnData,
          keys: Union[str, list[str]],
          basis="X_umap",
          dimensions_embedding=2,
          scale=True
          ):
    """
    This function applies the sckde function onto the adata using specific 
    adata keys. Essentially it weights the umap coordinates according to 
    specific genes.

    Parameters
    ----------
    basis : 
        
    dimensions_embedding : 
        
    scale : 
        
    adata : AnnData
        
    keys : Union[str, list[str]]
        

    Returns
    -------
    None     

    """
    if isinstance(keys, str):
        keys = [keys]
    elif isinstance(keys, list):
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
    # Grab the first n dimension UMAP axes
    um = adata.obsm[basis][:, 0:dimensions_embedding]

    # Compute a n dimensional grid over the embedding space

    kdes = {}
    for key in keys:
        axes, z = wkde_fft(um, adata[:, key])
        density = get_density(um, axes, z)
        kdes[key] = density

    # Find product of all keys.
    prod = np.ones(len(adata))

    for key,val in kdes.items():
        prod *= val
        
    return prod

def wkde_fft(data:np.ndarray, w:np.ndarray, n:int=1000) -> np.ndarray:
    """
    This calculates the weighted kde
    given a numpy array (that is either one
    or two dimensional and a weight array.
    
    Parameters
    ----------
    data : np.ndarray

    w : np.ndarray

    n : int = 100

    bw : str = 'ISJ'

    Returns:
    np.ndarray

    """


    # Coerce data and w to be numpy arrays.
    if not isinstance(data, np.ndarray):
        data = np.ndarray(data)

    if not isinstance(w, np.ndarray):
        w = w.X.todense().flatten()
        print(w)
        w = np.array(w)

    # Ensure that data is either one or two dimensions
    if len(data.shape) > 2:
        raise TooManyDimensions("There are too many dimesions!")

    n_dims = len(data.shape)

    # Normalize weights to match Nebulosa
    w = w / w.sum() * len(w)

    # Fit weighted FFT KDE
    kde = FFTKDE(kernel="gaussian")
    grid, z = kde.fit(data, weights=w).evaluate(n)
   
    # Extract out unique grid values per axis
    axes = [np.unique(grid[:, i]) for i in range(n_dims)]

    # Reshape grid  = depending on dimension
    if n_dims < 2:
        z = z.reshape(n)
    else:
        z = z.reshape(n, n)
     
    return axes, z

def get_density(data:np.ndarray, axes:np.ndarray, z:np.ndarray) -> np.ndarray:
    """
    Function to extract out density from the results of wkde_fft.

    Parameters
    ----------
    data : np.ndarray
        Data

    axes : np.ndarray
        Axes output from wkde_fft.

    z : np.ndarray
        Z output from wkde_fft.
    
    Returns
    -------
    np.ndarray
        Returns a numpy array that is of dimension (n_points, ). 
        This returns the density value per point given the results of
        wkde_fft.

    """
    ndim = len(axes)
    
    idxs = []
    for dim in range(ndim):
        indices = (np.searchsorted(axes[dim], data[:, dim], side="right") - 1)
        indices = np.clip(indices, 0, len(axes[dim]) - 1)
        idxs.append(indices)
    
    return z[tuple(idxs)]

