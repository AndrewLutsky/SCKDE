"""
Tests to test nebulosa equivalent functions.

"""

import os
import pandas as pd
import pytest
from sckde.utils import nebulosa_wkde2d, nebulosa_get_density
import numpy as np
import glob

test_dir = os.path.join(os.path.dirname(__file__), "nebulosa_tests")
tests = list(glob.glob(os.path.join(test_dir, "*/")))

@pytest.mark.parametrize("test_name", tests)
def test_wkde2d(test_name):
    """ 
    Tests the wkde2d fxn using an assert. The tests
    were generated in Nebulosa, and we aim to get within 1e-6.
    Importantly we already generate the plug in bandwidth in the test
    as it is difficult to replicate as that lies in the ks::hpi module
    in R. Instead we read in the bandwidth. 
    
    Parameters
    ----------
    test_name : str
        This is the folder directory inside of nebulosa_tests.

    Returns
    -------
    None
    """

    inp = pd.read_csv(os.path.join(test_name, "input.csv"))
    out = pd.read_csv(os.path.join(test_name, "z.csv")).values
    h = pd.read_csv(os.path.join(test_name, "h.csv")).values
    _, _, density = nebulosa_wkde2d(x = inp['x'], y = inp['y'], w = inp['w'],h = h)
    assert(np.allclose(density, out))

@pytest.mark.parametrize("test_name", tests)
def test_get_density(test_name):
    
    """
    This function tests the equivalent get density function in python 
    for weighted kde.

    Parameters
    ----------
    test_name : str
        This is the folder directory inside of nebulosa_tests.
        
    Returns
    -------
    None
    """
    inp = pd.read_csv(os.path.join(test_name, "input.csv"))
    out = pd.read_csv(os.path.join(test_name, "z.csv")).values
    h = pd.read_csv(os.path.join(test_name, "h.csv")).values
    dens_neb = pd.read_csv(os.path.join(test_name, "get_dens_z.csv")).values.flatten()
    x_int, y_int, density = nebulosa_wkde2d(x = inp['x'], y = inp['y'], w = inp['w'],h = h)
    dens = nebulosa_get_density(inp['x'], inp['y'], x_int, y_int, density)
    dens = np.array(dens)
    dens_neb = np.array(dens_neb)
    print("max diff:", np.max(np.abs(dens_neb - dens)))
    print("failing indices:", np.where(~np.isclose(dens_neb, dens)))
    print("dens_neb at failing:", dens_neb[~np.isclose(dens_neb, dens)][:5])
    print("dens at failing:", dens[~np.isclose(dens_neb, dens)][:5])
    assert(np.allclose(dens_neb, dens))



