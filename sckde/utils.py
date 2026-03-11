import numpy as np
from KDEpy import FFTKDE
from sckde.exceptions import UnequalArrayLength, InvalidBandwidthArray
from scipy.stats import norm
from collections.abc import Iterable
import functools

def nebulosa_wkde2d(x, y, w, h=None, n=100):
    """
    Create the 2d weighted kde over the points (x, y).

    This function emulates nebulosa::wkde2d and should match
    that function exactly.

    The weighted kernel density estimation is calculated as

    ```
    f(x; H) = \frac{1}{n} \Sum_{i=1}^{n}\left w_i K_h \left x - X_i \right \right)
    ```

    where at each x value (point in d-dimensional space), the kernel is evaluated
    using plug in bandwidth.

    So we calculate the kernel


    Parameters:
    -----------
    x : list[float]
        Coordinates in the x dimension.

    y : list[float]
        Coordinates in the y dimension.
    w : list[float]
        Weights of each point.
    h : float or tuple[float], optional
        Bandwidth(s) for the x and y directions. If it is a single
        float then the bandwidth is applied in both the x and y directions.
        len(h) <= 2.
    n : int, optional
        Number of grids/voxels in each embedding dimension. In this
        case two dimensions.

    Returns
    -------

    dict with keys:
        x : np.array
            Array of x coordinates.
        y : np.array
            Array of y coordinates.

        z : np.array
            2D array of shape (len(x), len(y)) that
            contains the weighted estimated density.
            The rows correspond to various x-values
            and the columns correspond to various
            y-values.

    Examples
    --------
    # TODO
    """
    # Check to make sure the arrays are the same lengths.
    if len(x) != len(y):
        raise UnequalArrayLengths("X and Y are of unequal length!")
    elif len(x) != len(w):
        raise UnequalArrayLengths("X and W are of unequal length!")
    # We don't need to check Y and W

    # Store length.
    length = len(x)

    # Attempt to cast to numpy array.
    x = np.array(x)
    y = np.array(y)
    w = np.array(w)
    x_fit = x.reshape(-1, 1) if x.ndim == 1 else x
    y_fit = y.reshape(-1, 1) if y.ndim == 1 else y
    # Now let us compute the bandwidths for x and y if h does not exist
    # using improved sheather-jones similar to ks::hpi.
    if h is None:
        h = [FFTKDE(bw="ISJ").fit(x).bw, FFTKDE(bw="ISJ").fit(y).bw]
    elif isinstance(h, float):
        h = [h, h]
    elif isinstance(h, Iterable) and len(h) <= 2:
        h = h.reshape(-1, 1)
        pass
    else:
        raise InvalidBandwidthArray(
            "The bandwidth matrix passed was either non-iterable and not a float, or had length > 2!"
        )

    # Unlike nebulosa, we don't pass in an adjust to the bandwidth as we
    # can manually change that in h.

    # Find the range over the data
    range_data = ((x.min(), x.max()), (y.min(), y.max()))

    # Compute the intervals for the range of the data.
    x_int = np.linspace(range_data[0][0], range_data[0][1], n)
    y_int = np.linspace(range_data[1][0], range_data[1][1], n)

    # Compute the distance matrix between x_int and x as well as
    # y_int and y normalized by the axes bandwidth.
    x_sub = np.subtract.outer(x_int, x) / h[0]
    y_sub = np.subtract.outer(y_int, y) / h[1]

    # This computes an outer matrix of size (n, length). where each row
    # is the vector w.
    w_mat = np.outer(np.ones(n), w)

    # Compute the density
    x_mat = norm.pdf(x_sub) * w
    y_mat = norm.pdf(y_sub) * w
    z = (x_mat @ y_mat.T) / (n * np.sum(w) * h[0] * h[1])

    return x_int, y_int, z


def nebulosa_get_density(x, y, x_int, y_int, z) -> np.array:
    """
    Find where in the intervals all the values of x and y
    and plug those into the z matrix.

    Parameters
    ----------
    x : np.array
        The x-coordinates of all the points in the embedding.

    y : np.array
        The y-coordinates of all the points in the embedding.

    x_int : np.array
        The interval array for x.

    y_int : np.array
        The interval array for y.

    z : np.array
        The density matrix.


    Returns
    -------
    np.array
        1D array of density values for each (x,y) point.

    """
    # Find where x values lie along the x intervals.
    x_idx = (np.searchsorted(x_int, x, side="right") - 1).astype(int)

    # Find where y values lie along the y intervals.
    y_idx = (np.searchsorted(y_int, y, side="right") - 1).astype(int)

    x_idx = np.clip(x_idx, 0, len(x_int) - 1)
    y_idx = np.clip(y_idx, 0, len(y_int) - 1)
    # Return those indices of the density
    print(type(x_int), type(y_int))
    return z[x_idx, y_idx]

def _todo(func):
    """
    Wrapper/decorator for todos.
    Parameters
    ----------
    func : function
        

    Returns
    -------
    Nothing
        
    """

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        raise NotImplementedError(f"{func.__name__} is not yet implemented.")
    return wrapper
