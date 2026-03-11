"""A list of exceptions to be used in the package."""


class KeysNotFound(Exception):
    """
    An Exception thrown if specific keys are not available
    in the given AnnData layer.

    Parameters
    ----------
    key : str
        The key that was not found in the AnnData object

    Returns
    -------

    None

    Example
    -------

    >>> raise KeysNotFound(f"{key}")

    """

    def __init__(self, key):
        super().__init__(f"{key} not found in the anndata!")


class UnequalArrayLength(Exception):
    """
    An Exception thrown if the arrays are not the correct
    length.

    Parameters
    ----------
    message : str
        The accompanying message to explain to the user
        the arrays that are not matching.

    Returns
    -------

    None

    Example
    -------

    >>> Raise UnequalArrayLength("A and B don't have the same length!")

    """

    def __init__(self, message):
        super().__init__(message)


class InvalidBandwidthArray(Exception):
    """
    An Exception thrown if the bandwidth array does not match
    the expected format.

    Parameters
    ----------
    message : str
        A message that states that the bandwidth array
        does not match the expected format.

    Returns
    -------

    None

    Example
    -------
    >>> Raise InvalidBandwidthArray("H has unexpected format!")
    """

    def __init__(self, message):
        super().__init__(message)

class TooManyDimensions(Exception):
    """
    An Exception thrown if the data has too 
    many dimensions as input in the embedding
    space.

    Parameters
    ----------
    message : str
        A message that states that there are too
        many dimensions for a given input.
    Returns
    -------

    None

    Example
    -------
    >>> Raise TooManyDimensions("There are too many dimensions!")
    """

    def __init__(self, message):
        super().__init__(message)
