"""
Author: Andrew Lutsky
Date Created: Jan 2026


This is my package that was created to replicate graphs
from Nebulosa. It has the ability to generate density values
that can be attached to the observations.
"""

# Import and make available sckde
from .sckde import sckde
from .utils import nebulosa_wkde2d
from . import plotting as pl
from . import exceptions as execeptions

__version__ == "0.1.0"
__author__ == "Andrew Lutsky"


