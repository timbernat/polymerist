"""A unified set of tools for setting up general organic polymer systems for molecular dynamics"""

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

# from ._version import __version__
from importlib.metadata import version
__version__ = version(__name__)

# used to check whether package is installed and importable
from .genutils import importutils

from .genutils.importutils.pkgiter import module_hierarchy
from .maths.combinatorics.tables import pascal