"""A unified set of tools for setting up general organic polymer systems for MD within the OpenFF framework"""

# Add imports here
from ._version import __version__
from .genutils import importutils
# from .genutils.importutils import register_submodules, module_by_pkg_str

# _MODULE_SELF = module_by_pkg_str(__package__) # keep reference to own module
# register_submodules(_MODULE_SELF, recursive=True, blacklist=['analysis'])